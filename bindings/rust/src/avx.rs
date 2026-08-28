// Copyright Supranational LLC
// Licensed under the Apache License, Version 2.0, see LICENSE for details.
// SPDX-License-Identifier: Apache-2.0

//! Expert AVX ABI bindings.
//!
//! The `b` fields below are public for raw FFI compatibility, not because the
//! byte representation is a supported serialization format. Library-produced
//! values and exact copies may be passed back to CatBLST. Inspecting or
//! modifying the bytes requires understanding the selected AVX representation;
//! arbitrary contents need not satisfy the bounds or algebraic preconditions
//! of any operation other than an explicitly documented validator.
#![doc = include_str!("../docs/avx.md")]

use crate::*;
use alloc::alloc::{alloc, dealloc};
use alloc::vec::Vec;
use core::alloc::Layout;
use core::mem::{size_of, MaybeUninit};
use core::ptr::{self, NonNull};
use core::sync::atomic::{compiler_fence, Ordering};

#[repr(C, align(64))]
#[derive(Copy, Clone)]
pub struct blst_avx_p1 {
    /// Opaque expert-ABI storage. Do not inspect or modify casually.
    pub b: [blst_byte; 384],
}

#[repr(C, align(64))]
#[derive(Copy, Clone)]
pub struct blst_avx_p1_affine {
    /// Opaque expert-ABI storage. Do not inspect or modify casually.
    pub b: [blst_byte; 256],
}

#[repr(C, align(64))]
#[derive(Copy, Clone)]
pub struct blst_avx_p2 {
    /// Opaque expert-ABI storage. Do not inspect or modify casually.
    pub b: [blst_byte; 768],
}

#[repr(C, align(64))]
#[derive(Copy, Clone)]
pub struct blst_avx_p2_affine {
    /// Opaque expert-ABI storage. Do not inspect or modify casually.
    pub b: [blst_byte; 512],
}

#[repr(C, align(64))]
#[derive(Copy, Clone)]
pub struct blst_avx_fp12 {
    /// Opaque expert-ABI storage. Do not inspect or modify casually.
    pub b: [blst_byte; 1536],
}

/// Semantic alias for AVX Fp12 storage known to contain a GT element.
pub type blst_avx_gt = blst_avx_fp12;

#[repr(C)]
pub struct blst_avx_prepared_g2_batch {
    _private: [u8; 0],
}

#[repr(C)]
pub struct blst_avx_prepared_g2_batch_workspace {
    _private: [u8; 0],
}

/// Owned, immutable native G2 Miller-line table.
///
/// Allocation and line preparation happen once during construction. Pairing
/// evaluation reuses the caller-owned native table without allocating.
pub struct NativeAvxPreparedG2Lines {
    storage: NonNull<u8>,
    layout: Layout,
}

/// Construction failure for a raw native prepared-G2 table.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NativeAvxPreparedG2LinesError {
    /// This build or host cannot execute the CatBLST AVX backend.
    Unavailable,
    /// The backend reported an invalid storage size or alignment.
    InvalidLayout,
    /// Allocation of the immutable native table failed.
    AllocationFailed,
    /// Native line preparation rejected the validated input or storage.
    InitializationFailed,
}

impl NativeAvxPreparedG2Lines {
    pub fn new(
        g2: &blst_avx_p2,
    ) -> Result<Self, NativeAvxPreparedG2LinesError> {
        if !unsafe { blst_avx_enabled() } {
            return Err(NativeAvxPreparedG2LinesError::Unavailable);
        }
        let size = unsafe { blst_avx_prepared_g2_lines_sizeof() };
        let align = unsafe { blst_avx_prepared_g2_lines_alignof() };
        let layout = Layout::from_size_align(size, align)
            .map_err(|_| NativeAvxPreparedG2LinesError::InvalidLayout)?;
        if size == 0 {
            return Err(NativeAvxPreparedG2LinesError::InvalidLayout);
        }
        let storage = NonNull::new(unsafe { alloc(layout) })
            .ok_or(NativeAvxPreparedG2LinesError::AllocationFailed)?;
        if !unsafe {
            blst_avx_prepared_g2_lines_init(
                storage.as_ptr().cast::<core::ffi::c_void>(),
                layout.size(),
                g2,
            )
        } {
            unsafe { dealloc(storage.as_ptr(), layout) };
            return Err(NativeAvxPreparedG2LinesError::InitializationFailed);
        }
        Ok(Self { storage, layout })
    }

    pub fn pairing_into(
        &self,
        out: &mut blst_avx_fp12,
        g1: &blst_avx_p1,
    ) -> Result<(), Error> {
        if unsafe {
            blst_avx_pairing_prepared_g2(
                out,
                self.storage.as_ptr().cast::<core::ffi::c_void>(),
                self.layout.size(),
                g1,
            )
        } {
            Ok(())
        } else {
            Err(Error::BackendFailure)
        }
    }

    pub fn storage_bytes(&self) -> usize {
        self.layout.size()
    }
}

impl core::fmt::Debug for NativeAvxPreparedG2Lines {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("NativeAvxPreparedG2Lines")
            .field("storage_bytes", &self.storage_bytes())
            .finish()
    }
}

impl Drop for NativeAvxPreparedG2Lines {
    fn drop(&mut self) {
        unsafe { dealloc(self.storage.as_ptr(), self.layout) };
    }
}

// SAFETY: initialization completes before publication. Evaluation receives a
// const table pointer and CatBLST does not mutate prepared line storage.
unsafe impl Send for NativeAvxPreparedG2Lines {}
unsafe impl Sync for NativeAvxPreparedG2Lines {}

pub type blst_avx_random_bytes_fn = Option<
    unsafe extern "C" fn(
        context: *mut core::ffi::c_void,
        out: *mut blst_byte,
        out_len: usize,
    ) -> bool,
>;

unsafe extern "C" fn copy_batch_membership_seed(
    context: *mut core::ffi::c_void,
    out: *mut blst_byte,
    out_len: usize,
) -> bool {
    if context.is_null() || out.is_null() || out_len != 32 {
        return false;
    }
    ptr::copy_nonoverlapping(context.cast::<blst_byte>(), out, out_len);
    true
}

/// Tests whether every AVX point belongs to G1 using CatBLST-owned challenges.
///
/// `random_seed` must be freshly sampled from a cryptographically secure RNG
/// after the batch is fixed. Under the SHA-256 PRF assumption, the randomized
/// check targets a false-accept probability below 2^-128. The seed is erased
/// from this function's local storage before return. Empty input succeeds
/// without consulting the seed when the AVX backend is available. Pointer-list
/// allocation failure returns `false` rather than panicking.
pub fn avx_p1s_in_g1(
    points: &[&blst_avx_p1],
    mut random_seed: [blst_byte; 32],
) -> bool {
    let mut pointers = Vec::new();
    if pointers.try_reserve_exact(points.len()).is_err() {
        random_seed.zeroize();
        return false;
    }
    pointers.extend(points.iter().map(|point| *point as *const blst_avx_p1));
    let accepted = unsafe {
        blst_avx_p1s_in_g1(
            pointers.as_ptr(),
            pointers.len(),
            Some(copy_batch_membership_seed),
            random_seed.as_mut_ptr().cast(),
        )
    };
    random_seed.zeroize();
    accepted
}

/// Tests whether every AVX point belongs to the prime-order BLS12-381 G2
/// subgroup.
///
/// `random_seed` must be freshly sampled after the batch is fixed. Under the
/// SHA-256 PRF assumption, a fixed invalid batch passes with probability at
/// most `13^-22` plus the PRF distinguishing advantage. `workspace` must have
/// at least [`blst_avx_p2s_in_g2_workspace_sizeof`] bytes and must be disjoint
/// from the points. Empty input succeeds without consulting either buffer when
/// the AVX backend is available. Pointer-list allocation failure returns
/// `false` rather than panicking.
pub fn avx_p2s_in_g2(
    points: &[&blst_avx_p2],
    mut random_seed: [blst_byte; 32],
    workspace: &mut [u8],
) -> bool {
    let mut pointers = Vec::new();
    if pointers.try_reserve_exact(points.len()).is_err() {
        random_seed.zeroize();
        return false;
    }
    pointers.extend(points.iter().map(|point| *point as *const blst_avx_p2));
    let accepted = unsafe {
        blst_avx_p2s_in_g2(
            pointers.as_ptr(),
            pointers.len(),
            Some(copy_batch_membership_seed),
            random_seed.as_mut_ptr().cast(),
            workspace.as_mut_ptr().cast(),
            workspace.len(),
        )
    };
    random_seed.zeroize();
    accepted
}

impl Default for blst_avx_p1 {
    fn default() -> Self {
        Self { b: [0; 384] }
    }
}

impl Default for blst_avx_p1_affine {
    fn default() -> Self {
        Self { b: [0; 256] }
    }
}

impl Default for blst_avx_p2 {
    fn default() -> Self {
        Self { b: [0; 768] }
    }
}

impl Default for blst_avx_p2_affine {
    fn default() -> Self {
        Self { b: [0; 512] }
    }
}

impl Default for blst_avx_fp12 {
    fn default() -> Self {
        Self { b: [0; 1536] }
    }
}

impl core::fmt::Debug for blst_avx_p1 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("blst_avx_p1").finish_non_exhaustive()
    }
}

impl core::fmt::Debug for blst_avx_p1_affine {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("blst_avx_p1_affine").finish_non_exhaustive()
    }
}

impl core::fmt::Debug for blst_avx_p2 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("blst_avx_p2").finish_non_exhaustive()
    }
}

impl core::fmt::Debug for blst_avx_p2_affine {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("blst_avx_p2_affine").finish_non_exhaustive()
    }
}

impl core::fmt::Debug for blst_avx_fp12 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("blst_avx_fp12").finish_non_exhaustive()
    }
}

impl PartialEq for blst_avx_p1 {
    fn eq(&self, other: &Self) -> bool {
        self.b == other.b
    }
}

impl PartialEq for blst_avx_p1_affine {
    fn eq(&self, other: &Self) -> bool {
        self.b == other.b
    }
}

impl PartialEq for blst_avx_p2 {
    fn eq(&self, other: &Self) -> bool {
        self.b == other.b
    }
}

impl PartialEq for blst_avx_p2_affine {
    fn eq(&self, other: &Self) -> bool {
        self.b == other.b
    }
}

impl PartialEq for blst_avx_fp12 {
    fn eq(&self, other: &Self) -> bool {
        self.b == other.b
    }
}

const SCALAR_BITS: usize = 255;
const SCRATCH_ALIGNMENT: usize = 64;

#[inline]
fn wipe_scalar_bytes(bytes: &mut [u8; 32]) {
    for byte in bytes {
        unsafe { ptr::write_volatile(byte, 0) };
    }
    compiler_fence(Ordering::SeqCst);
}

/// Failure from the safe, backend-native AVX curve API.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Error {
    /// This build or host cannot execute the CatBLST AVX backend.
    Unavailable,
    /// A canonical point or target-group encoding was malformed.
    InvalidEncoding,
    /// A decoded curve point was not in the prime-order subgroup.
    InvalidSubgroup,
    /// Point and scalar counts differed.
    MsmLengthMismatch { points: usize, scalars: usize },
    /// The typed MSM API does not accept an empty input.
    EmptyMsm,
    /// Prepared-G2 batch construction requires at least one operand.
    EmptyPreparedBatch,
    /// Sparse prepared-pairing index and output counts differed.
    PairingLengthMismatch { indices: usize, outputs: usize },
    /// A sparse prepared-pairing index exceeded the immutable batch.
    PreparedIndexOutOfRange { index: usize, len: usize },
    /// Caller-owned scratch was too small for the requested operation.
    ScratchTooSmall { required: usize, available: usize },
    /// A fixed-capacity pairing workspace cannot accept another pair.
    PairingCapacityExceeded { required: usize, capacity: usize },
    /// A required allocation failed while constructing reusable storage.
    AllocationFailed,
    /// A native CatBLST operation rejected its inputs or failed.
    BackendFailure,
}

#[allow(non_upper_case_globals)]
impl Error {
    /// Compatibility name for the former Monad-local generic kernel failure.
    #[deprecated(note = "use Error::BackendFailure")]
    pub const KernelFailure: Self = Self::BackendFailure;

    /// Compatibility name for the former Monad-local middle-product failure.
    #[deprecated(note = "use Error::BackendFailure")]
    pub const MiddleProductFailure: Self = Self::BackendFailure;
}

impl core::fmt::Display for Error {
    fn fmt(
        &self,
        formatter: &mut core::fmt::Formatter<'_>,
    ) -> core::fmt::Result {
        match self {
            Self::Unavailable => formatter.write_str("the CatBLST AVX backend is unavailable"),
            Self::InvalidEncoding => formatter.write_str("invalid canonical group encoding"),
            Self::InvalidSubgroup => formatter.write_str("point is outside the prime-order subgroup"),
            Self::MsmLengthMismatch { points, scalars } => {
                write!(formatter, "MSM length mismatch: points={points}, scalars={scalars}")
            }
            Self::EmptyMsm => formatter.write_str("MSM input must not be empty"),
            Self::EmptyPreparedBatch => {
                formatter.write_str("prepared G2 batch must not be empty")
            }
            Self::PairingLengthMismatch { indices, outputs } => write!(
                formatter,
                "prepared pairing length mismatch: indices={indices}, outputs={outputs}"
            ),
            Self::PreparedIndexOutOfRange { index, len } => write!(
                formatter,
                "prepared pairing index {index} is outside batch length {len}"
            ),
            Self::ScratchTooSmall { required, available } => write!(
                formatter,
                "MSM scratch is too small: required={required}, available={available}"
            ),
            Self::PairingCapacityExceeded { required, capacity } => write!(
                formatter,
                "pairing capacity exceeded: required={required}, capacity={capacity}"
            ),
            Self::AllocationFailed => formatter.write_str("reusable AVX storage allocation failed"),
            Self::BackendFailure => formatter.write_str("CatBLST AVX operation failed"),
        }
    }
}

#[cfg(blst_std)]
impl std::error::Error for Error {}

/// Returns whether the native AVX implementation is available.
#[inline]
pub fn is_available() -> bool {
    unsafe { blst_avx_enabled() }
}

#[inline]
fn ensure_available() -> Result<(), Error> {
    if is_available() {
        Ok(())
    } else {
        Err(Error::Unavailable)
    }
}

/// Backend-native G1 projective point.
///
/// The representation is opaque, 64-byte aligned, and incompatible with
/// [`crate::scalar::G1Projective`]. Canonical encoding is the only supported
/// interchange boundary between backends.
#[repr(transparent)]
#[derive(Clone, Copy)]
pub struct G1Projective(blst_avx_p1);

const _: () = {
    assert!(size_of::<G1Projective>() == size_of::<blst_avx_p1>());
    assert!(
        core::mem::align_of::<G1Projective>()
            == core::mem::align_of::<blst_avx_p1>()
    );
};

impl G1Projective {
    /// Returns the G1 identity without allocating.
    pub fn identity() -> Result<Self, Error> {
        ensure_available()?;
        let mut point = MaybeUninit::<blst_avx_p1>::uninit();
        if unsafe { blst_avx_p1_identity(point.as_mut_ptr()) } {
            Ok(Self(unsafe { point.assume_init() }))
        } else {
            Err(Error::BackendFailure)
        }
    }

    /// Returns the canonical G1 generator.
    pub fn generator() -> Result<Self, Error> {
        ensure_available()?;
        let generator = unsafe { blst_avx_p1_generator() };
        if generator.is_null() {
            return Err(Error::Unavailable);
        }
        let mut point = MaybeUninit::<blst_avx_p1>::uninit();
        // Copy bytes because portable AVX builds may expose static generator
        // storage with weaker alignment than the opaque public shell.
        unsafe {
            ptr::copy_nonoverlapping(
                generator.cast::<u8>(),
                point.as_mut_ptr().cast::<u8>(),
                size_of::<blst_avx_p1>(),
            );
        }
        Ok(Self(unsafe { point.assume_init() }))
    }

    /// Strictly decodes a compressed point and checks G1 subgroup membership.
    pub fn from_compressed(encoded: &[u8; 48]) -> Result<Self, Error> {
        ensure_available()?;
        let mut point = MaybeUninit::<blst_avx_p1>::uninit();
        if unsafe {
            blst_avx_p1_uncompress(point.as_mut_ptr(), encoded.as_ptr())
        } != BLST_ERROR::BLST_SUCCESS
        {
            return Err(Error::InvalidEncoding);
        }
        let point = unsafe { point.assume_init() };
        if !unsafe { blst_avx_p1_in_g1(&point) } {
            return Err(Error::InvalidSubgroup);
        }
        Ok(Self(point))
    }

    /// Returns the canonical 48-byte compressed encoding.
    pub fn to_compressed(&self) -> Result<[u8; 48], Error> {
        let mut encoded = MaybeUninit::<[u8; 48]>::uninit();
        if unsafe { blst_avx_p1_compress(encoded.as_mut_ptr().cast(), &self.0) }
        {
            Ok(unsafe { encoded.assume_init() })
        } else {
            Err(Error::BackendFailure)
        }
    }

    /// Multiplies by a typed canonical scalar in constant time.
    ///
    /// The operation allocates nothing. Runtime is independent of the scalar
    /// contents for the fixed 255-bit width.
    pub fn multiply(&self, scalar: &Scalar) -> Result<Self, Error> {
        let mut scalar_le = scalar.to_bytes_le();
        let mut out = MaybeUninit::<blst_avx_p1>::uninit();
        let succeeded = unsafe {
            blst_avx_p1_mult(
                out.as_mut_ptr(),
                &self.0,
                scalar_le.as_ptr(),
                SCALAR_BITS,
            )
        };
        wipe_scalar_bytes(&mut scalar_le);
        if succeeded {
            Ok(Self(unsafe { out.assume_init() }))
        } else {
            Err(Error::BackendFailure)
        }
    }

    /// Adds two validated G1 points without allocating.
    pub fn add(&self, other: &Self) -> Result<Self, Error> {
        let mut out = MaybeUninit::<blst_avx_p1>::uninit();
        if unsafe { blst_avx_p1_add(out.as_mut_ptr(), &self.0, &other.0) } {
            Ok(Self(unsafe { out.assume_init() }))
        } else {
            Err(Error::BackendFailure)
        }
    }

    /// Returns whether this point is the identity.
    #[inline]
    pub fn is_identity(&self) -> bool {
        unsafe { blst_avx_p1_is_inf(&self.0) }
    }

    fn as_raw_affine(&self) -> Result<blst_avx_p1_affine, Error> {
        let mut affine = MaybeUninit::<blst_avx_p1_affine>::uninit();
        if unsafe { blst_avx_p1_to_affine(affine.as_mut_ptr(), &self.0) } {
            Ok(unsafe { affine.assume_init() })
        } else {
            Err(Error::BackendFailure)
        }
    }

    #[inline]
    pub(crate) fn as_raw(&self) -> &blst_avx_p1 {
        &self.0
    }

    #[inline]
    pub(crate) fn as_raw_slice(points: &[Self]) -> &[blst_avx_p1] {
        unsafe {
            core::slice::from_raw_parts(points.as_ptr().cast(), points.len())
        }
    }
}

impl core::fmt::Debug for G1Projective {
    fn fmt(
        &self,
        formatter: &mut core::fmt::Formatter<'_>,
    ) -> core::fmt::Result {
        formatter
            .debug_tuple("G1Projective")
            .field(&"<native-avx>")
            .finish()
    }
}

impl PartialEq for G1Projective {
    fn eq(&self, other: &Self) -> bool {
        unsafe { blst_avx_p1_is_equal(&self.0, &other.0) }
    }
}

impl Eq for G1Projective {}

/// Backend-native G2 projective point.
///
/// The representation is opaque, 64-byte aligned, and incompatible with
/// [`crate::scalar::G2Projective`].
#[repr(transparent)]
#[derive(Clone, Copy)]
pub struct G2Projective(blst_avx_p2);

const _: () = {
    assert!(size_of::<G2Projective>() == size_of::<blst_avx_p2>());
    assert!(
        core::mem::align_of::<G2Projective>()
            == core::mem::align_of::<blst_avx_p2>()
    );
};

impl G2Projective {
    /// Returns the G2 identity without allocating.
    pub fn identity() -> Result<Self, Error> {
        ensure_available()?;
        let mut point = MaybeUninit::<blst_avx_p2>::uninit();
        if unsafe { blst_avx_p2_identity(point.as_mut_ptr()) } {
            Ok(Self(unsafe { point.assume_init() }))
        } else {
            Err(Error::BackendFailure)
        }
    }

    /// Returns the canonical G2 generator.
    pub fn generator() -> Result<Self, Error> {
        ensure_available()?;
        let generator = unsafe { blst_avx_p2_generator() };
        if generator.is_null() {
            return Err(Error::Unavailable);
        }
        let mut point = MaybeUninit::<blst_avx_p2>::uninit();
        unsafe {
            ptr::copy_nonoverlapping(
                generator.cast::<u8>(),
                point.as_mut_ptr().cast::<u8>(),
                size_of::<blst_avx_p2>(),
            );
        }
        Ok(Self(unsafe { point.assume_init() }))
    }

    /// Strictly decodes a compressed point and checks G2 subgroup membership.
    pub fn from_compressed(encoded: &[u8; 96]) -> Result<Self, Error> {
        ensure_available()?;
        let mut point = MaybeUninit::<blst_avx_p2>::uninit();
        if unsafe {
            blst_avx_p2_uncompress(point.as_mut_ptr(), encoded.as_ptr())
        } != BLST_ERROR::BLST_SUCCESS
        {
            return Err(Error::InvalidEncoding);
        }
        let point = unsafe { point.assume_init() };
        if !unsafe { blst_avx_p2_in_g2(&point) } {
            return Err(Error::InvalidSubgroup);
        }
        Ok(Self(point))
    }

    /// Returns the canonical 96-byte compressed encoding.
    pub fn to_compressed(&self) -> Result<[u8; 96], Error> {
        let mut encoded = MaybeUninit::<[u8; 96]>::uninit();
        if unsafe { blst_avx_p2_compress(encoded.as_mut_ptr().cast(), &self.0) }
        {
            Ok(unsafe { encoded.assume_init() })
        } else {
            Err(Error::BackendFailure)
        }
    }

    /// Multiplies by a typed canonical scalar in constant time.
    pub fn multiply(&self, scalar: &Scalar) -> Result<Self, Error> {
        let mut scalar_le = scalar.to_bytes_le();
        let mut out = MaybeUninit::<blst_avx_p2>::uninit();
        let succeeded = unsafe {
            blst_avx_p2_mult(
                out.as_mut_ptr(),
                &self.0,
                scalar_le.as_ptr(),
                SCALAR_BITS,
            )
        };
        wipe_scalar_bytes(&mut scalar_le);
        if succeeded {
            Ok(Self(unsafe { out.assume_init() }))
        } else {
            Err(Error::BackendFailure)
        }
    }

    /// Adds two validated G2 points without allocating.
    pub fn add(&self, other: &Self) -> Result<Self, Error> {
        let mut out = MaybeUninit::<blst_avx_p2>::uninit();
        if unsafe { blst_avx_p2_add(out.as_mut_ptr(), &self.0, &other.0) } {
            Ok(Self(unsafe { out.assume_init() }))
        } else {
            Err(Error::BackendFailure)
        }
    }

    /// Returns whether this point is the identity.
    #[inline]
    pub fn is_identity(&self) -> bool {
        unsafe { blst_avx_p2_is_inf(&self.0) }
    }

    fn as_raw_affine(&self) -> Result<blst_avx_p2_affine, Error> {
        let mut affine = MaybeUninit::<blst_avx_p2_affine>::uninit();
        if unsafe { blst_avx_p2_to_affine(affine.as_mut_ptr(), &self.0) } {
            Ok(unsafe { affine.assume_init() })
        } else {
            Err(Error::BackendFailure)
        }
    }

    #[inline]
    pub(crate) fn as_raw(&self) -> &blst_avx_p2 {
        &self.0
    }

    #[inline]
    pub(crate) fn as_raw_slice(points: &[Self]) -> &[blst_avx_p2] {
        unsafe {
            core::slice::from_raw_parts(points.as_ptr().cast(), points.len())
        }
    }
}

impl core::fmt::Debug for G2Projective {
    fn fmt(
        &self,
        formatter: &mut core::fmt::Formatter<'_>,
    ) -> core::fmt::Result {
        formatter
            .debug_tuple("G2Projective")
            .field(&"<native-avx>")
            .finish()
    }
}

impl PartialEq for G2Projective {
    fn eq(&self, other: &Self) -> bool {
        unsafe { blst_avx_p2_is_equal(&self.0, &other.0) }
    }
}

impl Eq for G2Projective {}

/// Backend-native target-group element.
///
/// The representation is opaque and incompatible with [`crate::scalar::Gt`].
/// Serialization uses CatBLST's shared canonical twelve-coefficient encoding.
#[repr(transparent)]
#[derive(Clone, Copy)]
pub struct Gt(blst_avx_gt);

const _: () = {
    assert!(size_of::<Gt>() == size_of::<blst_avx_gt>());
    assert!(
        core::mem::align_of::<Gt>() == core::mem::align_of::<blst_avx_gt>()
    );
};

impl Gt {
    /// Returns the target-group identity.
    pub fn identity() -> Result<Self, Error> {
        ensure_available()?;
        let identity = unsafe { blst_avx_fp12_one() };
        if identity.is_null() {
            return Err(Error::Unavailable);
        }
        let mut value = MaybeUninit::<blst_avx_gt>::uninit();
        unsafe {
            ptr::copy_nonoverlapping(
                identity.cast::<u8>(),
                value.as_mut_ptr().cast::<u8>(),
                size_of::<blst_avx_gt>(),
            );
        }
        Ok(Self(unsafe { value.assume_init() }))
    }

    /// Strictly decodes CatBLST's canonical 576-byte GT encoding.
    pub fn from_bytes_be(encoded: &[u8; 576]) -> Result<Self, Error> {
        ensure_available()?;
        let mut value = MaybeUninit::<blst_avx_gt>::uninit();
        if unsafe {
            blst_avx_gt_deserialize(value.as_mut_ptr(), encoded.as_ptr())
        } {
            Ok(Self(unsafe { value.assume_init() }))
        } else {
            Err(Error::InvalidEncoding)
        }
    }

    /// Returns CatBLST's canonical 576-byte GT encoding.
    pub fn to_bytes_be(&self) -> Result<[u8; 576], Error> {
        let mut encoded = MaybeUninit::<[u8; 576]>::uninit();
        if unsafe {
            blst_avx_gt_serialize(encoded.as_mut_ptr().cast(), &self.0)
        } {
            Ok(unsafe { encoded.assume_init() })
        } else {
            Err(Error::BackendFailure)
        }
    }

    /// Computes one pairing into caller-owned native GT storage.
    pub fn pairing_into(
        out: &mut Self,
        g1: &G1Projective,
        g2: &G2Projective,
    ) -> Result<(), Error> {
        if unsafe { blst_avx_pairing(&mut out.0, g1.as_raw(), g2.as_raw()) } {
            Ok(())
        } else {
            Err(Error::BackendFailure)
        }
    }

    /// Computes one pairing directly into native GT storage.
    pub fn pairing(
        g1: &G1Projective,
        g2: &G2Projective,
    ) -> Result<Self, Error> {
        let mut out = MaybeUninit::<blst_avx_gt>::uninit();
        if unsafe {
            blst_avx_pairing(out.as_mut_ptr(), g1.as_raw(), g2.as_raw())
        } {
            Ok(Self(unsafe { out.assume_init() }))
        } else {
            Err(Error::BackendFailure)
        }
    }

    /// Multiplies two target-group elements without allocating.
    pub fn multiply(&self, other: &Self) -> Result<Self, Error> {
        let mut out = MaybeUninit::<blst_avx_gt>::uninit();
        if unsafe { blst_avx_fp12_mul(out.as_mut_ptr(), &self.0, &other.0) } {
            Ok(Self(unsafe { out.assume_init() }))
        } else {
            Err(Error::BackendFailure)
        }
    }

    /// Squares a target-group element without allocating.
    pub fn square(&self) -> Result<Self, Error> {
        let mut out = MaybeUninit::<blst_avx_gt>::uninit();
        if unsafe { blst_avx_fp12_sqr(out.as_mut_ptr(), &self.0) } {
            Ok(Self(unsafe { out.assume_init() }))
        } else {
            Err(Error::BackendFailure)
        }
    }

    /// Returns the target-group inverse without allocating.
    pub fn inverse(&self) -> Result<Self, Error> {
        let mut out = MaybeUninit::<blst_avx_gt>::uninit();
        if unsafe { blst_avx_gt_inverse(out.as_mut_ptr(), &self.0) } {
            Ok(Self(unsafe { out.assume_init() }))
        } else {
            Err(Error::BackendFailure)
        }
    }

    /// Raises this element to a typed canonical scalar in constant time.
    ///
    /// Runtime and memory access are independent of the scalar contents at the
    /// fixed 255-bit width. The operation allocates nothing.
    pub fn pow(&self, scalar: &Scalar) -> Result<Self, Error> {
        let mut scalar_le = scalar.to_bytes_le();
        let mut out = MaybeUninit::<blst_avx_gt>::uninit();
        let succeeded = unsafe {
            blst_avx_gt_mult(
                out.as_mut_ptr(),
                &self.0,
                scalar_le.as_ptr(),
                SCALAR_BITS,
            )
        };
        wipe_scalar_bytes(&mut scalar_le);
        if succeeded {
            Ok(Self(unsafe { out.assume_init() }))
        } else {
            Err(Error::BackendFailure)
        }
    }

    #[inline]
    pub(crate) fn as_raw_slice_mut(values: &mut [Self]) -> &mut [blst_avx_gt] {
        unsafe {
            core::slice::from_raw_parts_mut(
                values.as_mut_ptr().cast(),
                values.len(),
            )
        }
    }
}

impl core::fmt::Debug for Gt {
    fn fmt(
        &self,
        formatter: &mut core::fmt::Formatter<'_>,
    ) -> core::fmt::Result {
        formatter.debug_tuple("Gt").field(&"<native-avx>").finish()
    }
}

impl PartialEq for Gt {
    fn eq(&self, other: &Self) -> bool {
        unsafe { blst_avx_fp12_is_equal(&self.0, &other.0) }
    }
}

impl Eq for Gt {}

/// Immutable native Miller lines for one validated G2 point.
#[derive(Debug)]
pub struct PreparedG2Lines {
    inner: NativeAvxPreparedG2Lines,
}

impl PreparedG2Lines {
    /// Allocates and prepares one immutable native line table.
    pub fn new(g2: &G2Projective) -> Result<Self, Error> {
        let inner =
            NativeAvxPreparedG2Lines::new(g2.as_raw()).map_err(|error| {
                match error {
                    NativeAvxPreparedG2LinesError::Unavailable => {
                        Error::Unavailable
                    }
                    NativeAvxPreparedG2LinesError::AllocationFailed => {
                        Error::AllocationFailed
                    }
                    NativeAvxPreparedG2LinesError::InvalidLayout
                    | NativeAvxPreparedG2LinesError::InitializationFailed => {
                        Error::BackendFailure
                    }
                }
            })?;
        Ok(Self { inner })
    }

    /// Evaluates a pairing with the prepared G2 operand without allocating.
    pub fn pairing(&self, g1: &G1Projective) -> Result<Gt, Error> {
        let mut out = MaybeUninit::<blst_avx_gt>::uninit();
        if unsafe {
            blst_avx_pairing_prepared_g2(
                out.as_mut_ptr(),
                self.inner.storage.as_ptr().cast::<core::ffi::c_void>(),
                self.inner.layout.size(),
                g1.as_raw(),
            )
        } {
            Ok(Gt(unsafe { out.assume_init() }))
        } else {
            Err(Error::BackendFailure)
        }
    }

    /// Returns the owned native table size in bytes.
    pub fn storage_bytes(&self) -> usize {
        self.inner.storage_bytes()
    }
}

/// Immutable native Miller-line tables for a sequence of validated G2 points.
///
/// Construction allocates and prepares every operand. Evaluations use a
/// distinct [`PreparedG2BatchWorkspace`], so this object may be shared by
/// concurrent workers.
pub struct PreparedG2Batch {
    inner: NonNull<blst_avx_prepared_g2_batch>,
    len: usize,
    storage_bytes: usize,
}

impl PreparedG2Batch {
    /// Allocates and prepares one native line table per input point.
    pub fn new(g2: &[G2Projective]) -> Result<Self, Error> {
        ensure_available()?;
        if g2.is_empty() {
            return Err(Error::EmptyPreparedBatch);
        }
        let inner = NonNull::new(unsafe {
            blst_avx_prepared_g2_batch_new(
                G2Projective::as_raw_slice(g2).as_ptr(),
                g2.len(),
            )
        })
        .ok_or(Error::AllocationFailed)?;
        let len = unsafe { blst_avx_prepared_g2_batch_len(inner.as_ptr()) };
        if len != g2.len() {
            unsafe { blst_avx_prepared_g2_batch_free(inner.as_ptr()) };
            return Err(Error::BackendFailure);
        }
        let storage_bytes =
            unsafe { blst_avx_prepared_g2_batch_storage_bytes(inner.as_ptr()) };
        Ok(Self {
            inner,
            len,
            storage_bytes,
        })
    }

    /// Returns the number of prepared operands.
    #[inline]
    pub fn len(&self) -> usize {
        self.len
    }

    /// Returns whether the prepared sequence is empty. Construction rejects
    /// empty input, so this is always false for a live object.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Returns CatBLST's accounting estimate for immutable owned storage.
    #[inline]
    pub fn storage_bytes(&self) -> usize {
        self.storage_bytes
    }
}

impl core::fmt::Debug for PreparedG2Batch {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("PreparedG2Batch")
            .field("len", &self.len)
            .field("storage_bytes", &self.storage_bytes)
            .finish()
    }
}

impl Drop for PreparedG2Batch {
    fn drop(&mut self) {
        unsafe { blst_avx_prepared_g2_batch_free(self.inner.as_ptr()) };
    }
}

// SAFETY: the native object is immutable after successful construction.
unsafe impl Send for PreparedG2Batch {}
unsafe impl Sync for PreparedG2Batch {}

/// Caller-owned scratch for allocation-free prepared-G2 batch evaluation.
///
/// A workspace is mutable and must be used by only one evaluation at a time.
pub struct PreparedG2BatchWorkspace {
    inner: NonNull<blst_avx_prepared_g2_batch_workspace>,
    capacity: usize,
}

impl PreparedG2BatchWorkspace {
    /// Allocates all Miller and batch-final-exponentiation scratch required by
    /// `batch`. No subsequent prefix evaluation allocates.
    pub fn new(batch: &PreparedG2Batch) -> Result<Self, Error> {
        let inner = NonNull::new(unsafe {
            blst_avx_prepared_g2_batch_workspace_new(batch.inner.as_ptr())
        })
        .ok_or(Error::AllocationFailed)?;
        Ok(Self {
            inner,
            capacity: batch.len,
        })
    }

    /// Evaluates the prepared prefix selected by `out.len()`.
    ///
    /// The shared G1 operand is converted to affine once. The output slice may
    /// be empty and must not overlap the inputs or owned native objects.
    pub fn pairing_prefix_into(
        &mut self,
        g1: &G1Projective,
        batch: &PreparedG2Batch,
        out: &mut [Gt],
    ) -> Result<(), Error> {
        if out.len() > batch.len || out.len() > self.capacity {
            return Err(Error::PairingCapacityExceeded {
                required: out.len(),
                capacity: core::cmp::min(batch.len, self.capacity),
            });
        }
        if unsafe {
            blst_avx_prepared_g2_batch_pairing_prefix(
                Gt::as_raw_slice_mut(out).as_mut_ptr(),
                out.len(),
                g1.as_raw(),
                batch.inner.as_ptr(),
                self.inner.as_ptr(),
            )
        } {
            Ok(())
        } else {
            Err(Error::BackendFailure)
        }
    }

    /// Evaluates selected prepared operands in `indices` order.
    ///
    /// The output length must equal the index count. Repeated indices are
    /// permitted. Bounds are checked before the native call publishes output.
    pub fn pairing_indices_into(
        &mut self,
        g1: &G1Projective,
        batch: &PreparedG2Batch,
        indices: &[usize],
        out: &mut [Gt],
    ) -> Result<(), Error> {
        if indices.len() != out.len() {
            return Err(Error::PairingLengthMismatch {
                indices: indices.len(),
                outputs: out.len(),
            });
        }
        if indices.len() > self.capacity {
            return Err(Error::PairingCapacityExceeded {
                required: indices.len(),
                capacity: self.capacity,
            });
        }
        if let Some(&index) = indices.iter().find(|&&index| index >= batch.len)
        {
            return Err(Error::PreparedIndexOutOfRange {
                index,
                len: batch.len,
            });
        }
        if unsafe {
            blst_avx_prepared_g2_batch_pairing_indices(
                Gt::as_raw_slice_mut(out).as_mut_ptr(),
                indices.as_ptr(),
                indices.len(),
                g1.as_raw(),
                batch.inner.as_ptr(),
                self.inner.as_ptr(),
            )
        } {
            Ok(())
        } else {
            Err(Error::BackendFailure)
        }
    }

    /// Returns the maximum prefix length this workspace can evaluate.
    #[inline]
    pub fn capacity(&self) -> usize {
        self.capacity
    }
}

impl core::fmt::Debug for PreparedG2BatchWorkspace {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("PreparedG2BatchWorkspace")
            .field("capacity", &self.capacity)
            .finish()
    }
}

impl Drop for PreparedG2BatchWorkspace {
    fn drop(&mut self) {
        unsafe {
            blst_avx_prepared_g2_batch_workspace_free(self.inner.as_ptr())
        };
    }
}

// SAFETY: ownership is unique and the workspace is not Sync.
unsafe impl Send for PreparedG2BatchWorkspace {}

/// Caller-owned affine staging for allocation-free batched pairings.
#[derive(Debug)]
pub struct PairingWorkspace {
    g1: Vec<blst_avx_p1_affine>,
    g2: Vec<blst_avx_p2_affine>,
    max_pairs: usize,
}

impl PairingWorkspace {
    /// Allocates fixed capacity outside the evaluation hot path.
    pub fn for_pairs(max_pairs: usize) -> Result<Self, Error> {
        ensure_available()?;
        let mut g1 = Vec::new();
        let mut g2 = Vec::new();
        g1.try_reserve_exact(max_pairs)
            .map_err(|_| Error::AllocationFailed)?;
        g2.try_reserve_exact(max_pairs)
            .map_err(|_| Error::AllocationFailed)?;
        Ok(Self { g1, g2, max_pairs })
    }

    /// Removes staged pairs while retaining all allocation.
    pub fn clear(&mut self) {
        self.g1.clear();
        self.g2.clear();
    }

    /// Converts and stages one pair without allocating.
    pub fn push(
        &mut self,
        g1: &G1Projective,
        g2: &G2Projective,
    ) -> Result<(), Error> {
        let required = self.g1.len() + 1;
        if required > self.max_pairs {
            return Err(Error::PairingCapacityExceeded {
                required,
                capacity: self.max_pairs,
            });
        }
        let g1 = g1.as_raw_affine()?;
        let g2 = g2.as_raw_affine()?;
        self.g1.push(g1);
        self.g2.push(g2);
        Ok(())
    }

    /// Evaluates all staged pairs without allocating; an empty product is one.
    pub fn evaluate(&self) -> Result<Gt, Error> {
        if self.g1.len() != self.g2.len() {
            return Err(Error::BackendFailure);
        }
        let mut out = MaybeUninit::<blst_avx_gt>::uninit();
        let (g1, g2) = if self.g1.is_empty() {
            (ptr::null(), ptr::null())
        } else {
            (self.g1.as_ptr(), self.g2.as_ptr())
        };
        if unsafe {
            blst_avx_pairing_n(out.as_mut_ptr(), g1, g2, self.g1.len())
        } {
            Ok(Gt(unsafe { out.assume_init() }))
        } else {
            Err(Error::BackendFailure)
        }
    }
}

#[repr(C, align(64))]
#[derive(Clone, Copy, Debug)]
struct AlignedScratchBlock([u8; SCRATCH_ALIGNMENT]);

/// Caller-owned, reusable scratch for native projective G1 Pippenger MSM.
#[derive(Debug)]
pub struct G1MsmWorkspace {
    scratch: Vec<MaybeUninit<AlignedScratchBlock>>,
}

impl G1MsmWorkspace {
    /// Allocates enough aligned scratch for at most `max_points` points.
    pub fn for_points(max_points: usize) -> Result<Self, Error> {
        ensure_available()?;
        let pippenger_required =
            unsafe { blst_avx_p1s_mult_pippenger_scratch_sizeof(max_points) };
        let secret_required = if max_points >= 8 {
            unsafe { blst_avx_p1s_mult_ct_scratch_sizeof() }
        } else {
            0
        };
        let required = pippenger_required.max(secret_required);
        let blocks = required.div_ceil(SCRATCH_ALIGNMENT);
        let mut scratch = Vec::new();
        scratch
            .try_reserve_exact(blocks)
            .map_err(|_| Error::AllocationFailed)?;
        scratch.resize_with(blocks, MaybeUninit::uninit);
        Ok(Self { scratch })
    }

    /// Returns the reusable scratch capacity in bytes.
    pub fn capacity_bytes(&self) -> usize {
        self.scratch.len() * SCRATCH_ALIGNMENT
    }

    /// Computes a constant-time MSM over secret typed scalars.
    ///
    /// The points must already be validated G1 subgroup points. The scalar
    /// slice retains CatBLST's native Montgomery representation across the
    /// Rust/C++ boundary; the kernel performs fixed-time GLV preparation
    /// internally. The operation allocates nothing and uses fixed 32-, 16-,
    /// and 8-point chunks before a constant-time scalar tail. Large public
    /// point tables occupy this workspace rather than the thread stack.
    pub fn secret_msm_into(
        &mut self,
        out: &mut G1Projective,
        points: &[G1Projective],
        scalars: &[Scalar],
    ) -> Result<(), Error> {
        if points.len() != scalars.len() {
            return Err(Error::MsmLengthMismatch {
                points: points.len(),
                scalars: scalars.len(),
            });
        }
        let required = if points.len() >= 8 {
            unsafe { blst_avx_p1s_mult_ct_scratch_sizeof() }
        } else {
            0
        };
        let available = self.capacity_bytes();
        if available < required {
            return Err(Error::ScratchTooSmall {
                required,
                available,
            });
        }
        let scratch = if self.scratch.is_empty() {
            ptr::null_mut()
        } else {
            self.scratch.as_mut_ptr().cast::<core::ffi::c_void>()
        };
        if unsafe {
            blst_avx_p1s_mult_ct(
                &mut out.0,
                G1Projective::as_raw_slice(points).as_ptr(),
                scalars.as_ptr().cast::<blst_fr>(),
                points.len(),
                scratch,
                available,
            )
        } {
            Ok(())
        } else {
            Err(Error::BackendFailure)
        }
    }

    /// Computes one MSM into caller-owned output without allocating.
    ///
    /// The Pippenger kernel is variable-time, so every scalar must be public.
    pub fn msm_into(
        &mut self,
        out: &mut G1Projective,
        points: &[G1Projective],
        scalars_le: &[[u8; 32]],
    ) -> Result<(), Error> {
        if points.len() != scalars_le.len() {
            return Err(Error::MsmLengthMismatch {
                points: points.len(),
                scalars: scalars_le.len(),
            });
        }
        if points.is_empty() {
            return Err(Error::EmptyMsm);
        }
        let required =
            unsafe { blst_avx_p1s_mult_pippenger_scratch_sizeof(points.len()) };
        let available = self.capacity_bytes();
        if available < required {
            return Err(Error::ScratchTooSmall {
                required,
                available,
            });
        }
        let scratch = if self.scratch.is_empty() {
            ptr::null_mut()
        } else {
            self.scratch.as_mut_ptr().cast::<blst_limb_t>()
        };
        if unsafe {
            blst_avx_p1s_mult_pippenger(
                &mut out.0,
                G1Projective::as_raw_slice(points).as_ptr(),
                points.len(),
                scalars_le.as_ptr().cast::<u8>(),
                SCALAR_BITS,
                scratch,
            )
        } {
            Ok(())
        } else {
            Err(Error::BackendFailure)
        }
    }
}

extern "C" {
    pub fn blst_avx_enabled() -> bool;
    pub fn blst_avx_p1_sizeof() -> usize;
    pub fn blst_avx_p1_alignof() -> usize;
    pub fn blst_avx_p1_affine_sizeof() -> usize;
    pub fn blst_avx_p1_affine_alignof() -> usize;
    pub fn blst_avx_p2_sizeof() -> usize;
    pub fn blst_avx_p2_alignof() -> usize;
    pub fn blst_avx_p2_affine_sizeof() -> usize;
    pub fn blst_avx_p2_affine_alignof() -> usize;
    pub fn blst_avx_fp12_sizeof() -> usize;
    pub fn blst_avx_fp12_alignof() -> usize;

    pub fn blst_avx_p1_add(
        out: *mut blst_avx_p1,
        a: *const blst_avx_p1,
        b: *const blst_avx_p1,
    ) -> bool;
    pub fn blst_avx_p1_from_scalar(
        out: *mut blst_avx_p1,
        input: *const blst_p1,
    ) -> bool;
    pub fn blst_p1_from_avx(
        out: *mut blst_p1,
        input: *const blst_avx_p1,
    ) -> bool;
    pub fn blst_avx_p1_mult(
        out: *mut blst_avx_p1,
        p: *const blst_avx_p1,
        scalar: *const blst_byte,
        nbits: usize,
    ) -> bool;
    pub fn blst_avx_p1s_mult_ct(
        out: *mut blst_avx_p1,
        points: *const blst_avx_p1,
        scalars: *const blst_fr,
        npoints: usize,
        scratch: *mut core::ffi::c_void,
        scratch_bytes: usize,
    ) -> bool;
    pub fn blst_avx_p1s_mult_ct_scratch_sizeof() -> usize;
    pub fn blst_avx_p1_to_affine(
        out: *mut blst_avx_p1_affine,
        input: *const blst_avx_p1,
    ) -> bool;
    pub fn blst_avx_p1_affine_from_scalar(
        out: *mut blst_avx_p1_affine,
        input: *const blst_p1_affine,
    ) -> bool;
    pub fn blst_p1_affine_from_avx(
        out: *mut blst_p1_affine,
        input: *const blst_avx_p1_affine,
    ) -> bool;
    pub fn blst_avx_p1_from_affine(
        out: *mut blst_avx_p1,
        input: *const blst_avx_p1_affine,
    ) -> bool;
    pub fn blst_avx_p1_in_g1(p: *const blst_avx_p1) -> bool;
    pub fn blst_avx_p1_is_inf(p: *const blst_avx_p1) -> bool;
    pub fn blst_avx_p1_is_equal(
        a: *const blst_avx_p1,
        b: *const blst_avx_p1,
    ) -> bool;
    pub fn blst_avx_p1_compress(
        out: *mut blst_byte,
        input: *const blst_avx_p1,
    ) -> bool;
    pub fn blst_avx_p1_serialize(
        out: *mut blst_byte,
        input: *const blst_avx_p1,
    ) -> bool;
    pub fn blst_avx_p1_deserialize(
        out: *mut blst_avx_p1,
        input: *const blst_byte,
    ) -> BLST_ERROR;
    pub fn blst_avx_p1_uncompress(
        out: *mut blst_avx_p1,
        input: *const blst_byte,
    ) -> BLST_ERROR;
    pub fn blst_avx_p1_mult_pippenger_scratch_sizeof(npoints: usize) -> usize;
    pub fn blst_avx_p1_mult_pippenger_bucket_used_scratch_sizeof(
        npoints: usize,
    ) -> usize;
    pub fn blst_avx_p1_mult_pippenger_glv_points_scratch_sizeof(
        npoints: usize,
    ) -> usize;
    pub fn blst_avx_p1_mult_pippenger_glv_scalars_scratch_sizeof(
        npoints: usize,
    ) -> usize;
    pub fn blst_avx_p1_mult_pippenger(
        out: *mut blst_avx_p1,
        points: *const blst_avx_p1_affine,
        scalars: *const blst_byte,
        npoints: usize,
        nbits: usize,
        scratch: *mut core::ffi::c_void,
        scratch_bytes: usize,
        bucket_used: *mut blst_byte,
        bucket_used_bytes: usize,
    ) -> bool;
    pub fn blst_avx_p1_mult_pippenger_glv_with_scratch(
        out: *mut blst_avx_p1,
        points: *const blst_avx_p1_affine,
        scalars: *const blst_byte,
        npoints: usize,
        scratch: *mut core::ffi::c_void,
        scratch_bytes: usize,
        bucket_used: *mut blst_byte,
        bucket_used_bytes: usize,
        split_points: *mut blst_avx_p1_affine,
        split_points_bytes: usize,
        split_scalars: *mut blst_byte,
        split_scalars_bytes: usize,
    ) -> bool;
    pub fn blst_avx_p1_mult_pippenger_glv_window_scratch_sizeof(
        npoints: usize,
        window: usize,
    ) -> usize;
    pub fn blst_avx_p1_mult_pippenger_glv_window_bucket_used_scratch_sizeof(
        npoints: usize,
        window: usize,
    ) -> usize;
    pub fn blst_avx_p1_mult_pippenger_glv_window_with_scratch(
        out: *mut blst_avx_p1,
        points: *const blst_avx_p1_affine,
        scalars: *const blst_byte,
        npoints: usize,
        window: usize,
        scratch: *mut core::ffi::c_void,
        scratch_bytes: usize,
        bucket_used: *mut blst_byte,
        bucket_used_bytes: usize,
        split_points: *mut blst_avx_p1_affine,
        split_points_bytes: usize,
        split_scalars: *mut blst_byte,
        split_scalars_bytes: usize,
    ) -> bool;
    pub fn blst_avx_p1s_mult_pippenger_scratch_sizeof(npoints: usize) -> usize;
    pub fn blst_avx_p1s_mult_pippenger(
        out: *mut blst_avx_p1,
        points: *const blst_avx_p1,
        npoints: usize,
        scalars: *const blst_byte,
        nbits: usize,
        scratch: *mut blst_limb_t,
    ) -> bool;
    pub fn blst_avx_p1_cneg(p: *mut blst_avx_p1, cbit: bool) -> bool;
    pub fn blst_avx_p1_batch_subgroup_check(
        points: *const *const blst_avx_p1,
        count: usize,
        coefficients: *const blst_byte,
        width: usize,
        shift: *const blst_avx_p1,
    ) -> bool;
    pub fn blst_avx_p1s_in_g1(
        points: *const *const blst_avx_p1,
        count: usize,
        fill_random: blst_avx_random_bytes_fn,
        random_context: *mut core::ffi::c_void,
    ) -> bool;
    pub fn blst_avx_p1_identity(out: *mut blst_avx_p1) -> bool;
    pub fn blst_avx_p1_generator() -> *const blst_avx_p1;

    pub fn blst_avx_p2_add(
        out: *mut blst_avx_p2,
        a: *const blst_avx_p2,
        b: *const blst_avx_p2,
    ) -> bool;
    pub fn blst_avx_p2_from_scalar(
        out: *mut blst_avx_p2,
        input: *const blst_p2,
    ) -> bool;
    pub fn blst_p2_from_avx(
        out: *mut blst_p2,
        input: *const blst_avx_p2,
    ) -> bool;
    pub fn blst_avx_p2_affine_from_scalar(
        out: *mut blst_avx_p2_affine,
        input: *const blst_p2_affine,
    ) -> bool;
    pub fn blst_p2_affine_from_avx(
        out: *mut blst_p2_affine,
        input: *const blst_avx_p2_affine,
    ) -> bool;
    pub fn blst_avx_p2_from_affine(
        out: *mut blst_avx_p2,
        input: *const blst_avx_p2_affine,
    ) -> bool;
    pub fn blst_avx_p2_cneg(p: *mut blst_avx_p2, cbit: bool) -> bool;
    pub fn blst_avx_p2_batch_subgroup_check(
        points: *const *const blst_avx_p2,
        count: usize,
        coefficients: *const blst_byte,
        width: usize,
        shift: *const blst_avx_p2,
    ) -> bool;
    pub fn blst_avx_p2s_in_g2(
        points: *const *const blst_avx_p2,
        count: usize,
        fill_random: blst_avx_random_bytes_fn,
        random_context: *mut core::ffi::c_void,
        workspace: *mut core::ffi::c_void,
        workspace_bytes: usize,
    ) -> bool;
    pub fn blst_avx_p2s_in_g2_workspace_sizeof() -> usize;
    pub fn blst_avx_p2_mult(
        out: *mut blst_avx_p2,
        p: *const blst_avx_p2,
        scalar: *const blst_byte,
        nbits: usize,
    ) -> bool;
    pub fn blst_avx_p2_mult_pippenger_scratch_sizeof(npoints: usize) -> usize;
    pub fn blst_avx_p2_mult_pippenger_bucket_used_scratch_sizeof(
        npoints: usize,
    ) -> usize;
    pub fn blst_avx_p2_mult_pippenger(
        out: *mut blst_avx_p2,
        points: *const blst_avx_p2_affine,
        scalars: *const blst_byte,
        npoints: usize,
        nbits: usize,
        scratch: *mut core::ffi::c_void,
        scratch_bytes: usize,
        bucket_used: *mut blst_byte,
        bucket_used_bytes: usize,
    ) -> bool;
    pub fn blst_avx_p2_to_affine(
        out: *mut blst_avx_p2_affine,
        input: *const blst_avx_p2,
    ) -> bool;
    pub fn blst_avx_p2_in_g2(p: *const blst_avx_p2) -> bool;
    pub fn blst_avx_p2_is_inf(p: *const blst_avx_p2) -> bool;
    pub fn blst_avx_p2_is_equal(
        a: *const blst_avx_p2,
        b: *const blst_avx_p2,
    ) -> bool;
    pub fn blst_avx_p2_compress(
        out: *mut blst_byte,
        input: *const blst_avx_p2,
    ) -> bool;
    pub fn blst_avx_p2_serialize(
        out: *mut blst_byte,
        input: *const blst_avx_p2,
    ) -> bool;
    pub fn blst_avx_p2_uncompress(
        out: *mut blst_avx_p2,
        input: *const blst_byte,
    ) -> BLST_ERROR;
    pub fn blst_avx_p2_deserialize(
        out: *mut blst_avx_p2,
        input: *const blst_byte,
    ) -> BLST_ERROR;
    pub fn blst_avx_p2s_mult_pippenger_scratch_sizeof(npoints: usize) -> usize;
    pub fn blst_avx_p2s_mult_pippenger(
        out: *mut blst_avx_p2,
        points: *const blst_avx_p2,
        npoints: usize,
        scalars: *const blst_byte,
        nbits: usize,
        scratch: *mut blst_limb_t,
    ) -> bool;
    pub fn blst_avx_p2_identity(out: *mut blst_avx_p2) -> bool;
    pub fn blst_avx_p2_generator() -> *const blst_avx_p2;

    pub fn blst_avx_fp12_is_equal(
        a: *const blst_avx_fp12,
        b: *const blst_avx_fp12,
    ) -> bool;
    pub fn blst_avx_fp12_from_scalar(
        out: *mut blst_avx_fp12,
        input: *const blst_fp12,
    ) -> bool;
    pub fn blst_fp12_from_avx(
        out: *mut blst_fp12,
        input: *const blst_avx_fp12,
    ) -> bool;
    /// Serializes a validated AVX GT value in the same canonical coefficient
    /// order as [`crate::scalar::Gt::to_bytes_be`].
    ///
    /// The 576-byte output contains twelve canonical 48-byte big-endian Fp
    /// coefficients ordered by Fp2 index, then Fp6 half, then Fp component.
    /// Both pointers must be non-null. The operation allocates nothing.
    pub fn blst_avx_gt_serialize(
        out: *mut blst_byte,
        input: *const blst_avx_gt,
    ) -> bool;
    /// Strictly decodes the canonical CatBLST GT coefficient order.
    ///
    /// Every Fp coefficient must be canonical and the assembled value must lie
    /// in GT; the identity is accepted. Failure returns `false` and leaves
    /// `out` unchanged. Input may overlap output. Validation is variable-time
    /// with respect to the public bytes, and the operation allocates nothing.
    pub fn blst_avx_gt_deserialize(
        out: *mut blst_avx_gt,
        input: *const blst_byte,
    ) -> bool;
    /// Multiplies a validated AVX GT element by a little-endian scalar.
    ///
    /// Runtime and memory access are independent of a canonical scalar's
    /// contents for fixed public `nbits`. `nbits` must not exceed 256, and the
    /// operation allocates nothing. Exact output/input aliasing is supported.
    pub fn blst_avx_gt_mult(
        out: *mut blst_avx_gt,
        input: *const blst_avx_gt,
        scalar: *const blst_byte,
        nbits: usize,
    ) -> bool;
    pub fn blst_avx_fp12_mul(
        out: *mut blst_avx_fp12,
        a: *const blst_avx_fp12,
        b: *const blst_avx_fp12,
    ) -> bool;
    pub fn blst_avx_fp12_sqr(
        out: *mut blst_avx_fp12,
        input: *const blst_avx_fp12,
    ) -> bool;
    pub fn blst_avx_gt_inverse(
        out: *mut blst_avx_gt,
        input: *const blst_avx_gt,
    ) -> bool;
    pub fn blst_avx_fp12_final_exp(
        out: *mut blst_avx_fp12,
        input: *const blst_avx_fp12,
    ) -> bool;
    pub fn blst_avx_fp12_in_group(input: *const blst_avx_fp12) -> bool;
    pub fn blst_avx_fp12_one() -> *const blst_avx_fp12;
    pub fn blst_avx_gt_mult_pippenger_scratch_sizeof(npoints: usize) -> usize;
    pub fn blst_avx_gt_mult_pippenger_bucket_used_scratch_sizeof(
        npoints: usize,
    ) -> usize;
    pub fn blst_avx_gt_mult_pippenger(
        out: *mut blst_avx_gt,
        points: *const blst_avx_gt,
        scalars: *const blst_byte,
        npoints: usize,
        scratch: *mut core::ffi::c_void,
        scratch_bytes: usize,
        bucket_used: *mut blst_byte,
        bucket_used_bytes: usize,
    ) -> bool;
    pub fn blst_avx_pairing(
        out: *mut blst_avx_fp12,
        g1: *const blst_avx_p1,
        g2: *const blst_avx_p2,
    ) -> bool;
    pub fn blst_avx_pairing_n(
        out: *mut blst_avx_fp12,
        g1: *const blst_avx_p1_affine,
        g2: *const blst_avx_p2_affine,
        n: usize,
    ) -> bool;
    pub fn blst_avx_prepared_g2_lines_sizeof() -> usize;
    pub fn blst_avx_prepared_g2_lines_alignof() -> usize;
    pub fn blst_avx_prepared_g2_lines_init(
        storage: *mut core::ffi::c_void,
        storage_size: usize,
        g2: *const blst_avx_p2,
    ) -> bool;
    pub fn blst_avx_pairing_prepared_g2(
        out: *mut blst_avx_fp12,
        storage: *const core::ffi::c_void,
        storage_size: usize,
        g1: *const blst_avx_p1,
    ) -> bool;
    pub fn blst_avx_prepared_g2_batch_new(
        g2: *const blst_avx_p2,
        count: usize,
    ) -> *mut blst_avx_prepared_g2_batch;
    pub fn blst_avx_prepared_g2_batch_len(
        batch: *const blst_avx_prepared_g2_batch,
    ) -> usize;
    pub fn blst_avx_prepared_g2_batch_storage_bytes(
        batch: *const blst_avx_prepared_g2_batch,
    ) -> usize;
    pub fn blst_avx_prepared_g2_batch_free(
        batch: *mut blst_avx_prepared_g2_batch,
    );
    pub fn blst_avx_prepared_g2_batch_workspace_new(
        batch: *const blst_avx_prepared_g2_batch,
    ) -> *mut blst_avx_prepared_g2_batch_workspace;
    pub fn blst_avx_prepared_g2_batch_workspace_free(
        workspace: *mut blst_avx_prepared_g2_batch_workspace,
    );
    pub fn blst_avx_prepared_g2_batch_pairing_prefix(
        out: *mut blst_avx_gt,
        count: usize,
        g1: *const blst_avx_p1,
        batch: *const blst_avx_prepared_g2_batch,
        workspace: *mut blst_avx_prepared_g2_batch_workspace,
    ) -> bool;
    pub fn blst_avx_prepared_g2_batch_pairing_indices(
        out: *mut blst_avx_gt,
        indices: *const usize,
        count: usize,
        g1: *const blst_avx_p1,
        batch: *const blst_avx_prepared_g2_batch,
        workspace: *mut blst_avx_prepared_g2_batch_workspace,
    ) -> bool;
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AvxMsmError {
    AvxUnavailable,
    ScalarLengthMismatch,
    ScratchSizeOverflow,
    AllocationFailed,
    ConversionFailed,
    BackendFailure,
}

fn avx_fp12_storage_len(bytes: usize) -> usize {
    let element_size = size_of::<blst_avx_fp12>();
    let whole_elements = bytes / element_size;
    whole_elements + usize::from(whole_elements * element_size != bytes)
}

/// Computes a variable-time AVX GLV4 Pippenger MSM over validated GT elements.
///
/// Every point must already be a GT element; this function does not validate
/// membership. `scalars` must contain exactly one 256-bit little-endian scalar
/// per point. Scalar contents are public and are interpreted modulo r. Debug
/// builds assert GT membership before entering the native backend. Temporary
/// allocation failure is reported as [`AvxMsmError::AllocationFailed`].
pub fn gt_msm(
    points: &[blst_fp12],
    scalars: &[u8],
) -> Result<blst_fp12, AvxMsmError> {
    let scalar_bytes = points
        .len()
        .checked_mul(32)
        .ok_or(AvxMsmError::ScratchSizeOverflow)?;
    if scalars.len() != scalar_bytes {
        return Err(AvxMsmError::ScalarLengthMismatch);
    }

    if !unsafe { blst_avx_enabled() } {
        return Err(AvxMsmError::AvxUnavailable);
    }

    #[cfg(debug_assertions)]
    for point in points {
        debug_assert!(
            unsafe { blst_fp12_in_group(point) },
            "gt_msm requires every input to be a validated GT element"
        );
    }

    if points.is_empty() {
        return Ok(unsafe { *blst_gt_identity() });
    }

    let mut avx_points = Vec::new();
    avx_points
        .try_reserve_exact(points.len())
        .map_err(|_| AvxMsmError::AllocationFailed)?;
    for point in points {
        let mut avx_point = MaybeUninit::<blst_avx_fp12>::uninit();
        let ok =
            unsafe { blst_avx_fp12_from_scalar(avx_point.as_mut_ptr(), point) };
        if !ok {
            return Err(AvxMsmError::ConversionFailed);
        }
        unsafe { avx_points.push(avx_point.assume_init()) };
    }

    let scratch_bytes =
        unsafe { blst_avx_gt_mult_pippenger_scratch_sizeof(points.len()) };
    let bucket_used_bytes = unsafe {
        blst_avx_gt_mult_pippenger_bucket_used_scratch_sizeof(points.len())
    };
    if scratch_bytes == 0 || bucket_used_bytes == 0 {
        return Err(AvxMsmError::ScratchSizeOverflow);
    }

    let scratch_len = avx_fp12_storage_len(scratch_bytes);
    // The foreign call owns these temporary writes. Keep the Rust lengths at
    // zero so no uninitialized element is ever claimed as an initialized value.
    let mut scratch = Vec::<MaybeUninit<blst_avx_fp12>>::new();
    scratch
        .try_reserve_exact(scratch_len)
        .map_err(|_| AvxMsmError::AllocationFailed)?;
    let mut bucket_used = Vec::<MaybeUninit<u8>>::new();
    bucket_used
        .try_reserve_exact(bucket_used_bytes)
        .map_err(|_| AvxMsmError::AllocationFailed)?;

    let mut avx_out = MaybeUninit::<blst_avx_fp12>::uninit();
    let ok = unsafe {
        blst_avx_gt_mult_pippenger(
            avx_out.as_mut_ptr(),
            avx_points.as_ptr(),
            scalars.as_ptr(),
            points.len(),
            scratch.as_mut_ptr().cast(),
            scratch_bytes,
            bucket_used.as_mut_ptr().cast(),
            bucket_used_bytes,
        )
    };
    if !ok {
        return Err(AvxMsmError::BackendFailure);
    }

    let avx_out = unsafe { avx_out.assume_init() };
    let mut out = MaybeUninit::<blst_fp12>::uninit();
    let ok = unsafe { blst_fp12_from_avx(out.as_mut_ptr(), &avx_out) };
    if ok {
        Ok(unsafe { out.assume_init() })
    } else {
        Err(AvxMsmError::ConversionFailed)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn avx_fp12_storage_len_handles_boundaries_without_overflow() {
        let element_size = size_of::<blst_avx_fp12>();
        assert!(element_size > 1);
        assert_eq!(avx_fp12_storage_len(0), 0);
        assert_eq!(avx_fp12_storage_len(1), 1);
        assert_eq!(avx_fp12_storage_len(element_size), 1);
        assert_eq!(avx_fp12_storage_len(element_size + 1), 2);

        let maximum = usize::MAX;
        let expected =
            maximum / element_size + usize::from(maximum % element_size != 0);
        assert_eq!(avx_fp12_storage_len(maximum), expected);
    }

    #[test]
    fn rust_layout_matches_avx_abi() {
        assert_eq!(core::mem::size_of::<blst_avx_p1>(), unsafe {
            blst_avx_p1_sizeof()
        });
        assert!(
            core::mem::align_of::<blst_avx_p1>()
                >= unsafe { blst_avx_p1_alignof() }
        );
        assert_eq!(core::mem::size_of::<blst_avx_p1_affine>(), unsafe {
            blst_avx_p1_affine_sizeof()
        });
        assert!(
            core::mem::align_of::<blst_avx_p1_affine>()
                >= unsafe { blst_avx_p1_affine_alignof() }
        );
        assert_eq!(core::mem::size_of::<blst_avx_p2>(), unsafe {
            blst_avx_p2_sizeof()
        });
        assert!(
            core::mem::align_of::<blst_avx_p2>()
                >= unsafe { blst_avx_p2_alignof() }
        );
        assert_eq!(core::mem::size_of::<blst_avx_p2_affine>(), unsafe {
            blst_avx_p2_affine_sizeof()
        });
        assert!(
            core::mem::align_of::<blst_avx_p2_affine>()
                >= unsafe { blst_avx_p2_affine_alignof() }
        );
        assert_eq!(core::mem::size_of::<blst_avx_fp12>(), unsafe {
            blst_avx_fp12_sizeof()
        });
        assert!(
            core::mem::align_of::<blst_avx_fp12>()
                >= unsafe { blst_avx_fp12_alignof() }
        );
        if unsafe { blst_avx_enabled() } {
            assert_eq!(
                unsafe { blst_avx_prepared_g2_lines_sizeof() },
                68 * 1024
            );
            assert_eq!(unsafe { blst_avx_prepared_g2_lines_alignof() }, 64);
        }
    }

    #[test]
    fn native_gt_canonical_roundtrip() {
        if !unsafe { blst_avx_enabled() } {
            return;
        }

        let mut value = MaybeUninit::<blst_avx_gt>::uninit();
        assert!(unsafe {
            blst_avx_pairing(
                value.as_mut_ptr(),
                blst_avx_p1_generator(),
                blst_avx_p2_generator(),
            )
        });
        let value = unsafe { value.assume_init() };
        assert!(unsafe { blst_avx_fp12_is_equal(&value, &value) });
        assert!(unsafe { blst_avx_fp12_in_group(&value) });
        let mut encoded = [0u8; 576];
        assert!(unsafe { blst_avx_gt_serialize(encoded.as_mut_ptr(), &value) });
        let mut decoded = MaybeUninit::<blst_avx_gt>::uninit();
        assert!(unsafe {
            blst_avx_gt_deserialize(decoded.as_mut_ptr(), encoded.as_ptr())
        });
        let decoded = unsafe { decoded.assume_init() };
        assert!(unsafe { blst_avx_fp12_is_equal(&value, &decoded) });

        let mut inverse = core::mem::MaybeUninit::<blst_avx_gt>::uninit();
        assert!(unsafe { blst_avx_gt_inverse(inverse.as_mut_ptr(), &value) });
        let inverse = unsafe { inverse.assume_init() };
        let mut product = core::mem::MaybeUninit::<blst_avx_gt>::uninit();
        assert!(unsafe {
            blst_avx_fp12_mul(product.as_mut_ptr(), &value, &inverse)
        });
        let product = unsafe { product.assume_init() };
        assert!(unsafe {
            blst_avx_fp12_is_equal(&product, blst_avx_fp12_one())
        });

        let invalid = [0u8; 576];
        let mut invalid_out = MaybeUninit::<blst_avx_gt>::uninit();
        assert!(!unsafe {
            blst_avx_gt_deserialize(invalid_out.as_mut_ptr(), invalid.as_ptr())
        });
    }

    #[test]
    fn native_batched_pairing_matches_individual_product() {
        if !unsafe { blst_avx_enabled() } {
            return;
        }
        let mut g1_affine = [blst_avx_p1_affine::default(); 2];
        let mut g2_affine = [blst_avx_p2_affine::default(); 2];
        for point in &mut g1_affine {
            assert!(unsafe {
                blst_avx_p1_to_affine(point, blst_avx_p1_generator())
            });
        }
        for point in &mut g2_affine {
            assert!(unsafe {
                blst_avx_p2_to_affine(point, blst_avx_p2_generator())
            });
        }
        let mut batched = MaybeUninit::<blst_avx_gt>::uninit();
        assert!(unsafe {
            blst_avx_pairing_n(
                batched.as_mut_ptr(),
                g1_affine.as_ptr(),
                g2_affine.as_ptr(),
                2,
            )
        });
        let batched = unsafe { batched.assume_init() };
        let mut single = MaybeUninit::<blst_avx_gt>::uninit();
        assert!(unsafe {
            blst_avx_pairing(
                single.as_mut_ptr(),
                blst_avx_p1_generator(),
                blst_avx_p2_generator(),
            )
        });
        let single = unsafe { single.assume_init() };
        let mut expected = MaybeUninit::<blst_avx_gt>::uninit();
        assert!(unsafe {
            blst_avx_fp12_mul(expected.as_mut_ptr(), &single, &single)
        });
        assert!(unsafe {
            blst_avx_fp12_is_equal(&batched, &expected.assume_init())
        });

        let mut empty = MaybeUninit::<blst_avx_gt>::uninit();
        assert!(unsafe {
            blst_avx_pairing_n(empty.as_mut_ptr(), ptr::null(), ptr::null(), 0)
        });
        assert!(unsafe {
            blst_avx_fp12_is_equal(&empty.assume_init(), blst_avx_fp12_one())
        });
    }

    #[test]
    fn native_prepared_pairing_matches_direct_pairing() {
        if !unsafe { blst_avx_enabled() } {
            return;
        }

        let g1 = unsafe { &*blst_avx_p1_generator() };
        let g2 = unsafe { &*blst_avx_p2_generator() };
        let prepared = NativeAvxPreparedG2Lines::new(g2).unwrap();
        let mut actual = blst_avx_fp12::default();
        prepared.pairing_into(&mut actual, g1).unwrap();

        let mut expected = blst_avx_fp12::default();
        assert!(unsafe { blst_avx_pairing(&mut expected, g1, g2) });
        assert!(unsafe { blst_avx_fp12_is_equal(&actual, &expected) });
    }

    #[test]
    fn prepared_g2_batch_matches_direct_pairings_for_prefixes() {
        if !is_available() {
            return;
        }

        let g1 = G1Projective::generator().unwrap();
        let generator = G2Projective::generator().unwrap();
        let doubled = generator.multiply(&Scalar::from(2u64)).unwrap();
        let identity = G2Projective::identity().unwrap();
        let points = [generator, doubled, identity];
        let batch = PreparedG2Batch::new(&points).unwrap();
        assert_eq!(batch.len(), points.len());
        assert!(batch.storage_bytes() >= points.len() * 68 * 1024);

        let mut workspace = PreparedG2BatchWorkspace::new(&batch).unwrap();
        let gt_identity = Gt::identity().unwrap();
        let mut out = [gt_identity; 3];
        workspace
            .pairing_prefix_into(&g1, &batch, &mut out[..2])
            .unwrap();
        assert_eq!(out[0], Gt::pairing(&g1, &points[0]).unwrap());
        assert_eq!(out[1], Gt::pairing(&g1, &points[1]).unwrap());

        workspace
            .pairing_prefix_into(&g1, &batch, &mut out)
            .unwrap();
        assert_eq!(out[2], gt_identity);
        workspace
            .pairing_indices_into(&g1, &batch, &[1, 0], &mut out[..2])
            .unwrap();
        assert_eq!(out[0], Gt::pairing(&g1, &points[1]).unwrap());
        assert_eq!(out[1], Gt::pairing(&g1, &points[0]).unwrap());
        workspace.pairing_prefix_into(&g1, &batch, &mut []).unwrap();

        let g1_identity = G1Projective::identity().unwrap();
        workspace
            .pairing_prefix_into(&g1_identity, &batch, &mut out)
            .unwrap();
        assert_eq!(out, [gt_identity; 3]);
    }

    #[test]
    fn prepared_g2_batch_workspace_rejects_an_oversized_prefix() {
        if !is_available() {
            return;
        }

        let g1 = G1Projective::generator().unwrap();
        let batch = PreparedG2Batch::new(&[G2Projective::generator().unwrap()])
            .unwrap();
        let mut workspace = PreparedG2BatchWorkspace::new(&batch).unwrap();
        let identity = Gt::identity().unwrap();
        let mut out = [identity; 2];
        assert_eq!(
            workspace.pairing_prefix_into(&g1, &batch, &mut out),
            Err(Error::PairingCapacityExceeded {
                required: 2,
                capacity: 1,
            })
        );
        assert_eq!(
            workspace.pairing_indices_into(&g1, &batch, &[0], &mut out),
            Err(Error::PairingLengthMismatch {
                indices: 1,
                outputs: 2,
            })
        );
        assert_eq!(
            workspace.pairing_indices_into(&g1, &batch, &[1], &mut out[..1]),
            Err(Error::PreparedIndexOutOfRange { index: 1, len: 1 })
        );
        assert!(matches!(
            PreparedG2Batch::new(&[]),
            Err(Error::EmptyPreparedBatch)
        ));
    }

    #[test]
    fn gt_msm_rejects_bad_scalar_length() {
        let points = [blst_fp12::default()];
        assert_eq!(
            gt_msm(&points, &[]),
            Err(AvxMsmError::ScalarLengthMismatch)
        );
    }

    #[test]
    fn gt_msm_identity_inputs_return_identity() {
        if unsafe { !blst_avx_enabled() } {
            return;
        }

        let identity = unsafe { *blst_gt_identity() };
        let points = [identity, identity];
        let mut scalars = [0u8; 64];
        scalars[0] = 3;
        scalars[32] = 5;

        assert_eq!(gt_msm(&points, &scalars).unwrap(), identity);
    }

    #[test]
    fn gt_msm_empty_returns_identity() {
        if unsafe { !blst_avx_enabled() } {
            return;
        }

        assert_eq!(gt_msm(&[], &[]).unwrap(), unsafe { *blst_gt_identity() });
    }

    #[test]
    fn opaque_equality_is_backend_independent_storage_equality() {
        let p1 = blst_avx_p1::default();
        let p2 = blst_avx_p2::default();
        let fp12 = blst_avx_fp12::default();
        assert_eq!(p1, p1);
        assert_eq!(p2, p2);
        assert_eq!(fp12, fp12);

        let mut different = fp12;
        different.b[0] = 1;
        assert_ne!(fp12, different);
    }

    #[test]
    fn typed_curve_roundtrips_and_matches_scalar_wire_encodings() {
        if !is_available() {
            return;
        }

        let g1 = G1Projective::generator().unwrap();
        let g2 = G2Projective::generator().unwrap();
        let g1_bytes = g1.to_compressed().unwrap();
        let g2_bytes = g2.to_compressed().unwrap();
        assert_eq!(G1Projective::from_compressed(&g1_bytes).unwrap(), g1);
        assert_eq!(G2Projective::from_compressed(&g2_bytes).unwrap(), g2);
        assert_eq!(
            g1_bytes,
            crate::scalar::G1Affine::generator().to_compressed()
        );
        assert_eq!(
            g2_bytes,
            crate::scalar::G2Affine::generator().to_compressed()
        );

        let gt = Gt::pairing(&g1, &g2).unwrap();
        let gt_bytes = gt.to_bytes_be().unwrap();
        assert_eq!(Gt::from_bytes_be(&gt_bytes).unwrap(), gt);
        assert_eq!(
            gt_bytes,
            crate::scalar::pairing(
                &crate::scalar::G1Affine::generator(),
                &crate::scalar::G2Affine::generator(),
            )
            .to_bytes_be()
        );
    }

    #[test]
    fn typed_gt_pow_matches_constant_time_scalar_backend() {
        if !is_available() {
            return;
        }

        let scalar = crate::Scalar::from(7u64);
        let native = Gt::pairing(
            &G1Projective::generator().unwrap(),
            &G2Projective::generator().unwrap(),
        )
        .unwrap()
        .pow(&scalar)
        .unwrap();
        let portable = crate::scalar::pairing(
            &crate::scalar::G1Affine::generator(),
            &crate::scalar::G2Affine::generator(),
        ) * scalar;
        assert_eq!(native.to_bytes_be().unwrap(), portable.to_bytes_be());
    }

    #[test]
    fn typed_workspaces_reuse_fixed_capacity() {
        if !is_available() {
            return;
        }

        let g1 = G1Projective::generator().unwrap();
        let g2 = G2Projective::generator().unwrap();
        let mut pairing = PairingWorkspace::for_pairs(1).unwrap();
        pairing.push(&g1, &g2).unwrap();
        assert_eq!(pairing.evaluate().unwrap(), Gt::pairing(&g1, &g2).unwrap());
        assert_eq!(
            pairing.push(&g1, &g2),
            Err(Error::PairingCapacityExceeded {
                required: 2,
                capacity: 1,
            })
        );
        pairing.clear();
        pairing.push(&g1, &g2).unwrap();

        let mut msm = G1MsmWorkspace::for_points(1).unwrap();
        let capacity = msm.capacity_bytes();
        let mut out = G1Projective::identity().unwrap();
        let one = crate::Scalar::from(1u64).to_bytes_le();
        msm.msm_into(&mut out, &[g1], &[one]).unwrap();
        assert_eq!(out, g1);
        assert_eq!(msm.capacity_bytes(), capacity);
    }

    #[test]
    fn secret_g1_msm_matches_independent_constant_time_products() {
        if !is_available() {
            return;
        }

        let generator = G1Projective::generator().unwrap();
        let points = (1u64..=65)
            .map(|factor| generator.multiply(&Scalar::from(factor)).unwrap())
            .collect::<Vec<_>>();
        let scalars = (19u64..=83).map(Scalar::from).collect::<Vec<_>>();
        let mut workspace = G1MsmWorkspace::for_points(points.len()).unwrap();

        for count in [0usize, 1, 7, 8, 15, 16, 31, 32, 33, 63, 64, 65] {
            let mut expected = G1Projective::identity().unwrap();
            for (point, scalar) in points[..count].iter().zip(&scalars[..count])
            {
                expected =
                    expected.add(&point.multiply(scalar).unwrap()).unwrap();
            }

            let mut actual = G1Projective::identity().unwrap();
            workspace
                .secret_msm_into(
                    &mut actual,
                    &points[..count],
                    &scalars[..count],
                )
                .unwrap();
            assert_eq!(actual, expected, "count={count}");
        }

        let mut out = G1Projective::identity().unwrap();
        assert_eq!(
            workspace.secret_msm_into(&mut out, &points[..1], &scalars[..0]),
            Err(Error::MsmLengthMismatch {
                points: 1,
                scalars: 0,
            })
        );

        let mut undersized = G1MsmWorkspace::for_points(7).unwrap();
        let available = undersized.capacity_bytes();
        let required = unsafe { blst_avx_p1s_mult_ct_scratch_sizeof() };
        assert!(available < required);
        assert_eq!(
            undersized.secret_msm_into(&mut out, &points[..8], &scalars[..8]),
            Err(Error::ScratchTooSmall {
                required,
                available,
            })
        );
    }

    #[cfg(debug_assertions)]
    #[test]
    fn gt_msm_debug_diagnoses_non_gt_input() {
        if unsafe { !blst_avx_enabled() } {
            return;
        }

        let invalid: blst_fp12 = unsafe { core::mem::zeroed() };
        let result = std::panic::catch_unwind(|| gt_msm(&[invalid], &[0; 32]));
        assert!(result.is_err());
    }
}
