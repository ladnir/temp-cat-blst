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

use crate::*;
use alloc::vec::Vec;
use core::mem::{size_of, MaybeUninit};
use core::ptr;

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

/// Tests whether every AVX point belongs to G2 using CatBLST-owned challenges.
///
/// The random-seed requirements and empty-input behavior match
/// [`avx_p1s_in_g1`].
pub fn avx_p2s_in_g2(
    points: &[&blst_avx_p2],
    mut random_seed: [blst_byte; 32],
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
    ) -> bool;
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
    pub fn blst_avx_gt_serialize(
        out: *mut blst_byte,
        input: *const blst_avx_gt,
    ) -> bool;
    pub fn blst_avx_gt_deserialize(
        out: *mut blst_avx_gt,
        input: *const blst_byte,
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
        assert!(unsafe { blst_avx_fp12_mul(product.as_mut_ptr(), &value, &inverse) });
        let product = unsafe { product.assume_init() };
        assert!(unsafe { blst_avx_fp12_is_equal(&product, blst_avx_fp12_one()) });

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
            assert!(unsafe { blst_avx_p1_to_affine(point, blst_avx_p1_generator()) });
        }
        for point in &mut g2_affine {
            assert!(unsafe { blst_avx_p2_to_affine(point, blst_avx_p2_generator()) });
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
        assert!(unsafe { blst_avx_fp12_mul(expected.as_mut_ptr(), &single, &single) });
        assert!(unsafe { blst_avx_fp12_is_equal(&batched, &expected.assume_init()) });

        let mut empty = MaybeUninit::<blst_avx_gt>::uninit();
        assert!(unsafe { blst_avx_pairing_n(empty.as_mut_ptr(), ptr::null(), ptr::null(), 0) });
        assert!(unsafe { blst_avx_fp12_is_equal(&empty.assume_init(), blst_avx_fp12_one()) });
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
