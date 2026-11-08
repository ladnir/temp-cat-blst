// Copyright Category Labs
// Licensed under the Apache License, Version 2.0, see LICENSE for details.
// SPDX-License-Identifier: Apache-2.0

//! Lean, `blstrs`-inspired Rust types over category-blst's scalar backend.
//!
//! This module provides the typed curve surface needed by protocol code without
//! importing the `ff`, `group`, or `pairing` trait ecosystems. Every wrapper is
//! transparent over the corresponding native scalar-ABI object. Operations are
//! synchronous, allocation-free, and expose raw borrows for zero-conversion
//! kernel interoperation.
//!
//! # Safe-type invariants
//!
//! - [`Scalar`] contains a reduced element of the BLS12-381 scalar field in
//!   Montgomery form. Zero is valid.
//! - G1 and G2 wrappers contain an on-curve point in the prime-order subgroup.
//!   The group identity is valid.
//! - [`crate::scalar::MillerLoopResult`] contains a nonzero product of
//!   Miller-loop outputs from validated subgroup points. It is not a
//!   [`crate::scalar::Gt`] value before final exponentiation.
//! - [`crate::scalar::Gt`] contains an element of the BLS12-381 pairing target
//!   group. The target-group identity is valid.
//!
//! Safe constructors establish these invariants, and safe arithmetic preserves
//! them. Raw mutable access and unchecked raw construction transfer the stated
//! invariant to the caller. A caller must restore the invariant before a raw
//! mutable borrow ends. Passing an invalid wrapper to a safe method can violate
//! the preconditions of the underlying native operation.
//!
//! # Timing and storage
//!
//! Scalar-field addition, subtraction, multiplication, and negation, plus the
//! unsuffixed point and target-group multiplication methods, are intended to be
//! constant-time with respect to valid scalar and group values. Scalar
//! multiplication uses a fixed 256-bit width. Methods with a `_vartime` suffix
//! require public scalars. No constant-time guarantee is made for validation,
//! decoding, encoding, equality, zero or identity tests, affine conversion,
//! pairing, or other methods unless their documentation says otherwise.
//!
//! Arithmetic operations use stack or caller-owned storage and allocate
//! nothing. [`PreparedG2Batch`] is the deliberate exception: construction owns
//! immutable native line tables on the heap, while every evaluation remains
//! allocation-free. Independent operations on disjoint mutable values are safe
//! to run concurrently. Immutable values, including a prepared batch, may be
//! shared under Rust's ordinary aliasing rules. Safe methods enforce their
//! aliasing requirements through Rust borrows; an unsafe caller passing raw
//! pointers to a native kernel remains responsible for that kernel's aliasing
//! contract.
//!
//! # Encoding boundary
//!
//! Scalar decoders distinguish strict canonical decoding from wide modular
//! reduction. Point decoders check the encoding, curve equation, and prime-order
//! subgroup; they accept the identity. Target-group decoding checks all twelve
//! base-field coefficients and target-group membership.
//!
//! The types in this module are not representation-compatible with the native
//! types in [`crate::avx`]. A protocol instance should retain one backend's
//! native point and target-group types. Canonical encodings are the compatibility
//! boundary between scalar and AVX implementations.

use crate::*;
use core::fmt;
use core::iter::Sum;
use core::mem::MaybeUninit;
use core::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use core::ptr::NonNull;

#[path = "scalar_msm.rs"]
mod scalar_msm;
pub use scalar_msm::{
    CtMsmError, G1CtMsmWorkspace, G1FixedCtMsmTable, G2CtMsmWorkspace,
};

/// A reduced scalar-field element in category-blst's native Montgomery form.
///
/// Let `r` denote the BLS12-381 scalar-field order. The represented value lies
/// in `Fr = Z/rZ`; its native limbs contain the reduced Montgomery encoding of
/// that value. This type permits zero and is distinct from [`blst_scalar`], the
/// canonical non-Montgomery byte representation used by the raw secret-key API.
///
/// Addition, subtraction, multiplication, and negation are intended to be
/// constant-time with respect to field values. These operations allocate
/// nothing. [`Default`] returns [`Scalar::ZERO`], and sums of empty iterators
/// return zero. [`Debug`](fmt::Debug) deliberately redacts the scalar value.
#[repr(transparent)]
#[derive(Copy, Clone, PartialEq, Eq)]
pub struct Scalar(blst_fr);

impl Scalar {
    /// The additive identity `0` in `Fr`.
    pub const ZERO: Self = Self(blst_fr { l: [0; 4] });

    /// The multiplicative identity `1` in `Fr`.
    ///
    /// The stored limbs are the category-blst Montgomery representation, not the
    /// conventional integer limbs of one.
    pub const ONE: Self = Self(blst_fr {
        l: [
            0x0000_0001_ffff_fffe,
            0x5884_b7fa_0003_4802,
            0x998c_4fef_ecbc_4ff5,
            0x1824_b159_acc5_056f,
        ],
    });

    /// Borrows the reduced native Montgomery-form value without conversion.
    ///
    /// The returned reference is valid for the lifetime of `self`. Reading it
    /// does not alter the scalar invariant and allocates nothing.
    #[inline]
    pub const fn as_raw(&self) -> &blst_fr {
        &self.0
    }

    /// Mutably borrows the native Montgomery-form value without conversion.
    ///
    /// This accessor is intended for native kernels that write one reduced
    /// `blst_fr` result into caller-owned storage. It allocates and converts
    /// nothing.
    ///
    /// # Safety
    ///
    /// Before the borrow ends, the caller must write a reduced Montgomery-form
    /// element of `Fr`. Uninitialized, unreduced, or conventional-form limbs do
    /// not satisfy [`Scalar`]'s invariant.
    #[inline]
    pub unsafe fn as_raw_mut(&mut self) -> &mut blst_fr {
        &mut self.0
    }

    /// Wraps a reduced native Montgomery-form scalar without conversion.
    ///
    /// This constructor performs no validation, reduction, or allocation.
    ///
    /// # Safety
    ///
    /// `raw` must contain a reduced Montgomery-form element of `Fr`. A canonical
    /// integer or [`blst_scalar`] representation is not a valid input.
    #[inline]
    pub const unsafe fn from_raw_unchecked(raw: blst_fr) -> Self {
        Self(raw)
    }

    /// Consumes this wrapper and returns its reduced native Montgomery-form value.
    ///
    /// This operation copies no additional representation and allocates nothing.
    #[inline]
    pub const fn into_raw(self) -> blst_fr {
        self.0
    }

    /// Returns `true` exactly when this scalar represents `0` in `Fr`.
    ///
    /// This observation allocates nothing. No constant-time guarantee is made.
    #[inline]
    pub fn is_zero(&self) -> bool {
        self.0 == Self::ZERO.0
    }

    /// Strictly decodes a canonical 32-byte big-endian scalar.
    ///
    /// The byte string represents an unsigned integer. This method returns
    /// `Some` exactly when that integer lies in `[0, r)`, including zero. It does
    /// not reduce a noncanonical integer. Success converts the value to native
    /// Montgomery form; failure returns `None` and exposes no partial value.
    /// The operation allocates nothing.
    #[inline]
    pub fn from_bytes_be(bytes: &[u8; 32]) -> Option<Self> {
        unsafe {
            let mut encoded = MaybeUninit::<blst_scalar>::uninit();
            blst_scalar_from_bendian(encoded.as_mut_ptr(), bytes.as_ptr());
            let encoded = encoded.assume_init();
            if !blst_scalar_fr_check(&encoded) {
                return None;
            }
            let mut out = MaybeUninit::<blst_fr>::uninit();
            blst_fr_from_scalar(out.as_mut_ptr(), &encoded);
            Some(Self(out.assume_init()))
        }
    }

    /// Strictly decodes a canonical 32-byte little-endian scalar.
    ///
    /// The byte string represents an unsigned integer. This method returns
    /// `Some` exactly when that integer lies in `[0, r)`, including zero. It does
    /// not reduce a noncanonical integer. Success converts the value to native
    /// Montgomery form; failure returns `None` and exposes no partial value.
    /// The operation allocates nothing.
    #[inline]
    pub fn from_bytes_le(bytes: &[u8; 32]) -> Option<Self> {
        unsafe {
            let mut encoded = MaybeUninit::<blst_scalar>::uninit();
            blst_scalar_from_lendian(encoded.as_mut_ptr(), bytes.as_ptr());
            let encoded = encoded.assume_init();
            if !blst_scalar_fr_check(&encoded) {
                return None;
            }
            let mut out = MaybeUninit::<blst_fr>::uninit();
            blst_fr_from_scalar(out.as_mut_ptr(), &encoded);
            Some(Self(out.assume_init()))
        }
    }

    /// Reduces a 64-byte big-endian unsigned integer modulo `r`.
    ///
    /// Every 64-byte input succeeds, and zero is a valid result. This method is
    /// deterministic modular reduction, not strict decoding or rejection
    /// sampling. A protocol should use it only when its scalar-derivation rule
    /// specifies wide reduction. The operation allocates nothing.
    #[inline]
    pub fn from_bytes_wide_be(bytes: &[u8; 64]) -> Self {
        unsafe {
            let mut encoded = MaybeUninit::<blst_scalar>::uninit();
            let _ = blst_scalar_from_be_bytes(
                encoded.as_mut_ptr(),
                bytes.as_ptr(),
                bytes.len(),
            );
            let encoded = encoded.assume_init();
            let mut out = MaybeUninit::<blst_fr>::uninit();
            blst_fr_from_scalar(out.as_mut_ptr(), &encoded);
            Self(out.assume_init())
        }
    }

    /// Reduces a 64-byte little-endian unsigned integer modulo `r`.
    ///
    /// Every 64-byte input succeeds, and zero is a valid result. This method is
    /// deterministic modular reduction, not strict decoding or rejection
    /// sampling. A protocol should use it only when its scalar-derivation rule
    /// specifies wide reduction. The operation allocates nothing.
    #[inline]
    pub fn from_bytes_wide_le(bytes: &[u8; 64]) -> Self {
        unsafe {
            let mut encoded = MaybeUninit::<blst_scalar>::uninit();
            let _ = blst_scalar_from_le_bytes(
                encoded.as_mut_ptr(),
                bytes.as_ptr(),
                bytes.len(),
            );
            let encoded = encoded.assume_init();
            let mut out = MaybeUninit::<blst_fr>::uninit();
            blst_fr_from_scalar(out.as_mut_ptr(), &encoded);
            Self(out.assume_init())
        }
    }

    /// Encodes this scalar as the unique 32-byte big-endian integer in `[0, r)`.
    ///
    /// This operation converts from Montgomery form and allocates nothing.
    #[inline]
    pub fn to_bytes_be(&self) -> [u8; 32] {
        unsafe {
            let mut encoded = MaybeUninit::<blst_scalar>::uninit();
            blst_scalar_from_fr(encoded.as_mut_ptr(), &self.0);
            let encoded = encoded.assume_init();
            let mut out = MaybeUninit::<[u8; 32]>::uninit();
            blst_bendian_from_scalar(out.as_mut_ptr().cast::<u8>(), &encoded);
            out.assume_init()
        }
    }

    /// Encodes this scalar as the unique 32-byte little-endian integer in `[0, r)`.
    ///
    /// This operation converts from Montgomery form and allocates nothing.
    #[inline]
    pub fn to_bytes_le(&self) -> [u8; 32] {
        unsafe {
            let mut encoded = MaybeUninit::<blst_scalar>::uninit();
            blst_scalar_from_fr(encoded.as_mut_ptr(), &self.0);
            let encoded = encoded.assume_init();
            let mut out = MaybeUninit::<[u8; 32]>::uninit();
            blst_lendian_from_scalar(out.as_mut_ptr().cast::<u8>(), &encoded);
            out.assume_init()
        }
    }

    /// Computes the multiplicative inverse in `Fr`.
    ///
    /// For nonzero `self`, the returned value is the unique scalar `x` such that
    /// `self * x == Scalar::ONE`. Zero returns `None` before calling the native
    /// inverse, whose domain excludes zero. The native inversion on the nonzero
    /// path is intended to be constant-time, but this method necessarily
    /// exposes whether the input is zero through its return value and initial
    /// zero test. It allocates nothing.
    #[inline]
    pub fn invert(&self) -> Option<Self> {
        if self.is_zero() {
            return None;
        }
        unsafe {
            let mut out = MaybeUninit::<blst_fr>::uninit();
            blst_fr_inverse(out.as_mut_ptr(), &self.0);
            Some(Self(out.assume_init()))
        }
    }

    #[inline]
    fn add_inner(lhs: &Self, rhs: &Self) -> Self {
        unsafe {
            let mut out = MaybeUninit::<blst_fr>::uninit();
            blst_fr_add(out.as_mut_ptr(), &lhs.0, &rhs.0);
            Self(out.assume_init())
        }
    }

    #[inline]
    fn sub_inner(lhs: &Self, rhs: &Self) -> Self {
        unsafe {
            let mut out = MaybeUninit::<blst_fr>::uninit();
            blst_fr_sub(out.as_mut_ptr(), &lhs.0, &rhs.0);
            Self(out.assume_init())
        }
    }

    #[inline]
    fn mul_inner(lhs: &Self, rhs: &Self) -> Self {
        unsafe {
            let mut out = MaybeUninit::<blst_fr>::uninit();
            blst_fr_mul(out.as_mut_ptr(), &lhs.0, &rhs.0);
            Self(out.assume_init())
        }
    }
}

impl Default for Scalar {
    #[inline]
    fn default() -> Self {
        Self::ZERO
    }
}

impl fmt::Debug for Scalar {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str("Scalar(..)")
    }
}

impl From<u64> for Scalar {
    #[inline]
    fn from(value: u64) -> Self {
        unsafe {
            let words = [value, 0, 0, 0];
            let mut out = MaybeUninit::<blst_fr>::uninit();
            blst_fr_from_uint64(out.as_mut_ptr(), words.as_ptr());
            Self(out.assume_init())
        }
    }
}

macro_rules! impl_scalar_binary_op {
    ($trait:ident, $method:ident, $inner:ident) => {
        impl $trait for Scalar {
            type Output = Scalar;
            #[inline]
            fn $method(self, rhs: Scalar) -> Scalar {
                Scalar::$inner(&self, &rhs)
            }
        }

        impl<'a> $trait<&'a Scalar> for Scalar {
            type Output = Scalar;
            #[inline]
            fn $method(self, rhs: &'a Scalar) -> Scalar {
                Scalar::$inner(&self, rhs)
            }
        }

        impl<'a> $trait<Scalar> for &'a Scalar {
            type Output = Scalar;
            #[inline]
            fn $method(self, rhs: Scalar) -> Scalar {
                Scalar::$inner(self, &rhs)
            }
        }

        impl<'a, 'b> $trait<&'b Scalar> for &'a Scalar {
            type Output = Scalar;
            #[inline]
            fn $method(self, rhs: &'b Scalar) -> Scalar {
                Scalar::$inner(self, rhs)
            }
        }
    };
}

impl_scalar_binary_op!(Add, add, add_inner);
impl_scalar_binary_op!(Sub, sub, sub_inner);
impl_scalar_binary_op!(Mul, mul, mul_inner);

impl Neg for Scalar {
    type Output = Scalar;
    #[inline]
    fn neg(self) -> Scalar {
        unsafe {
            let mut out = MaybeUninit::<blst_fr>::uninit();
            blst_fr_cneg(out.as_mut_ptr(), &self.0, true);
            Scalar(out.assume_init())
        }
    }
}

impl Neg for &Scalar {
    type Output = Scalar;
    #[inline]
    fn neg(self) -> Scalar {
        -*self
    }
}

macro_rules! impl_scalar_assign_op {
    ($trait:ident, $method:ident, $op:tt) => {
        impl $trait<Scalar> for Scalar {
            #[inline]
            fn $method(&mut self, rhs: Scalar) {
                *self = *self $op rhs;
            }
        }

        impl<'a> $trait<&'a Scalar> for Scalar {
            #[inline]
            fn $method(&mut self, rhs: &'a Scalar) {
                *self = *self $op rhs;
            }
        }
    };
}

impl_scalar_assign_op!(AddAssign, add_assign, +);
impl_scalar_assign_op!(SubAssign, sub_assign, -);
impl_scalar_assign_op!(MulAssign, mul_assign, *);

impl Sum for Scalar {
    fn sum<I: Iterator<Item = Scalar>>(iter: I) -> Self {
        iter.fold(Self::ZERO, |sum, value| sum + value)
    }
}

impl<'a> Sum<&'a Scalar> for Scalar {
    fn sum<I: Iterator<Item = &'a Scalar>>(iter: I) -> Self {
        iter.fold(Self::ZERO, |sum, value| sum + value)
    }
}

/// A G1 point in the prime-order subgroup, stored in scalar-backend projective
/// coordinates.
///
/// The identity is valid. The coordinates are native arithmetic storage: they
/// are not unique and are not a wire encoding. Equality compares mathematical
/// points rather than coordinate limbs. [`Default`] and an empty iterator sum
/// return the identity. Group addition and subtraction are complete for safe
/// values, and multiplication by [`Scalar`] uses the fixed-width
/// constant-time entry point. Operations allocate nothing.
#[repr(transparent)]
#[derive(Copy, Clone, Debug)]
pub struct G1Projective(blst_p1);

/// A G1 point in the prime-order subgroup, stored in scalar-backend affine
/// coordinates.
///
/// The identity is valid. This native arithmetic representation is not a wire
/// encoding. Safe decoding checks both the curve equation and subgroup
/// membership. Equality compares mathematical points. Operations allocate
/// nothing.
#[repr(transparent)]
#[derive(Copy, Clone, Debug)]
pub struct G1Affine(blst_p1_affine);

/// A G2 point in the prime-order subgroup, stored in scalar-backend projective
/// coordinates.
///
/// The identity is valid. The coordinates are native arithmetic storage: they
/// are not unique and are not a wire encoding. Equality compares mathematical
/// points rather than coordinate limbs. [`Default`] and an empty iterator sum
/// return the identity. Group addition and subtraction are complete for safe
/// values, and multiplication by [`Scalar`] uses the fixed-width
/// constant-time entry point. Operations allocate nothing.
#[repr(transparent)]
#[derive(Copy, Clone, Debug)]
pub struct G2Projective(blst_p2);

/// A G2 point in the prime-order subgroup, stored in scalar-backend affine
/// coordinates.
///
/// The identity is valid. This native arithmetic representation is not a wire
/// encoding. Safe decoding checks both the curve equation and subgroup
/// membership. Equality compares mathematical points. Operations allocate
/// nothing.
#[repr(transparent)]
#[derive(Copy, Clone, Debug)]
pub struct G2Affine(blst_p2_affine);

/// An affine batch was given an output slice of the wrong length.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AffineBatchLengthError {
    pub expected: usize,
    pub actual: usize,
}

impl G1Projective {
    /// Converts a secret-derived G1 point without a coordinate-dependent branch.
    /// The identity is supported. This operation allocates nothing.
    #[inline]
    pub fn to_affine_ct(&self) -> G1Affine {
        unsafe {
            let mut out = MaybeUninit::<blst_p1_affine>::uninit();
            blst_p1_to_affine_ct(out.as_mut_ptr(), &self.0);
            G1Affine(out.assume_init())
        }
    }

    /// Batch-converts secret-derived G1 points using one fixed-chain inverse.
    /// Input and output lengths must match. Empty input is valid. This operation
    /// allocates nothing, and execution does not depend on point coordinates.
    #[inline]
    pub fn batch_to_affine_ct(
        input: &[Self],
        output: &mut [G1Affine],
    ) -> Result<(), AffineBatchLengthError> {
        if input.len() != output.len() {
            return Err(AffineBatchLengthError {
                expected: input.len(),
                actual: output.len(),
            });
        }
        if input.is_empty() {
            return Ok(());
        }
        let descriptors = [input.as_ptr().cast::<blst_p1>(), core::ptr::null()];
        unsafe {
            blst_p1s_to_affine_ct(
                output.as_mut_ptr().cast::<blst_p1_affine>(),
                descriptors.as_ptr(),
                input.len(),
            );
        }
        Ok(())
    }
}

macro_rules! impl_scalar_curve_group {
    (
        $projective:ident,
        $affine:ident,
        $raw_projective:ty,
        $raw_affine:ty,
        $compressed_size:expr,
        $uncompressed_size:expr,
        $generator:ident,
        $affine_generator:ident,
        $is_equal:ident,
        $affine_is_equal:ident,
        $is_inf:ident,
        $affine_is_inf:ident,
        $on_curve:ident,
        $affine_on_curve:ident,
        $in_group:ident,
        $affine_in_group:ident,
        $add_or_double:ident,
        $sub:ident,
        $double:ident,
        $cneg:ident,
        $mult:ident,
        $mult_vartime:ident,
        $to_affine:ident,
        $from_affine:ident,
        $compress:ident,
        $affine_compress:ident,
        $serialize:ident,
        $affine_serialize:ident,
        $uncompress:ident,
        $deserialize:ident
    ) => {
        impl $projective {
            /// Returns the group identity in native projective form.
            ///
            /// This operation allocates nothing.
            #[inline]
            pub fn identity() -> Self {
                Self(<$raw_projective>::default())
            }

            /// Returns the category-blst standard generator.
            ///
            /// The returned point is nonidentity and generates the prime-order
            /// subgroup. This operation allocates nothing.
            #[inline]
            pub fn generator() -> Self {
                unsafe { Self(*$generator()) }
            }

            /// Borrows the native scalar-backend projective representation.
            ///
            /// This is a zero-conversion borrow. The coordinates are not a
            /// canonical encoding and must not be serialized directly.
            #[inline]
            pub const fn as_raw(&self) -> &$raw_projective {
                &self.0
            }

            /// Mutably borrows the native projective representation.
            ///
            /// This accessor is intended for a native kernel that writes a
            /// projective subgroup point into caller-owned storage. It performs
            /// no validation, conversion, or allocation.
            ///
            /// # Safety
            ///
            /// Before the borrow ends, the caller must leave a valid native
            /// projective representation of an on-curve point in the
            /// prime-order subgroup. The identity is permitted.
            #[inline]
            pub unsafe fn as_raw_mut(&mut self) -> &mut $raw_projective {
                &mut self.0
            }

            /// Wraps a native projective subgroup point without conversion.
            ///
            /// This constructor performs no validation or allocation.
            ///
            /// # Safety
            ///
            /// `raw` must be a valid native projective representation of an
            /// on-curve point in the prime-order subgroup. The identity is
            /// permitted.
            #[inline]
            pub const unsafe fn from_raw_unchecked(
                raw: $raw_projective,
            ) -> Self {
                Self(raw)
            }

            /// Consumes the wrapper and returns its native projective storage.
            ///
            /// The returned coordinates retain the subgroup invariant but are
            /// not a canonical encoding. This operation allocates nothing.
            #[inline]
            pub const fn into_raw(self) -> $raw_projective {
                self.0
            }

            /// Returns `true` exactly when this point is the group identity.
            ///
            /// This observation performs no validation and allocates nothing.
            #[inline]
            pub fn is_identity(&self) -> bool {
                unsafe { $is_inf(&self.0) }
            }

            /// Returns `true` when the native representation lies on the curve.
            ///
            /// Safe values already satisfy this property. This method is useful
            /// as a defensive assertion after unsafe native interoperation; it
            /// does not establish the complete wrapper invariant by itself.
            #[inline]
            pub fn is_on_curve(&self) -> bool {
                unsafe { $on_curve(&self.0) }
            }

            /// Returns `true` when the native point lies in the prime-order subgroup.
            ///
            /// Safe values already satisfy this property. This observation does
            /// not alter or re-wrap the point and allocates nothing.
            #[inline]
            pub fn is_in_subgroup(&self) -> bool {
                unsafe { $in_group(&self.0) }
            }

            /// `blstrs`-inspired spelling for [`Self::is_in_subgroup`].
            ///
            /// This is an exact alias, including its preconditions and timing
            /// behavior.
            #[inline]
            pub fn is_torsion_free(&self) -> bool {
                self.is_in_subgroup()
            }

            /// Returns `self + self`.
            ///
            /// Safe inputs require no curve or subgroup validation. The result
            /// remains in the prime-order subgroup. This operation allocates
            /// nothing.
            #[inline]
            pub fn double(&self) -> Self {
                unsafe {
                    let mut out = MaybeUninit::<$raw_projective>::uninit();
                    $double(out.as_mut_ptr(), &self.0);
                    Self(out.assume_init())
                }
            }

            /// Multiplies by a possibly secret scalar using the fixed-width
            /// constant-time scalar multiplication entry point.
            ///
            /// This computes `[scalar] self`. The scalar is encoded canonically
            /// into a 256-bit stack buffer before the native call. The result
            /// remains in the prime-order subgroup, including when either input
            /// is zero or the identity. The operation allocates nothing.
            #[inline]
            pub fn multiply(&self, scalar: &Scalar) -> Self {
                let encoded = scalar.to_bytes_le();
                unsafe {
                    let mut out = MaybeUninit::<$raw_projective>::uninit();
                    $mult(out.as_mut_ptr(), &self.0, encoded.as_ptr(), 256);
                    Self(out.assume_init())
                }
            }

            /// Multiplies by a public scalar using the variable-time entry point.
            ///
            /// This computes `[scalar] self` using a fixed 256-bit canonical
            /// scalar encoding. `scalar` must be public because execution time
            /// and memory access may depend on it. The result remains in the
            /// prime-order subgroup. The operation allocates nothing.
            #[inline]
            pub fn multiply_vartime(&self, scalar: &Scalar) -> Self {
                let encoded = scalar.to_bytes_le();
                unsafe {
                    let mut out = MaybeUninit::<$raw_projective>::uninit();
                    $mult_vartime(
                        out.as_mut_ptr(),
                        &self.0,
                        encoded.as_ptr(),
                        256,
                    );
                    Self(out.assume_init())
                }
            }

            /// Converts this point to the equivalent affine representation.
            ///
            /// The identity remains the identity. This operation performs no
            /// validation and allocates nothing. No constant-time guarantee is
            /// made for this conversion.
            #[inline]
            pub fn to_affine(&self) -> $affine {
                unsafe {
                    let mut out = MaybeUninit::<$raw_affine>::uninit();
                    $to_affine(out.as_mut_ptr(), &self.0);
                    $affine(out.assume_init())
                }
            }

            /// Encodes this point in canonical compressed form.
            ///
            /// The fixed-size output is the unique category-blst compressed
            /// encoding of the mathematical point. The identity is encodable.
            /// This operation allocates nothing.
            #[inline]
            pub fn to_compressed(&self) -> [u8; $compressed_size] {
                unsafe {
                    let mut out =
                        MaybeUninit::<[u8; $compressed_size]>::uninit();
                    $compress(out.as_mut_ptr().cast::<u8>(), &self.0);
                    out.assume_init()
                }
            }

            /// Encodes this point in canonical uncompressed form.
            ///
            /// The fixed-size output is the unique category-blst uncompressed
            /// encoding of the mathematical point. The identity is encodable.
            /// This operation allocates nothing.
            #[inline]
            pub fn to_uncompressed(&self) -> [u8; $uncompressed_size] {
                unsafe {
                    let mut out =
                        MaybeUninit::<[u8; $uncompressed_size]>::uninit();
                    $serialize(out.as_mut_ptr().cast::<u8>(), &self.0);
                    out.assume_init()
                }
            }

            /// Strictly decodes and subgroup-validates a compressed point.
            ///
            /// This has the same accepted encodings and errors as the
            /// corresponding affine decoder, followed by an allocation-free
            /// affine-to-projective conversion. The identity is accepted.
            #[inline]
            pub fn from_compressed(
                encoded: &[u8; $compressed_size],
            ) -> Result<Self, BLST_ERROR> {
                $affine::from_compressed(encoded).map(Self::from)
            }

            /// Strictly decodes and subgroup-validates an uncompressed point.
            ///
            /// This has the same accepted encodings and errors as the
            /// corresponding affine decoder, followed by an allocation-free
            /// affine-to-projective conversion. The identity is accepted.
            #[inline]
            pub fn from_uncompressed(
                encoded: &[u8; $uncompressed_size],
            ) -> Result<Self, BLST_ERROR> {
                $affine::from_uncompressed(encoded).map(Self::from)
            }

            #[inline]
            fn add_inner(lhs: &Self, rhs: &Self) -> Self {
                unsafe {
                    let mut out = MaybeUninit::<$raw_projective>::uninit();
                    $add_or_double(out.as_mut_ptr(), &lhs.0, &rhs.0);
                    Self(out.assume_init())
                }
            }

            #[inline]
            fn sub_inner(lhs: &Self, rhs: &Self) -> Self {
                unsafe {
                    let mut out = MaybeUninit::<$raw_projective>::uninit();
                    $sub(out.as_mut_ptr(), &lhs.0, &rhs.0);
                    Self(out.assume_init())
                }
            }
        }

        impl $affine {
            /// Returns the group identity in native affine form.
            ///
            /// This operation allocates nothing.
            #[inline]
            pub fn identity() -> Self {
                Self(<$raw_affine>::default())
            }

            /// Returns the category-blst standard generator in affine form.
            ///
            /// The returned point is nonidentity and generates the prime-order
            /// subgroup. This operation allocates nothing.
            #[inline]
            pub fn generator() -> Self {
                unsafe { Self(*$affine_generator()) }
            }

            /// Borrows the native scalar-backend affine representation.
            ///
            /// This is a zero-conversion borrow. The native coordinates are not
            /// a canonical encoding and must not be serialized directly.
            #[inline]
            pub const fn as_raw(&self) -> &$raw_affine {
                &self.0
            }

            /// Mutably borrows the native affine representation.
            ///
            /// This accessor is intended for a native kernel that writes an
            /// affine subgroup point into caller-owned storage. It performs no
            /// validation, conversion, or allocation.
            ///
            /// # Safety
            ///
            /// Before the borrow ends, the caller must leave a valid native
            /// affine representation of an on-curve point in the prime-order
            /// subgroup. The identity is permitted.
            #[inline]
            pub unsafe fn as_raw_mut(&mut self) -> &mut $raw_affine {
                &mut self.0
            }

            /// Wraps a native affine subgroup point without conversion.
            ///
            /// This constructor performs no validation or allocation.
            ///
            /// # Safety
            ///
            /// `raw` must be a valid native affine representation of an
            /// on-curve point in the prime-order subgroup. The identity is
            /// permitted.
            #[inline]
            pub const unsafe fn from_raw_unchecked(raw: $raw_affine) -> Self {
                Self(raw)
            }

            /// Consumes the wrapper and returns its native affine storage.
            ///
            /// The returned coordinates retain the subgroup invariant but are
            /// not a canonical encoding. This operation allocates nothing.
            #[inline]
            pub const fn into_raw(self) -> $raw_affine {
                self.0
            }

            /// Returns `true` exactly when this point is the group identity.
            ///
            /// This observation performs no validation and allocates nothing.
            #[inline]
            pub fn is_identity(&self) -> bool {
                unsafe { $affine_is_inf(&self.0) }
            }

            /// Returns `true` when the native representation lies on the curve.
            ///
            /// Safe values already satisfy this property. This method is useful
            /// as a defensive assertion after unsafe native interoperation; it
            /// does not establish the complete wrapper invariant by itself.
            #[inline]
            pub fn is_on_curve(&self) -> bool {
                unsafe { $affine_on_curve(&self.0) }
            }

            /// Returns `true` when the native point lies in the prime-order subgroup.
            ///
            /// Safe values already satisfy this property. This observation does
            /// not alter or re-wrap the point and allocates nothing.
            #[inline]
            pub fn is_in_subgroup(&self) -> bool {
                unsafe { $affine_in_group(&self.0) }
            }

            /// `blstrs`-inspired spelling for [`Self::is_in_subgroup`].
            ///
            /// This is an exact alias, including its preconditions and timing
            /// behavior.
            #[inline]
            pub fn is_torsion_free(&self) -> bool {
                self.is_in_subgroup()
            }

            /// Encodes this point in canonical compressed form.
            ///
            /// The fixed-size output is the unique category-blst compressed
            /// encoding of the mathematical point. The identity is encodable.
            /// This operation allocates nothing.
            #[inline]
            pub fn to_compressed(&self) -> [u8; $compressed_size] {
                unsafe {
                    let mut out =
                        MaybeUninit::<[u8; $compressed_size]>::uninit();
                    $affine_compress(out.as_mut_ptr().cast::<u8>(), &self.0);
                    out.assume_init()
                }
            }

            /// Encodes this point in canonical uncompressed form.
            ///
            /// The fixed-size output is the unique category-blst uncompressed
            /// encoding of the mathematical point. The identity is encodable.
            /// This operation allocates nothing.
            #[inline]
            pub fn to_uncompressed(&self) -> [u8; $uncompressed_size] {
                unsafe {
                    let mut out =
                        MaybeUninit::<[u8; $uncompressed_size]>::uninit();
                    $affine_serialize(out.as_mut_ptr().cast::<u8>(), &self.0);
                    out.assume_init()
                }
            }

            /// Strictly decodes a canonical compressed point and establishes
            /// prime-order subgroup membership.
            ///
            /// The native decoder checks compressed-format flags, canonical
            /// field elements, and the curve equation. Its exact [`BLST_ERROR`]
            /// is returned on failure. A decoded point outside the prime-order
            /// subgroup returns [`BLST_ERROR::BLST_POINT_NOT_IN_GROUP`]. The
            /// canonical identity encoding is accepted. Failure exposes no
            /// partially decoded point, and the operation allocates nothing.
            #[inline]
            pub fn from_compressed(
                encoded: &[u8; $compressed_size],
            ) -> Result<Self, BLST_ERROR> {
                unsafe {
                    let mut out = MaybeUninit::<$raw_affine>::uninit();
                    let error = $uncompress(out.as_mut_ptr(), encoded.as_ptr());
                    if error != BLST_ERROR::BLST_SUCCESS {
                        return Err(error);
                    }
                    let out = out.assume_init();
                    if !$affine_in_group(&out) {
                        return Err(BLST_ERROR::BLST_POINT_NOT_IN_GROUP);
                    }
                    Ok(Self(out))
                }
            }

            /// Strictly decodes a canonical uncompressed point and establishes
            /// prime-order subgroup membership.
            ///
            /// The native decoder checks uncompressed-format flags, canonical
            /// field elements, and the curve equation. Its exact [`BLST_ERROR`]
            /// is returned on failure. A decoded point outside the prime-order
            /// subgroup returns [`BLST_ERROR::BLST_POINT_NOT_IN_GROUP`]. The
            /// canonical identity encoding is accepted. Failure exposes no
            /// partially decoded point, and the operation allocates nothing.
            #[inline]
            pub fn from_uncompressed(
                encoded: &[u8; $uncompressed_size],
            ) -> Result<Self, BLST_ERROR> {
                unsafe {
                    let mut out = MaybeUninit::<$raw_affine>::uninit();
                    let error =
                        $deserialize(out.as_mut_ptr(), encoded.as_ptr());
                    if error != BLST_ERROR::BLST_SUCCESS {
                        return Err(error);
                    }
                    let out = out.assume_init();
                    if !$affine_in_group(&out) {
                        return Err(BLST_ERROR::BLST_POINT_NOT_IN_GROUP);
                    }
                    Ok(Self(out))
                }
            }
        }

        impl Default for $projective {
            #[inline]
            fn default() -> Self {
                Self::identity()
            }
        }

        impl Default for $affine {
            #[inline]
            fn default() -> Self {
                Self::identity()
            }
        }

        impl PartialEq for $projective {
            #[inline]
            fn eq(&self, other: &Self) -> bool {
                unsafe { $is_equal(&self.0, &other.0) }
            }
        }

        impl Eq for $projective {}

        impl PartialEq for $affine {
            #[inline]
            fn eq(&self, other: &Self) -> bool {
                unsafe { $affine_is_equal(&self.0, &other.0) }
            }
        }

        impl Eq for $affine {}

        impl From<$projective> for $affine {
            #[inline]
            fn from(point: $projective) -> Self {
                point.to_affine()
            }
        }

        impl From<&$projective> for $affine {
            #[inline]
            fn from(point: &$projective) -> Self {
                point.to_affine()
            }
        }

        impl From<$affine> for $projective {
            #[inline]
            fn from(point: $affine) -> Self {
                unsafe {
                    let mut out = MaybeUninit::<$raw_projective>::uninit();
                    $from_affine(out.as_mut_ptr(), &point.0);
                    Self(out.assume_init())
                }
            }
        }

        impl From<&$affine> for $projective {
            #[inline]
            fn from(point: &$affine) -> Self {
                Self::from(*point)
            }
        }

        impl Add for $projective {
            type Output = Self;
            #[inline]
            fn add(self, rhs: Self) -> Self {
                Self::add_inner(&self, &rhs)
            }
        }

        impl Add<&$projective> for $projective {
            type Output = Self;
            #[inline]
            fn add(self, rhs: &$projective) -> Self {
                Self::add_inner(&self, rhs)
            }
        }

        impl Add<$projective> for &$projective {
            type Output = $projective;
            #[inline]
            fn add(self, rhs: $projective) -> $projective {
                $projective::add_inner(self, &rhs)
            }
        }

        impl Add<&$projective> for &$projective {
            type Output = $projective;
            #[inline]
            fn add(self, rhs: &$projective) -> $projective {
                $projective::add_inner(self, rhs)
            }
        }

        impl Sub for $projective {
            type Output = Self;
            #[inline]
            fn sub(self, rhs: Self) -> Self {
                Self::sub_inner(&self, &rhs)
            }
        }

        impl Sub<&$projective> for $projective {
            type Output = Self;
            #[inline]
            fn sub(self, rhs: &$projective) -> Self {
                Self::sub_inner(&self, rhs)
            }
        }

        impl Sub<$projective> for &$projective {
            type Output = $projective;
            #[inline]
            fn sub(self, rhs: $projective) -> $projective {
                $projective::sub_inner(self, &rhs)
            }
        }

        impl Sub<&$projective> for &$projective {
            type Output = $projective;
            #[inline]
            fn sub(self, rhs: &$projective) -> $projective {
                $projective::sub_inner(self, rhs)
            }
        }

        impl Neg for $projective {
            type Output = Self;
            #[inline]
            fn neg(mut self) -> Self {
                unsafe { $cneg(&mut self.0, true) };
                self
            }
        }

        impl Neg for &$projective {
            type Output = $projective;
            #[inline]
            fn neg(self) -> $projective {
                -*self
            }
        }

        impl AddAssign<$projective> for $projective {
            #[inline]
            fn add_assign(&mut self, rhs: $projective) {
                *self = *self + rhs;
            }
        }

        impl AddAssign<&$projective> for $projective {
            #[inline]
            fn add_assign(&mut self, rhs: &$projective) {
                *self = *self + rhs;
            }
        }

        impl SubAssign<$projective> for $projective {
            #[inline]
            fn sub_assign(&mut self, rhs: $projective) {
                *self = *self - rhs;
            }
        }

        impl SubAssign<&$projective> for $projective {
            #[inline]
            fn sub_assign(&mut self, rhs: &$projective) {
                *self = *self - rhs;
            }
        }

        impl Mul<Scalar> for $projective {
            type Output = Self;
            #[inline]
            fn mul(self, rhs: Scalar) -> Self {
                self.multiply(&rhs)
            }
        }

        impl Mul<&Scalar> for $projective {
            type Output = Self;
            #[inline]
            fn mul(self, rhs: &Scalar) -> Self {
                self.multiply(rhs)
            }
        }

        impl Mul<Scalar> for &$projective {
            type Output = $projective;
            #[inline]
            fn mul(self, rhs: Scalar) -> $projective {
                self.multiply(&rhs)
            }
        }

        impl Mul<&Scalar> for &$projective {
            type Output = $projective;
            #[inline]
            fn mul(self, rhs: &Scalar) -> $projective {
                self.multiply(rhs)
            }
        }

        impl Sum for $projective {
            fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
                iter.fold(Self::identity(), |sum, value| sum + value)
            }
        }

        impl<'a> Sum<&'a $projective> for $projective {
            fn sum<I: Iterator<Item = &'a $projective>>(iter: I) -> Self {
                iter.fold(Self::identity(), |sum, value| sum + value)
            }
        }
    };
}

impl_scalar_curve_group!(
    G1Projective,
    G1Affine,
    blst_p1,
    blst_p1_affine,
    48,
    96,
    blst_p1_generator,
    blst_p1_affine_generator,
    blst_p1_is_equal,
    blst_p1_affine_is_equal,
    blst_p1_is_inf,
    blst_p1_affine_is_inf,
    blst_p1_on_curve,
    blst_p1_affine_on_curve,
    blst_p1_in_g1,
    blst_p1_affine_in_g1,
    blst_p1_add_or_double,
    blst_p1_sub,
    blst_p1_double,
    blst_p1_cneg,
    blst_p1_mult,
    blst_p1_mult_vartime,
    blst_p1_to_affine,
    blst_p1_from_affine,
    blst_p1_compress,
    blst_p1_affine_compress,
    blst_p1_serialize,
    blst_p1_affine_serialize,
    blst_p1_uncompress,
    blst_p1_deserialize
);

impl_scalar_curve_group!(
    G2Projective,
    G2Affine,
    blst_p2,
    blst_p2_affine,
    96,
    192,
    blst_p2_generator,
    blst_p2_affine_generator,
    blst_p2_is_equal,
    blst_p2_affine_is_equal,
    blst_p2_is_inf,
    blst_p2_affine_is_inf,
    blst_p2_on_curve,
    blst_p2_affine_on_curve,
    blst_p2_in_g2,
    blst_p2_affine_in_g2,
    blst_p2_add_or_double,
    blst_p2_sub,
    blst_p2_double,
    blst_p2_cneg,
    blst_p2_mult,
    blst_p2_mult_vartime,
    blst_p2_to_affine,
    blst_p2_from_affine,
    blst_p2_compress,
    blst_p2_affine_compress,
    blst_p2_serialize,
    blst_p2_affine_serialize,
    blst_p2_uncompress,
    blst_p2_deserialize
);

/// A product of Miller loops that has not yet undergone final exponentiation.
///
/// This is deliberately distinct from [`Gt`]. The shared `blst_fp12` storage
/// does not by itself establish target-group membership. Safe values are
/// nonzero products of Miller-loop contributions evaluated on G1 and G2
/// prime-order subgroup points. Multiplication combines products, [`Default`]
/// is the empty product, and all operations allocate nothing.
#[repr(transparent)]
#[derive(Copy, Clone, Debug)]
pub struct MillerLoopResult(blst_fp12);

impl MillerLoopResult {
    /// Returns the empty Miller-loop product, the multiplicative Fp12 identity.
    ///
    /// This nonzero value is a valid input to final exponentiation. This
    /// operation allocates nothing.
    #[inline]
    pub fn identity() -> Self {
        unsafe { Self(*blst_fp12_one()) }
    }

    /// Computes one Miller-loop contribution from validated subgroup points.
    ///
    /// This computes the unfinalized contribution for `e(g1, g2)`. The type
    /// invariants of `g1` and `g2` satisfy the native routine's curve and
    /// subgroup preconditions; this method performs no additional validation.
    /// If either point is the identity, the contribution is the empty product.
    /// The operation does not perform final exponentiation and allocates
    /// nothing.
    #[inline]
    pub fn from_pair(g1: &G1Affine, g2: &G2Affine) -> Self {
        unsafe {
            let mut out = MaybeUninit::<blst_fp12>::uninit();
            blst_miller_loop(out.as_mut_ptr(), g2.as_raw(), g1.as_raw());
            Self(out.assume_init())
        }
    }

    /// Computes the product of Miller loops over corresponding affine slices.
    ///
    /// Element `i` of `g1` is paired with element `i` of `g2`. Both slices must
    /// have the same length; the output is the product of their unfinalized
    /// pairing contributions. A mismatch returns
    /// [`PairingError::LengthMismatch`] before any native call. Empty slices
    /// return the empty product.
    ///
    /// The point wrappers satisfy the native routine's curve and subgroup
    /// preconditions, so this method performs no per-point validation. One
    /// batched native call processes every pair. The two fixed-size descriptor
    /// arrays live on the stack; the operation does not allocate or create
    /// pointer arrays proportional to the batch size. It does not perform final
    /// exponentiation.
    #[inline]
    pub fn from_pairs(
        g1: &[G1Affine],
        g2: &[G2Affine],
    ) -> Result<Self, PairingError> {
        if g1.len() != g2.len() {
            return Err(PairingError::LengthMismatch {
                g1: g1.len(),
                g2: g2.len(),
            });
        }
        if g1.is_empty() {
            return Ok(Self::identity());
        }

        // Transparent wrappers preserve the native element stride, so the C
        // batch API can continue contiguously after each first descriptor.
        let g1_descriptors =
            [g1.as_ptr().cast::<blst_p1_affine>(), core::ptr::null()];
        let g2_descriptors =
            [g2.as_ptr().cast::<blst_p2_affine>(), core::ptr::null()];
        unsafe {
            let mut out = MaybeUninit::<blst_fp12>::uninit();
            blst_miller_loop_n(
                out.as_mut_ptr(),
                g2_descriptors.as_ptr(),
                g1_descriptors.as_ptr(),
                g1.len(),
            );
            Ok(Self(out.assume_init()))
        }
    }

    /// Applies final exponentiation and returns the corresponding [`Gt`] value.
    ///
    /// A safe [`MillerLoopResult`] is nonzero and therefore satisfies the native
    /// final-exponentiation precondition. The operation performs no validation
    /// and allocates nothing. No constant-time guarantee is made.
    #[inline]
    pub fn final_exponentiation(&self) -> Gt {
        unsafe {
            let mut out = MaybeUninit::<blst_fp12>::uninit();
            blst_final_exp(out.as_mut_ptr(), &self.0);
            Gt(out.assume_init())
        }
    }

    /// Borrows the native Fp12 Miller-loop storage without conversion.
    ///
    /// The returned Fp12 value is not necessarily in [`Gt`] before final
    /// exponentiation. It must not be reinterpreted as a target-group value.
    #[inline]
    pub const fn as_raw(&self) -> &blst_fp12 {
        &self.0
    }

    /// Wraps a trusted native Miller-loop result without conversion.
    ///
    /// This constructor performs no nonzero, provenance, or subgroup check.
    ///
    /// # Safety
    ///
    /// `raw` must be nonzero and must equal a product of Miller-loop
    /// contributions computed from G1 and G2 prime-order subgroup points. A
    /// general Fp12 value does not satisfy this contract merely because it is
    /// nonzero.
    #[inline]
    pub const unsafe fn from_raw_unchecked(raw: blst_fp12) -> Self {
        Self(raw)
    }

    #[inline]
    fn multiply_inner(lhs: &Self, rhs: &Self) -> Self {
        unsafe {
            let mut out = MaybeUninit::<blst_fp12>::uninit();
            blst_fp12_mul(out.as_mut_ptr(), &lhs.0, &rhs.0);
            Self(out.assume_init())
        }
    }
}

/// Construction or shape errors for scalar prepared-G2 batch evaluation.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum PreparedG2BatchError {
    /// A prepared batch must own at least one operand.
    EmptyBatch,
    /// Native immutable line-table allocation failed.
    AllocationFailed,
    /// Per-pair G1 and prepared-index counts differed.
    InputLengthMismatch {
        /// Number of G1 inputs.
        g1: usize,
        /// Number of prepared indices.
        indices: usize,
    },
    /// A selected prepared index exceeded the immutable batch.
    IndexOutOfRange {
        /// Rejected index.
        index: usize,
        /// Prepared batch length.
        len: usize,
    },
    /// The native boundary rejected an input after Rust validation.
    BackendFailure,
}

/// Immutable scalar Miller-line tables for validated G2 subgroup points.
///
/// Construction allocates and prepares one table per operand. The object owns
/// no mutable evaluation state, so independent workers may share it. Evaluation
/// accepts native scalar G1 projective points and performs no representation
/// conversion through another backend. Identity operands contribute the
/// Miller identity.
pub struct PreparedG2Batch {
    inner: NonNull<blst_scalar_prepared_g2_batch>,
    len: usize,
    storage_bytes: usize,
}

impl PreparedG2Batch {
    /// Allocates immutable prepared lines for every supplied G2 point.
    ///
    /// Each point must already satisfy [`G2Projective`]'s prime-order subgroup
    /// invariant. Construction is variable-time and rejects an empty slice.
    pub fn new(g2: &[G2Projective]) -> Result<Self, PreparedG2BatchError> {
        if g2.is_empty() {
            return Err(PreparedG2BatchError::EmptyBatch);
        }
        let inner = NonNull::new(unsafe {
            blst_scalar_prepared_g2_batch_new(
                g2.as_ptr().cast::<blst_p2>(),
                g2.len(),
            )
        })
        .ok_or(PreparedG2BatchError::AllocationFailed)?;
        let len = unsafe { blst_scalar_prepared_g2_batch_len(inner.as_ptr()) };
        if len != g2.len() {
            unsafe { blst_scalar_prepared_g2_batch_free(inner.as_ptr()) };
            return Err(PreparedG2BatchError::BackendFailure);
        }
        let storage_bytes = unsafe {
            blst_scalar_prepared_g2_batch_storage_bytes(inner.as_ptr())
        };
        if storage_bytes == 0 {
            unsafe { blst_scalar_prepared_g2_batch_free(inner.as_ptr()) };
            return Err(PreparedG2BatchError::BackendFailure);
        }
        Ok(Self {
            inner,
            len,
            storage_bytes,
        })
    }

    /// Returns the number of prepared operands.
    #[inline]
    pub const fn len(&self) -> usize {
        self.len
    }

    /// Returns whether the prepared sequence is empty.
    ///
    /// Construction rejects empty input, so this is always false for a live
    /// value and exists only for ordinary collection ergonomics.
    #[inline]
    pub const fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Returns native object and line-table bytes owned by this batch.
    #[inline]
    pub const fn storage_bytes(&self) -> usize {
        self.storage_bytes
    }

    /// Multiplies selected prepared Miller-loop contributions.
    ///
    /// Element `i` of `g1` is paired with prepared operand `indices[i]`. The
    /// slices must have equal length and every index must be in range. The
    /// result remains outside the target group until
    /// [`MillerLoopResult::final_exponentiation`]. The call allocates nothing,
    /// retains no state, and may run concurrently against the same immutable
    /// batch. Timing may depend on public slice length, indices, and point
    /// identities.
    pub fn miller_loop_product_indices(
        &self,
        g1: &[G1Projective],
        indices: &[usize],
    ) -> Result<MillerLoopResult, PreparedG2BatchError> {
        if g1.len() != indices.len() {
            return Err(PreparedG2BatchError::InputLengthMismatch {
                g1: g1.len(),
                indices: indices.len(),
            });
        }
        if let Some(&index) = indices.iter().find(|&&index| index >= self.len) {
            return Err(PreparedG2BatchError::IndexOutOfRange {
                index,
                len: self.len,
            });
        }
        let mut out = MaybeUninit::<blst_fp12>::uninit();
        if unsafe {
            blst_scalar_prepared_g2_batch_miller_loop_product_indices(
                out.as_mut_ptr(),
                g1.as_ptr().cast::<blst_p1>(),
                indices.as_ptr(),
                g1.len(),
                self.inner.as_ptr(),
            )
        } {
            Ok(unsafe {
                MillerLoopResult::from_raw_unchecked(out.assume_init())
            })
        } else {
            Err(PreparedG2BatchError::BackendFailure)
        }
    }

    /// Evaluates a contiguous range against one shared G1 point.
    ///
    /// Element `i` of `out` receives the pairing with prepared operand
    /// `begin + i`. The range must fit inside the immutable batch. Independent
    /// calls may share the batch when their output slices do not overlap. The
    /// operation allocates no memory and converts the shared G1 point to affine
    /// once per range.
    pub fn pairing_range_into(
        &self,
        g1: &G1Projective,
        begin: usize,
        out: &mut [Gt],
    ) -> Result<(), PreparedG2BatchError> {
        // SAFETY: `MaybeUninit<Gt>` has the same layout as `Gt`. The delegated
        // call initializes every element on success. On failure, `out` retains
        // the same unspecified-output contract as the native operation.
        let uninit = unsafe {
            core::slice::from_raw_parts_mut(
                out.as_mut_ptr().cast::<MaybeUninit<Gt>>(),
                out.len(),
            )
        };
        self.pairing_range_uninit_into(g1, begin, uninit)
    }

    /// Evaluates a contiguous range into uninitialized output storage.
    ///
    /// On success every output element is initialized and satisfies the `Gt`
    /// invariant. On failure the caller must continue to treat every element
    /// as uninitialized. This avoids initializing a task-owned result buffer
    /// that the native range will overwrite completely.
    pub fn pairing_range_uninit_into(
        &self,
        g1: &G1Projective,
        begin: usize,
        out: &mut [MaybeUninit<Gt>],
    ) -> Result<(), PreparedG2BatchError> {
        if begin > self.len {
            return Err(PreparedG2BatchError::IndexOutOfRange {
                index: begin,
                len: self.len,
            });
        }
        let end = begin.checked_add(out.len()).ok_or(
            PreparedG2BatchError::IndexOutOfRange {
                index: usize::MAX,
                len: self.len,
            },
        )?;
        if end > self.len {
            return Err(PreparedG2BatchError::IndexOutOfRange {
                index: end - 1,
                len: self.len,
            });
        }
        if unsafe {
            blst_scalar_prepared_g2_batch_pairing_range(
                out.as_mut_ptr().cast::<blst_fp12>(),
                begin,
                out.len(),
                g1.as_raw(),
                self.inner.as_ptr(),
            )
        } {
            Ok(())
        } else {
            Err(PreparedG2BatchError::BackendFailure)
        }
    }
}

impl core::fmt::Debug for PreparedG2Batch {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("PreparedG2Batch")
            .field("len", &self.len)
            .field("storage_bytes", &self.storage_bytes)
            .finish_non_exhaustive()
    }
}

impl Drop for PreparedG2Batch {
    fn drop(&mut self) {
        unsafe { blst_scalar_prepared_g2_batch_free(self.inner.as_ptr()) };
    }
}

// The native object owns immutable line tables and exposes only const
// evaluation. Destruction still requires exclusive Rust ownership.
unsafe impl Send for PreparedG2Batch {}
unsafe impl Sync for PreparedG2Batch {}

impl Default for MillerLoopResult {
    #[inline]
    fn default() -> Self {
        Self::identity()
    }
}

impl PartialEq for MillerLoopResult {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        unsafe { blst_fp12_is_equal(&self.0, &other.0) }
    }
}

impl Eq for MillerLoopResult {}

impl Mul for MillerLoopResult {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: Self) -> Self {
        Self::multiply_inner(&self, &rhs)
    }
}

impl MulAssign for MillerLoopResult {
    #[inline]
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

/// A validated element of the order-`r` BLS12-381 pairing target group.
///
/// Group operations use the additive conventions familiar from `blstrs` even
/// though the underlying Fp12 group law is multiplication. Thus `a + b`
/// multiplies the Fp12 values, `-a` computes the group inverse, and `a * s`
/// exponentiates by scalar `s`. [`Default`] and an empty iterator sum return the
/// target-group identity. Equality compares Fp12 values representing group
/// elements. Operations allocate nothing.
#[repr(transparent)]
#[derive(Copy, Clone, Debug)]
pub struct Gt(blst_fp12);

impl Gt {
    /// Returns the target-group identity.
    ///
    /// Under the underlying Fp12 representation this is the multiplicative
    /// identity. This operation allocates nothing.
    #[inline]
    pub fn identity() -> Self {
        unsafe { Self(*blst_gt_identity()) }
    }

    /// Returns the pairing of the standard G1 and G2 generators.
    ///
    /// The returned value is nonidentity and generates the order-`r` target
    /// group. This operation allocates nothing.
    #[inline]
    pub fn generator() -> Self {
        unsafe { Self(*blst_gt_generator()) }
    }

    /// Computes a pairing from validated subgroup points.
    ///
    /// This is final exponentiation applied to one Miller-loop contribution.
    /// The input type invariants satisfy the native curve and subgroup
    /// preconditions; no additional validation is performed. If either input is
    /// the identity, the result is [`Gt::identity`]. The operation allocates
    /// nothing. No constant-time guarantee is made.
    #[inline]
    pub fn pairing(g1: &G1Affine, g2: &G2Affine) -> Self {
        MillerLoopResult::from_pair(g1, g2).final_exponentiation()
    }

    /// Borrows the native target-group representation without conversion.
    ///
    /// The returned Fp12 value retains the target-group invariant. It is native
    /// arithmetic storage, not the canonical 576-byte encoding.
    #[inline]
    pub const fn as_raw(&self) -> &blst_fp12 {
        &self.0
    }

    /// Mutably borrows the native target-group representation without conversion.
    ///
    /// This accessor is intended for a native kernel that writes a target-group
    /// value into caller-owned storage. It performs no validation, conversion,
    /// or allocation.
    ///
    /// # Safety
    ///
    /// Before the borrow ends, the caller must leave a valid Fp12 element in the
    /// order-`r` target subgroup. The identity is permitted.
    #[inline]
    pub unsafe fn as_raw_mut(&mut self) -> &mut blst_fp12 {
        &mut self.0
    }

    /// Wraps a native target-group value without conversion.
    ///
    /// This constructor performs no target-group check or allocation.
    ///
    /// # Safety
    ///
    /// `raw` must lie in the order-`r` BLS12-381 target subgroup. The identity
    /// is permitted. A general Fp12 or unfinalized Miller-loop result does not
    /// satisfy this precondition.
    #[inline]
    pub const unsafe fn from_raw_unchecked(raw: blst_fp12) -> Self {
        Self(raw)
    }

    /// Consumes the wrapper and returns its native Fp12 storage.
    ///
    /// The returned value retains the target-group invariant but is not the
    /// canonical wire encoding. This operation allocates nothing.
    #[inline]
    pub const fn into_raw(self) -> blst_fp12 {
        self.0
    }

    /// Returns `true` exactly when this is the target-group identity.
    ///
    /// This observation performs no validation and allocates nothing.
    #[inline]
    pub fn is_identity(&self) -> bool {
        unsafe { blst_gt_is_identity(&self.0) }
    }

    /// Returns `true` when the underlying Fp12 value lies in the target group.
    ///
    /// Safe values already satisfy this property. This method is useful as a
    /// defensive assertion after unsafe native interoperation; it does not
    /// alter or re-wrap the value and allocates nothing.
    #[inline]
    pub fn is_in_subgroup(&self) -> bool {
        unsafe { blst_fp12_in_group(&self.0) }
    }

    /// Returns `self + self` under additive notation.
    ///
    /// Equivalently, this squares the underlying Fp12 element. The result
    /// remains in the target group and the operation allocates nothing.
    #[inline]
    pub fn double(&self) -> Self {
        unsafe {
            let mut out = MaybeUninit::<blst_fp12>::uninit();
            blst_gt_double(out.as_mut_ptr(), &self.0);
            Self(out.assume_init())
        }
    }

    /// Multiplies by a possibly secret scalar using the fixed-width
    /// constant-time target-group entry point.
    ///
    /// Under the underlying multiplicative convention this computes `self`
    /// raised to `scalar`. The scalar is encoded canonically into a 256-bit
    /// stack buffer before the native call. The operation allocates nothing.
    #[inline]
    pub fn multiply(&self, scalar: &Scalar) -> Self {
        let encoded = scalar.to_bytes_le();
        unsafe {
            let mut out = MaybeUninit::<blst_fp12>::uninit();
            blst_gt_mult(out.as_mut_ptr(), &self.0, encoded.as_ptr(), 256);
            Self(out.assume_init())
        }
    }

    /// Multiplies by a public scalar using the variable-time entry point.
    ///
    /// Under the underlying multiplicative convention this computes `self`
    /// raised to `scalar`, using a fixed 256-bit canonical scalar encoding.
    /// `scalar` must be public because execution time and memory access may
    /// depend on it. The operation allocates nothing.
    #[inline]
    pub fn multiply_vartime(&self, scalar: &Scalar) -> Self {
        let encoded = scalar.to_bytes_le();
        unsafe {
            let mut out = MaybeUninit::<blst_fp12>::uninit();
            blst_gt_mult_vartime(
                out.as_mut_ptr(),
                &self.0,
                encoded.as_ptr(),
                256,
            );
            Self(out.assume_init())
        }
    }

    /// Encodes this target-group element as twelve canonical 48-byte
    /// big-endian base-field elements.
    ///
    /// Category-blst orders the coefficients by Fp2 index, then Fp6 half, then
    /// base-field component. This order is part of the CatBLST wire contract;
    /// it need not match an in-memory `blst_fp12` traversal. The identity is
    /// encodable. The operation allocates nothing.
    #[inline]
    pub fn to_bytes_be(&self) -> [u8; 576] {
        unsafe {
            let mut out = MaybeUninit::<[u8; 576]>::uninit();
            blst_bendian_from_fp12(out.as_mut_ptr().cast::<u8>(), &self.0);
            out.assume_init()
        }
    }

    /// Strictly decodes the canonical 576-byte target-group representation.
    ///
    /// The coefficient order is Fp2 index, then Fp6 half, then base-field
    /// component, matching [`Gt::to_bytes_be`]. Each coefficient is re-encoded
    /// to reject integers at least the base-field modulus. A noncanonical
    /// coefficient returns [`BLST_ERROR::BLST_BAD_ENCODING`]. An assembled Fp12
    /// value outside the target group returns
    /// [`BLST_ERROR::BLST_POINT_NOT_IN_GROUP`]. The identity is accepted.
    /// Failure exposes no partially decoded value, and the operation allocates
    /// nothing.
    pub fn from_bytes_be(encoded: &[u8; 576]) -> Result<Self, BLST_ERROR> {
        let mut raw = blst_fp12 {
            fp6: [blst_fp6::default(), blst_fp6::default()],
        };
        let mut offset = 0;
        // Category-blst's established GT wire order interleaves the two Fp6
        // halves for each Fp2 coefficient: fp2, fp6, then base-field component.
        for fp2_index in 0..3 {
            for fp6_index in 0..2 {
                for fp_index in 0..2 {
                    let fp =
                        &mut raw.fp6[fp6_index].fp2[fp2_index].fp[fp_index];
                    unsafe {
                        blst_fp_from_bendian(fp, encoded[offset..].as_ptr());
                    }
                    offset += 48;
                }
            }
        }

        let candidate = Self(raw);
        if candidate.to_bytes_be() != *encoded {
            return Err(BLST_ERROR::BLST_BAD_ENCODING);
        }
        if !candidate.is_in_subgroup() {
            return Err(BLST_ERROR::BLST_POINT_NOT_IN_GROUP);
        }
        Ok(candidate)
    }

    #[inline]
    fn add_inner(lhs: &Self, rhs: &Self) -> Self {
        unsafe {
            let mut out = MaybeUninit::<blst_fp12>::uninit();
            blst_gt_add(out.as_mut_ptr(), &lhs.0, &rhs.0);
            Self(out.assume_init())
        }
    }

    #[inline]
    fn sub_inner(lhs: &Self, rhs: &Self) -> Self {
        unsafe {
            let mut out = MaybeUninit::<blst_fp12>::uninit();
            blst_gt_sub(out.as_mut_ptr(), &lhs.0, &rhs.0);
            Self(out.assume_init())
        }
    }
}

impl Default for Gt {
    #[inline]
    fn default() -> Self {
        Self::identity()
    }
}

impl PartialEq for Gt {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        unsafe { blst_gt_is_equal(&self.0, &other.0) }
    }
}

impl Eq for Gt {}

macro_rules! impl_gt_binary_op {
    ($trait:ident, $method:ident, $inner:ident) => {
        impl $trait for Gt {
            type Output = Gt;
            #[inline]
            fn $method(self, rhs: Gt) -> Gt {
                Gt::$inner(&self, &rhs)
            }
        }
        impl $trait<&Gt> for Gt {
            type Output = Gt;
            #[inline]
            fn $method(self, rhs: &Gt) -> Gt {
                Gt::$inner(&self, rhs)
            }
        }
        impl $trait<Gt> for &Gt {
            type Output = Gt;
            #[inline]
            fn $method(self, rhs: Gt) -> Gt {
                Gt::$inner(self, &rhs)
            }
        }
        impl $trait<&Gt> for &Gt {
            type Output = Gt;
            #[inline]
            fn $method(self, rhs: &Gt) -> Gt {
                Gt::$inner(self, rhs)
            }
        }
    };
}

impl_gt_binary_op!(Add, add, add_inner);
impl_gt_binary_op!(Sub, sub, sub_inner);

impl Neg for Gt {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        unsafe {
            let mut out = MaybeUninit::<blst_fp12>::uninit();
            blst_gt_neg(out.as_mut_ptr(), &self.0);
            Self(out.assume_init())
        }
    }
}

impl Neg for &Gt {
    type Output = Gt;
    #[inline]
    fn neg(self) -> Gt {
        -*self
    }
}

impl AddAssign for Gt {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl AddAssign<&Gt> for Gt {
    #[inline]
    fn add_assign(&mut self, rhs: &Gt) {
        *self = *self + rhs;
    }
}

impl SubAssign for Gt {
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl SubAssign<&Gt> for Gt {
    #[inline]
    fn sub_assign(&mut self, rhs: &Gt) {
        *self = *self - rhs;
    }
}

impl Mul<Scalar> for Gt {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: Scalar) -> Self {
        self.multiply(&rhs)
    }
}

impl Mul<&Scalar> for Gt {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: &Scalar) -> Self {
        self.multiply(rhs)
    }
}

impl Mul<Scalar> for &Gt {
    type Output = Gt;
    #[inline]
    fn mul(self, rhs: Scalar) -> Gt {
        self.multiply(&rhs)
    }
}

impl Mul<&Scalar> for &Gt {
    type Output = Gt;
    #[inline]
    fn mul(self, rhs: &Scalar) -> Gt {
        self.multiply(rhs)
    }
}

impl Sum for Gt {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::identity(), |sum, value| sum + value)
    }
}

impl<'a> Sum<&'a Gt> for Gt {
    fn sum<I: Iterator<Item = &'a Gt>>(iter: I) -> Self {
        iter.fold(Self::identity(), |sum, value| sum + value)
    }
}

/// Computes the pairing of validated affine subgroup points.
///
/// This is equivalent to [`Gt::pairing`], including its input invariants,
/// identity behavior, timing scope, and allocation behavior.
#[inline]
pub fn pairing(g1: &G1Affine, g2: &G2Affine) -> Gt {
    Gt::pairing(g1, g2)
}

/// Shape errors reported before an allocation-free batched pairing call.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum PairingError {
    /// The G1 and G2 slices contain different numbers of points.
    LengthMismatch {
        /// Number of G1 inputs supplied by the caller.
        g1: usize,
        /// Number of G2 inputs supplied by the caller.
        g2: usize,
    },
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::mem::{align_of, size_of};

    #[test]
    fn wrappers_are_native_transparent_storage() {
        assert_eq!(size_of::<Scalar>(), size_of::<blst_fr>());
        assert_eq!(align_of::<Scalar>(), align_of::<blst_fr>());
        assert_eq!(size_of::<G1Projective>(), size_of::<blst_p1>());
        assert_eq!(align_of::<G1Projective>(), align_of::<blst_p1>());
        assert_eq!(size_of::<G1Affine>(), size_of::<blst_p1_affine>());
        assert_eq!(align_of::<G1Affine>(), align_of::<blst_p1_affine>());
        assert_eq!(size_of::<G2Projective>(), size_of::<blst_p2>());
        assert_eq!(align_of::<G2Projective>(), align_of::<blst_p2>());
        assert_eq!(size_of::<G2Affine>(), size_of::<blst_p2_affine>());
        assert_eq!(align_of::<G2Affine>(), align_of::<blst_p2_affine>());
        assert_eq!(size_of::<MillerLoopResult>(), size_of::<blst_fp12>());
        assert_eq!(align_of::<MillerLoopResult>(), align_of::<blst_fp12>());
        assert_eq!(size_of::<Gt>(), size_of::<blst_fp12>());
        assert_eq!(align_of::<Gt>(), align_of::<blst_fp12>());

        let scalar = Scalar::from(7);
        assert_eq!(
            scalar.as_raw() as *const _,
            &scalar as *const _ as *const _
        );
        let g1 = G1Projective::generator();
        assert_eq!(g1.as_raw() as *const _, &g1 as *const _ as *const _);
        let gt = Gt::generator();
        assert_eq!(gt.as_raw() as *const _, &gt as *const _ as *const _);
    }

    #[test]
    fn constant_time_g1_normalization_matches_public_conversion() {
        let points = [
            G1Projective::identity(),
            G1Projective::generator(),
            G1Projective::generator() * Scalar::from(7),
            G1Projective::generator() * Scalar::from(65537),
        ];
        let expected = points.map(|point| point.to_affine());
        for (point, affine) in points.iter().zip(expected.iter()) {
            assert_eq!(point.to_affine_ct(), *affine);
        }
        let mut output = [G1Affine::identity(); 4];
        G1Projective::batch_to_affine_ct(&points, &mut output).unwrap();
        assert_eq!(output, expected);
        assert_eq!(G1Projective::batch_to_affine_ct(&[], &mut []), Ok(()));
        assert_eq!(
            G1Projective::batch_to_affine_ct(&points, &mut output[..3]),
            Err(AffineBatchLengthError {
                expected: 4,
                actual: 3,
            })
        );
    }

    #[test]
    fn scalar_constants_and_canonical_encodings_are_correct() {
        assert!(Scalar::ZERO.is_zero());
        assert_eq!(Scalar::ONE, Scalar::from(1));
        assert_eq!(Scalar::ZERO.to_bytes_be(), [0; 32]);
        assert_eq!(Scalar::ZERO.to_bytes_le(), [0; 32]);

        let scalar = Scalar::from(0x0123_4567_89ab_cdef);
        let mut expected_be = [0u8; 32];
        expected_be[24..]
            .copy_from_slice(&0x0123_4567_89ab_cdefu64.to_be_bytes());
        let mut expected_le = [0u8; 32];
        expected_le[..8]
            .copy_from_slice(&0x0123_4567_89ab_cdefu64.to_le_bytes());
        assert_eq!(scalar.to_bytes_be(), expected_be);
        assert_eq!(scalar.to_bytes_le(), expected_le);
        assert_eq!(Scalar::from_bytes_be(&expected_be), Some(scalar));
        assert_eq!(Scalar::from_bytes_le(&expected_le), Some(scalar));

        let noncanonical = [0xff; 32];
        assert_eq!(Scalar::from_bytes_be(&noncanonical), None);
        assert_eq!(Scalar::from_bytes_le(&noncanonical), None);
    }

    #[test]
    fn wide_scalar_reduction_and_field_arithmetic_agree() {
        let mut wide_be = [0u8; 64];
        wide_be[56..].copy_from_slice(&19u64.to_be_bytes());
        let mut wide_le = [0u8; 64];
        wide_le[..8].copy_from_slice(&19u64.to_le_bytes());
        assert_eq!(Scalar::from_bytes_wide_be(&wide_be), Scalar::from(19));
        assert_eq!(Scalar::from_bytes_wide_le(&wide_le), Scalar::from(19));

        let seven = Scalar::from(7);
        let eleven = Scalar::from(11);
        assert_eq!(seven + eleven, Scalar::from(18));
        assert_eq!(eleven - seven, Scalar::from(4));
        assert_eq!(seven * eleven, Scalar::from(77));
        assert_eq!(seven * seven.invert().unwrap(), Scalar::ONE);
        assert_eq!(Scalar::ZERO.invert(), None);
        assert_eq!(seven + (-seven), Scalar::ZERO);
    }

    #[test]
    fn scalar_group_identities_and_roundtrips_are_validated() {
        let g1_identity = G1Projective::identity();
        assert!(g1_identity.is_identity());
        assert!(g1_identity.is_on_curve());
        assert!(g1_identity.is_in_subgroup());
        assert!(G1Affine::from(g1_identity).is_identity());

        let g2_identity = G2Projective::identity();
        assert!(g2_identity.is_identity());
        assert!(g2_identity.is_on_curve());
        assert!(g2_identity.is_in_subgroup());
        assert!(G2Affine::from(g2_identity).is_identity());

        let g1 = G1Projective::generator() * Scalar::from(17);
        let g2 = G2Projective::generator() * Scalar::from(29);
        assert_eq!(G1Projective::from_compressed(&g1.to_compressed()), Ok(g1));
        assert_eq!(
            G1Projective::from_uncompressed(&g1.to_uncompressed()),
            Ok(g1)
        );
        assert_eq!(G2Projective::from_compressed(&g2.to_compressed()), Ok(g2));
        assert_eq!(
            G2Projective::from_uncompressed(&g2.to_uncompressed()),
            Ok(g2)
        );

        assert_eq!(g1 + g1, g1.double());
        assert_eq!(g2 + g2, g2.double());
        assert_eq!(g1 - g1, G1Projective::identity());
        assert_eq!(g2 - g2, G2Projective::identity());
        assert_eq!(g1 + (-g1), G1Projective::identity());
        assert_eq!(g2 + (-g2), G2Projective::identity());
    }

    #[test]
    fn malformed_point_encodings_are_rejected() {
        assert_ne!(
            G1Affine::from_compressed(&[0xff; 48]),
            Ok(G1Affine::identity())
        );
        assert_ne!(
            G2Affine::from_compressed(&[0xff; 96]),
            Ok(G2Affine::identity())
        );
    }

    #[test]
    fn pairing_and_batched_miller_loop_match_group_arithmetic() {
        let g1 = G1Affine::generator();
        let g2 = G2Affine::generator();
        let pairing_value = pairing(&g1, &g2);
        assert_eq!(pairing_value, Gt::generator());
        assert!(pairing_value.is_in_subgroup());

        let scalar = Scalar::from(13);
        let scaled_g1 = G1Affine::from(G1Projective::from(g1) * scalar);
        assert_eq!(pairing(&scaled_g1, &g2), pairing_value * scalar);

        let g1s = [g1, scaled_g1];
        let g2s = [g2, g2];
        let batch = MillerLoopResult::from_pairs(&g1s, &g2s)
            .unwrap()
            .final_exponentiation();
        assert_eq!(
            batch,
            pairing(&g1s[0], &g2s[0]) + pairing(&g1s[1], &g2s[1])
        );
        assert_eq!(
            MillerLoopResult::from_pairs(&g1s, &g2s[..1]),
            Err(PairingError::LengthMismatch { g1: 2, g2: 1 })
        );
    }

    #[test]
    fn prepared_g2_miller_products_compose_before_final_exponentiation() {
        let g1 = G1Projective::generator();
        let g1_doubled = g1 * Scalar::from(2);
        let g1_identity = G1Projective::identity();
        let g2 = G2Projective::generator();
        let g2_tripled = g2 * Scalar::from(3);
        let g2_identity = G2Projective::identity();
        let batch =
            PreparedG2Batch::new(&[g2, g2_tripled, g2_identity]).unwrap();
        let g1_points = [g1, g1_doubled, g1_identity];
        let indices = [0usize, 1, 2];

        let complete = batch
            .miller_loop_product_indices(&g1_points, &indices)
            .unwrap()
            .final_exponentiation();
        let left = batch
            .miller_loop_product_indices(&g1_points[..1], &indices[..1])
            .unwrap();
        let right = batch
            .miller_loop_product_indices(&g1_points[1..], &indices[1..])
            .unwrap();
        let split = (left * right).final_exponentiation();
        let expected = Gt::pairing(&g1.to_affine(), &g2.to_affine())
            + Gt::pairing(&g1_doubled.to_affine(), &g2_tripled.to_affine());

        assert_eq!(complete, expected);
        assert_eq!(split, expected);
        assert_eq!(
            batch
                .miller_loop_product_indices(&[], &[])
                .unwrap()
                .final_exponentiation(),
            Gt::identity()
        );
        assert_eq!(batch.len(), 3);
        assert!(!batch.is_empty());
        assert!(batch.storage_bytes() > 0);
    }

    #[test]
    fn prepared_g2_pairing_ranges_match_direct_pairings() {
        let g1 = G1Projective::generator() * Scalar::from(5);
        let generator = G2Projective::generator();
        let operands = [
            generator,
            generator * Scalar::from(2),
            G2Projective::identity(),
            generator * Scalar::from(7),
        ];
        let batch = PreparedG2Batch::new(&operands).unwrap();
        let mut out = [Gt::identity(); 3];
        batch.pairing_range_into(&g1, 1, &mut out).unwrap();
        for (offset, actual) in out.iter().enumerate() {
            assert_eq!(
                *actual,
                Gt::pairing(&g1.to_affine(), &operands[1 + offset].to_affine())
            );
        }
        batch.pairing_range_into(&g1, batch.len(), &mut []).unwrap();
    }

    #[test]
    #[cfg_attr(
        target_arch = "wasm32",
        ignore = "requires native threads; sequential partition/range tests remain enabled"
    )]
    fn prepared_g2_miller_products_match_partitions_and_concurrent_ranges() {
        const COUNT: usize = 37;
        let g1_generator = G1Projective::generator();
        let g2_generator = G2Projective::generator();
        let g1_identity = G1Projective::identity();
        let g2_identity = G2Projective::identity();
        let g1: Vec<_> = (0..COUNT)
            .map(|i| {
                if i % 11 == 0 {
                    g1_identity
                } else {
                    g1_generator * Scalar::from((i + 1) as u64)
                }
            })
            .collect();
        let g2: Vec<_> = (0..COUNT)
            .map(|i| {
                if i % 13 == 0 {
                    g2_identity
                } else {
                    g2_generator * Scalar::from((2 * i + 1) as u64)
                }
            })
            .collect();
        let indices: Vec<_> = (0..COUNT).map(|i| (7 * i + 3) % COUNT).collect();
        let batch = PreparedG2Batch::new(&g2).unwrap();

        let expected = batch
            .miller_loop_product_indices(&g1, &indices)
            .unwrap()
            .final_exponentiation();
        for grain in [1usize, 2, 5, 10, 16, 31, COUNT] {
            let product = g1
                .chunks(grain)
                .zip(indices.chunks(grain))
                .map(|(points, prepared_indices)| {
                    batch
                        .miller_loop_product_indices(points, prepared_indices)
                        .unwrap()
                })
                .fold(MillerLoopResult::identity(), |left, right| left * right);
            assert_eq!(product.final_exponentiation(), expected);
        }

        let concurrent = std::thread::scope(|scope| {
            let batch_ref = &batch;
            let g1_ref = &g1;
            let indices_ref = &indices;
            let handles: Vec<_> = [(0, 10), (10, 20), (20, 30), (30, COUNT)]
                .iter()
                .copied()
                .map(|(begin, end)| {
                    scope.spawn(move || {
                        batch_ref
                            .miller_loop_product_indices(
                                &g1_ref[begin..end],
                                &indices_ref[begin..end],
                            )
                            .unwrap()
                    })
                })
                .collect();
            handles
                .into_iter()
                .rev()
                .map(|handle| handle.join().unwrap())
                .fold(MillerLoopResult::identity(), |left, right| left * right)
                .final_exponentiation()
        });
        assert_eq!(concurrent, expected);
    }

    #[test]
    fn prepared_g2_batch_rejects_invalid_shapes() {
        let g1 = G1Projective::generator();
        let batch = PreparedG2Batch::new(&[G2Projective::generator()]).unwrap();

        assert!(matches!(
            PreparedG2Batch::new(&[]),
            Err(PreparedG2BatchError::EmptyBatch)
        ));
        assert_eq!(
            batch.miller_loop_product_indices(&[g1], &[]),
            Err(PreparedG2BatchError::InputLengthMismatch {
                g1: 1,
                indices: 0,
            })
        );
        assert_eq!(
            batch.miller_loop_product_indices(&[g1], &[1]),
            Err(PreparedG2BatchError::IndexOutOfRange { index: 1, len: 1 })
        );
        let mut out = [Gt::identity(); 1];
        assert_eq!(
            batch.pairing_range_into(&g1, 1, &mut out),
            Err(PreparedG2BatchError::IndexOutOfRange { index: 1, len: 1 })
        );
    }

    #[test]
    fn target_group_encoding_is_canonical_and_validated() {
        let gt = Gt::generator() * Scalar::from(7919);
        let encoded = gt.to_bytes_be();
        assert_eq!(Gt::from_bytes_be(&encoded), Ok(gt));

        assert_eq!(
            Gt::from_bytes_be(&[0xff; 576]),
            Err(BLST_ERROR::BLST_BAD_ENCODING)
        );
    }

    #[cfg(feature = "avx")]
    #[test]
    fn canonical_encodings_match_native_avx_representations() {
        if !unsafe { crate::avx::blst_avx_enabled() } {
            return;
        }

        let scalar = Scalar::from(37);

        let scalar_g1 = G1Projective::generator() * scalar;
        let mut avx_g1 = crate::avx::blst_avx_p1::default();
        assert!(unsafe {
            crate::avx::blst_avx_p1_from_scalar(&mut avx_g1, scalar_g1.as_raw())
        });
        let mut avx_g1_encoded = [0u8; 48];
        assert!(unsafe {
            crate::avx::blst_avx_p1_compress(
                avx_g1_encoded.as_mut_ptr(),
                &avx_g1,
            )
        });
        assert_eq!(avx_g1_encoded, scalar_g1.to_compressed());

        let scalar_g2 = G2Projective::generator() * scalar;
        let mut avx_g2 = crate::avx::blst_avx_p2::default();
        assert!(unsafe {
            crate::avx::blst_avx_p2_from_scalar(&mut avx_g2, scalar_g2.as_raw())
        });
        let mut avx_g2_encoded = [0u8; 96];
        assert!(unsafe {
            crate::avx::blst_avx_p2_compress(
                avx_g2_encoded.as_mut_ptr(),
                &avx_g2,
            )
        });
        assert_eq!(avx_g2_encoded, scalar_g2.to_compressed());

        let scalar_gt = Gt::generator() * scalar;
        let mut avx_gt = crate::avx::blst_avx_gt::default();
        assert!(unsafe {
            crate::avx::blst_avx_fp12_from_scalar(
                &mut avx_gt,
                scalar_gt.as_raw(),
            )
        });
        let mut avx_gt_encoded = [0u8; 576];
        assert!(unsafe {
            crate::avx::blst_avx_gt_serialize(
                avx_gt_encoded.as_mut_ptr(),
                &avx_gt,
            )
        });
        assert_eq!(avx_gt_encoded, scalar_gt.to_bytes_be());
    }
}
