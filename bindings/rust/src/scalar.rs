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
//! The types in this module are not representation-compatible with the native
//! types in [`crate::avx`]. A protocol instance should retain one backend's
//! native point and target-group types. Canonical encodings are the compatibility
//! boundary between scalar and AVX implementations.

use crate::*;
use core::fmt;
use core::iter::Sum;
use core::mem::MaybeUninit;
use core::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

/// A reduced scalar-field element in category-blst's native Montgomery form.
///
/// This type permits zero and is distinct from [`blst_scalar`], the canonical
/// non-Montgomery byte representation used by the raw secret-key API.
#[repr(transparent)]
#[derive(Copy, Clone, PartialEq, Eq)]
pub struct Scalar(blst_fr);

impl Scalar {
    /// The additive identity.
    pub const ZERO: Self = Self(blst_fr { l: [0; 4] });

    /// The multiplicative identity in category-blst's Montgomery representation.
    pub const ONE: Self = Self(blst_fr {
        l: [
            0x0000_0001_ffff_fffe,
            0x5884_b7fa_0003_4802,
            0x998c_4fef_ecbc_4ff5,
            0x1824_b159_acc5_056f,
        ],
    });

    /// Borrows the native Montgomery-form value without conversion.
    #[inline]
    pub const fn as_raw(&self) -> &blst_fr {
        &self.0
    }

    /// Mutably borrows the native Montgomery-form value without conversion.
    ///
    /// # Safety
    ///
    /// The caller must leave a reduced field element in the object.
    #[inline]
    pub unsafe fn as_raw_mut(&mut self) -> &mut blst_fr {
        &mut self.0
    }

    /// Wraps a reduced native Montgomery-form scalar without conversion.
    ///
    /// # Safety
    ///
    /// `raw` must be reduced modulo the BLS12-381 scalar-field order.
    #[inline]
    pub const unsafe fn from_raw_unchecked(raw: blst_fr) -> Self {
        Self(raw)
    }

    /// Returns the native Montgomery-form value without conversion.
    #[inline]
    pub const fn into_raw(self) -> blst_fr {
        self.0
    }

    /// Returns `true` when this scalar is zero.
    #[inline]
    pub fn is_zero(&self) -> bool {
        self.0 == Self::ZERO.0
    }

    /// Strictly decodes a canonical 32-byte big-endian scalar.
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

    /// Reduces a 64-byte big-endian integer modulo the scalar-field order.
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

    /// Reduces a 64-byte little-endian integer modulo the scalar-field order.
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

    /// Encodes this scalar canonically in big-endian order.
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

    /// Encodes this scalar canonically in little-endian order.
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

    /// Computes the multiplicative inverse, or returns `None` for zero.
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

/// A prime-order G1 point in the scalar backend's projective representation.
#[repr(transparent)]
#[derive(Copy, Clone, Debug)]
pub struct G1Projective(blst_p1);

/// A prime-order G1 point in the scalar backend's affine representation.
#[repr(transparent)]
#[derive(Copy, Clone, Debug)]
pub struct G1Affine(blst_p1_affine);

/// A prime-order G2 point in the scalar backend's projective representation.
#[repr(transparent)]
#[derive(Copy, Clone, Debug)]
pub struct G2Projective(blst_p2);

/// A prime-order G2 point in the scalar backend's affine representation.
#[repr(transparent)]
#[derive(Copy, Clone, Debug)]
pub struct G2Affine(blst_p2_affine);

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
            /// Returns the group identity.
            #[inline]
            pub fn identity() -> Self {
                Self(<$raw_projective>::default())
            }

            /// Returns the standard group generator.
            #[inline]
            pub fn generator() -> Self {
                unsafe { Self(*$generator()) }
            }

            /// Borrows the native scalar-backend representation without conversion.
            #[inline]
            pub const fn as_raw(&self) -> &$raw_projective {
                &self.0
            }

            /// Mutably borrows the native representation without conversion.
            ///
            /// # Safety
            ///
            /// The caller must preserve the on-curve, prime-order subgroup invariant.
            #[inline]
            pub unsafe fn as_raw_mut(&mut self) -> &mut $raw_projective {
                &mut self.0
            }

            /// Wraps a native projective subgroup point without conversion.
            ///
            /// # Safety
            ///
            /// `raw` must be on-curve and in the prime-order subgroup.
            #[inline]
            pub const unsafe fn from_raw_unchecked(
                raw: $raw_projective,
            ) -> Self {
                Self(raw)
            }

            /// Returns the native projective value without conversion.
            #[inline]
            pub const fn into_raw(self) -> $raw_projective {
                self.0
            }

            /// Returns `true` when this point is the group identity.
            #[inline]
            pub fn is_identity(&self) -> bool {
                unsafe { $is_inf(&self.0) }
            }

            /// Returns `true` when the underlying value lies on the curve.
            #[inline]
            pub fn is_on_curve(&self) -> bool {
                unsafe { $on_curve(&self.0) }
            }

            /// Returns `true` when the underlying value lies in the prime-order subgroup.
            #[inline]
            pub fn is_in_subgroup(&self) -> bool {
                unsafe { $in_group(&self.0) }
            }

            /// `blstrs`-inspired spelling for [`Self::is_in_subgroup`].
            #[inline]
            pub fn is_torsion_free(&self) -> bool {
                self.is_in_subgroup()
            }

            /// Doubles this point.
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

            /// Converts this point to affine form.
            #[inline]
            pub fn to_affine(&self) -> $affine {
                unsafe {
                    let mut out = MaybeUninit::<$raw_affine>::uninit();
                    $to_affine(out.as_mut_ptr(), &self.0);
                    $affine(out.assume_init())
                }
            }

            /// Encodes this point in canonical compressed form.
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
            #[inline]
            pub fn from_compressed(
                encoded: &[u8; $compressed_size],
            ) -> Result<Self, BLST_ERROR> {
                $affine::from_compressed(encoded).map(Self::from)
            }

            /// Strictly decodes and subgroup-validates an uncompressed point.
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
            /// Returns the affine encoding of the group identity.
            #[inline]
            pub fn identity() -> Self {
                Self(<$raw_affine>::default())
            }

            /// Returns the standard affine group generator.
            #[inline]
            pub fn generator() -> Self {
                unsafe { Self(*$affine_generator()) }
            }

            /// Borrows the native scalar-backend representation without conversion.
            #[inline]
            pub const fn as_raw(&self) -> &$raw_affine {
                &self.0
            }

            /// Mutably borrows the native representation without conversion.
            ///
            /// # Safety
            ///
            /// The caller must preserve the on-curve, prime-order subgroup invariant.
            #[inline]
            pub unsafe fn as_raw_mut(&mut self) -> &mut $raw_affine {
                &mut self.0
            }

            /// Wraps a native affine subgroup point without conversion.
            ///
            /// # Safety
            ///
            /// `raw` must be on-curve and in the prime-order subgroup.
            #[inline]
            pub const unsafe fn from_raw_unchecked(raw: $raw_affine) -> Self {
                Self(raw)
            }

            /// Returns the native affine value without conversion.
            #[inline]
            pub const fn into_raw(self) -> $raw_affine {
                self.0
            }

            /// Returns `true` when this point is the group identity.
            #[inline]
            pub fn is_identity(&self) -> bool {
                unsafe { $affine_is_inf(&self.0) }
            }

            /// Returns `true` when the underlying value lies on the curve.
            #[inline]
            pub fn is_on_curve(&self) -> bool {
                unsafe { $affine_on_curve(&self.0) }
            }

            /// Returns `true` when the underlying value lies in the prime-order subgroup.
            #[inline]
            pub fn is_in_subgroup(&self) -> bool {
                unsafe { $affine_in_group(&self.0) }
            }

            /// `blstrs`-inspired spelling for [`Self::is_in_subgroup`].
            #[inline]
            pub fn is_torsion_free(&self) -> bool {
                self.is_in_subgroup()
            }

            /// Encodes this point in canonical compressed form.
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
/// does not by itself establish target-group membership.
#[repr(transparent)]
#[derive(Copy, Clone, Debug)]
pub struct MillerLoopResult(blst_fp12);

impl MillerLoopResult {
    /// Returns the empty Miller-loop product.
    #[inline]
    pub fn identity() -> Self {
        unsafe { Self(*blst_fp12_one()) }
    }

    /// Computes one Miller loop from validated subgroup points.
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
    /// This operation uses fixed-size descriptor storage on the stack and does
    /// not allocate. It returns an error when the slice lengths differ.
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

    /// Applies final exponentiation and returns a target-group element.
    #[inline]
    pub fn final_exponentiation(&self) -> Gt {
        unsafe {
            let mut out = MaybeUninit::<blst_fp12>::uninit();
            blst_final_exp(out.as_mut_ptr(), &self.0);
            Gt(out.assume_init())
        }
    }

    /// Borrows the native Fp12 Miller-loop storage without conversion.
    #[inline]
    pub const fn as_raw(&self) -> &blst_fp12 {
        &self.0
    }

    /// Wraps a trusted native Miller-loop result without conversion.
    ///
    /// # Safety
    ///
    /// `raw` must be a nonzero product of valid Miller-loop contributions.
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

/// A validated element of the BLS12-381 pairing target group.
///
/// Group operations use the additive conventions familiar from `blstrs` even
/// though the underlying Fp12 group law is multiplication.
#[repr(transparent)]
#[derive(Copy, Clone, Debug)]
pub struct Gt(blst_fp12);

impl Gt {
    /// Returns the target-group identity.
    #[inline]
    pub fn identity() -> Self {
        unsafe { Self(*blst_gt_identity()) }
    }

    /// Returns the pairing of the standard G1 and G2 generators.
    #[inline]
    pub fn generator() -> Self {
        unsafe { Self(*blst_gt_generator()) }
    }

    /// Computes a pairing from validated subgroup points.
    #[inline]
    pub fn pairing(g1: &G1Affine, g2: &G2Affine) -> Self {
        MillerLoopResult::from_pair(g1, g2).final_exponentiation()
    }

    /// Borrows the native target-group representation without conversion.
    #[inline]
    pub const fn as_raw(&self) -> &blst_fp12 {
        &self.0
    }

    /// Mutably borrows the native target-group representation without conversion.
    ///
    /// # Safety
    ///
    /// The caller must preserve target-group membership.
    #[inline]
    pub unsafe fn as_raw_mut(&mut self) -> &mut blst_fp12 {
        &mut self.0
    }

    /// Wraps a native target-group value without conversion.
    ///
    /// # Safety
    ///
    /// `raw` must lie in the BLS12-381 target group.
    #[inline]
    pub const unsafe fn from_raw_unchecked(raw: blst_fp12) -> Self {
        Self(raw)
    }

    /// Returns the native target-group value without conversion.
    #[inline]
    pub const fn into_raw(self) -> blst_fp12 {
        self.0
    }

    /// Returns `true` when this is the target-group identity.
    #[inline]
    pub fn is_identity(&self) -> bool {
        unsafe { blst_gt_is_identity(&self.0) }
    }

    /// Returns `true` when the underlying Fp12 value lies in the target group.
    #[inline]
    pub fn is_in_subgroup(&self) -> bool {
        unsafe { blst_fp12_in_group(&self.0) }
    }

    /// Doubles this target-group element under additive notation.
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
    /// Each Fp coefficient is re-encoded to reject values at least the base
    /// field modulus, and the assembled Fp12 value is checked for target-group
    /// membership.
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
#[inline]
pub fn pairing(g1: &G1Affine, g2: &G2Affine) -> Gt {
    Gt::pairing(g1, g2)
}

/// Shape errors reported by allocation-free batched pairing helpers.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum PairingError {
    /// G1 and G2 slices must contain the same number of points.
    LengthMismatch { g1: usize, g2: usize },
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
