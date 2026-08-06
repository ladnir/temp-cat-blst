// Copyright Supranational LLC
// Licensed under the Apache License, Version 2.0, see LICENSE for details.
// SPDX-License-Identifier: Apache-2.0

//! Low-level scalar-ABI Rust conveniences.
//!
//! This module adds thin, allocation-free Rust helpers on top of the raw C ABI
//! types:
//! - equality for point and Fp12 types
//! - Fp12 arithmetic helpers
//! - scalar hashing helper
//! - small Rust-side smoke tests of the raw boundary types

use crate::*;
use core::mem::MaybeUninit;
use core::ptr;

fn checked_scalar_bits(scalar: &[u8], nbits: usize) {
    assert!(nbits <= 256, "scalar bit length exceeds 256");
    assert!(
        scalar.len() >= (nbits + 7) / 8,
        "scalar byte slice is too short"
    );
}

impl blst_p1 {
    /// Multiplies this previously validated G1 element by `scalar mod r`.
    ///
    /// This method does not validate that `self` is on-curve and in the
    /// prime-order subgroup. Canonical scalars are handled in constant time for
    /// fixed public `nbits`; noncanonical inputs are accepted, but timing may
    /// reveal that reduction was required.
    ///
    /// # Panics
    ///
    /// Panics when `nbits > 256` or `scalar` is shorter than
    /// `ceil(nbits / 8)` bytes.
    pub fn mult(&self, scalar: &[u8], nbits: usize) -> Self {
        checked_scalar_bits(scalar, nbits);
        let mut out = MaybeUninit::<Self>::uninit();
        unsafe {
            blst_p1_mult(out.as_mut_ptr(), self, scalar.as_ptr(), nbits);
            out.assume_init()
        }
    }

    /// Optimized variable-time subgroup multiplication for public scalars.
    ///
    /// `self` must already be a validated prime-order G1 element. This method
    /// does not validate membership and must not be used with a secret scalar.
    ///
    /// # Panics
    ///
    /// Panics when `nbits > 256` or `scalar` is too short.
    pub fn mult_vartime(&self, scalar: &[u8], nbits: usize) -> Self {
        checked_scalar_bits(scalar, nbits);
        let mut out = MaybeUninit::<Self>::uninit();
        unsafe {
            blst_p1_mult_vartime(
                out.as_mut_ptr(),
                self,
                scalar.as_ptr(),
                nbits,
            );
            out.assume_init()
        }
    }

    /// Legacy upstream `unchecked_mult`: constant-time exact multiplication of
    /// an on-curve point that need not belong to G1.
    ///
    /// This method performs no curve validation and applies the exact unsigned
    /// scalar rather than reducing it as part of a subgroup contract.
    ///
    /// # Panics
    ///
    /// Panics when `nbits > 256` or `scalar` is too short.
    pub fn unchecked_mult(&self, scalar: &[u8], nbits: usize) -> Self {
        checked_scalar_bits(scalar, nbits);
        let mut out = MaybeUninit::<Self>::uninit();
        unsafe {
            blst_p1_unchecked_mult(
                out.as_mut_ptr(),
                self,
                scalar.as_ptr(),
                nbits,
            );
            out.assume_init()
        }
    }
}

impl blst_p2 {
    /// Multiplies this previously validated G2 element by `scalar mod r`.
    ///
    /// This method does not validate that `self` is on-curve and in the
    /// prime-order subgroup. Canonical scalars are handled in constant time for
    /// fixed public `nbits`; noncanonical inputs are accepted, but timing may
    /// reveal that reduction was required.
    ///
    /// # Panics
    ///
    /// Panics when `nbits > 256` or `scalar` is shorter than
    /// `ceil(nbits / 8)` bytes.
    pub fn mult(&self, scalar: &[u8], nbits: usize) -> Self {
        checked_scalar_bits(scalar, nbits);
        let mut out = MaybeUninit::<Self>::uninit();
        unsafe {
            blst_p2_mult(out.as_mut_ptr(), self, scalar.as_ptr(), nbits);
            out.assume_init()
        }
    }

    /// Optimized variable-time subgroup multiplication for public scalars.
    ///
    /// `self` must already be a validated prime-order G2 element. This method
    /// does not validate membership and must not be used with a secret scalar.
    ///
    /// # Panics
    ///
    /// Panics when `nbits > 256` or `scalar` is too short.
    pub fn mult_vartime(&self, scalar: &[u8], nbits: usize) -> Self {
        checked_scalar_bits(scalar, nbits);
        let mut out = MaybeUninit::<Self>::uninit();
        unsafe {
            blst_p2_mult_vartime(
                out.as_mut_ptr(),
                self,
                scalar.as_ptr(),
                nbits,
            );
            out.assume_init()
        }
    }

    /// Legacy upstream `unchecked_mult`: constant-time exact multiplication of
    /// an on-curve point that need not belong to G2.
    ///
    /// This method performs no curve validation and applies the exact unsigned
    /// scalar rather than reducing it as part of a subgroup contract.
    ///
    /// # Panics
    ///
    /// Panics when `nbits > 256` or `scalar` is too short.
    pub fn unchecked_mult(&self, scalar: &[u8], nbits: usize) -> Self {
        checked_scalar_bits(scalar, nbits);
        let mut out = MaybeUninit::<Self>::uninit();
        unsafe {
            blst_p2_unchecked_mult(
                out.as_mut_ptr(),
                self,
                scalar.as_ptr(),
                nbits,
            );
            out.assume_init()
        }
    }
}

/// Equality of on-curve values. This implementation does not validate either operand.
impl PartialEq for blst_p1 {
    fn eq(&self, other: &Self) -> bool {
        unsafe { blst_p1_is_equal(self, other) }
    }
}

/// Equality of on-curve values. This implementation does not validate either operand.
impl PartialEq for blst_p1_affine {
    fn eq(&self, other: &Self) -> bool {
        unsafe { blst_p1_affine_is_equal(self, other) }
    }
}

/// Equality of on-curve values. This implementation does not validate either operand.
impl PartialEq for blst_p2 {
    fn eq(&self, other: &Self) -> bool {
        unsafe { blst_p2_is_equal(self, other) }
    }
}

/// Equality of on-curve values. This implementation does not validate either operand.
impl PartialEq for blst_p2_affine {
    fn eq(&self, other: &Self) -> bool {
        unsafe { blst_p2_affine_is_equal(self, other) }
    }
}

impl Default for blst_fp12 {
    fn default() -> Self {
        unsafe { *blst_fp12_one() }
    }
}

impl PartialEq for blst_fp12 {
    fn eq(&self, other: &Self) -> bool {
        unsafe { blst_fp12_is_equal(self, other) }
    }
}

impl core::ops::Mul for blst_fp12 {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        let mut out = MaybeUninit::<blst_fp12>::uninit();
        unsafe {
            blst_fp12_mul(out.as_mut_ptr(), &self, &other);
            out.assume_init()
        }
    }
}

impl core::ops::MulAssign for blst_fp12 {
    fn mul_assign(&mut self, other: Self) {
        unsafe { blst_fp12_mul(self, self, &other) }
    }
}

impl blst_fp12 {
    /// Multiplies this GT element by `scalar mod r`.
    ///
    /// `self` must already lie in GT; storage and this method do not validate
    /// membership. Canonical scalars are handled in constant time for fixed
    /// public `nbits`. Noncanonical inputs are accepted, but timing may reveal
    /// that reduction was required.
    ///
    /// # Panics
    ///
    /// Panics when `nbits > 256` or `scalar` is too short.
    pub fn gt_mult(&self, scalar: &[u8], nbits: usize) -> Self {
        checked_scalar_bits(scalar, nbits);
        let mut out = MaybeUninit::<Self>::uninit();
        unsafe {
            blst_gt_mult(out.as_mut_ptr(), self, scalar.as_ptr(), nbits);
            out.assume_init()
        }
    }

    /// Optimized variable-time GT multiplication for public scalars only.
    ///
    /// `self` must already lie in GT. This method does not validate membership
    /// and must not be used with a secret scalar.
    ///
    /// # Panics
    ///
    /// Panics when `nbits > 256` or `scalar` is too short.
    pub fn gt_mult_vartime(&self, scalar: &[u8], nbits: usize) -> Self {
        checked_scalar_bits(scalar, nbits);
        let mut out = MaybeUninit::<Self>::uninit();
        unsafe {
            blst_gt_mult_vartime(
                out.as_mut_ptr(),
                self,
                scalar.as_ptr(),
                nbits,
            );
            out.assume_init()
        }
    }

    /// Computes one Miller loop from valid on-curve affine inputs.
    ///
    /// The inputs need not be subgroup elements for the formula itself, but
    /// this method performs no curve or subgroup validation. The result is a
    /// general non-final-exponentiated Fp12 value, not yet a GT element.
    pub fn miller_loop(q: &blst_p2_affine, p: &blst_p1_affine) -> Self {
        let mut out = MaybeUninit::<blst_fp12>::uninit();
        unsafe {
            blst_miller_loop(out.as_mut_ptr(), q, p);
            out.assume_init()
        }
    }

    /// Computes the product of Miller loops over paired affine slices.
    ///
    /// Every input must be a valid on-curve affine point. No curve or subgroup
    /// validation is performed. Empty slices return the Fp12 identity.
    ///
    /// # Panics
    ///
    /// Panics when `q.len() != p.len()`.
    pub fn miller_loop_n(q: &[blst_p2_affine], p: &[blst_p1_affine]) -> Self {
        let n_elems = q.len();
        if n_elems != p.len() {
            panic!("inputs' lengths mismatch");
        }
        if n_elems == 0 {
            return Self::default();
        }
        let qs: [*const _; 2] = [&q[0], ptr::null()];
        let ps: [*const _; 2] = [&p[0], ptr::null()];
        let mut out = MaybeUninit::<blst_fp12>::uninit();
        unsafe {
            blst_miller_loop_n(out.as_mut_ptr(), &qs[0], &ps[0], n_elems);
            out.assume_init()
        }
    }

    /// Applies final exponentiation to a nonzero general Fp12/Miller-loop value.
    ///
    /// The result lies in GT. Zero is outside the supported domain and is not
    /// reported through a `Result`.
    pub fn final_exp(&self) -> Self {
        let mut out = MaybeUninit::<blst_fp12>::uninit();
        unsafe {
            blst_final_exp(out.as_mut_ptr(), self);
            out.assume_init()
        }
    }

    /// Returns `true` when this general Fp12 value lies in GT.
    pub fn in_group(&self) -> bool {
        unsafe { blst_fp12_in_group(self) }
    }

    /// Returns `true` when two trusted Miller-loop contributions have equal
    /// images after final exponentiation.
    ///
    /// The raw `blst_fp12` storage does not encode that domain, and this method
    /// performs no independent input validation. Callers should normally pass
    /// values returned by [`Self::miller_loop`] or [`Self::miller_loop_n`].
    pub fn finalverify(a: &Self, b: &Self) -> bool {
        unsafe { blst_fp12_finalverify(a, b) }
    }

    /// Serializes this Fp12 value into its big-endian 576-byte encoding.
    pub fn to_bendian(&self) -> [u8; 48 * 12] {
        let mut out = MaybeUninit::<[u8; 48 * 12]>::uninit();
        unsafe {
            blst_bendian_from_fp12(out.as_mut_ptr() as *mut u8, self);
            out.assume_init()
        }
    }
}

impl blst_scalar {
    /// Hashes arbitrary bytes into a scalar when the expanded output is valid.
    pub fn hash_to(msg: &[u8], dst: &[u8]) -> Option<Self> {
        unsafe {
            let mut out = <Self>::default();
            let mut elem = [0u8; 48];
            blst_expand_message_xmd(
                elem.as_mut_ptr(),
                elem.len(),
                msg.as_ptr(),
                msg.len(),
                dst.as_ptr(),
                dst.len(),
            );
            if blst_scalar_from_be_bytes(&mut out, elem.as_ptr(), elem.len()) {
                Some(out)
            } else {
                None
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        blst_fp, blst_fr, blst_p1, blst_p1_affine, blst_p2, blst_p2_affine,
    };

    use super::*;

    fn scalar_bytes(value: u64) -> [u8; 32] {
        let mut out = [0u8; 32];
        for (i, byte) in value.to_le_bytes().iter().enumerate() {
            out[i] = *byte;
        }
        out
    }

    fn fr_to_words(value: &blst_fr) -> [u64; 4] {
        let mut out = [0u64; 4];
        unsafe {
            blst_uint64_from_fr(out.as_mut_ptr(), value);
        }
        out
    }

    fn fp_to_words(value: &blst_fp) -> [u64; 6] {
        let mut out = [0u64; 6];
        unsafe {
            blst_uint64_from_fp(out.as_mut_ptr(), value);
        }
        out
    }

    #[test]
    fn c_api_primitives_accept_nonzero_points() {
        unsafe {
            let p_scalar = scalar_bytes(3);
            let q_scalar = scalar_bytes(11);

            let mut p = blst_p1::default();
            let mut q = blst_p2::default();
            blst_p1_mult(&mut p, blst_p1_generator(), p_scalar.as_ptr(), 256);
            blst_p2_mult(&mut q, blst_p2_generator(), q_scalar.as_ptr(), 256);

            assert!(blst_p1_in_g1(&p));
            assert!(blst_p2_in_g2(&q));

            let mut p_affine = blst_p1_affine::default();
            let mut q_affine = blst_p2_affine::default();
            blst_p1_to_affine(&mut p_affine, &p);
            blst_p2_to_affine(&mut q_affine, &q);

            let gt = blst_fp12::miller_loop(&q_affine, &p_affine).final_exp();
            assert!(gt.in_group());
            assert_eq!(gt * blst_fp12::default(), gt);
            assert!(blst_fp12::finalverify(&gt, &gt));
        }
    }

    #[test]
    fn constant_time_and_variable_time_group_multiplication_agree() {
        let scalar = scalar_bytes(0x0123_4567_89ab_cdef);
        unsafe {
            let g1 = *blst_p1_generator();
            assert_eq!(g1.mult(&scalar, 256), g1.mult_vartime(&scalar, 256));

            let g2 = *blst_p2_generator();
            assert_eq!(g2.mult(&scalar, 256), g2.mult_vartime(&scalar, 256));

            let gt = *blst_gt_generator();
            assert_eq!(
                gt.gt_mult(&scalar, 256),
                gt.gt_mult_vartime(&scalar, 256)
            );
        }
    }

    #[test]
    fn multi_miller_loop_matches_ordered_product() {
        let scalars = [2u64, 3, 5, 7];
        let mut ps = [blst_p1_affine::default(); 4];
        let mut qs = [blst_p2_affine::default(); 4];

        for (i, scalar) in scalars.iter().copied().enumerate() {
            let scalar = scalar_bytes(scalar);
            let mut p = blst_p1::default();
            let mut q = blst_p2::default();
            unsafe {
                blst_p1_mult(&mut p, blst_p1_generator(), scalar.as_ptr(), 256);
                blst_p2_mult(&mut q, blst_p2_generator(), scalar.as_ptr(), 256);
                blst_p1_to_affine(&mut ps[i], &p);
                blst_p2_to_affine(&mut qs[i], &q);
            }
        }

        let mut expected = blst_fp12::default();
        for i in 0..ps.len() {
            expected *= blst_fp12::miller_loop(&qs[i], &ps[i]);
        }

        assert_eq!(blst_fp12::miller_loop_n(&qs, &ps), expected);
        assert_eq!(blst_fp12::miller_loop_n(&[], &[]), blst_fp12::default());
    }

    #[test]
    fn group_c_api_extensions_link_and_agree() {
        unsafe {
            let two = scalar_bytes(2);

            let g1 = *blst_p1_generator();
            let mut g1_twice = blst_p1::default();
            let mut g1_difference = blst_p1::default();
            blst_p1_mult(&mut g1_twice, &g1, two.as_ptr(), 256);
            blst_p1_sub(&mut g1_difference, &g1_twice, &g1);
            assert!(blst_p1_is_equal(&g1_difference, &g1));

            let mut g1_affine = blst_p1_affine::default();
            blst_p1_to_affine(&mut g1_affine, &g1);
            blst_p1_affine_cneg(&mut g1_affine, true);
            blst_p1_sub_affine(&mut g1_difference, &g1, &g1_affine);
            assert!(blst_p1_is_equal(&g1_difference, &g1_twice));

            let g2 = *blst_p2_generator();
            let mut g2_twice = blst_p2::default();
            let mut g2_difference = blst_p2::default();
            blst_p2_mult(&mut g2_twice, &g2, two.as_ptr(), 256);
            blst_p2_sub(&mut g2_difference, &g2_twice, &g2);
            assert!(blst_p2_is_equal(&g2_difference, &g2));

            let mut g2_affine = blst_p2_affine::default();
            blst_p2_to_affine(&mut g2_affine, &g2);
            blst_p2_affine_cneg(&mut g2_affine, true);
            blst_p2_sub_affine(&mut g2_difference, &g2, &g2_affine);
            assert!(blst_p2_is_equal(&g2_difference, &g2_twice));

            let mut g1_pairing_affine = blst_p1_affine::default();
            let mut g2_pairing_affine = blst_p2_affine::default();
            blst_p1_to_affine(&mut g1_pairing_affine, &g1);
            blst_p2_to_affine(&mut g2_pairing_affine, &g2);
            let gt =
                blst_fp12::miller_loop(&g2_pairing_affine, &g1_pairing_affine)
                    .final_exp();
            assert!(blst_gt_is_equal(blst_gt_generator(), &gt));
            assert!(blst_gt_is_identity(blst_gt_identity()));

            let mut gt_added = blst_fp12::default();
            let mut gt_doubled = blst_fp12::default();
            let mut gt_negated = blst_fp12::default();
            let mut gt_cancelled = blst_fp12::default();
            let mut gt_difference = blst_fp12::default();
            blst_gt_add(&mut gt_added, &gt, &gt);
            blst_gt_double(&mut gt_doubled, &gt);
            blst_gt_neg(&mut gt_negated, &gt);
            blst_gt_add(&mut gt_cancelled, &gt, &gt_negated);
            blst_gt_sub(&mut gt_difference, &gt_doubled, &gt);
            assert!(blst_gt_is_equal(&gt_added, &gt_doubled));
            assert!(blst_gt_is_identity(&gt_cancelled));
            assert!(blst_gt_is_equal(&gt_difference, &gt));

            let mut gt_twice = blst_fp12::default();
            blst_gt_mult(&mut gt_twice, &gt, two.as_ptr(), 256);
            assert_eq!(gt_twice, gt * gt);
        }
    }

    #[test]
    fn field_c_api_primitives_round_trip() {
        unsafe {
            let five_words = [5u64, 0, 0, 0];
            let one_fr_words = [1u64, 0, 0, 0];
            let mut five = blst_fr::default();
            let mut one_fr = blst_fr::default();
            blst_fr_from_uint64(&mut five, five_words.as_ptr());
            blst_fr_from_uint64(&mut one_fr, one_fr_words.as_ptr());

            let mut inverse = blst_fr::default();
            let mut product = blst_fr::default();
            let mut power = blst_fr::default();
            let mut power_zero = blst_fr::default();
            blst_fr_inverse(&mut inverse, &five);
            blst_fr_mul(&mut product, &five, &inverse);
            let exponent_three = [3u8];
            blst_fr_pow(&mut power, &five, exponent_three.as_ptr(), 2);
            blst_fr_pow(&mut power_zero, &five, core::ptr::null(), 0);
            assert_eq!(fr_to_words(&product), one_fr_words);
            assert_eq!(fr_to_words(&power), [125, 0, 0, 0]);
            assert_eq!(fr_to_words(&power_zero), one_fr_words);

            let fields = [one_fr, five];
            let mut bulk = [blst_scalar::default(), blst_scalar::default()];
            let mut individual =
                [blst_scalar::default(), blst_scalar::default()];
            assert!(blst_scalars_from_fr(
                bulk.as_mut_ptr(),
                fields.as_ptr(),
                fields.len()
            ));
            blst_scalar_from_fr(&mut individual[0], &fields[0]);
            blst_scalar_from_fr(&mut individual[1], &fields[1]);
            assert_eq!(bulk, individual);

            let three_words = [3u64, 0, 0, 0, 0, 0];
            let one_fp_words = [1u64, 0, 0, 0, 0, 0];
            let mut three = blst_fp::default();
            let mut one_fp = blst_fp::default();
            blst_fp_from_uint64(&mut three, three_words.as_ptr());
            blst_fp_from_uint64(&mut one_fp, one_fp_words.as_ptr());

            let mut fp_inverse = blst_fp::default();
            let mut fp_product = blst_fp::default();
            blst_fp_inverse(&mut fp_inverse, &three);
            blst_fp_mul(&mut fp_product, &three, &fp_inverse);
            assert_eq!(fp_to_words(&fp_product), one_fp_words);

            let mut square = blst_fp::default();
            let mut sqrt = blst_fp::default();
            blst_fp_sqr(&mut square, &three);
            assert!(blst_fp_is_square(&square));
            assert!(blst_fp_sqrt(&mut sqrt, &square));
            blst_fp_sqr(&mut fp_product, &sqrt);
            assert_eq!(fp_to_words(&fp_product), fp_to_words(&square));
        }
    }

    #[test]
    fn variable_length_fp_decoders_accept_empty_and_link() {
        unsafe {
            let seven_words = [7u64, 0, 0, 0, 0, 0];
            let mut empty_le = blst_fp::default();
            let mut empty_be = blst_fp::default();
            blst_fp_from_uint64(&mut empty_le, seven_words.as_ptr());
            blst_fp_from_uint64(&mut empty_be, seven_words.as_ptr());
            blst_fp_from_le_bytes(&mut empty_le, core::ptr::null(), 0);
            blst_fp_from_be_bytes(&mut empty_be, core::ptr::null(), 0);
            assert_eq!(empty_le, blst_fp::default());
            assert_eq!(empty_be, blst_fp::default());

            let input = core::array::from_fn::<_, 49, _>(|i| {
                (i.wrapping_mul(131).wrapping_add(0x5b) & 0xff) as u8
            });
            let mut reversed = input;
            reversed.reverse();
            let mut little = blst_fp::default();
            let mut big_reversed = blst_fp::default();
            blst_fp_from_le_bytes(&mut little, input.as_ptr(), input.len());
            blst_fp_from_be_bytes(
                &mut big_reversed,
                reversed.as_ptr(),
                reversed.len(),
            );
            assert_eq!(little, big_reversed);
        }
    }
}
