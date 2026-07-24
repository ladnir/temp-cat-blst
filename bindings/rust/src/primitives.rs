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

#[cfg(feature = "std")]
use std::sync::{
    atomic::{AtomicUsize, Ordering},
    mpsc::channel,
    Arc,
};

impl PartialEq for blst_p1 {
    fn eq(&self, other: &Self) -> bool {
        unsafe { blst_p1_is_equal(self, other) }
    }
}

impl PartialEq for blst_p1_affine {
    fn eq(&self, other: &Self) -> bool {
        unsafe { blst_p1_affine_is_equal(self, other) }
    }
}

impl PartialEq for blst_p2 {
    fn eq(&self, other: &Self) -> bool {
        unsafe { blst_p2_is_equal(self, other) }
    }
}

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
    /// Computes one Miller loop from affine inputs.
    pub fn miller_loop(q: &blst_p2_affine, p: &blst_p1_affine) -> Self {
        let mut out = MaybeUninit::<blst_fp12>::uninit();
        unsafe {
            blst_miller_loop(out.as_mut_ptr(), q, p);
            out.assume_init()
        }
    }

    #[cfg(not(feature = "std"))]
    /// Computes the product of Miller loops over paired affine slices.
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

    #[cfg(feature = "std")]
    /// Computes the product of Miller loops over paired affine slices.
    pub fn miller_loop_n(q: &[blst_p2_affine], p: &[blst_p1_affine]) -> Self {
        let n_elems = q.len();
        if n_elems != p.len() {
            panic!("inputs' lengths mismatch");
        }
        if n_elems == 0 {
            return Self::default();
        }

        let pool = crate::mt::da_pool();

        let mut n_workers = pool.max_count();
        if n_workers == 1 {
            let qs: [*const _; 2] = [&q[0], ptr::null()];
            let ps: [*const _; 2] = [&p[0], ptr::null()];
            let mut out = MaybeUninit::<blst_fp12>::uninit();
            unsafe {
                blst_miller_loop_n(out.as_mut_ptr(), &qs[0], &ps[0], n_elems);
                return out.assume_init();
            }
        }

        let (tx, rx) = channel();
        let counter = Arc::new(AtomicUsize::new(0));

        let stride = core::cmp::min((n_elems + n_workers - 1) / n_workers, 16);
        n_workers = core::cmp::min((n_elems + stride - 1) / stride, n_workers);
        for _ in 0..n_workers {
            let tx = tx.clone();
            let counter = counter.clone();

            pool.joined_execute(move || {
                let mut acc = blst_fp12::default();
                let mut tmp = MaybeUninit::<blst_fp12>::uninit();
                let mut qs: [*const _; 2] = [ptr::null(), ptr::null()];
                let mut ps: [*const _; 2] = [ptr::null(), ptr::null()];

                loop {
                    let work = counter.fetch_add(stride, Ordering::Relaxed);
                    if work >= n_elems {
                        break;
                    }
                    let n = core::cmp::min(n_elems - work, stride);
                    qs[0] = &q[work];
                    ps[0] = &p[work];
                    unsafe {
                        blst_miller_loop_n(tmp.as_mut_ptr(), &qs[0], &ps[0], n);
                        acc *= tmp.assume_init();
                    }
                }

                tx.send(acc).expect("disaster");
            });
        }

        let mut acc = rx.recv().unwrap();
        for _ in 1..n_workers {
            acc *= rx.recv().unwrap();
        }

        acc
    }

    /// Applies the final exponentiation to this Fp12 value.
    pub fn final_exp(&self) -> Self {
        let mut out = MaybeUninit::<blst_fp12>::uninit();
        unsafe {
            blst_final_exp(out.as_mut_ptr(), self);
            out.assume_init()
        }
    }

    /// Returns `true` when this value lies in the pairing target subgroup.
    pub fn in_group(&self) -> bool {
        unsafe { blst_fp12_in_group(self) }
    }

    /// Returns `true` when `a` and `b` pass the final verification check.
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
}
