// Copyright Supranational LLC
// Licensed under the Apache License, Version 2.0, see LICENSE for details.
// SPDX-License-Identifier: Apache-2.0

//! Expert AVX ABI bindings.

use crate::*;

#[repr(C, align(64))]
#[derive(Copy, Clone)]
pub struct blst_avx_p1 {
    pub b: [byte; 384],
}

#[repr(C, align(64))]
#[derive(Copy, Clone)]
pub struct blst_avx_p1_affine {
    pub b: [byte; 256],
}

#[repr(C, align(64))]
#[derive(Copy, Clone)]
pub struct blst_avx_p2 {
    pub b: [byte; 768],
}

#[repr(C, align(64))]
#[derive(Copy, Clone)]
pub struct blst_avx_p2_affine {
    pub b: [byte; 512],
}

#[repr(C, align(64))]
#[derive(Copy, Clone)]
pub struct blst_avx_fp12 {
    pub b: [byte; 1536],
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
        unsafe { blst_avx_p1_is_equal(self, other) }
    }
}

impl Eq for blst_avx_p1 {}

impl PartialEq for blst_avx_p1_affine {
    fn eq(&self, other: &Self) -> bool {
        self.b == other.b
    }
}

impl Eq for blst_avx_p1_affine {}

impl PartialEq for blst_avx_p2 {
    fn eq(&self, other: &Self) -> bool {
        unsafe { blst_avx_p2_is_equal(self, other) }
    }
}

impl Eq for blst_avx_p2 {}

impl PartialEq for blst_avx_p2_affine {
    fn eq(&self, other: &Self) -> bool {
        self.b == other.b
    }
}

impl Eq for blst_avx_p2_affine {}

impl PartialEq for blst_avx_fp12 {
    fn eq(&self, other: &Self) -> bool {
        unsafe { blst_avx_fp12_is_equal(self, other) }
    }
}

impl Eq for blst_avx_fp12 {}

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
    pub fn blst_avx_p1_mult(
        out: *mut blst_avx_p1,
        p: *const blst_avx_p1,
        scalar: *const byte,
        nbits: usize,
    ) -> bool;
    pub fn blst_avx_p1_to_affine(
        out: *mut blst_avx_p1_affine,
        input: *const blst_avx_p1,
    ) -> bool;
    pub fn blst_avx_p1_in_g1(p: *const blst_avx_p1) -> bool;
    pub fn blst_avx_p1_is_inf(p: *const blst_avx_p1) -> bool;
    pub fn blst_avx_p1_is_equal(a: *const blst_avx_p1, b: *const blst_avx_p1) -> bool;
    pub fn blst_avx_p1_compress(out: *mut byte, input: *const blst_avx_p1) -> bool;
    pub fn blst_avx_p1_serialize(out: *mut byte, input: *const blst_avx_p1) -> bool;
    pub fn blst_avx_p1_deserialize(out: *mut blst_avx_p1, input: *const byte) -> BLST_ERROR;
    pub fn blst_avx_p1_uncompress(out: *mut blst_avx_p1, input: *const byte) -> BLST_ERROR;
    pub fn blst_avx_p1s_mult_pippenger_scratch_sizeof(npoints: usize) -> usize;
    pub fn blst_avx_p1s_mult_pippenger(
        out: *mut blst_avx_p1,
        points: *const blst_avx_p1,
        npoints: usize,
        scalars: *const byte,
        nbits: usize,
        scratch: *mut limb_t,
    ) -> bool;
    pub fn blst_avx_p1_generator() -> *const blst_avx_p1;

    pub fn blst_avx_p2_add(
        out: *mut blst_avx_p2,
        a: *const blst_avx_p2,
        b: *const blst_avx_p2,
    ) -> bool;
    pub fn blst_avx_p2_mult(
        out: *mut blst_avx_p2,
        p: *const blst_avx_p2,
        scalar: *const byte,
        nbits: usize,
    ) -> bool;
    pub fn blst_avx_p2_to_affine(
        out: *mut blst_avx_p2_affine,
        input: *const blst_avx_p2,
    ) -> bool;
    pub fn blst_avx_p2_in_g2(p: *const blst_avx_p2) -> bool;
    pub fn blst_avx_p2_is_inf(p: *const blst_avx_p2) -> bool;
    pub fn blst_avx_p2_is_equal(a: *const blst_avx_p2, b: *const blst_avx_p2) -> bool;
    pub fn blst_avx_p2_compress(out: *mut byte, input: *const blst_avx_p2) -> bool;
    pub fn blst_avx_p2_uncompress(out: *mut blst_avx_p2, input: *const byte) -> BLST_ERROR;
    pub fn blst_avx_p2s_mult_pippenger_scratch_sizeof(npoints: usize) -> usize;
    pub fn blst_avx_p2s_mult_pippenger(
        out: *mut blst_avx_p2,
        points: *const blst_avx_p2,
        npoints: usize,
        scalars: *const byte,
        nbits: usize,
        scratch: *mut limb_t,
    ) -> bool;
    pub fn blst_avx_p2_generator() -> *const blst_avx_p2;

    pub fn blst_avx_fp12_is_equal(
        a: *const blst_avx_fp12,
        b: *const blst_avx_fp12,
    ) -> bool;
    pub fn blst_avx_pairing(
        out: *mut blst_avx_fp12,
        g1: *const blst_avx_p1,
        g2: *const blst_avx_p2,
    ) -> bool;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rust_layout_matches_avx_abi() {
        if unsafe { !blst_avx_enabled() } {
            return;
        }

        assert_eq!(core::mem::size_of::<blst_avx_p1>(), unsafe {
            blst_avx_p1_sizeof()
        });
        assert_eq!(core::mem::align_of::<blst_avx_p1>(), unsafe {
            blst_avx_p1_alignof()
        });
        assert_eq!(core::mem::size_of::<blst_avx_p2>(), unsafe {
            blst_avx_p2_sizeof()
        });
        assert_eq!(core::mem::align_of::<blst_avx_p2>(), unsafe {
            blst_avx_p2_alignof()
        });
        assert_eq!(core::mem::size_of::<blst_avx_fp12>(), unsafe {
            blst_avx_fp12_sizeof()
        });
        assert_eq!(core::mem::align_of::<blst_avx_fp12>(), unsafe {
            blst_avx_fp12_alignof()
        });
    }
}
