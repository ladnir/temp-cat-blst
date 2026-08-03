// Copyright Supranational LLC
// Licensed under the Apache License, Version 2.0, see LICENSE for details.
// SPDX-License-Identifier: Apache-2.0

//! Expert AVX ABI bindings.

use crate::*;
use alloc::vec::Vec;
use core::mem::{size_of, MaybeUninit};

#[repr(C, align(64))]
#[derive(Copy, Clone)]
pub struct blst_avx_p1 {
    pub b: [blst_byte; 384],
}

#[repr(C, align(64))]
#[derive(Copy, Clone)]
pub struct blst_avx_p1_affine {
    pub b: [blst_byte; 256],
}

#[repr(C, align(64))]
#[derive(Copy, Clone)]
pub struct blst_avx_p2 {
    pub b: [blst_byte; 768],
}

#[repr(C, align(64))]
#[derive(Copy, Clone)]
pub struct blst_avx_p2_affine {
    pub b: [blst_byte; 512],
}

#[repr(C, align(64))]
#[derive(Copy, Clone)]
pub struct blst_avx_fp12 {
    pub b: [blst_byte; 1536],
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
    pub fn blst_avx_p1_generator() -> *const blst_avx_p1;

    pub fn blst_avx_p2_add(
        out: *mut blst_avx_p2,
        a: *const blst_avx_p2,
        b: *const blst_avx_p2,
    ) -> bool;
    pub fn blst_avx_p2_cneg(p: *mut blst_avx_p2, cbit: bool) -> bool;
    pub fn blst_avx_p2_batch_subgroup_check(
        points: *const *const blst_avx_p2,
        count: usize,
        coefficients: *const blst_byte,
        width: usize,
        shift: *const blst_avx_p2,
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
    pub fn blst_avx_fp12_mult_pippenger_scratch_sizeof(npoints: usize)
        -> usize;
    pub fn blst_avx_fp12_mult_pippenger_bucket_used_scratch_sizeof(
        npoints: usize,
    ) -> usize;
    pub fn blst_avx_fp12_mult_pippenger(
        out: *mut blst_avx_fp12,
        points: *const blst_avx_fp12,
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
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AvxMsmError {
    AvxUnavailable,
    ScalarLengthMismatch,
    ScratchSizeOverflow,
    ConversionFailed,
    BackendFailure,
}

fn avx_fp12_storage_len(bytes: usize) -> usize {
    let element_size = size_of::<blst_avx_fp12>();
    let whole_elements = bytes / element_size;
    whole_elements + usize::from(whole_elements * element_size != bytes)
}

/// Computes a variable-time AVX GLV4 Pippenger MSM over GT/Fp12 elements.
///
/// `scalars` must contain exactly one 256-bit little-endian scalar per point.
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

    if points.is_empty() {
        return Ok(blst_fp12::default());
    }

    let mut avx_points = Vec::with_capacity(points.len());
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
        unsafe { blst_avx_fp12_mult_pippenger_scratch_sizeof(points.len()) };
    let bucket_used_bytes = unsafe {
        blst_avx_fp12_mult_pippenger_bucket_used_scratch_sizeof(points.len())
    };
    if scratch_bytes == 0 || bucket_used_bytes == 0 {
        return Err(AvxMsmError::ScratchSizeOverflow);
    }

    let scratch_len = avx_fp12_storage_len(scratch_bytes);
    // The foreign call owns these temporary writes. Keep the Rust lengths at
    // zero so no uninitialized element is ever claimed as an initialized value.
    let mut scratch =
        Vec::<MaybeUninit<blst_avx_fp12>>::with_capacity(scratch_len);
    let mut bucket_used =
        Vec::<MaybeUninit<u8>>::with_capacity(bucket_used_bytes);

    let mut avx_out = MaybeUninit::<blst_avx_fp12>::uninit();
    let ok = unsafe {
        blst_avx_fp12_mult_pippenger(
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
        if unsafe { !blst_avx_enabled() } {
            return;
        }

        assert_eq!(core::mem::size_of::<blst_avx_p1>(), unsafe {
            blst_avx_p1_sizeof()
        });
        assert!(
            core::mem::align_of::<blst_avx_p1>()
                >= unsafe { blst_avx_p1_alignof() }
        );
        assert_eq!(core::mem::size_of::<blst_avx_p2>(), unsafe {
            blst_avx_p2_sizeof()
        });
        assert!(
            core::mem::align_of::<blst_avx_p2>()
                >= unsafe { blst_avx_p2_alignof() }
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

        let points = [blst_fp12::default(), blst_fp12::default()];
        let mut scalars = [0u8; 64];
        scalars[0] = 3;
        scalars[32] = 5;

        assert_eq!(gt_msm(&points, &scalars).unwrap(), blst_fp12::default());
    }
}
