// Copyright Supranational LLC
// Licensed under the Apache License, Version 2.0, see LICENSE for details.
// SPDX-License-Identifier: Apache-2.0

//! Middle-product bindings.
//!
//! This module exposes the middle-product API as a small safe Rust layer over
//! the plain C ABI:
//! - precompute `b` once
//! - reuse the returned handle
//! - write outputs into caller-provided `&mut [blst_fp12]`

use crate::*;

#[repr(C)]
pub struct blst_middle_product_precomputed_b {
    _private: [u8; 0],
}

extern "C" {
    fn blst_middle_product_precompute_b(
        b: *const blst_p2,
        n: usize,
    ) -> *mut blst_middle_product_precomputed_b;
    fn blst_middle_product(
        a: *const blst_p1,
        n: usize,
        precomputed_b: *const blst_middle_product_precomputed_b,
        out: *mut blst_fp12,
    ) -> bool;
    fn blst_middle_product_precomputed_b_free(
        precomputed_b: *mut blst_middle_product_precomputed_b,
    );
}

/// Errors returned by the safe Rust middle-product wrapper.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MiddleProductError {
    EmptyInput,
    InvalidBLength,
    LengthMismatch,
    OutputTooSmall,
    PrecomputeFailed,
    BackendFailure,
}

/// Opaque owner for a native precomputed-`b` handle.
///
/// The underlying native object is freed automatically on drop.
#[derive(Debug)]
pub struct MiddleProductPrecomputedB {
    raw: *mut blst_middle_product_precomputed_b,
    n: usize,
}

impl MiddleProductPrecomputedB {
    /// Precomputes `b`, where `b.len()` must be exactly `2 * n`.
    pub fn new(b: &[blst_p2]) -> Result<Self, MiddleProductError> {
        precompute_middle_product_b(b)
    }

    /// Returns the logical output length associated with this precompute.
    pub fn len(&self) -> usize {
        self.n
    }

    /// Returns `true` when the precompute has zero logical length.
    pub fn is_empty(&self) -> bool {
        self.n == 0
    }
}

impl Drop for MiddleProductPrecomputedB {
    fn drop(&mut self) {
        if !self.raw.is_null() {
            unsafe { blst_middle_product_precomputed_b_free(self.raw) };
        }
    }
}

/// Precomputes the `b` operand for repeated middle-product evaluations.
pub fn precompute_middle_product_b(
    b: &[blst_p2],
) -> Result<MiddleProductPrecomputedB, MiddleProductError> {
    if b.is_empty() {
        return Err(MiddleProductError::EmptyInput);
    }
    if b.len() % 2 != 0 {
        return Err(MiddleProductError::InvalidBLength);
    }

    let n = b.len() / 2;
    let raw = unsafe { blst_middle_product_precompute_b(b.as_ptr(), n) };
    if raw.is_null() {
        return Err(MiddleProductError::PrecomputeFailed);
    }

    Ok(MiddleProductPrecomputedB { raw, n })
}

/// Computes the middle product for `a` using a precomputed `b`.
///
/// `out` must have length at least `a.len()`, and `a.len()` must match the
/// logical length of `precomputed_b`.
pub fn middle_product(
    a: &[blst_p1],
    precomputed_b: &MiddleProductPrecomputedB,
    out: &mut [blst_fp12],
) -> Result<(), MiddleProductError> {
    if a.len() != precomputed_b.n {
        return Err(MiddleProductError::LengthMismatch);
    }
    if out.len() < precomputed_b.n {
        return Err(MiddleProductError::OutputTooSmall);
    }

    let ok = unsafe {
        blst_middle_product(
            a.as_ptr(),
            precomputed_b.n,
            precomputed_b.raw,
            out.as_mut_ptr(),
        )
    };
    if ok {
        Ok(())
    } else {
        Err(MiddleProductError::BackendFailure)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn zero_fp12() -> blst_fp12 {
        blst_fp12 {
            fp6: [blst_fp6::default(), blst_fp6::default()],
        }
    }

    #[test]
    fn precompute_rejects_empty_input() {
        let b: [blst_p2; 0] = [];
        assert!(matches!(
            precompute_middle_product_b(&b),
            Err(MiddleProductError::EmptyInput)
        ));
    }

    #[test]
    fn precompute_rejects_odd_length() {
        let b = [blst_p2::default()];
        assert!(matches!(
            precompute_middle_product_b(&b),
            Err(MiddleProductError::InvalidBLength)
        ));
    }

    #[test]
    fn middle_product_rejects_length_mismatch() {
        let b = [blst_p2::default(), blst_p2::default()];
        let precomputed = precompute_middle_product_b(&b).unwrap();
        let a = [blst_p1::default(), blst_p1::default()];
        let mut out = [zero_fp12()];

        assert_eq!(
            middle_product(&a, &precomputed, &mut out),
            Err(MiddleProductError::LengthMismatch)
        );
    }

    #[test]
    fn middle_product_rejects_small_output_buffer() {
        let b = [blst_p2::default(), blst_p2::default()];
        let precomputed = precompute_middle_product_b(&b).unwrap();
        let a = [blst_p1::default()];
        let mut out: [blst_fp12; 0] = [];

        assert_eq!(
            middle_product(&a, &precomputed, &mut out),
            Err(MiddleProductError::OutputTooSmall)
        );
    }

    #[test]
    fn middle_product_accepts_zero_boundary_points() {
        let b = [blst_p2::default(), blst_p2::default()];
        let precomputed = MiddleProductPrecomputedB::new(&b).unwrap();
        assert_eq!(precomputed.len(), 1);
        assert!(!precomputed.is_empty());

        let a = [blst_p1::default()];
        let mut out = [zero_fp12()];

        middle_product(&a, &precomputed, &mut out).unwrap();
    }
}
