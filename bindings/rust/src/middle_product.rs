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
//! - recover the constant-term contribution with `update`
//! - combine GT outputs with `combine`
//!
//! For `a.len() == n` and a precomputed `b.len() == 2*n`, output `k` is
//! `product_i pairing(a[i], b[n + k - i])`. All G1/G2 values must already be
//! reduced, on-curve members of their prime-order subgroups; infinity is the
//! identity. The native implementation does not validate those properties.
//! The returned `blst_fp12` values semantically belong to GT.

use crate::*;

#[repr(C)]
pub struct blst_middle_product_precomputed_b {
    _private: [u8; 0],
}

#[repr(C)]
pub struct blst_middle_product_precomputed_b_scalar {
    _private: [u8; 0],
}

#[repr(C)]
pub struct blst_middle_product_precomputed_b_avx {
    _private: [u8; 0],
}

#[repr(C)]
pub struct blst_middle_product_workspace_scalar {
    _private: [u8; 0],
}

#[repr(C)]
pub struct blst_middle_product_workspace_avx {
    _private: [u8; 0],
}

extern "C" {
    fn blst_middle_product_precompute_b_scalar(
        b: *const blst_p2,
        n: usize,
    ) -> *mut blst_middle_product_precomputed_b_scalar;
    fn blst_middle_product_scalar(
        a: *const blst_p1,
        n: usize,
        precomputed_b: *const blst_middle_product_precomputed_b_scalar,
        out: *mut blst_fp12,
    ) -> bool;
    fn blst_middle_product_workspace_scalar_new(
        precomputed_b: *const blst_middle_product_precomputed_b_scalar,
    ) -> *mut blst_middle_product_workspace_scalar;
    fn blst_middle_product_scalar_with_workspace(
        a: *const blst_p1,
        n: usize,
        precomputed_b: *const blst_middle_product_precomputed_b_scalar,
        workspace: *mut blst_middle_product_workspace_scalar,
        out: *mut blst_fp12,
    ) -> bool;
    fn blst_middle_product_workspace_scalar_free(
        workspace: *mut blst_middle_product_workspace_scalar,
    );
    fn blst_middle_product_update_scalar(
        a0: *const blst_p1,
        precomputed_b: *const blst_middle_product_precomputed_b_scalar,
        out: *mut blst_fp12,
    ) -> bool;
    fn blst_middle_product_precomputed_b_scalar_free(
        precomputed_b: *mut blst_middle_product_precomputed_b_scalar,
    );
    fn blst_middle_product_precompute_b_avx(
        b: *const blst_p2,
        n: usize,
    ) -> *mut blst_middle_product_precomputed_b_avx;
    fn blst_middle_product_avx(
        a: *const blst_p1,
        n: usize,
        precomputed_b: *const blst_middle_product_precomputed_b_avx,
        out: *mut blst_fp12,
    ) -> bool;
    fn blst_middle_product_workspace_avx_new(
        precomputed_b: *const blst_middle_product_precomputed_b_avx,
    ) -> *mut blst_middle_product_workspace_avx;
    fn blst_middle_product_avx_with_workspace(
        a: *const blst_p1,
        n: usize,
        precomputed_b: *const blst_middle_product_precomputed_b_avx,
        workspace: *mut blst_middle_product_workspace_avx,
        out: *mut blst_fp12,
    ) -> bool;
    fn blst_middle_product_workspace_avx_free(
        workspace: *mut blst_middle_product_workspace_avx,
    );
    fn blst_middle_product_update_avx(
        a0: *const blst_p1,
        precomputed_b: *const blst_middle_product_precomputed_b_avx,
        out: *mut blst_fp12,
    ) -> bool;
    fn blst_middle_product_precomputed_b_avx_free(
        precomputed_b: *mut blst_middle_product_precomputed_b_avx,
    );
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
    fn blst_middle_product_update(
        a0: *const blst_p1,
        precomputed_b: *const blst_middle_product_precomputed_b,
        out: *mut blst_fp12,
    ) -> bool;
    fn blst_middle_product_combine(
        mid: *const blst_fp12,
        cnst: *const blst_fp12,
        n: usize,
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
    WorkspaceCreationFailed,
    BackendFailure,
}

/// Opaque owner for a native precomputed-`b` handle.
///
/// The underlying native object is immutable after construction and freed
/// automatically on drop. The wrapper does not expose concurrent mutation or
/// manual lifetime management.
#[derive(Debug)]
pub struct MiddleProductPrecomputedB {
    raw: *mut blst_middle_product_precomputed_b,
    n: usize,
}

impl MiddleProductPrecomputedB {
    /// Precomputes `b`, where `b.len()` must be exactly `2 * n` and every
    /// element must already be a reduced G2 subgroup point.
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

/// Opaque precomputation whose evaluations always use the scalar backend.
#[derive(Debug)]
pub struct ScalarMiddleProductPrecomputedB {
    raw: *mut blst_middle_product_precomputed_b_scalar,
    n: usize,
}

impl ScalarMiddleProductPrecomputedB {
    pub fn new(b: &[blst_p2]) -> Result<Self, MiddleProductError> {
        validate_precompute_input(b)?;
        let n = b.len() / 2;
        let raw =
            unsafe { blst_middle_product_precompute_b_scalar(b.as_ptr(), n) };
        if raw.is_null() {
            return Err(MiddleProductError::PrecomputeFailed);
        }
        Ok(Self { raw, n })
    }

    pub fn len(&self) -> usize {
        self.n
    }

    pub fn is_empty(&self) -> bool {
        self.n == 0
    }
}

impl Drop for ScalarMiddleProductPrecomputedB {
    fn drop(&mut self) {
        if !self.raw.is_null() {
            unsafe { blst_middle_product_precomputed_b_scalar_free(self.raw) };
        }
    }
}

/// Caller-owned reusable scalar FFT workspace.
#[derive(Debug)]
pub struct ScalarMiddleProductWorkspace {
    raw: *mut blst_middle_product_workspace_scalar,
    n: usize,
}

impl ScalarMiddleProductWorkspace {
    pub fn new(
        precomputed_b: &ScalarMiddleProductPrecomputedB,
    ) -> Result<Self, MiddleProductError> {
        let raw = unsafe {
            blst_middle_product_workspace_scalar_new(precomputed_b.raw)
        };
        if raw.is_null() {
            return Err(MiddleProductError::WorkspaceCreationFailed);
        }
        Ok(Self {
            raw,
            n: precomputed_b.n,
        })
    }
}

impl Drop for ScalarMiddleProductWorkspace {
    fn drop(&mut self) {
        if !self.raw.is_null() {
            unsafe { blst_middle_product_workspace_scalar_free(self.raw) };
        }
    }
}

/// Opaque precomputation whose evaluations always use the AVX backend.
///
/// Construction returns `PrecomputeFailed` when the linked CatBLST build does
/// not contain AVX support.
#[derive(Debug)]
pub struct AvxMiddleProductPrecomputedB {
    raw: *mut blst_middle_product_precomputed_b_avx,
    n: usize,
}

impl AvxMiddleProductPrecomputedB {
    pub fn new(b: &[blst_p2]) -> Result<Self, MiddleProductError> {
        validate_precompute_input(b)?;
        let n = b.len() / 2;
        let raw =
            unsafe { blst_middle_product_precompute_b_avx(b.as_ptr(), n) };
        if raw.is_null() {
            return Err(MiddleProductError::PrecomputeFailed);
        }
        Ok(Self { raw, n })
    }

    pub fn len(&self) -> usize {
        self.n
    }

    pub fn is_empty(&self) -> bool {
        self.n == 0
    }
}

impl Drop for AvxMiddleProductPrecomputedB {
    fn drop(&mut self) {
        if !self.raw.is_null() {
            unsafe { blst_middle_product_precomputed_b_avx_free(self.raw) };
        }
    }
}

/// Caller-owned reusable AVX FFT and scalar-to-AVX conversion workspace.
#[derive(Debug)]
pub struct AvxMiddleProductWorkspace {
    raw: *mut blst_middle_product_workspace_avx,
    n: usize,
}

impl AvxMiddleProductWorkspace {
    pub fn new(
        precomputed_b: &AvxMiddleProductPrecomputedB,
    ) -> Result<Self, MiddleProductError> {
        let raw =
            unsafe { blst_middle_product_workspace_avx_new(precomputed_b.raw) };
        if raw.is_null() {
            return Err(MiddleProductError::WorkspaceCreationFailed);
        }
        Ok(Self {
            raw,
            n: precomputed_b.n,
        })
    }
}

impl Drop for AvxMiddleProductWorkspace {
    fn drop(&mut self) {
        if !self.raw.is_null() {
            unsafe { blst_middle_product_workspace_avx_free(self.raw) };
        }
    }
}

fn validate_precompute_input(b: &[blst_p2]) -> Result<(), MiddleProductError> {
    if b.is_empty() {
        return Err(MiddleProductError::EmptyInput);
    }
    if b.len() % 2 != 0 {
        return Err(MiddleProductError::InvalidBLength);
    }
    Ok(())
}

/// Precomputes the `b` operand for repeated middle-product evaluations.
///
/// Every value in `b` must be a reduced, on-curve G2 subgroup point. Infinity
/// is accepted as the identity. The native operation is variable-time and
/// allocates; allocation, unsupported size, and backend failures are reported
/// as `MiddleProductError::PrecomputeFailed`.
pub fn precompute_middle_product_b(
    b: &[blst_p2],
) -> Result<MiddleProductPrecomputedB, MiddleProductError> {
    validate_precompute_input(b)?;

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
///
/// Every value in `a` must be a reduced, on-curve G1 subgroup point. Infinity
/// is accepted as the identity. Results are GT elements stored as `blst_fp12`.
/// The native operation is variable-time and can allocate on each call. On
/// error, the first `precomputed_b.len()` elements of `out` are unchanged.
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

/// Computes a middle product through an explicitly scalar precompute.
pub fn middle_product_scalar(
    a: &[blst_p1],
    precomputed_b: &ScalarMiddleProductPrecomputedB,
    out: &mut [blst_fp12],
) -> Result<(), MiddleProductError> {
    validate_evaluation_buffers(a.len(), precomputed_b.n, out.len())?;
    let ok = unsafe {
        blst_middle_product_scalar(
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

/// Computes a scalar middle product without allocating in the native call.
///
/// `workspace` must have been constructed for a scalar precompute with the
/// same logical length. It may be reused by subsequent calls, but is mutably
/// borrowed so concurrent evaluations require distinct workspaces.
pub fn middle_product_scalar_with_workspace(
    a: &[blst_p1],
    precomputed_b: &ScalarMiddleProductPrecomputedB,
    workspace: &mut ScalarMiddleProductWorkspace,
    out: &mut [blst_fp12],
) -> Result<(), MiddleProductError> {
    validate_evaluation_buffers(a.len(), precomputed_b.n, out.len())?;
    if workspace.n != precomputed_b.n {
        return Err(MiddleProductError::LengthMismatch);
    }
    if unsafe {
        blst_middle_product_scalar_with_workspace(
            a.as_ptr(),
            precomputed_b.n,
            precomputed_b.raw,
            workspace.raw,
            out.as_mut_ptr(),
        )
    } {
        Ok(())
    } else {
        Err(MiddleProductError::BackendFailure)
    }
}

/// Computes a middle product through an explicitly AVX precompute.
pub fn middle_product_avx(
    a: &[blst_p1],
    precomputed_b: &AvxMiddleProductPrecomputedB,
    out: &mut [blst_fp12],
) -> Result<(), MiddleProductError> {
    validate_evaluation_buffers(a.len(), precomputed_b.n, out.len())?;
    let ok = unsafe {
        blst_middle_product_avx(
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

/// Computes an AVX middle product without allocating in the native call.
pub fn middle_product_avx_with_workspace(
    a: &[blst_p1],
    precomputed_b: &AvxMiddleProductPrecomputedB,
    workspace: &mut AvxMiddleProductWorkspace,
    out: &mut [blst_fp12],
) -> Result<(), MiddleProductError> {
    validate_evaluation_buffers(a.len(), precomputed_b.n, out.len())?;
    if workspace.n != precomputed_b.n {
        return Err(MiddleProductError::LengthMismatch);
    }
    if unsafe {
        blst_middle_product_avx_with_workspace(
            a.as_ptr(),
            precomputed_b.n,
            precomputed_b.raw,
            workspace.raw,
            out.as_mut_ptr(),
        )
    } {
        Ok(())
    } else {
        Err(MiddleProductError::BackendFailure)
    }
}

fn validate_evaluation_buffers(
    input_len: usize,
    expected: usize,
    output_len: usize,
) -> Result<(), MiddleProductError> {
    if input_len != expected {
        return Err(MiddleProductError::LengthMismatch);
    }
    if output_len < expected {
        return Err(MiddleProductError::OutputTooSmall);
    }
    Ok(())
}

/// Computes the constant-term contribution for one G1 point against a
/// precomputed `b` tail.
///
/// `a0` must be a reduced, on-curve G1 subgroup point; infinity is accepted as
/// the identity. Results are GT elements stored as `blst_fp12`. The native
/// operation is variable-time and can allocate on each call. On error, the
/// first `precomputed_b.len()` elements of `out` are unchanged.
pub fn update(
    a0: &blst_p1,
    precomputed_b: &MiddleProductPrecomputedB,
    out: &mut [blst_fp12],
) -> Result<(), MiddleProductError> {
    if out.len() < precomputed_b.n {
        return Err(MiddleProductError::OutputTooSmall);
    }

    let ok = unsafe {
        blst_middle_product_update(a0, precomputed_b.raw, out.as_mut_ptr())
    };
    if ok {
        Ok(())
    } else {
        Err(MiddleProductError::BackendFailure)
    }
}

/// Computes the constant-term contribution through the scalar backend.
pub fn update_scalar(
    a0: &blst_p1,
    precomputed_b: &ScalarMiddleProductPrecomputedB,
    out: &mut [blst_fp12],
) -> Result<(), MiddleProductError> {
    if out.len() < precomputed_b.n {
        return Err(MiddleProductError::OutputTooSmall);
    }
    if unsafe {
        blst_middle_product_update_scalar(
            a0,
            precomputed_b.raw,
            out.as_mut_ptr(),
        )
    } {
        Ok(())
    } else {
        Err(MiddleProductError::BackendFailure)
    }
}

/// Computes the constant-term contribution through the AVX backend.
pub fn update_avx(
    a0: &blst_p1,
    precomputed_b: &AvxMiddleProductPrecomputedB,
    out: &mut [blst_fp12],
) -> Result<(), MiddleProductError> {
    if out.len() < precomputed_b.n {
        return Err(MiddleProductError::OutputTooSmall);
    }
    if unsafe {
        blst_middle_product_update_avx(a0, precomputed_b.raw, out.as_mut_ptr())
    } {
        Ok(())
    } else {
        Err(MiddleProductError::BackendFailure)
    }
}

/// Combines two GT output vectors elementwise without allocating.
///
/// Every input must be a reduced GT element stored as `blst_fp12`; arbitrary
/// Fp12 values are outside this function's contract. This safe wrapper uses
/// separate Rust borrows, while the underlying C operation also permits exact
/// output aliasing and forbids partial overlap.
pub fn combine(
    mid: &[blst_fp12],
    cnst: &[blst_fp12],
    out: &mut [blst_fp12],
) -> Result<(), MiddleProductError> {
    if mid.len() != cnst.len() {
        return Err(MiddleProductError::LengthMismatch);
    }
    if out.len() < mid.len() {
        return Err(MiddleProductError::OutputTooSmall);
    }
    if mid.is_empty() {
        return Ok(());
    }

    let ok = unsafe {
        blst_middle_product_combine(
            mid.as_ptr(),
            cnst.as_ptr(),
            mid.len(),
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

    #[test]
    fn update_accepts_zero_boundary_points() {
        let b = [blst_p2::default(), blst_p2::default()];
        let precomputed = MiddleProductPrecomputedB::new(&b).unwrap();
        let a0 = blst_p1::default();
        let mut out = [zero_fp12()];

        update(&a0, &precomputed, &mut out).unwrap();
        assert_eq!(out[0], blst_fp12::default());
    }

    #[test]
    fn explicit_scalar_backend_accepts_zero_boundary_points() {
        let b = [blst_p2::default(), blst_p2::default()];
        let precomputed = ScalarMiddleProductPrecomputedB::new(&b).unwrap();
        let mut workspace =
            ScalarMiddleProductWorkspace::new(&precomputed).unwrap();
        let a = [blst_p1::default()];
        let mut middle = [zero_fp12()];
        let mut middle_reused = [zero_fp12()];
        let mut constant = [zero_fp12()];

        middle_product_scalar(&a, &precomputed, &mut middle).unwrap();
        middle_product_scalar_with_workspace(
            &a,
            &precomputed,
            &mut workspace,
            &mut middle_reused,
        )
        .unwrap();
        update_scalar(&a[0], &precomputed, &mut constant).unwrap();
        assert_eq!(middle_reused, middle);
        assert_eq!(middle[0], blst_fp12::default());
        assert_eq!(constant[0], blst_fp12::default());
    }

    #[test]
    fn explicit_avx_backend_matches_scalar_or_reports_unavailable() {
        let b = unsafe { [*blst_p2_generator(), *blst_p2_generator()] };
        let avx = AvxMiddleProductPrecomputedB::new(&b);
        if !unsafe { blst_avx_enabled() } {
            assert!(matches!(avx, Err(MiddleProductError::PrecomputeFailed)));
            return;
        }

        let scalar = ScalarMiddleProductPrecomputedB::new(&b).unwrap();
        let avx = avx.unwrap();
        let mut scalar_workspace =
            ScalarMiddleProductWorkspace::new(&scalar).unwrap();
        let mut avx_workspace = AvxMiddleProductWorkspace::new(&avx).unwrap();
        let a = unsafe { [*blst_p1_generator()] };
        let mut scalar_middle = [zero_fp12()];
        let mut avx_middle = [zero_fp12()];
        let mut scalar_reused = [zero_fp12()];
        let mut avx_reused = [zero_fp12()];
        let mut scalar_constant = [zero_fp12()];
        let mut avx_constant = [zero_fp12()];

        middle_product_scalar(&a, &scalar, &mut scalar_middle).unwrap();
        middle_product_avx(&a, &avx, &mut avx_middle).unwrap();
        middle_product_scalar_with_workspace(
            &a,
            &scalar,
            &mut scalar_workspace,
            &mut scalar_reused,
        )
        .unwrap();
        middle_product_avx_with_workspace(
            &a,
            &avx,
            &mut avx_workspace,
            &mut avx_reused,
        )
        .unwrap();
        update_scalar(&a[0], &scalar, &mut scalar_constant).unwrap();
        update_avx(&a[0], &avx, &mut avx_constant).unwrap();

        assert_eq!(avx_middle, scalar_middle);
        assert_eq!(scalar_reused, scalar_middle);
        assert_eq!(avx_reused, scalar_middle);
        assert_eq!(avx_constant, scalar_constant);
    }

    #[test]
    fn combine_rejects_length_mismatch() {
        let mid = [zero_fp12()];
        let cnst = [zero_fp12(), zero_fp12()];
        let mut out = [zero_fp12(); 2];

        assert_eq!(
            combine(&mid, &cnst, &mut out),
            Err(MiddleProductError::LengthMismatch)
        );
    }

    #[test]
    fn combine_accepts_empty_vectors() {
        let mid: [blst_fp12; 0] = [];
        let cnst: [blst_fp12; 0] = [];
        let mut out: [blst_fp12; 0] = [];

        combine(&mid, &cnst, &mut out).unwrap();
    }
}
