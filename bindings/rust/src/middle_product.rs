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
//! - keep explicit AVX inputs and outputs native with
//!   `middle_product_avx_native_with_workspace`
//! - split scalar or native-AVX evaluation into caller-scheduled disjoint ranges
//!
//! For `a.len() == n` and a precomputed `b.len() == 2*n`, output `k` is
//! `product_i pairing(a[i], b[n + k - i])`. All G1/G2 values must already be
//! reduced, on-curve members of their prime-order subgroups; infinity is the
//! identity. The native implementation does not validate those properties.
//! The returned `blst_fp12` values semantically belong to GT.

use core::ops::Range;

use crate::*;

/// Stateless phase identifier for the ranged middle-product interface.
///
/// Input loading and output extraction use their typed operations. Forward-J0,
/// forward-nonzero, and inverse-butterfly phases use a zero-based FFT stage;
/// every other phase requires stage zero.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(C)]
pub enum MiddleProductRangePhase {
    InputLoad = 0,
    ForwardBitReverse = 1,
    ForwardJ0 = 2,
    ForwardNonzero = 3,
    AffineConversion = 4,
    PointwisePairing = 5,
    InverseBitReverse = 6,
    InverseButterflies = 7,
    Normalization = 8,
    Extraction = 9,
}

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

#[repr(C)]
pub struct blst_middle_product_workspace_avx_native {
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
    fn blst_middle_product_scalar_range_work_count(
        precomputed_b: *const blst_middle_product_precomputed_b_scalar,
        phase: MiddleProductRangePhase,
        stage: usize,
        count: *mut usize,
    ) -> bool;
    fn blst_middle_product_scalar_input_range(
        a: *const blst_p1,
        n: usize,
        precomputed_b: *const blst_middle_product_precomputed_b_scalar,
        workspace: *mut blst_middle_product_workspace_scalar,
        begin: usize,
        end: usize,
    ) -> bool;
    fn blst_middle_product_scalar_execute_range(
        precomputed_b: *const blst_middle_product_precomputed_b_scalar,
        workspace: *mut blst_middle_product_workspace_scalar,
        phase: MiddleProductRangePhase,
        stage: usize,
        begin: usize,
        end: usize,
    ) -> bool;
    fn blst_middle_product_scalar_output_range(
        precomputed_b: *const blst_middle_product_precomputed_b_scalar,
        workspace: *mut blst_middle_product_workspace_scalar,
        out: *mut blst_fp12,
        n: usize,
        begin: usize,
        end: usize,
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
    fn blst_middle_product_precompute_b_avx_native(
        b: *const blst_avx_p2,
        n: usize,
    ) -> *mut blst_middle_product_precomputed_b_avx;
    fn blst_middle_product_workspace_avx_native_new(
        precomputed_b: *const blst_middle_product_precomputed_b_avx,
    ) -> *mut blst_middle_product_workspace_avx_native;
    fn blst_middle_product_avx_native_with_workspace(
        a: *const blst_avx_p1,
        n: usize,
        precomputed_b: *const blst_middle_product_precomputed_b_avx,
        workspace: *mut blst_middle_product_workspace_avx_native,
        out: *mut blst_avx_gt,
    ) -> bool;
    fn blst_middle_product_avx_native_range_work_count(
        precomputed_b: *const blst_middle_product_precomputed_b_avx,
        phase: MiddleProductRangePhase,
        stage: usize,
        count: *mut usize,
    ) -> bool;
    fn blst_middle_product_avx_native_input_range(
        a: *const blst_avx_p1,
        n: usize,
        precomputed_b: *const blst_middle_product_precomputed_b_avx,
        workspace: *mut blst_middle_product_workspace_avx_native,
        begin: usize,
        end: usize,
    ) -> bool;
    fn blst_middle_product_avx_native_execute_range(
        precomputed_b: *const blst_middle_product_precomputed_b_avx,
        workspace: *mut blst_middle_product_workspace_avx_native,
        phase: MiddleProductRangePhase,
        stage: usize,
        begin: usize,
        end: usize,
    ) -> bool;
    fn blst_middle_product_avx_native_output_range(
        precomputed_b: *const blst_middle_product_precomputed_b_avx,
        workspace: *mut blst_middle_product_workspace_avx_native,
        out: *mut blst_avx_gt,
        n: usize,
        begin: usize,
        end: usize,
    ) -> bool;
    fn blst_middle_product_workspace_avx_native_free(
        workspace: *mut blst_middle_product_workspace_avx_native,
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

// SAFETY: construction publishes a fully initialized native allocation.
// Evaluation receives the allocation through a const pointer, the native
// implementation does not mutate it, and destruction requires exclusive Rust
// ownership. Concurrent evaluations therefore share immutable storage only.
unsafe impl Send for MiddleProductPrecomputedB {}
unsafe impl Sync for MiddleProductPrecomputedB {}

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

// SAFETY: construction publishes a fully initialized native allocation.
// Scalar evaluation receives it through a const pointer and does not mutate
// the precomputed table. Drop requires exclusive Rust ownership.
unsafe impl Send for ScalarMiddleProductPrecomputedB {}
unsafe impl Sync for ScalarMiddleProductPrecomputedB {}

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

    /// Returns the work-item count for one ranged phase and stage.
    #[inline]
    pub fn range_work_count(
        &self,
        phase: MiddleProductRangePhase,
        stage: usize,
    ) -> Result<usize, MiddleProductError> {
        let mut count = 0;
        if unsafe {
            blst_middle_product_scalar_range_work_count(
                self.raw, phase, stage, &mut count,
            )
        } {
            Ok(count)
        } else {
            Err(MiddleProductError::BackendFailure)
        }
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
    capacity: usize,
}

// SAFETY: the native allocation has no thread affinity. Concurrent access is
// exposed only through unsafe range methods whose caller must enforce disjoint
// same-phase ranges and phase/stage barriers. Safe one-shot evaluation retains
// an exclusive mutable borrow.
unsafe impl Send for ScalarMiddleProductWorkspace {}
unsafe impl Sync for ScalarMiddleProductWorkspace {}

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
            capacity: precomputed_b.n,
        })
    }

    /// Returns the largest logical precompute accepted by this workspace.
    pub fn capacity(&self) -> usize {
        self.capacity
    }

    /// Loads one disjoint range into shared workspace storage.
    ///
    /// # Safety
    ///
    /// Concurrent callers must use disjoint ranges from this same input-load
    /// phase. The caller must keep `a`, `precomputed_b`, and this workspace
    /// alive through the call and must establish a barrier before advancing to
    /// the next phase. No call may overlap workspace destruction or safe
    /// one-shot evaluation.
    #[inline]
    pub unsafe fn input_range(
        &self,
        a: &[blst_p1],
        precomputed_b: &ScalarMiddleProductPrecomputedB,
        range: Range<usize>,
    ) -> Result<(), MiddleProductError> {
        validate_range_inputs(a.len(), precomputed_b.n, self.capacity)?;
        if unsafe {
            blst_middle_product_scalar_input_range(
                a.as_ptr(),
                precomputed_b.n,
                precomputed_b.raw,
                self.raw,
                range.start,
                range.end,
            )
        } {
            Ok(())
        } else {
            Err(MiddleProductError::BackendFailure)
        }
    }

    /// Executes one disjoint phase/stage range against shared workspace.
    ///
    /// # Safety
    ///
    /// Concurrent calls must use disjoint ranges of the same phase and stage.
    /// The caller supplies every stage and phase barrier and must exclude safe
    /// one-shot evaluation and destruction for the duration of the call.
    #[inline]
    pub unsafe fn execute_range(
        &self,
        precomputed_b: &ScalarMiddleProductPrecomputedB,
        phase: MiddleProductRangePhase,
        stage: usize,
        range: Range<usize>,
    ) -> Result<(), MiddleProductError> {
        validate_workspace_capacity(self.capacity, precomputed_b.n)?;
        if unsafe {
            blst_middle_product_scalar_execute_range(
                precomputed_b.raw,
                self.raw,
                phase,
                stage,
                range.start,
                range.end,
            )
        } {
            Ok(())
        } else {
            Err(MiddleProductError::BackendFailure)
        }
    }

    /// Extracts one disjoint range into a caller-owned output allocation.
    ///
    /// # Safety
    ///
    /// `out` must be valid and initialized for `out_len` target-group elements.
    /// Concurrent output calls must address disjoint ranges of the same output
    /// allocation. The caller must keep the allocation, precompute, and
    /// workspace alive, must have completed all earlier phases, and must
    /// exclude destruction and safe one-shot evaluation.
    #[inline]
    pub unsafe fn output_range(
        &self,
        precomputed_b: &ScalarMiddleProductPrecomputedB,
        out: *mut blst_fp12,
        out_len: usize,
        range: Range<usize>,
    ) -> Result<(), MiddleProductError> {
        validate_output_range(
            out.cast_const(),
            out_len,
            precomputed_b.n,
            self.capacity,
        )?;
        if unsafe {
            blst_middle_product_scalar_output_range(
                precomputed_b.raw,
                self.raw,
                out,
                precomputed_b.n,
                range.start,
                range.end,
            )
        } {
            Ok(())
        } else {
            Err(MiddleProductError::BackendFailure)
        }
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

// SAFETY: construction publishes a fully initialized native allocation. AVX
// evaluation reads the table through a const pointer and uses caller-owned
// workspace for mutable scratch. Drop requires exclusive Rust ownership.
unsafe impl Send for AvxMiddleProductPrecomputedB {}
unsafe impl Sync for AvxMiddleProductPrecomputedB {}

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
    capacity: usize,
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
            capacity: precomputed_b.n,
        })
    }

    /// Returns the largest logical precompute accepted by this workspace.
    pub fn capacity(&self) -> usize {
        self.capacity
    }
}

impl Drop for AvxMiddleProductWorkspace {
    fn drop(&mut self) {
        if !self.raw.is_null() {
            unsafe { blst_middle_product_workspace_avx_free(self.raw) };
        }
    }
}

/// Opaque owner for a precomputed native-AVX G2 operand.
///
/// Unlike [`AvxMiddleProductPrecomputedB`], construction accepts AVX-native G2
/// values directly and performs no scalar/native conversion.
#[derive(Debug)]
pub struct NativeAvxMiddleProductPrecomputedB {
    raw: *mut blst_middle_product_precomputed_b_avx,
    n: usize,
}

// The native allocation is immutable after construction. Evaluation accepts
// it through a const pointer, and destruction requires exclusive ownership.
unsafe impl Send for NativeAvxMiddleProductPrecomputedB {}
unsafe impl Sync for NativeAvxMiddleProductPrecomputedB {}

impl NativeAvxMiddleProductPrecomputedB {
    pub fn new(b: &[blst_avx_p2]) -> Result<Self, MiddleProductError> {
        validate_precompute_input(b)?;
        let n = b.len() / 2;
        let raw = unsafe {
            blst_middle_product_precompute_b_avx_native(b.as_ptr(), n)
        };
        if raw.is_null() {
            return Err(MiddleProductError::PrecomputeFailed);
        }
        Ok(Self { raw, n })
    }

    /// Precomputes directly from typed AVX G2 values without representation
    /// conversion.
    pub fn new_typed(
        b: &[crate::avx::G2Projective],
    ) -> Result<Self, MiddleProductError> {
        Self::new(crate::avx::G2Projective::as_raw_slice(b))
    }

    pub fn len(&self) -> usize {
        self.n
    }

    pub fn is_empty(&self) -> bool {
        self.n == 0
    }

    /// Returns the work-item count for one native-AVX ranged phase and stage.
    #[inline]
    pub fn range_work_count(
        &self,
        phase: MiddleProductRangePhase,
        stage: usize,
    ) -> Result<usize, MiddleProductError> {
        let mut count = 0;
        if unsafe {
            blst_middle_product_avx_native_range_work_count(
                self.raw, phase, stage, &mut count,
            )
        } {
            Ok(count)
        } else {
            Err(MiddleProductError::BackendFailure)
        }
    }
}

impl Drop for NativeAvxMiddleProductPrecomputedB {
    fn drop(&mut self) {
        if !self.raw.is_null() {
            unsafe { blst_middle_product_precomputed_b_avx_free(self.raw) };
        }
    }
}

/// Caller-owned native AVX FFT workspace.
///
/// This workspace contains only backend-native FFT scratch. It has no scalar
/// input-conversion or scalar output-staging buffers.
#[derive(Debug)]
pub struct NativeAvxMiddleProductWorkspace {
    raw: *mut blst_middle_product_workspace_avx_native,
    capacity: usize,
}

// SAFETY: as for the scalar workspace, safe one-shot evaluation is exclusive.
// Shared access is available only through unsafe range methods with an
// explicit disjointness and barrier contract.
unsafe impl Send for NativeAvxMiddleProductWorkspace {}
unsafe impl Sync for NativeAvxMiddleProductWorkspace {}

impl NativeAvxMiddleProductWorkspace {
    pub fn new(
        precomputed_b: &NativeAvxMiddleProductPrecomputedB,
    ) -> Result<Self, MiddleProductError> {
        let raw = unsafe {
            blst_middle_product_workspace_avx_native_new(precomputed_b.raw)
        };
        if raw.is_null() {
            return Err(MiddleProductError::WorkspaceCreationFailed);
        }
        Ok(Self {
            raw,
            capacity: precomputed_b.n,
        })
    }

    /// Returns the largest logical precompute accepted by this workspace.
    pub fn capacity(&self) -> usize {
        self.capacity
    }

    /// Loads one disjoint native-AVX input range into shared workspace.
    ///
    /// # Safety
    ///
    /// The caller must enforce the same disjoint-range, lifetime, and barrier
    /// contract as [`ScalarMiddleProductWorkspace::input_range`].
    #[inline]
    pub unsafe fn input_range(
        &self,
        a: &[crate::avx::G1Projective],
        precomputed_b: &NativeAvxMiddleProductPrecomputedB,
        range: Range<usize>,
    ) -> Result<(), MiddleProductError> {
        validate_range_inputs(a.len(), precomputed_b.n, self.capacity)?;
        if unsafe {
            blst_middle_product_avx_native_input_range(
                crate::avx::G1Projective::as_raw_slice(a).as_ptr(),
                precomputed_b.n,
                precomputed_b.raw,
                self.raw,
                range.start,
                range.end,
            )
        } {
            Ok(())
        } else {
            Err(MiddleProductError::BackendFailure)
        }
    }

    /// Executes one disjoint native-AVX phase/stage range.
    ///
    /// # Safety
    ///
    /// Concurrent calls must use disjoint ranges of the same phase and stage.
    /// The caller supplies every barrier and excludes one-shot evaluation and
    /// destruction while any range call is active.
    #[inline]
    pub unsafe fn execute_range(
        &self,
        precomputed_b: &NativeAvxMiddleProductPrecomputedB,
        phase: MiddleProductRangePhase,
        stage: usize,
        range: Range<usize>,
    ) -> Result<(), MiddleProductError> {
        validate_workspace_capacity(self.capacity, precomputed_b.n)?;
        if unsafe {
            blst_middle_product_avx_native_execute_range(
                precomputed_b.raw,
                self.raw,
                phase,
                stage,
                range.start,
                range.end,
            )
        } {
            Ok(())
        } else {
            Err(MiddleProductError::BackendFailure)
        }
    }

    /// Extracts one disjoint range into native-AVX output storage.
    ///
    /// # Safety
    ///
    /// `out` must be valid and initialized for `out_len` `crate::avx::Gt`
    /// values. Concurrent calls must address disjoint ranges. All earlier
    /// phases must be complete, and no one-shot evaluation or destruction may
    /// overlap this call.
    #[inline]
    pub unsafe fn output_range(
        &self,
        precomputed_b: &NativeAvxMiddleProductPrecomputedB,
        out: *mut crate::avx::Gt,
        out_len: usize,
        range: Range<usize>,
    ) -> Result<(), MiddleProductError> {
        validate_output_range(
            out.cast_const(),
            out_len,
            precomputed_b.n,
            self.capacity,
        )?;
        if unsafe {
            blst_middle_product_avx_native_output_range(
                precomputed_b.raw,
                self.raw,
                out.cast(),
                precomputed_b.n,
                range.start,
                range.end,
            )
        } {
            Ok(())
        } else {
            Err(MiddleProductError::BackendFailure)
        }
    }

    /// Evaluates directly over typed native AVX inputs and outputs.
    ///
    /// The call performs no allocation and does not cross the scalar
    /// representation. `out` must contain at least `a.len()` elements.
    pub fn evaluate_into(
        &mut self,
        a: &[crate::avx::G1Projective],
        precomputed_b: &NativeAvxMiddleProductPrecomputedB,
        out: &mut [crate::avx::Gt],
    ) -> Result<(), MiddleProductError> {
        middle_product_avx_native_with_workspace(
            crate::avx::G1Projective::as_raw_slice(a),
            precomputed_b,
            self,
            crate::avx::Gt::as_raw_slice_mut(out),
        )
    }
}

impl Drop for NativeAvxMiddleProductWorkspace {
    fn drop(&mut self) {
        if !self.raw.is_null() {
            unsafe { blst_middle_product_workspace_avx_native_free(self.raw) };
        }
    }
}

fn validate_precompute_input<T>(b: &[T]) -> Result<(), MiddleProductError> {
    if b.is_empty() {
        return Err(MiddleProductError::EmptyInput);
    }
    if b.len() % 2 != 0 {
        return Err(MiddleProductError::InvalidBLength);
    }
    Ok(())
}

#[inline]
fn validate_workspace_capacity(
    capacity: usize,
    required: usize,
) -> Result<(), MiddleProductError> {
    if capacity < required {
        Err(MiddleProductError::LengthMismatch)
    } else {
        Ok(())
    }
}

#[inline]
fn validate_range_inputs(
    input_len: usize,
    expected: usize,
    capacity: usize,
) -> Result<(), MiddleProductError> {
    if input_len != expected {
        return Err(MiddleProductError::LengthMismatch);
    }
    validate_workspace_capacity(capacity, expected)
}

#[inline]
fn validate_output_range<T>(
    out: *const T,
    out_len: usize,
    expected: usize,
    capacity: usize,
) -> Result<(), MiddleProductError> {
    if out_len < expected {
        return Err(MiddleProductError::OutputTooSmall);
    }
    if out.is_null() {
        return Err(MiddleProductError::LengthMismatch);
    }
    validate_workspace_capacity(capacity, expected)
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
/// `workspace` must have been constructed for a scalar precompute at least as
/// large as `precomputed_b`. It may be reused by subsequent calls and smaller
/// precomputes, but is mutably borrowed so concurrent one-shot evaluations
/// require distinct workspaces.
pub fn middle_product_scalar_with_workspace(
    a: &[blst_p1],
    precomputed_b: &ScalarMiddleProductPrecomputedB,
    workspace: &mut ScalarMiddleProductWorkspace,
    out: &mut [blst_fp12],
) -> Result<(), MiddleProductError> {
    validate_evaluation_buffers(a.len(), precomputed_b.n, out.len())?;
    if workspace.capacity < precomputed_b.n {
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
/// A workspace created for a larger precompute may evaluate a smaller one.
pub fn middle_product_avx_with_workspace(
    a: &[blst_p1],
    precomputed_b: &AvxMiddleProductPrecomputedB,
    workspace: &mut AvxMiddleProductWorkspace,
    out: &mut [blst_fp12],
) -> Result<(), MiddleProductError> {
    validate_evaluation_buffers(a.len(), precomputed_b.n, out.len())?;
    if workspace.capacity < precomputed_b.n {
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

/// Computes an AVX middle product using native input and output types.
///
/// Construction and evaluation never cross the scalar representation. The
/// native call allocates nothing after `workspace` is constructed, and a
/// workspace created for a larger precompute may evaluate a smaller one.
pub fn middle_product_avx_native_with_workspace(
    a: &[blst_avx_p1],
    precomputed_b: &NativeAvxMiddleProductPrecomputedB,
    workspace: &mut NativeAvxMiddleProductWorkspace,
    out: &mut [blst_avx_gt],
) -> Result<(), MiddleProductError> {
    validate_evaluation_buffers(a.len(), precomputed_b.n, out.len())?;
    if workspace.capacity < precomputed_b.n {
        return Err(MiddleProductError::LengthMismatch);
    }
    if unsafe {
        blst_middle_product_avx_native_with_workspace(
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

    fn assert_send_sync<T: Send + Sync>() {}

    #[test]
    fn immutable_precomputations_are_send_and_sync() {
        assert_send_sync::<MiddleProductPrecomputedB>();
        assert_send_sync::<ScalarMiddleProductPrecomputedB>();
        assert_send_sync::<AvxMiddleProductPrecomputedB>();
        assert_send_sync::<NativeAvxMiddleProductPrecomputedB>();
        assert_send_sync::<ScalarMiddleProductWorkspace>();
        assert_send_sync::<NativeAvxMiddleProductWorkspace>();
    }

    fn scalar_ranged_serial(
        a: &[blst_p1],
        precomputed: &ScalarMiddleProductPrecomputedB,
        workspace: &ScalarMiddleProductWorkspace,
        out: &mut [blst_fp12],
    ) {
        let execute = |phase, stage| {
            let count = precomputed.range_work_count(phase, stage).unwrap();
            unsafe {
                workspace
                    .execute_range(precomputed, phase, stage, 0..count)
                    .unwrap();
            }
        };

        let input_count = precomputed
            .range_work_count(MiddleProductRangePhase::InputLoad, 0)
            .unwrap();
        unsafe {
            workspace
                .input_range(a, precomputed, 0..input_count)
                .unwrap();
        }
        execute(MiddleProductRangePhase::ForwardBitReverse, 0);
        for stage in 0usize.. {
            let Ok(j0_count) = precomputed
                .range_work_count(MiddleProductRangePhase::ForwardJ0, stage)
            else {
                break;
            };
            unsafe {
                workspace
                    .execute_range(
                        precomputed,
                        MiddleProductRangePhase::ForwardJ0,
                        stage,
                        0..j0_count,
                    )
                    .unwrap();
            }
            let nonzero_count = precomputed
                .range_work_count(
                    MiddleProductRangePhase::ForwardNonzero,
                    stage,
                )
                .unwrap();
            unsafe {
                workspace
                    .execute_range(
                        precomputed,
                        MiddleProductRangePhase::ForwardNonzero,
                        stage,
                        0..nonzero_count,
                    )
                    .unwrap();
            }
        }
        execute(MiddleProductRangePhase::AffineConversion, 0);
        execute(MiddleProductRangePhase::PointwisePairing, 0);
        execute(MiddleProductRangePhase::InverseBitReverse, 0);
        for stage in 0usize.. {
            let Ok(count) = precomputed.range_work_count(
                MiddleProductRangePhase::InverseButterflies,
                stage,
            ) else {
                break;
            };
            unsafe {
                workspace
                    .execute_range(
                        precomputed,
                        MiddleProductRangePhase::InverseButterflies,
                        stage,
                        0..count,
                    )
                    .unwrap();
            }
        }
        execute(MiddleProductRangePhase::Normalization, 0);
        unsafe {
            workspace
                .output_range(
                    precomputed,
                    out.as_mut_ptr(),
                    out.len(),
                    0..out.len(),
                )
                .unwrap();
        }
    }

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
    fn scalar_range_api_reuses_larger_workspace_and_matches_one_shot() {
        let g1 = unsafe { *blst_p1_generator() };
        let g2 = unsafe { *blst_p2_generator() };
        let large_b = [g2; 8];
        let small_b = [g2; 4];
        let large = ScalarMiddleProductPrecomputedB::new(&large_b).unwrap();
        let small = ScalarMiddleProductPrecomputedB::new(&small_b).unwrap();
        let workspace = ScalarMiddleProductWorkspace::new(&large).unwrap();
        let mut exclusive_workspace =
            ScalarMiddleProductWorkspace::new(&large).unwrap();
        let a = [g1; 2];
        let mut expected = [zero_fp12(); 2];
        let mut ranged = [zero_fp12(); 2];

        middle_product_scalar_with_workspace(
            &a,
            &small,
            &mut exclusive_workspace,
            &mut expected,
        )
        .unwrap();
        scalar_ranged_serial(&a, &small, &workspace, &mut ranged);

        assert_eq!(workspace.capacity(), large.len());
        assert_eq!(ranged, expected);
    }

    #[test]
    fn scalar_range_api_rejects_invalid_range_and_small_workspace() {
        let g2 = unsafe { *blst_p2_generator() };
        let small = ScalarMiddleProductPrecomputedB::new(&[g2; 2]).unwrap();
        let large = ScalarMiddleProductPrecomputedB::new(&[g2; 4]).unwrap();
        let workspace = ScalarMiddleProductWorkspace::new(&small).unwrap();

        assert_eq!(
            unsafe {
                workspace.execute_range(
                    &small,
                    MiddleProductRangePhase::ForwardBitReverse,
                    0,
                    1..0,
                )
            },
            Err(MiddleProductError::BackendFailure)
        );
        assert_eq!(
            unsafe {
                workspace.execute_range(
                    &large,
                    MiddleProductRangePhase::ForwardBitReverse,
                    0,
                    0..0,
                )
            },
            Err(MiddleProductError::LengthMismatch)
        );
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
    fn native_avx_backend_accepts_and_returns_only_native_values() {
        if !unsafe { blst_avx_enabled() } {
            return;
        }

        let scalar_b = unsafe { [*blst_p2_generator(), *blst_p2_generator()] };
        let scalar_a = unsafe { [*blst_p1_generator()] };
        let mut native_b = [blst_avx_p2::default(); 2];
        let mut native_a = [blst_avx_p1::default(); 1];
        for (native, scalar) in native_b.iter_mut().zip(scalar_b.iter()) {
            assert!(unsafe { blst_avx_p2_from_scalar(native, scalar) });
        }
        assert!(unsafe {
            blst_avx_p1_from_scalar(&mut native_a[0], &scalar_a[0])
        });

        let precomputed =
            NativeAvxMiddleProductPrecomputedB::new(&native_b).unwrap();
        let mut workspace =
            NativeAvxMiddleProductWorkspace::new(&precomputed).unwrap();
        let mut native_out = [blst_avx_gt::default(); 1];
        middle_product_avx_native_with_workspace(
            &native_a,
            &precomputed,
            &mut workspace,
            &mut native_out,
        )
        .unwrap();

        let scalar_precomputed =
            ScalarMiddleProductPrecomputedB::new(&scalar_b).unwrap();
        let mut scalar_workspace =
            ScalarMiddleProductWorkspace::new(&scalar_precomputed).unwrap();
        let mut expected = [zero_fp12(); 1];
        middle_product_scalar_with_workspace(
            &scalar_a,
            &scalar_precomputed,
            &mut scalar_workspace,
            &mut expected,
        )
        .unwrap();

        let mut converted = [zero_fp12(); 1];
        assert!(unsafe {
            blst_fp12_from_avx(&mut converted[0], &native_out[0])
        });
        assert_eq!(converted, expected);
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
