// Copyright Supranational LLC
// Licensed under the Apache License, Version 2.0, see LICENSE for details.
// SPDX-License-Identifier: Apache-2.0

//! Pairing-context and uniqueness helpers for the high-level Rust surface.

#[cfg(feature = "full-rust-bindings")]
use crate::*;

#[cfg(feature = "full-rust-bindings")]
use alloc::alloc::{alloc_zeroed, dealloc, handle_alloc_error};
#[cfg(feature = "full-rust-bindings")]
use alloc::vec;
#[cfg(feature = "full-rust-bindings")]
use alloc::vec::Vec;
#[cfg(feature = "full-rust-bindings")]
use core::alloc::Layout;
#[cfg(feature = "full-rust-bindings")]
use core::any::Any;
#[cfg(feature = "full-rust-bindings")]
use core::ptr::{self, NonNull};

#[cfg(feature = "full-rust-bindings")]
extern "C" {
    fn blst_pairing_alignof() -> usize;
}

/// Owned pairing context with correct native size and alignment.
#[cfg(feature = "full-rust-bindings")]
#[derive(Debug)]
pub struct Pairing {
    storage: NonNull<u8>,
    layout: Layout,
}

#[cfg(feature = "full-rust-bindings")]
impl Pairing {
    /// Allocates and initializes a pairing context for the given DST.
    pub fn new(hash_or_encode: bool, dst: &[u8]) -> Self {
        let size = unsafe { blst_pairing_sizeof() };
        let align = unsafe { blst_pairing_alignof() };
        let layout = Layout::from_size_align(size.max(1), align)
            .expect("invalid blst_pairing layout");
        let storage = unsafe {
            let ptr = alloc_zeroed(layout);
            if ptr.is_null() {
                handle_alloc_error(layout);
            }
            NonNull::new_unchecked(ptr)
        };
        let mut obj = Self { storage, layout };
        obj.init(hash_or_encode, dst);
        obj
    }

    /// Reinitializes an existing pairing context.
    pub fn init(&mut self, hash_or_encode: bool, dst: &[u8]) {
        unsafe {
            blst_pairing_init(
                self.ctx(),
                hash_or_encode,
                dst.as_ptr(),
                dst.len(),
            )
        }
    }

    fn ctx(&mut self) -> *mut blst_pairing {
        self.storage.as_ptr() as *mut blst_pairing
    }

    fn const_ctx(&self) -> *const blst_pairing {
        self.storage.as_ptr() as *const blst_pairing
    }

    /// Aggregates one verification item into the pairing context.
    pub fn aggregate(
        &mut self,
        pk: &dyn Any,
        pk_validate: bool,
        sig: &dyn Any,
        sig_groupcheck: bool,
        msg: &[u8],
        aug: &[u8],
    ) -> BLST_ERROR {
        if pk.is::<blst_p1_affine>() {
            unsafe {
                blst_pairing_chk_n_aggr_pk_in_g1(
                    self.ctx(),
                    match pk.downcast_ref::<blst_p1_affine>() {
                        Some(pk) => pk,
                        None => ptr::null(),
                    },
                    pk_validate,
                    match sig.downcast_ref::<blst_p2_affine>() {
                        Some(sig) => sig,
                        None => ptr::null(),
                    },
                    sig_groupcheck,
                    msg.as_ptr(),
                    msg.len(),
                    aug.as_ptr(),
                    aug.len(),
                )
            }
        } else if pk.is::<blst_p2_affine>() {
            unsafe {
                blst_pairing_chk_n_aggr_pk_in_g2(
                    self.ctx(),
                    match pk.downcast_ref::<blst_p2_affine>() {
                        Some(pk) => pk,
                        None => ptr::null(),
                    },
                    pk_validate,
                    match sig.downcast_ref::<blst_p1_affine>() {
                        Some(sig) => sig,
                        None => ptr::null(),
                    },
                    sig_groupcheck,
                    msg.as_ptr(),
                    msg.len(),
                    aug.as_ptr(),
                    aug.len(),
                )
            }
        } else {
            panic!("whaaaa?")
        }
    }

    /// Aggregates one weighted verification item into the pairing context.
    pub fn mul_n_aggregate(
        &mut self,
        pk: &dyn Any,
        pk_validate: bool,
        sig: &dyn Any,
        sig_groupcheck: bool,
        scalar: &[u8],
        nbits: usize,
        msg: &[u8],
        aug: &[u8],
    ) -> BLST_ERROR {
        if pk.is::<blst_p1_affine>() {
            unsafe {
                blst_pairing_chk_n_mul_n_aggr_pk_in_g1(
                    self.ctx(),
                    match pk.downcast_ref::<blst_p1_affine>() {
                        Some(pk) => pk,
                        None => ptr::null(),
                    },
                    pk_validate,
                    match sig.downcast_ref::<blst_p2_affine>() {
                        Some(sig) => sig,
                        None => ptr::null(),
                    },
                    sig_groupcheck,
                    scalar.as_ptr(),
                    nbits,
                    msg.as_ptr(),
                    msg.len(),
                    aug.as_ptr(),
                    aug.len(),
                )
            }
        } else if pk.is::<blst_p2_affine>() {
            unsafe {
                blst_pairing_chk_n_mul_n_aggr_pk_in_g2(
                    self.ctx(),
                    match pk.downcast_ref::<blst_p2_affine>() {
                        Some(pk) => pk,
                        None => ptr::null(),
                    },
                    pk_validate,
                    match sig.downcast_ref::<blst_p1_affine>() {
                        Some(sig) => sig,
                        None => ptr::null(),
                    },
                    sig_groupcheck,
                    scalar.as_ptr(),
                    nbits,
                    msg.as_ptr(),
                    msg.len(),
                    aug.as_ptr(),
                    aug.len(),
                )
            }
        } else {
            panic!("whaaaa?")
        }
    }

    /// Commits pending work in the pairing context.
    pub fn commit(&mut self) {
        unsafe { blst_pairing_commit(self.ctx()) }
    }

    /// Merges another compatible pairing context into this one.
    pub fn merge(&mut self, ctx1: &Self) -> BLST_ERROR {
        unsafe { blst_pairing_merge(self.ctx(), ctx1.const_ctx()) }
    }

    /// Aggregates one raw affine `(q, p)` pair into the context.
    pub fn raw_aggregate(&mut self, q: &blst_p2_affine, p: &blst_p1_affine) {
        unsafe { blst_pairing_raw_aggregate(self.ctx(), q, p) }
    }

    /// Runs the final verification step.
    pub fn finalverify(&self, sig: Option<&blst_fp12>) -> bool {
        unsafe {
            blst_pairing_finalverify(
                self.const_ctx(),
                sig.map_or(ptr::null(), |s| s as *const blst_fp12),
            )
        }
    }

    /// Computes the aggregate GT value for one signature.
    pub fn aggregated(gtsig: &mut blst_fp12, sig: &dyn Any) {
        unsafe {
            if sig.is::<blst_p1_affine>() {
                blst_aggregated_in_g1(
                    gtsig,
                    match sig.downcast_ref::<blst_p1_affine>() {
                        Some(sig) => sig,
                        None => ptr::null(),
                    },
                )
            } else {
                blst_aggregated_in_g2(
                    gtsig,
                    match sig.downcast_ref::<blst_p2_affine>() {
                        Some(sig) => sig,
                        None => ptr::null(),
                    },
                )
            }
        }
    }

    /// Returns the current pairing accumulator as Fp12.
    pub fn as_fp12(&mut self) -> blst_fp12 {
        unsafe { blst_pairing_as_fp12(self.ctx()) }
    }
}

#[cfg(feature = "full-rust-bindings")]
impl Drop for Pairing {
    fn drop(&mut self) {
        unsafe { dealloc(self.storage.as_ptr(), self.layout) };
    }
}

#[cfg(feature = "full-rust-bindings")]
unsafe impl Send for Pairing {}

/// Returns `true` when all message slices are distinct.
#[cfg(feature = "full-rust-bindings")]
pub fn uniq(msgs: &[&[u8]]) -> bool {
    let n_elems = msgs.len();

    if n_elems == 1 {
        return true;
    } else if n_elems == 2 {
        return msgs[0] != msgs[1];
    }

    let mut v: Vec<u64> = vec![0; unsafe { blst_uniq_sizeof(n_elems) } / 8];
    let ctx = v.as_mut_ptr() as *mut blst_uniq;

    unsafe { blst_uniq_init(ctx) };

    for msg in msgs.iter() {
        if !unsafe { blst_uniq_test(ctx, msg.as_ptr(), msg.len()) } {
            return false;
        }
    }

    true
}

/// Debug helper that prints a byte slice in hexadecimal.
#[cfg(feature = "std")]
pub fn print_bytes(bytes: &[u8], name: &str) {
    print!("{} ", name);
    for b in bytes.iter() {
        print!("{:02x}", b);
    }
    println!();
}
