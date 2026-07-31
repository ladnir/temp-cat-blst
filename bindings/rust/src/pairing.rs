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
use core::mem;
#[cfg(feature = "full-rust-bindings")]
use core::ptr::{self, NonNull};

#[cfg(feature = "full-rust-bindings")]
extern "C" {
    fn blst_pairing_alignof() -> usize;
}

#[cfg(feature = "full-rust-bindings")]
fn pairing_layout(
    context_size: usize,
    context_align: usize,
    dst_capacity: usize,
) -> Option<Layout> {
    context_size.checked_add(dst_capacity).and_then(|size| {
        Layout::from_size_align(size.max(1), context_align).ok()
    })
}

/// Owned pairing context with correct native size and alignment.
#[cfg(feature = "full-rust-bindings")]
#[derive(Debug)]
pub struct Pairing {
    storage: NonNull<u8>,
    layout: Layout,
    context_size: usize,
}

#[cfg(feature = "full-rust-bindings")]
impl Pairing {
    /// Allocates and initializes a pairing context with an owned copy of `dst`.
    ///
    /// The native context and DST occupy one allocation. The DST begins exactly
    /// at the end of the native context, which engages the native context's
    /// relocatable trailing-DST representation.
    pub fn new(hash_or_encode: bool, dst: &[u8]) -> Self {
        Self::allocate(hash_or_encode, dst)
    }

    fn allocate(hash_or_encode: bool, dst: &[u8]) -> Self {
        let context_size = unsafe { blst_pairing_sizeof() };
        let align = unsafe { blst_pairing_alignof() };
        let layout = pairing_layout(context_size, align, dst.len())
            .expect("invalid or overflowing blst_pairing layout");
        let storage = unsafe {
            let ptr = alloc_zeroed(layout);
            if ptr.is_null() {
                handle_alloc_error(layout);
            }
            NonNull::new_unchecked(ptr)
        };

        let mut obj = Self {
            storage,
            layout,
            context_size,
        };
        obj.initialize(hash_or_encode, dst);
        obj
    }

    /// Reinitializes the pairing context with an owned copy of `dst`.
    ///
    /// Existing storage is reused when its trailing capacity is sufficient.
    /// Growth initializes a replacement allocation before releasing the old
    /// context, so `dst` remains valid throughout reinitialization.
    pub fn init(&mut self, hash_or_encode: bool, dst: &[u8]) {
        if dst.len() <= self.dst_capacity() {
            self.initialize(hash_or_encode, dst);
        } else {
            let mut replacement = Self::allocate(hash_or_encode, dst);
            mem::swap(self, &mut replacement);
        }
    }

    fn initialize(&mut self, hash_or_encode: bool, dst: &[u8]) {
        let owned_dst = self.dst_storage();
        unsafe {
            if !dst.is_empty() {
                ptr::copy_nonoverlapping(dst.as_ptr(), owned_dst, dst.len());
            }
            blst_pairing_init(self.ctx(), hash_or_encode, owned_dst, dst.len())
        }
    }

    fn dst_storage(&self) -> *mut u8 {
        unsafe { self.storage.as_ptr().add(self.context_size) }
    }

    fn dst_capacity(&self) -> usize {
        self.layout.size() - self.context_size
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
    ///
    /// Returns [`BLST_ERROR::BLST_BAD_SCALAR`] without entering native code
    /// when `scalar` does not contain all `nbits` requested by the caller.
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
        let scalar_bytes = nbits / 8 + usize::from(nbits % 8 != 0);
        if scalar.len() < scalar_bytes {
            return BLST_ERROR::BLST_BAD_SCALAR;
        }

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

// SAFETY: `Pairing` exclusively owns its native context and trailing DST
// allocation. Moving it between threads does not change either address, and
// the native context has no thread affinity. Mutation still requires `&mut
// self`; `Pairing` deliberately does not implement `Sync` or `Clone`.
#[cfg(feature = "full-rust-bindings")]
unsafe impl Send for Pairing {}

#[cfg(all(test, feature = "full-rust-bindings"))]
mod tests {
    use super::*;
    use core::slice;

    fn assert_owned_dst(pairing: &Pairing, expected: &[u8]) {
        let native_dst = unsafe { blst_pairing_get_dst(pairing.const_ctx()) };
        assert_eq!(native_dst, pairing.dst_storage());
        assert_eq!(
            unsafe { slice::from_raw_parts(native_dst, expected.len()) },
            expected
        );
    }

    fn generators() -> (blst_p1_affine, blst_p2_affine) {
        unsafe { (*blst_p1_affine_generator(), *blst_p2_affine_generator()) }
    }

    fn aggregate_generator_pair(
        pairing: &mut Pairing,
        msg: &[u8],
    ) -> BLST_ERROR {
        let (g1, g2) = generators();
        pairing.aggregate(&g1, true, &g2, true, msg, &[])
    }

    fn assert_mul_n_result(scalar: &[u8], nbits: usize, expected: BLST_ERROR) {
        let (g1, g2) = generators();

        let mut g1_pairing = Pairing::new(true, b"scalar length DST");
        assert_eq!(
            g1_pairing.mul_n_aggregate(
                &g1,
                false,
                &g2,
                false,
                scalar,
                nbits,
                b"message",
                &[],
            ),
            expected
        );

        let mut g2_pairing = Pairing::new(true, b"scalar length DST");
        assert_eq!(
            g2_pairing.mul_n_aggregate(
                &g2,
                false,
                &g1,
                false,
                scalar,
                nbits,
                b"message",
                &[],
            ),
            expected
        );
    }

    #[cfg(feature = "std")]
    #[test]
    fn owns_temporary_dst_and_moves_safely_between_threads() {
        let pairing = {
            let temporary_dst = b"temporary pairing DST".to_vec();
            Pairing::new(true, &temporary_dst)
        };

        let pairing = std::thread::spawn(move || {
            assert_owned_dst(&pairing, b"temporary pairing DST");
            pairing
        })
        .join()
        .expect("pairing thread panicked");

        assert_owned_dst(&pairing, b"temporary pairing DST");
    }

    #[test]
    fn init_reuses_capacity_and_reallocates_before_growth() {
        let mut pairing = Pairing::new(true, b"initial-capacity");
        let initial_storage = pairing.storage;
        let initial_capacity = pairing.dst_capacity();

        pairing.init(false, b"short");
        assert_eq!(pairing.storage, initial_storage);
        assert_eq!(pairing.dst_capacity(), initial_capacity);
        assert_owned_dst(&pairing, b"short");

        let grown_dst = vec![0x5au8; initial_capacity + 1];
        pairing.init(true, &grown_dst);
        assert_ne!(pairing.storage, initial_storage);
        assert_eq!(pairing.dst_capacity(), grown_dst.len());
        assert_owned_dst(&pairing, &grown_dst);
    }

    #[test]
    fn empty_dst_uses_trailing_native_representation() {
        let mut pairing = Pairing::new(true, &[]);
        assert_eq!(
            unsafe { blst_pairing_get_dst(pairing.const_ctx()) },
            pairing.dst_storage()
        );
        assert_eq!(pairing.dst_capacity(), 0);
        assert_eq!(pairing.layout.size(), pairing.context_size);
        assert_eq!(pairing.layout.align(), unsafe { blst_pairing_alignof() });

        assert_eq!(
            aggregate_generator_pair(&mut pairing, b"empty DST message"),
            BLST_ERROR::BLST_SUCCESS
        );
        pairing.commit();
        assert!(!pairing.finalverify(None));
    }

    #[test]
    fn allocation_layout_rejects_size_overflow() {
        assert!(pairing_layout(1, 1, usize::MAX).is_none());
    }

    #[test]
    fn independently_owned_contexts_merge_without_borrowing_dst() {
        let mut left = {
            let dst = b"shared merge DST".to_vec();
            Pairing::new(true, &dst)
        };
        let mut right = {
            let dst = b"shared merge DST".to_vec();
            Pairing::new(true, &dst)
        };

        assert_eq!(
            aggregate_generator_pair(&mut left, b"left"),
            BLST_ERROR::BLST_SUCCESS
        );
        assert_eq!(
            aggregate_generator_pair(&mut right, b"right"),
            BLST_ERROR::BLST_SUCCESS
        );
        left.commit();
        right.commit();
        assert_eq!(left.merge(&right), BLST_ERROR::BLST_SUCCESS);
        drop(right);
        assert_owned_dst(&left, b"shared merge DST");
    }

    #[test]
    fn mul_n_aggregate_rejects_undersized_scalars_for_both_groups() {
        for (scalar, nbits) in [
            (&[][..], 1),
            (&[][..], 7),
            (&[][..], 8),
            (&[0u8][..], 9),
            (&[][..], usize::MAX),
        ] {
            assert_mul_n_result(scalar, nbits, BLST_ERROR::BLST_BAD_SCALAR);
        }
    }

    #[test]
    fn mul_n_aggregate_accepts_zero_and_exact_scalar_lengths() {
        for (scalar, nbits) in [
            (&[][..], 0),
            (&[1u8][..], 1),
            (&[1u8][..], 7),
            (&[1u8][..], 8),
            (&[1u8, 0u8][..], 9),
        ] {
            assert_mul_n_result(scalar, nbits, BLST_ERROR::BLST_SUCCESS);
        }
    }
}

/// Returns `true` when all message slices are distinct.
#[cfg(feature = "full-rust-bindings")]
pub fn uniq(msgs: &[&[u8]]) -> bool {
    let n_elems = msgs.len();

    if n_elems == 0 || n_elems == 1 {
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
