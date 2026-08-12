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
use core::mem;
#[cfg(feature = "full-rust-bindings")]
use core::ptr::{self, NonNull};

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

/// Owned pairing context with correct native size, alignment, and DST lifetime.
///
/// Each typed aggregate call selects G1-public-key/G2-signature or
/// G2-public-key/G1-signature orientation at compile time. The native context
/// records the first selected orientation and rejects later calls in the other
/// orientation. Aggregate zero or more items, call [`Pairing::commit`], then call
/// [`Pairing::finalverify`] or merge committed worker contexts. Aggregation may
/// resume after a commit. [`Pairing::raw_aggregate`] and
/// [`Pairing::miller_result`] are expert operations over trusted affine inputs
/// and pre-final-exponentiation values.
#[cfg(feature = "full-rust-bindings")]
#[derive(Debug)]
pub struct Pairing {
    storage: NonNull<u8>,
    layout: Layout,
    context_size: usize,
}

#[cfg(feature = "full-rust-bindings")]
mod sealed {
    pub trait Sealed {}
}

/// An affine signature representation supported by the pairing context.
///
/// This trait is sealed and implemented only for [`blst_p1_affine`] and
/// [`blst_p2_affine`]. It provides compile-time pairing orientation without
/// runtime type inspection.
#[cfg(feature = "full-rust-bindings")]
pub trait PairingSignature: sealed::Sealed + Sized {
    #[doc(hidden)]
    unsafe fn write_miller_result(out: *mut blst_fp12, signature: *const Self);
}

/// An affine public-key representation supported by the pairing context.
///
/// The associated signature type fixes the pairing orientation at compile
/// time: G1 public keys use G2 signatures and G2 public keys use G1
/// signatures. The trait is sealed and cannot be implemented externally.
#[cfg(feature = "full-rust-bindings")]
pub trait PairingPublicKey: sealed::Sealed + Sized {
    type Signature: PairingSignature;

    #[doc(hidden)]
    unsafe fn aggregate_into(
        ctx: *mut blst_pairing,
        pk: *const Self,
        pk_validate: bool,
        signature: *const Self::Signature,
        sig_groupcheck: bool,
        msg: &[u8],
        aug: &[u8],
    ) -> BLST_ERROR;

    #[doc(hidden)]
    unsafe fn mul_n_aggregate_into(
        ctx: *mut blst_pairing,
        pk: *const Self,
        pk_validate: bool,
        signature: *const Self::Signature,
        sig_groupcheck: bool,
        scalar: &[u8],
        nbits: usize,
        msg: &[u8],
        aug: &[u8],
    ) -> BLST_ERROR;
}

#[cfg(feature = "full-rust-bindings")]
impl sealed::Sealed for blst_p1_affine {}
#[cfg(feature = "full-rust-bindings")]
impl sealed::Sealed for blst_p2_affine {}

#[cfg(feature = "full-rust-bindings")]
impl PairingSignature for blst_p1_affine {
    unsafe fn write_miller_result(out: *mut blst_fp12, signature: *const Self) {
        unsafe { blst_aggregated_in_g1(out, signature) }
    }
}

#[cfg(feature = "full-rust-bindings")]
impl PairingSignature for blst_p2_affine {
    unsafe fn write_miller_result(out: *mut blst_fp12, signature: *const Self) {
        unsafe { blst_aggregated_in_g2(out, signature) }
    }
}

#[cfg(feature = "full-rust-bindings")]
impl PairingPublicKey for blst_p1_affine {
    type Signature = blst_p2_affine;

    unsafe fn aggregate_into(
        ctx: *mut blst_pairing,
        pk: *const Self,
        pk_validate: bool,
        signature: *const Self::Signature,
        sig_groupcheck: bool,
        msg: &[u8],
        aug: &[u8],
    ) -> BLST_ERROR {
        unsafe {
            blst_pairing_chk_n_aggr_pk_in_g1(
                ctx,
                pk,
                pk_validate,
                signature,
                sig_groupcheck,
                msg.as_ptr(),
                msg.len(),
                aug.as_ptr(),
                aug.len(),
            )
        }
    }

    unsafe fn mul_n_aggregate_into(
        ctx: *mut blst_pairing,
        pk: *const Self,
        pk_validate: bool,
        signature: *const Self::Signature,
        sig_groupcheck: bool,
        scalar: &[u8],
        nbits: usize,
        msg: &[u8],
        aug: &[u8],
    ) -> BLST_ERROR {
        unsafe {
            blst_pairing_chk_n_mul_n_aggr_pk_in_g1(
                ctx,
                pk,
                pk_validate,
                signature,
                sig_groupcheck,
                scalar.as_ptr(),
                nbits,
                msg.as_ptr(),
                msg.len(),
                aug.as_ptr(),
                aug.len(),
            )
        }
    }
}

#[cfg(feature = "full-rust-bindings")]
impl PairingPublicKey for blst_p2_affine {
    type Signature = blst_p1_affine;

    unsafe fn aggregate_into(
        ctx: *mut blst_pairing,
        pk: *const Self,
        pk_validate: bool,
        signature: *const Self::Signature,
        sig_groupcheck: bool,
        msg: &[u8],
        aug: &[u8],
    ) -> BLST_ERROR {
        unsafe {
            blst_pairing_chk_n_aggr_pk_in_g2(
                ctx,
                pk,
                pk_validate,
                signature,
                sig_groupcheck,
                msg.as_ptr(),
                msg.len(),
                aug.as_ptr(),
                aug.len(),
            )
        }
    }

    unsafe fn mul_n_aggregate_into(
        ctx: *mut blst_pairing,
        pk: *const Self,
        pk_validate: bool,
        signature: *const Self::Signature,
        sig_groupcheck: bool,
        scalar: &[u8],
        nbits: usize,
        msg: &[u8],
        aug: &[u8],
    ) -> BLST_ERROR {
        unsafe {
            blst_pairing_chk_n_mul_n_aggr_pk_in_g2(
                ctx,
                pk,
                pk_validate,
                signature,
                sig_groupcheck,
                scalar.as_ptr(),
                nbits,
                msg.as_ptr(),
                msg.len(),
                aug.as_ptr(),
                aug.len(),
            )
        }
    }
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

    /// Reinitializes the pairing context with an owned copy of `dst`,
    /// discarding all prior aggregate state.
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
    ///
    /// `P` selects the pairing orientation at compile time. `sig` is supplied
    /// once for an internally aggregated signature and is `None` for later
    /// public-key/message contributions. The signature is fully subgroup
    /// checked when `sig_groupcheck` is true. The public key is subgroup
    /// checked only when `pk_validate` is true and must otherwise have been
    /// validated by the caller.
    pub fn aggregate<P: PairingPublicKey>(
        &mut self,
        pk: &P,
        pk_validate: bool,
        sig: Option<&P::Signature>,
        sig_groupcheck: bool,
        msg: &[u8],
        aug: &[u8],
    ) -> BLST_ERROR {
        unsafe {
            P::aggregate_into(
                self.ctx(),
                pk,
                pk_validate,
                sig.map_or(ptr::null(), |signature| signature),
                sig_groupcheck,
                msg,
                aug,
            )
        }
    }

    /// Aggregates one weighted verification item into the pairing context.
    ///
    /// Returns [`BLST_ERROR::BLST_BAD_SCALAR`] without entering native code
    /// when `nbits` exceeds 256 or `scalar` does not contain all requested
    /// bits.
    pub fn mul_n_aggregate<P: PairingPublicKey>(
        &mut self,
        pk: &P,
        pk_validate: bool,
        sig: Option<&P::Signature>,
        sig_groupcheck: bool,
        scalar: &[u8],
        nbits: usize,
        msg: &[u8],
        aug: &[u8],
    ) -> BLST_ERROR {
        if nbits > 256 {
            return BLST_ERROR::BLST_BAD_SCALAR;
        }
        let scalar_bytes = nbits / 8 + usize::from(nbits % 8 != 0);
        if scalar.len() < scalar_bytes {
            return BLST_ERROR::BLST_BAD_SCALAR;
        }

        unsafe {
            P::mul_n_aggregate_into(
                self.ctx(),
                pk,
                pk_validate,
                sig.map_or(ptr::null(), |signature| signature),
                sig_groupcheck,
                scalar,
                nbits,
                msg,
                aug,
            )
        }
    }

    /// Flushes pending pairing terms into the Miller product.
    ///
    /// This is idempotent when no terms are pending and is required before
    /// `merge` or `finalverify`.
    pub fn commit(&mut self) {
        unsafe { blst_pairing_commit(self.ctx()) }
    }

    /// Merges another committed, same-orientation pairing context into this one.
    ///
    /// Both contexts must have no pending terms. The source is unchanged, and
    /// an error leaves this context unchanged.
    pub fn merge(&mut self, ctx1: &Self) -> BLST_ERROR {
        unsafe { blst_pairing_merge(self.ctx(), ctx1.const_ctx()) }
    }

    /// Aggregates one raw affine `(G2 q, G1 p)` pair into the context.
    ///
    /// No curve or subgroup validation is performed; inputs must be trusted
    /// on-curve affine points. Either infinity contributes the Miller identity.
    pub fn raw_aggregate(&mut self, q: &blst_p2_affine, p: &blst_p1_affine) {
        unsafe { blst_pairing_raw_aggregate(self.ctx(), q, p) }
    }

    /// Runs final verification on a committed context.
    ///
    /// `sig` is an optional trusted external Miller contribution, not a GT
    /// element. When absent, the internally accumulated signature is used.
    pub fn finalverify(&self, sig: Option<&blst_fp12>) -> bool {
        unsafe {
            blst_pairing_finalverify(
                self.const_ctx(),
                sig.map_or(ptr::null(), |s| s as *const blst_fp12),
            )
        }
    }

    /// Computes the pre-final-exponentiation Miller contribution for a fully
    /// validated affine G1 or G2 signature.
    pub fn aggregated<S: PairingSignature>(gtsig: &mut blst_fp12, sig: &S) {
        unsafe { S::write_miller_result(gtsig, sig) }
    }

    /// Commits pending pairs and returns their Miller product. This value has
    /// not undergone final exponentiation and therefore is not a GT element.
    pub fn miller_result(&mut self) -> blst_fp12 {
        let mut out = blst_fp12::default();
        unsafe { blst_pairing_get_miller_result(&mut out, self.ctx()) };
        out
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
        pairing.aggregate(&g1, true, Some(&g2), true, msg, &[])
    }

    fn assert_pairing_orientation<P, S>()
    where
        P: PairingPublicKey<Signature = S>,
        S: PairingSignature,
    {
    }

    fn assert_mul_n_result(scalar: &[u8], nbits: usize, expected: BLST_ERROR) {
        let (g1, g2) = generators();

        let mut g1_pairing = Pairing::new(true, b"scalar length DST");
        assert_eq!(
            g1_pairing.mul_n_aggregate(
                &g1,
                false,
                Some(&g2),
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
                Some(&g1),
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
    fn pairing_orientation_is_fixed_by_the_public_key_type() {
        assert_pairing_orientation::<blst_p1_affine, blst_p2_affine>();
        assert_pairing_orientation::<blst_p2_affine, blst_p1_affine>();
    }

    #[test]
    fn empty_pairing_exports_the_miller_identity() {
        let mut pairing = Pairing::new(true, b"miller identity DST");
        let result = pairing.miller_result();
        assert!(unsafe { blst_fp12_is_equal(&result, blst_fp12_one()) });
    }

    #[test]
    fn explicit_none_omits_the_signature_contribution() {
        let (g1, _) = generators();
        let mut pairing = Pairing::new(true, b"optional signature DST");
        assert_eq!(
            pairing.aggregate(&g1, true, None, false, b"message", &[],),
            BLST_ERROR::BLST_SUCCESS
        );
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
    fn mul_n_aggregate_rejects_scalars_wider_than_256_bits() {
        let scalar = [0u8; 33];
        assert_mul_n_result(&scalar, 257, BLST_ERROR::BLST_BAD_SCALAR);
        assert_mul_n_result(&scalar, usize::MAX, BLST_ERROR::BLST_BAD_SCALAR);
    }

    #[test]
    fn mul_n_aggregate_accepts_zero_and_exact_scalar_lengths() {
        let full_width = [1u8; 32];
        for (scalar, nbits) in [
            (&[][..], 0),
            (&[1u8][..], 1),
            (&[1u8][..], 7),
            (&[1u8][..], 8),
            (&[1u8, 0u8][..], 9),
            (&full_width[..], 256),
        ] {
            assert_mul_n_result(scalar, nbits, BLST_ERROR::BLST_SUCCESS);
        }
    }

    #[test]
    fn uniq_helper_distinguishes_messages_and_empty_slices() {
        assert!(uniq(&[]));
        assert!(uniq(&[b"one".as_ref()]));
        assert!(uniq(&[b"one".as_ref(), b"two".as_ref(), b"three".as_ref()]));
        assert!(!uniq(&[b"one".as_ref(), b"two".as_ref(), b"one".as_ref()]));
        assert!(!uniq(&[b"".as_ref(), b"empty".as_ref(), b"".as_ref()]));
    }

    #[test]
    fn uniq_ffi_reports_capacity_duplicate_and_invalid_input() {
        let capacity = 2;
        let size = unsafe { blst_uniq_sizeof(capacity) };
        let words = (size + mem::size_of::<u64>() - 1) / mem::size_of::<u64>();
        let mut storage = vec![0u64; words];
        let ctx = storage.as_mut_ptr() as *mut blst_uniq;

        assert!(unsafe { blst_uniq_init(ctx, capacity) });
        assert_eq!(
            unsafe { blst_uniq_insert(ctx, b"one".as_ptr(), 3) },
            BLST_UNIQ_RESULT::BLST_UNIQ_INSERTED
        );
        assert_eq!(
            unsafe { blst_uniq_insert(ctx, b"two".as_ptr(), 3) },
            BLST_UNIQ_RESULT::BLST_UNIQ_INSERTED
        );
        assert_eq!(
            unsafe { blst_uniq_insert(ctx, b"one".as_ptr(), 3) },
            BLST_UNIQ_RESULT::BLST_UNIQ_DUPLICATE
        );
        assert_eq!(
            unsafe { blst_uniq_insert(ctx, b"three".as_ptr(), 5) },
            BLST_UNIQ_RESULT::BLST_UNIQ_CAPACITY_EXCEEDED
        );
        assert_eq!(
            unsafe { blst_uniq_insert(ctx, ptr::null(), 1) },
            BLST_UNIQ_RESULT::BLST_UNIQ_INVALID_INPUT
        );
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

    let size = unsafe { blst_uniq_sizeof(n_elems) };
    let align = unsafe { blst_uniq_alignof() };
    let word_size = mem::size_of::<u64>();
    let words = match size.checked_add(word_size - 1) {
        Some(rounded) if size != 0 && align <= mem::align_of::<u64>() => {
            rounded / word_size
        }
        _ => return false,
    };
    let mut v: Vec<u64> = vec![0; words];
    let ctx = v.as_mut_ptr() as *mut blst_uniq;

    if !unsafe { blst_uniq_init(ctx, n_elems) } {
        return false;
    }

    for msg in msgs.iter() {
        match unsafe { blst_uniq_insert(ctx, msg.as_ptr(), msg.len()) } {
            BLST_UNIQ_RESULT::BLST_UNIQ_INSERTED => {}
            BLST_UNIQ_RESULT::BLST_UNIQ_DUPLICATE
            | BLST_UNIQ_RESULT::BLST_UNIQ_CAPACITY_EXCEEDED
            | BLST_UNIQ_RESULT::BLST_UNIQ_DEPTH_EXCEEDED
            | BLST_UNIQ_RESULT::BLST_UNIQ_INVALID_INPUT => return false,
            _ => return false,
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
