// Copyright Supranational LLC
// Licensed under the Apache License, Version 2.0, see LICENSE for details.
// SPDX-License-Identifier: Apache-2.0

#![doc = include_str!("../docs/overview.md")]
#![cfg_attr(not(blst_std), no_std)]
#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]

extern crate alloc;

#[cfg(all(test, feature = "full-rust-bindings"))]
use alloc::boxed::Box;
#[cfg(all(test, feature = "full-rust-bindings"))]
use alloc::vec;
#[cfg(feature = "full-rust-bindings")]
use alloc::vec::Vec;
#[cfg(feature = "full-rust-bindings")]
use core::mem::transmute;
#[cfg(feature = "full-rust-bindings")]
use core::mem::MaybeUninit;
#[cfg(feature = "full-rust-bindings")]
use core::ptr;
#[cfg(all(feature = "full-rust-bindings", feature = "serde"))]
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use zeroize::Zeroize;

include!("bindings.rs");
mod avx;
pub mod fft;
mod middle_product;
mod pairing;
mod primitives;

pub use avx::{
    avx_p1s_in_g1, avx_p2s_in_g2, blst_avx_enabled, blst_avx_fp12,
    blst_avx_fp12_alignof, blst_avx_fp12_final_exp, blst_avx_fp12_from_scalar,
    blst_avx_fp12_in_group, blst_avx_fp12_is_equal, blst_avx_fp12_mul,
    blst_avx_fp12_one, blst_avx_fp12_sizeof, blst_avx_fp12_sqr, blst_avx_gt,
    blst_avx_gt_mult_pippenger,
    blst_avx_gt_mult_pippenger_bucket_used_scratch_sizeof,
    blst_avx_gt_mult_pippenger_scratch_sizeof, blst_avx_p1, blst_avx_p1_add,
    blst_avx_p1_affine, blst_avx_p1_affine_alignof,
    blst_avx_p1_affine_from_scalar, blst_avx_p1_affine_sizeof,
    blst_avx_p1_alignof, blst_avx_p1_batch_subgroup_check, blst_avx_p1_cneg,
    blst_avx_p1_compress, blst_avx_p1_deserialize, blst_avx_p1_from_affine,
    blst_avx_p1_from_scalar, blst_avx_p1_generator, blst_avx_p1_in_g1,
    blst_avx_p1_is_equal, blst_avx_p1_is_inf, blst_avx_p1_mult,
    blst_avx_p1_mult_pippenger,
    blst_avx_p1_mult_pippenger_bucket_used_scratch_sizeof,
    blst_avx_p1_mult_pippenger_glv_points_scratch_sizeof,
    blst_avx_p1_mult_pippenger_glv_scalars_scratch_sizeof,
    blst_avx_p1_mult_pippenger_glv_window_bucket_used_scratch_sizeof,
    blst_avx_p1_mult_pippenger_glv_window_scratch_sizeof,
    blst_avx_p1_mult_pippenger_glv_window_with_scratch,
    blst_avx_p1_mult_pippenger_glv_with_scratch,
    blst_avx_p1_mult_pippenger_scratch_sizeof, blst_avx_p1_serialize,
    blst_avx_p1_sizeof, blst_avx_p1_to_affine, blst_avx_p1_uncompress,
    blst_avx_p1s_in_g1, blst_avx_p1s_mult_pippenger,
    blst_avx_p1s_mult_pippenger_scratch_sizeof, blst_avx_p2, blst_avx_p2_add,
    blst_avx_p2_affine, blst_avx_p2_affine_alignof,
    blst_avx_p2_affine_from_scalar, blst_avx_p2_affine_sizeof,
    blst_avx_p2_alignof, blst_avx_p2_batch_subgroup_check, blst_avx_p2_cneg,
    blst_avx_p2_compress, blst_avx_p2_deserialize, blst_avx_p2_from_affine,
    blst_avx_p2_from_scalar, blst_avx_p2_generator, blst_avx_p2_in_g2,
    blst_avx_p2_is_equal, blst_avx_p2_is_inf, blst_avx_p2_mult,
    blst_avx_p2_mult_pippenger,
    blst_avx_p2_mult_pippenger_bucket_used_scratch_sizeof,
    blst_avx_p2_mult_pippenger_scratch_sizeof, blst_avx_p2_serialize,
    blst_avx_p2_sizeof, blst_avx_p2_to_affine, blst_avx_p2_uncompress,
    blst_avx_p2s_in_g2, blst_avx_p2s_mult_pippenger,
    blst_avx_p2s_mult_pippenger_scratch_sizeof, blst_avx_pairing,
    blst_avx_random_bytes_fn, blst_fp12_from_avx, blst_p1_affine_from_avx,
    blst_p1_from_avx, blst_p2_affine_from_avx, blst_p2_from_avx, gt_msm,
    AvxMsmError,
};
pub use middle_product::{
    blst_middle_product_precomputed_b, combine, middle_product,
    precompute_middle_product_b, update, MiddleProductError,
    MiddleProductPrecomputedB,
};
#[cfg(blst_std)]
pub use pairing::print_bytes;
#[cfg(feature = "full-rust-bindings")]
pub use pairing::{uniq, Pairing, PairingPublicKey, PairingSignature};

#[cfg(feature = "full-rust-bindings")]
macro_rules! sig_variant_impl {
    (
        $name:expr,
        $pk:ty,
        $pk_aff:ty,
        $sig:ty,
        $sig_aff:ty,
        $sk_to_pk:ident,
        $hash_or_encode:expr,
        $hash_or_encode_to:ident,
        $sign:ident,
        $pk_eq:ident,
        $sig_eq:ident,
        $verify:ident,
        $pk_in_group:ident,
        $pk_to_aff:ident,
        $pk_from_aff:ident,
        $pk_ser:ident,
        $pk_comp:ident,
        $pk_deser:ident,
        $pk_uncomp:ident,
        $pk_comp_size:expr,
        $pk_ser_size:expr,
        $sig_in_group:ident,
        $sig_to_aff:ident,
        $sig_from_aff:ident,
        $sig_ser:ident,
        $sig_comp:ident,
        $sig_deser:ident,
        $sig_uncomp:ident,
        $sig_comp_size:expr,
        $sig_ser_size:expr,
        $pk_add_or_dbl:ident,
        $pk_add_or_dbl_aff:ident,
        $pk_cneg:ident,
        $sig_add_or_dbl:ident,
        $sig_add_or_dbl_aff:ident,
        $pk_is_inf:ident,
        $sig_is_inf:ident,
        $sig_aggr_in_group:ident,
    ) => {
        /// Secret Key
        #[repr(transparent)]
        #[derive(Default, Debug, Clone, Zeroize)]
        #[zeroize(drop)]
        pub struct SecretKey {
            value: blst_scalar,
        }

        impl SecretKey {
            /// Deterministically derives a canonical nonzero secret key with
            /// the upstream-compatible default-v4 keygen scheme.
            ///
            /// `key_info` may be empty. The native operation is allocation-free
            /// and the resulting secret key is erased when dropped.
            ///
            /// # Errors
            ///
            /// Returns [`BLST_ERROR::BLST_BAD_ENCODING`] when `ikm` contains
            /// fewer than 32 bytes.
            pub fn key_gen(
                ikm: &[u8],
                key_info: &[u8],
            ) -> Result<Self, BLST_ERROR> {
                if ikm.len() < 32 {
                    return Err(BLST_ERROR::BLST_BAD_ENCODING);
                }
                let mut sk = SecretKey::default();
                unsafe {
                    blst_keygen(
                        &mut sk.value,
                        ikm.as_ptr(),
                        ikm.len(),
                        key_info.as_ptr(),
                        key_info.len(),
                    );
                }
                Ok(sk)
            }

            /// Retains upstream BLST's legacy v3 derivation.
            ///
            /// Unlike v4 and later, v3 does not retry the negligible event in
            /// which reduction produces zero. This method is retained for
            /// compatibility; new protocols should use [`Self::key_gen`].
            ///
            /// # Errors
            ///
            /// Returns [`BLST_ERROR::BLST_BAD_ENCODING`] when `ikm` contains
            /// fewer than 32 bytes.
            pub fn key_gen_v3(
                ikm: &[u8],
                key_info: &[u8],
            ) -> Result<Self, BLST_ERROR> {
                if ikm.len() < 32 {
                    return Err(BLST_ERROR::BLST_BAD_ENCODING);
                }
                let mut sk = SecretKey::default();
                unsafe {
                    blst_keygen_v3(
                        &mut sk.value,
                        ikm.as_ptr(),
                        ikm.len(),
                        key_info.as_ptr(),
                        key_info.len(),
                    );
                }
                Ok(sk)
            }

            /// Derives a canonical nonzero secret key with the v4.5 interface.
            ///
            /// `salt` and `info` may be empty. Rust slices always select the
            /// explicit-salt form; use [`Self::key_gen`] for the default salt.
            ///
            /// # Errors
            ///
            /// Returns [`BLST_ERROR::BLST_BAD_ENCODING`] when `ikm` contains
            /// fewer than 32 bytes.
            pub fn key_gen_v4_5(
                ikm: &[u8],
                salt: &[u8],
                info: &[u8],
            ) -> Result<Self, BLST_ERROR> {
                if ikm.len() < 32 {
                    return Err(BLST_ERROR::BLST_BAD_ENCODING);
                }
                let mut sk = SecretKey::default();
                unsafe {
                    blst_keygen_v4_5(
                        &mut sk.value,
                        ikm.as_ptr(),
                        ikm.len(),
                        salt.as_ptr(),
                        salt.len(),
                        info.as_ptr(),
                        info.len(),
                    );
                }
                Ok(sk)
            }

            /// Derives a canonical nonzero secret key with the v5 interface.
            ///
            /// Empty `salt` and `info` slices are valid. The safe wrapper always
            /// supplies the non-null salt pointer required by the C operation.
            ///
            /// # Errors
            ///
            /// Returns [`BLST_ERROR::BLST_BAD_ENCODING`] when `ikm` contains
            /// fewer than 32 bytes.
            pub fn key_gen_v5(
                ikm: &[u8],
                salt: &[u8],
                info: &[u8],
            ) -> Result<Self, BLST_ERROR> {
                if ikm.len() < 32 {
                    return Err(BLST_ERROR::BLST_BAD_ENCODING);
                }
                let mut sk = SecretKey::default();
                unsafe {
                    blst_keygen_v5(
                        &mut sk.value,
                        ikm.as_ptr(),
                        ikm.len(),
                        salt.as_ptr(),
                        salt.len(),
                        info.as_ptr(),
                        info.len(),
                    );
                }
                Ok(sk)
            }

            /// Derives a canonical nonzero EIP-2333 master secret key.
            ///
            /// # Errors
            ///
            /// Returns [`BLST_ERROR::BLST_BAD_ENCODING`] when `ikm` contains
            /// fewer than 32 bytes.
            pub fn derive_master_eip2333(
                ikm: &[u8],
            ) -> Result<Self, BLST_ERROR> {
                if ikm.len() < 32 {
                    return Err(BLST_ERROR::BLST_BAD_ENCODING);
                }
                let mut sk = SecretKey::default();
                unsafe {
                    blst_derive_master_eip2333(
                        &mut sk.value,
                        ikm.as_ptr(),
                        ikm.len(),
                    );
                }
                Ok(sk)
            }

            /// Derives the canonical nonzero EIP-2333 child key at
            /// `child_index`. `self` is the validated parent secret key.
            #[must_use]
            pub fn derive_child_eip2333(&self, child_index: u32) -> Self {
                let mut sk = SecretKey::default();
                unsafe {
                    blst_derive_child_eip2333(
                        &mut sk.value,
                        &self.value,
                        child_index,
                    );
                }
                sk
            }

            /// Derives the prime-order public key corresponding to this
            /// canonical nonzero secret key.
            #[must_use]
            pub fn sk_to_pk(&self) -> PublicKey {
                // TODO - would the user like the serialized/compressed pk as well?
                let mut pk_aff = PublicKey::default();
                //let mut pk_ser = [0u8; $pk_ser_size];

                unsafe {
                    $sk_to_pk(
                        //pk_ser.as_mut_ptr(),
                        ptr::null_mut(),
                        &mut pk_aff.point,
                        &self.value,
                    );
                }
                pk_aff
            }

            /// Hashes `msg` to the signature subgroup with `dst` and `aug`,
            /// then returns its constant-time scalar product by this canonical
            /// nonzero secret key.
            ///
            /// The hash-to-curve step establishes the subgroup precondition
            /// required by the native signing operation.
            #[must_use]
            pub fn sign(
                &self,
                msg: &[u8],
                dst: &[u8],
                aug: &[u8],
            ) -> Signature {
                // TODO - would the user like the serialized/compressed sig as well?
                let mut q = <$sig>::default();
                let mut sig_aff = <$sig_aff>::default();
                //let mut sig_ser = [0u8; $sig_ser_size];
                unsafe {
                    $hash_or_encode_to(
                        &mut q,
                        msg.as_ptr(),
                        msg.len(),
                        dst.as_ptr(),
                        dst.len(),
                        aug.as_ptr(),
                        aug.len(),
                    );
                    $sign(ptr::null_mut(), &mut sig_aff, &q, &self.value);
                }
                Signature { point: sig_aff }
            }

            // TODO - formally speaking application is entitled to have
            // ultimate control over secret key storage, which means that
            // corresponding serialization/deserialization subroutines
            // should accept reference to where to store the result, as
            // opposite to returning one.

            // serialize
            pub fn serialize(&self) -> [u8; 32] {
                let mut sk_out = [0; 32];
                unsafe {
                    blst_bendian_from_scalar(sk_out.as_mut_ptr(), &self.value);
                }
                sk_out
            }

            // deserialize
            pub fn deserialize(sk_in: &[u8]) -> Result<Self, BLST_ERROR> {
                let mut sk = blst_scalar::default();
                if sk_in.len() != 32 {
                    return Err(BLST_ERROR::BLST_BAD_ENCODING);
                }
                unsafe {
                    blst_scalar_from_bendian(&mut sk, sk_in.as_ptr());
                    if !blst_sk_check(&sk) {
                        return Err(BLST_ERROR::BLST_BAD_ENCODING);
                    }
                }
                Ok(Self { value: sk })
            }

            pub fn to_bytes(&self) -> [u8; 32] {
                SecretKey::serialize(&self)
            }

            pub fn from_bytes(sk_in: &[u8]) -> Result<Self, BLST_ERROR> {
                SecretKey::deserialize(sk_in)
            }
        }

        #[cfg(feature = "serde-secret")]
        impl Serialize for SecretKey {
            fn serialize<S: Serializer>(
                &self,
                ser: S,
            ) -> Result<S::Ok, S::Error> {
                let bytes = zeroize::Zeroizing::new(self.serialize());
                ser.serialize_bytes(bytes.as_ref())
            }
        }

        #[cfg(feature = "serde-secret")]
        impl<'de> Deserialize<'de> for SecretKey {
            fn deserialize<D: Deserializer<'de>>(
                deser: D,
            ) -> Result<Self, D::Error> {
                let bytes: &[u8] = Deserialize::deserialize(deser)?;
                Self::deserialize(bytes).map_err(|e| {
                    <D::Error as serde::de::Error>::custom(alloc::format!(
                        "{:?}", e
                    ))
                })
            }
        }

        // From<by-value> traits are not provided to discourage duplication
        // of the secret key material.
        impl<'a> From<&'a SecretKey> for &'a blst_scalar {
            fn from(sk: &'a SecretKey) -> Self {
                unsafe { transmute::<&SecretKey, Self>(sk) }
            }
        }

        impl<'a> core::convert::TryFrom<&'a blst_scalar> for &'a SecretKey {
            type Error = BLST_ERROR;

            fn try_from(sk: &'a blst_scalar) -> Result<Self, Self::Error> {
                unsafe {
                    if !blst_sk_check(sk) {
                        return Err(BLST_ERROR::BLST_BAD_ENCODING);
                    }
                    Ok(transmute::<&blst_scalar, Self>(sk))
                }
            }
        }

        #[repr(transparent)]
        #[derive(Default, Debug, Clone, Copy)]
        /// An affine public-key representation.
        ///
        /// This type does not encode validation state. Decoding establishes a
        /// canonical on-curve point only; call [`PublicKey::validate`] (or use
        /// [`PublicKey::key_validate`]) before passing a value to operations
        /// whose contract requires a prime-order subgroup public key.
        pub struct PublicKey {
            point: $pk_aff,
        }

        impl PublicKey {
            // Core operations

            /// Validates prime-order subgroup membership and rejects infinity.
            ///
            /// Success establishes the public-key precondition required by
            /// signing-scheme verification and aggregation operations.
            pub fn validate(&self) -> Result<(), BLST_ERROR> {
                unsafe {
                    if $pk_is_inf(&self.point) {
                        return Err(BLST_ERROR::BLST_PK_IS_INFINITY);
                    }
                    if !$pk_in_group(&self.point) {
                        return Err(BLST_ERROR::BLST_POINT_NOT_IN_GROUP);
                    }
                }
                Ok(())
            }

            /// Decodes and fully validates a non-infinity subgroup public key.
            pub fn key_validate(key: &[u8]) -> Result<Self, BLST_ERROR> {
                let pk = PublicKey::from_bytes(key)?;
                pk.validate()?;
                Ok(pk)
            }

            pub fn from_aggregate(agg_pk: &AggregatePublicKey) -> Self {
                let mut pk_aff = <$pk_aff>::default();
                unsafe {
                    $pk_to_aff(&mut pk_aff, &agg_pk.point);
                }
                Self { point: pk_aff }
            }

            // Serdes

            /// Compresses an on-curve point. This does not validate subgroup membership.
            pub fn compress(&self) -> [u8; $pk_comp_size] {
                let mut pk_comp = [0u8; $pk_comp_size];
                unsafe {
                    $pk_comp(pk_comp.as_mut_ptr(), &self.point);
                }
                pk_comp
            }

            /// Serializes an on-curve point. This does not validate subgroup membership.
            pub fn serialize(&self) -> [u8; $pk_ser_size] {
                let mut pk_out = [0u8; $pk_ser_size];
                unsafe {
                    $pk_ser(pk_out.as_mut_ptr(), &self.point);
                }
                pk_out
            }

            /// Decodes a canonical compressed public-key-sized point.
            ///
            /// Success establishes on-curve status only, accepts infinity, and
            /// does not establish subgroup membership. Call [`Self::validate`]
            /// before subgroup-only operations unless membership is established
            /// by an approved batch-validation protocol.
            pub fn uncompress(pk_comp: &[u8]) -> Result<Self, BLST_ERROR> {
                if pk_comp.len() == $pk_comp_size && (pk_comp[0] & 0x80) != 0 {
                    let mut pk = <$pk_aff>::default();
                    let err = unsafe { $pk_uncomp(&mut pk, pk_comp.as_ptr()) };
                    if err != BLST_ERROR::BLST_SUCCESS {
                        return Err(err);
                    }
                    Ok(Self { point: pk })
                } else {
                    Err(BLST_ERROR::BLST_BAD_ENCODING)
                }
            }

            /// Decodes canonical compressed or uncompressed point input.
            ///
            /// Success establishes on-curve status only, accepts infinity, and
            /// does not establish subgroup membership. The returned Rust type
            /// does not record validation state; call [`Self::validate`] before
            /// subgroup-only operations unless an approved batch validator is used.
            pub fn deserialize(pk_in: &[u8]) -> Result<Self, BLST_ERROR> {
                if (pk_in.len() == $pk_ser_size && (pk_in[0] & 0x80) == 0)
                    || (pk_in.len() == $pk_comp_size && (pk_in[0] & 0x80) != 0)
                {
                    let mut pk = <$pk_aff>::default();
                    let err = unsafe { $pk_deser(&mut pk, pk_in.as_ptr()) };
                    if err != BLST_ERROR::BLST_SUCCESS {
                        return Err(err);
                    }
                    Ok(Self { point: pk })
                } else {
                    Err(BLST_ERROR::BLST_BAD_ENCODING)
                }
            }

            /// Alias for [`Self::deserialize`]; it does not perform subgroup validation.
            pub fn from_bytes(pk_in: &[u8]) -> Result<Self, BLST_ERROR> {
                PublicKey::deserialize(pk_in)
            }

            pub fn to_bytes(&self) -> [u8; $pk_comp_size] {
                self.compress()
            }
        }

        // Trait for equality comparisons which are equivalence relations.
        //
        // This means, that in addition to a == b and a != b being strict
        // inverses, the equality must be reflexive, symmetric and transitive.
        impl Eq for PublicKey {}

        impl PartialEq for PublicKey {
            fn eq(&self, other: &Self) -> bool {
                unsafe { $pk_eq(&self.point, &other.point) }
            }
        }

        #[cfg(feature = "serde")]
        impl Serialize for PublicKey {
            fn serialize<S: Serializer>(
                &self,
                ser: S,
            ) -> Result<S::Ok, S::Error> {
                ser.serialize_bytes(&self.serialize())
            }
        }

        #[cfg(feature = "serde")]
        /// Serde decoding has the same on-curve-only postcondition as
        /// [`PublicKey::deserialize`]; it does not establish subgroup membership.
        impl<'de> Deserialize<'de> for PublicKey {
            fn deserialize<D: Deserializer<'de>>(
                deser: D,
            ) -> Result<Self, D::Error> {
                let bytes: &[u8] = Deserialize::deserialize(deser)?;
                Self::deserialize(&bytes).map_err(|e| {
                    <D::Error as serde::de::Error>::custom(alloc::format!(
                        "{:?}", e
                    ))
                })
            }
        }

        impl From<PublicKey> for $pk_aff {
            fn from(pk: PublicKey) -> Self {
                pk.point
            }
        }

        impl<'a> From<&'a PublicKey> for &'a $pk_aff {
            fn from(pk: &'a PublicKey) -> Self {
                &pk.point
            }
        }

        impl From<$pk_aff> for PublicKey {
            fn from(point: $pk_aff) -> Self {
                Self { point }
            }
        }

        #[repr(transparent)]
        #[derive(Debug, Clone, Copy)]
        pub struct AggregatePublicKey {
            point: $pk,
        }

        impl AggregatePublicKey {
            pub fn from_public_key(pk: &PublicKey) -> Self {
                let mut agg_pk = <$pk>::default();
                unsafe {
                    $pk_from_aff(&mut agg_pk, &pk.point);
                }
                Self { point: agg_pk }
            }

            pub fn to_public_key(&self) -> PublicKey {
                let mut pk = <$pk_aff>::default();
                unsafe {
                    $pk_to_aff(&mut pk, &self.point);
                }
                PublicKey { point: pk }
            }

            // Aggregate
            pub fn aggregate(
                pks: &[&PublicKey],
                pks_validate: bool,
            ) -> Result<Self, BLST_ERROR> {
                if pks.len() == 0 {
                    return Err(BLST_ERROR::BLST_AGGR_TYPE_MISMATCH);
                }
                if pks_validate {
                    pks[0].validate()?;
                }
                let mut agg_pk = AggregatePublicKey::from_public_key(pks[0]);
                for s in pks.iter().skip(1) {
                    if pks_validate {
                        s.validate()?;
                    }
                    unsafe {
                        $pk_add_or_dbl_aff(
                            &mut agg_pk.point,
                            &agg_pk.point,
                            &s.point,
                        );
                    }
                }
                Ok(agg_pk)
            }

            pub fn aggregate_with_randomness(
                pks: &[PublicKey],
                randomness: &[u8],
                nbits: usize,
                pks_groupcheck: bool,
            ) -> Result<Self, BLST_ERROR> {
                if pks.len() == 0 {
                    return Err(BLST_ERROR::BLST_AGGR_TYPE_MISMATCH);
                }
                if pks_groupcheck {
                    pks.validate()?;
                }
                Ok(pks.mult(randomness, nbits))
            }

            pub fn aggregate_serialized(
                pks: &[&[u8]],
                pks_validate: bool,
            ) -> Result<Self, BLST_ERROR> {
                // TODO - threading
                if pks.len() == 0 {
                    return Err(BLST_ERROR::BLST_AGGR_TYPE_MISMATCH);
                }
                let mut pk = if pks_validate {
                    PublicKey::key_validate(pks[0])?
                } else {
                    PublicKey::from_bytes(pks[0])?
                };
                let mut agg_pk = AggregatePublicKey::from_public_key(&pk);
                for s in pks.iter().skip(1) {
                    pk = if pks_validate {
                        PublicKey::key_validate(s)?
                    } else {
                        PublicKey::from_bytes(s)?
                    };
                    unsafe {
                        $pk_add_or_dbl_aff(
                            &mut agg_pk.point,
                            &agg_pk.point,
                            &pk.point,
                        );
                    }
                }
                Ok(agg_pk)
            }

            pub fn add_aggregate(&mut self, agg_pk: &AggregatePublicKey) {
                unsafe {
                    $pk_add_or_dbl(&mut self.point, &self.point, &agg_pk.point);
                }
            }

            /// Subtracts `agg_pk` from this aggregate without modifying it.
            pub fn sub_aggregate(&mut self, agg_pk: &AggregatePublicKey) {
                unsafe {
                    let mut negated = agg_pk.clone();
                    $pk_cneg(&mut negated.point, true);
                    $pk_add_or_dbl(
                        &mut self.point,
                        &self.point,
                        &negated.point,
                    );
                }
            }

            pub fn add_public_key(
                &mut self,
                pk: &PublicKey,
                pk_validate: bool,
            ) -> Result<(), BLST_ERROR> {
                if pk_validate {
                    pk.validate()?;
                }
                unsafe {
                    $pk_add_or_dbl_aff(&mut self.point, &self.point, &pk.point);
                }
                Ok(())
            }
        }

        impl From<AggregatePublicKey> for $pk {
            fn from(pk: AggregatePublicKey) -> Self {
                pk.point
            }
        }

        impl<'a> From<&'a AggregatePublicKey> for &'a $pk {
            fn from(pk: &'a AggregatePublicKey) -> Self {
                &pk.point
            }
        }

        impl From<$pk> for AggregatePublicKey {
            fn from(point: $pk) -> Self {
                Self { point }
            }
        }

        #[repr(transparent)]
        #[derive(Debug, Clone, Copy)]
        /// An affine signature representation.
        ///
        /// This type does not encode validation state. Decoding establishes a
        /// canonical on-curve point only; call [`Signature::validate`] (or use
        /// [`Signature::sig_validate`]) before subgroup-only operations.
        pub struct Signature {
            point: $sig_aff,
        }

        impl Signature {
            // sig_infcheck, check for infinity, is a way to avoid going
            // into resource-consuming verification. Passing 'false' is
            // always cryptographically safe, but application might want
            // to guard against obviously bogus individual[!] signatures.
            /// Validates prime-order subgroup membership and, when requested,
            /// rejects infinity. Success establishes the signature precondition
            /// required by signing-scheme verification and aggregation.
            pub fn validate(
                &self,
                sig_infcheck: bool,
            ) -> Result<(), BLST_ERROR> {
                unsafe {
                    if sig_infcheck && $sig_is_inf(&self.point) {
                        return Err(BLST_ERROR::BLST_PK_IS_INFINITY);
                    }
                    if !$sig_in_group(&self.point) {
                        return Err(BLST_ERROR::BLST_POINT_NOT_IN_GROUP);
                    }
                }
                Ok(())
            }

            /// Decodes and validates a subgroup signature, optionally rejecting infinity.
            pub fn sig_validate(
                sig: &[u8],
                sig_infcheck: bool,
            ) -> Result<Self, BLST_ERROR> {
                let sig = Signature::from_bytes(sig)?;
                sig.validate(sig_infcheck)?;
                Ok(sig)
            }

            pub fn verify(
                &self,
                sig_groupcheck: bool,
                msg: &[u8],
                dst: &[u8],
                aug: &[u8],
                pk: &PublicKey,
                pk_validate: bool,
            ) -> BLST_ERROR {
                let aug_msg = [aug, msg].concat();
                self.aggregate_verify(
                    sig_groupcheck,
                    &[aug_msg.as_slice()],
                    dst,
                    &[pk],
                    pk_validate,
                )
            }

            pub fn aggregate_verify(
                &self,
                sig_groupcheck: bool,
                msgs: &[&[u8]],
                dst: &[u8],
                pks: &[&PublicKey],
                pks_validate: bool,
            ) -> BLST_ERROR {
                let n_elems = pks.len();
                if n_elems == 0 || msgs.len() != n_elems {
                    return BLST_ERROR::BLST_VERIFY_FAIL;
                }

                let mut pairing = Pairing::new($hash_or_encode, dst);

                let err = pairing.aggregate(
                    &pks[0].point,
                    pks_validate,
                    Some(&self.point),
                    sig_groupcheck,
                    &msgs[0],
                    &[],
                );
                if err != BLST_ERROR::BLST_SUCCESS {
                    return err;
                }

                for i in 1..n_elems {
                    let err = pairing.aggregate(
                        &pks[i].point,
                        pks_validate,
                        None,
                        false,
                        &msgs[i],
                        &[],
                    );
                    if err != BLST_ERROR::BLST_SUCCESS {
                        return err;
                    }
                }

                pairing.commit();

                if pairing.finalverify(None) {
                    BLST_ERROR::BLST_SUCCESS
                } else {
                    BLST_ERROR::BLST_VERIFY_FAIL
                }
            }

            // pks are assumed to be verified for proof of possession,
            // which implies that they are already group-checked
            pub fn fast_aggregate_verify(
                &self,
                sig_groupcheck: bool,
                msg: &[u8],
                dst: &[u8],
                pks: &[&PublicKey],
            ) -> BLST_ERROR {
                let agg_pk = match AggregatePublicKey::aggregate(pks, false) {
                    Ok(agg_sig) => agg_sig,
                    Err(err) => return err,
                };
                let pk = agg_pk.to_public_key();
                self.aggregate_verify(
                    sig_groupcheck,
                    &[msg],
                    dst,
                    &[&pk],
                    false,
                )
            }

            pub fn fast_aggregate_verify_pre_aggregated(
                &self,
                sig_groupcheck: bool,
                msg: &[u8],
                dst: &[u8],
                pk: &PublicKey,
            ) -> BLST_ERROR {
                self.aggregate_verify(sig_groupcheck, &[msg], dst, &[pk], false)
            }

            // https://ethresear.ch/t/fast-verification-of-multiple-bls-signatures/5407
            #[allow(clippy::too_many_arguments)]
            pub fn verify_multiple_aggregate_signatures(
                msgs: &[&[u8]],
                dst: &[u8],
                pks: &[&PublicKey],
                pks_validate: bool,
                sigs: &[&Signature],
                sigs_groupcheck: bool,
                rands: &[blst_scalar],
                rand_bits: usize,
            ) -> BLST_ERROR {
                let n_elems = pks.len();
                if n_elems == 0
                    || msgs.len() != n_elems
                    || sigs.len() != n_elems
                    || rands.len() != n_elems
                {
                    return BLST_ERROR::BLST_VERIFY_FAIL;
                }

                // TODO - check msg uniqueness?

                let mut pairing = Pairing::new($hash_or_encode, dst);

                for i in 0..n_elems {
                    let err = pairing.mul_n_aggregate(
                        &pks[i].point,
                        pks_validate,
                        Some(&sigs[i].point),
                        sigs_groupcheck,
                        &rands[i].b,
                        rand_bits,
                        msgs[i],
                        &[],
                    );
                    if err != BLST_ERROR::BLST_SUCCESS {
                        return err;
                    }
                }

                pairing.commit();

                if pairing.finalverify(None) {
                    BLST_ERROR::BLST_SUCCESS
                } else {
                    BLST_ERROR::BLST_VERIFY_FAIL
                }
            }

            pub fn from_aggregate(agg_sig: &AggregateSignature) -> Self {
                let mut sig_aff = <$sig_aff>::default();
                unsafe {
                    $sig_to_aff(&mut sig_aff, &agg_sig.point);
                }
                Self { point: sig_aff }
            }

            /// Compresses an on-curve point. This does not validate subgroup membership.
            pub fn compress(&self) -> [u8; $sig_comp_size] {
                let mut sig_comp = [0; $sig_comp_size];
                unsafe {
                    $sig_comp(sig_comp.as_mut_ptr(), &self.point);
                }
                sig_comp
            }

            /// Serializes an on-curve point. This does not validate subgroup membership.
            pub fn serialize(&self) -> [u8; $sig_ser_size] {
                let mut sig_out = [0; $sig_ser_size];
                unsafe {
                    $sig_ser(sig_out.as_mut_ptr(), &self.point);
                }
                sig_out
            }

            /// Decodes a canonical compressed signature-sized point.
            ///
            /// Success establishes on-curve status only, accepts infinity, and
            /// does not establish subgroup membership. Call [`Self::validate`]
            /// before subgroup-only operations unless membership is established
            /// by an approved batch-validation protocol.
            pub fn uncompress(sig_comp: &[u8]) -> Result<Self, BLST_ERROR> {
                if sig_comp.len() == $sig_comp_size && (sig_comp[0] & 0x80) != 0
                {
                    let mut sig = <$sig_aff>::default();
                    let err =
                        unsafe { $sig_uncomp(&mut sig, sig_comp.as_ptr()) };
                    if err != BLST_ERROR::BLST_SUCCESS {
                        return Err(err);
                    }
                    Ok(Self { point: sig })
                } else {
                    Err(BLST_ERROR::BLST_BAD_ENCODING)
                }
            }

            /// Decodes canonical compressed or uncompressed point input.
            ///
            /// Success establishes on-curve status only, accepts infinity, and
            /// does not establish subgroup membership. The returned Rust type
            /// does not record validation state; call [`Self::validate`] before
            /// subgroup-only operations unless an approved batch validator is used.
            pub fn deserialize(sig_in: &[u8]) -> Result<Self, BLST_ERROR> {
                if (sig_in.len() == $sig_ser_size && (sig_in[0] & 0x80) == 0)
                    || (sig_in.len() == $sig_comp_size
                        && (sig_in[0] & 0x80) != 0)
                {
                    let mut sig = <$sig_aff>::default();
                    let err = unsafe { $sig_deser(&mut sig, sig_in.as_ptr()) };
                    if err != BLST_ERROR::BLST_SUCCESS {
                        return Err(err);
                    }
                    Ok(Self { point: sig })
                } else {
                    Err(BLST_ERROR::BLST_BAD_ENCODING)
                }
            }

            /// Alias for [`Self::deserialize`]; it does not perform subgroup validation.
            pub fn from_bytes(sig_in: &[u8]) -> Result<Self, BLST_ERROR> {
                Signature::deserialize(sig_in)
            }

            pub fn to_bytes(&self) -> [u8; $sig_comp_size] {
                self.compress()
            }

            /// Returns whether this candidate is on-curve and in the signature subgroup.
            pub fn subgroup_check(&self) -> bool {
                unsafe { $sig_in_group(&self.point) }
            }
        }

        // Trait for equality comparisons which are equivalence relations.
        //
        // This means, that in addition to a == b and a != b being strict
        // inverses, the equality must be reflexive, symmetric and transitive.
        impl Eq for Signature {}

        impl PartialEq for Signature {
            fn eq(&self, other: &Self) -> bool {
                unsafe { $sig_eq(&self.point, &other.point) }
            }
        }

        #[cfg(feature = "serde")]
        impl Serialize for Signature {
            fn serialize<S: Serializer>(
                &self,
                ser: S,
            ) -> Result<S::Ok, S::Error> {
                ser.serialize_bytes(&self.serialize())
            }
        }

        #[cfg(feature = "serde")]
        /// Serde decoding has the same on-curve-only postcondition as
        /// [`Signature::deserialize`]; it does not establish subgroup membership.
        impl<'de> Deserialize<'de> for Signature {
            fn deserialize<D: Deserializer<'de>>(
                deser: D,
            ) -> Result<Self, D::Error> {
                let bytes: &[u8] = Deserialize::deserialize(deser)?;
                Self::deserialize(&bytes).map_err(|e| {
                    <D::Error as serde::de::Error>::custom(alloc::format!(
                        "{:?}", e
                    ))
                })
            }
        }

        impl From<Signature> for $sig_aff {
            fn from(sig: Signature) -> Self {
                sig.point
            }
        }

        impl<'a> From<&'a Signature> for &'a $sig_aff {
            fn from(sig: &'a Signature) -> Self {
                &sig.point
            }
        }

        impl From<$sig_aff> for Signature {
            fn from(point: $sig_aff) -> Self {
                Self { point }
            }
        }

        #[repr(transparent)]
        #[derive(Debug, Clone, Copy)]
        pub struct AggregateSignature {
            point: $sig,
        }

        impl AggregateSignature {
            pub fn validate(&self) -> Result<(), BLST_ERROR> {
                unsafe {
                    if !$sig_aggr_in_group(&self.point) {
                        return Err(BLST_ERROR::BLST_POINT_NOT_IN_GROUP);
                    }
                }
                Ok(())
            }

            pub fn from_signature(sig: &Signature) -> Self {
                let mut agg_sig = <$sig>::default();
                unsafe {
                    $sig_from_aff(&mut agg_sig, &sig.point);
                }
                Self { point: agg_sig }
            }

            pub fn to_signature(&self) -> Signature {
                let mut sig = <$sig_aff>::default();
                unsafe {
                    $sig_to_aff(&mut sig, &self.point);
                }
                Signature { point: sig }
            }

            // Aggregate
            pub fn aggregate(
                sigs: &[&Signature],
                sigs_groupcheck: bool,
            ) -> Result<Self, BLST_ERROR> {
                if sigs.len() == 0 {
                    return Err(BLST_ERROR::BLST_AGGR_TYPE_MISMATCH);
                }
                if sigs_groupcheck {
                    // We can't actually judge if input is individual or
                    // aggregated signature, so we can't enforce infinity
                    // check.
                    sigs[0].validate(false)?;
                }
                let mut agg_sig = AggregateSignature::from_signature(sigs[0]);
                for s in sigs.iter().skip(1) {
                    if sigs_groupcheck {
                        s.validate(false)?;
                    }
                    unsafe {
                        $sig_add_or_dbl_aff(
                            &mut agg_sig.point,
                            &agg_sig.point,
                            &s.point,
                        );
                    }
                }
                Ok(agg_sig)
            }

            pub fn aggregate_with_randomness(
                sigs: &[Signature],
                randomness: &[u8],
                nbits: usize,
                sigs_groupcheck: bool,
            ) -> Result<Self, BLST_ERROR> {
                if sigs.len() == 0 {
                    return Err(BLST_ERROR::BLST_AGGR_TYPE_MISMATCH);
                }
                if sigs_groupcheck {
                    sigs.validate()?;
                }
                Ok(sigs.mult(randomness, nbits))
            }

            pub fn aggregate_serialized(
                sigs: &[&[u8]],
                sigs_groupcheck: bool,
            ) -> Result<Self, BLST_ERROR> {
                // TODO - threading
                if sigs.len() == 0 {
                    return Err(BLST_ERROR::BLST_AGGR_TYPE_MISMATCH);
                }
                let mut sig = if sigs_groupcheck {
                    Signature::sig_validate(sigs[0], false)?
                } else {
                    Signature::from_bytes(sigs[0])?
                };
                let mut agg_sig = AggregateSignature::from_signature(&sig);
                for s in sigs.iter().skip(1) {
                    sig = if sigs_groupcheck {
                        Signature::sig_validate(s, false)?
                    } else {
                        Signature::from_bytes(s)?
                    };
                    unsafe {
                        $sig_add_or_dbl_aff(
                            &mut agg_sig.point,
                            &agg_sig.point,
                            &sig.point,
                        );
                    }
                }
                Ok(agg_sig)
            }

            pub fn add_aggregate(&mut self, agg_sig: &AggregateSignature) {
                unsafe {
                    $sig_add_or_dbl(
                        &mut self.point,
                        &self.point,
                        &agg_sig.point,
                    );
                }
            }

            pub fn add_signature(
                &mut self,
                sig: &Signature,
                sig_groupcheck: bool,
            ) -> Result<(), BLST_ERROR> {
                if sig_groupcheck {
                    sig.validate(false)?;
                }
                unsafe {
                    $sig_add_or_dbl_aff(
                        &mut self.point,
                        &self.point,
                        &sig.point,
                    );
                }
                Ok(())
            }

            pub fn subgroup_check(&self) -> bool {
                unsafe { $sig_aggr_in_group(&self.point) }
            }
        }

        impl From<AggregateSignature> for $sig {
            fn from(sig: AggregateSignature) -> Self {
                sig.point
            }
        }

        impl<'a> From<&'a AggregateSignature> for &'a $sig {
            fn from(sig: &'a AggregateSignature) -> Self {
                &sig.point
            }
        }

        impl From<$sig> for AggregateSignature {
            fn from(point: $sig) -> Self {
                Self { point }
            }
        }

        impl MultiPoint for [PublicKey] {
            type Output = AggregatePublicKey;

            fn mult(&self, scalars: &[u8], nbits: usize) -> Self::Output {
                Self::Output {
                    point: unsafe { transmute::<&[_], &[$pk_aff]>(self) }
                        .mult(scalars, nbits),
                }
            }

            fn add(&self) -> Self::Output {
                Self::Output {
                    point: unsafe { transmute::<&[_], &[$pk_aff]>(self) }.add(),
                }
            }

            fn validate(&self) -> Result<(), BLST_ERROR> {
                unsafe { transmute::<&[_], &[$pk_aff]>(self) }.validate()
            }
        }

        impl MultiPoint for [Signature] {
            type Output = AggregateSignature;

            fn mult(&self, scalars: &[u8], nbits: usize) -> Self::Output {
                Self::Output {
                    point: unsafe { transmute::<&[_], &[$sig_aff]>(self) }
                        .mult(scalars, nbits),
                }
            }

            fn add(&self) -> Self::Output {
                Self::Output {
                    point: unsafe { transmute::<&[_], &[$sig_aff]>(self) }
                        .add(),
                }
            }

            fn validate(&self) -> Result<(), BLST_ERROR> {
                unsafe { transmute::<&[_], &[$sig_aff]>(self) }.validate()
            }
        }

        #[cfg(test)]
        mod tests {
            use super::*;
            use rand_chacha::ChaCha20Rng;
            use rand_core::{RngCore, SeedableRng};

            // Testing only - do not use for production
            pub fn gen_random_key(
                rng: &mut rand_chacha::ChaCha20Rng,
            ) -> SecretKey {
                let mut ikm = [0u8; 32];
                rng.fill_bytes(&mut ikm);

                let mut sk = <blst_scalar>::default();
                unsafe {
                    blst_keygen(&mut sk, ikm.as_ptr(), 32, ptr::null(), 0);
                }
                SecretKey { value: sk }
            }

            #[test]
            fn test_sign_n_verify() {
                let ikm: [u8; 32] = [
                    0x93, 0xad, 0x7e, 0x65, 0xde, 0xad, 0x05, 0x2a, 0x08, 0x3a,
                    0x91, 0x0c, 0x8b, 0x72, 0x85, 0x91, 0x46, 0x4c, 0xca, 0x56,
                    0x60, 0x5b, 0xb0, 0x56, 0xed, 0xfe, 0x2b, 0x60, 0xa6, 0x3c,
                    0x48, 0x99,
                ];

                let sk = SecretKey::key_gen(&ikm, &[]).unwrap();
                let pk = sk.sk_to_pk();

                let dst = b"BLS_SIG_BLS12381G2_XMD:SHA-256_SSWU_RO_NUL_";
                let msg = b"hello foo";
                let sig = sk.sign(msg, dst, &[]);

                let err = sig.verify(true, msg, dst, &[], &pk, true);
                assert_eq!(err, BLST_ERROR::BLST_SUCCESS);
            }

            #[test]
            fn test_key_gen_v5_empty_salt() {
                let expected = [
                    0x19, 0xba, 0x3d, 0xca, 0x3f, 0x17, 0x4c, 0x27, 0x43, 0x2c,
                    0x1b, 0x8a, 0x5f, 0x91, 0xa1, 0xc9, 0xbe, 0xe5, 0x4f, 0x01,
                    0xec, 0xe4, 0x97, 0xd9, 0x40, 0x66, 0xf3, 0xce, 0xc2, 0x86,
                    0xc1, 0x10,
                ];
                let mut ikm = [0u8; 32];
                ikm[0] = 9;

                let first = SecretKey::key_gen_v5(&ikm, &[], &[]).unwrap();
                let second = SecretKey::key_gen_v5(&ikm, &[], &[]).unwrap();

                assert_eq!(first.to_bytes(), expected);
                assert_eq!(first.to_bytes(), second.to_bytes());
            }

            #[test]
            fn test_aggregate() {
                let num_msgs = 10;
                let dst = b"BLS_SIG_BLS12381G2_XMD:SHA-256_SSWU_RO_NUL_";

                let seed = [0u8; 32];
                let mut rng = ChaCha20Rng::from_seed(seed);

                let sks: Vec<_> =
                    (0..num_msgs).map(|_| gen_random_key(&mut rng)).collect();
                let pks =
                    sks.iter().map(|sk| sk.sk_to_pk()).collect::<Vec<_>>();
                let pks_refs: Vec<&PublicKey> =
                    pks.iter().map(|pk| pk).collect();
                let pks_rev: Vec<&PublicKey> =
                    pks.iter().rev().map(|pk| pk).collect();

                let pk_comp = pks[0].compress();
                let pk_uncomp = PublicKey::uncompress(&pk_comp);
                assert_eq!(pk_uncomp.is_ok(), true);

                let mut msgs: Vec<Vec<u8>> = vec![vec![]; num_msgs];
                for i in 0..num_msgs {
                    let msg_len = (rng.next_u64() & 0x3F) + 1;
                    msgs[i] = vec![0u8; msg_len as usize];
                    rng.fill_bytes(&mut msgs[i]);
                }

                let msgs_refs: Vec<&[u8]> =
                    msgs.iter().map(|m| m.as_slice()).collect();

                let sigs = sks
                    .iter()
                    .zip(msgs.iter())
                    .map(|(sk, m)| (sk.sign(m, dst, &[])))
                    .collect::<Vec<Signature>>();

                let mut errs = sigs
                    .iter()
                    .zip(msgs.iter())
                    .zip(pks.iter())
                    .map(|((s, m), pk)| (s.verify(true, m, dst, &[], pk, true)))
                    .collect::<Vec<BLST_ERROR>>();
                assert_eq!(errs, vec![BLST_ERROR::BLST_SUCCESS; num_msgs]);

                // Swap message/public key pairs to create bad signature
                errs = sigs
                    .iter()
                    .zip(msgs.iter())
                    .zip(pks.iter().rev())
                    .map(|((s, m), pk)| (s.verify(true, m, dst, &[], pk, true)))
                    .collect::<Vec<BLST_ERROR>>();
                assert_ne!(errs, vec![BLST_ERROR::BLST_SUCCESS; num_msgs]);

                let sig_refs =
                    sigs.iter().map(|s| s).collect::<Vec<&Signature>>();
                let agg = match AggregateSignature::aggregate(&sig_refs, true) {
                    Ok(agg) => agg,
                    Err(err) => panic!("aggregate failure: {:?}", err),
                };

                let agg_sig = agg.to_signature();
                let mut result = agg_sig
                    .aggregate_verify(false, &msgs_refs, dst, &pks_refs, false);
                assert_eq!(result, BLST_ERROR::BLST_SUCCESS);

                // Swap message/public key pairs to create bad signature
                result = agg_sig
                    .aggregate_verify(false, &msgs_refs, dst, &pks_rev, false);
                assert_ne!(result, BLST_ERROR::BLST_SUCCESS);
            }

            #[test]
            fn test_sub_aggregate() {
                let pks = (1u8..=4)
                    .map(|tag| {
                        let mut ikm = [tag; 32];
                        ikm[0] ^= 0xa5;
                        SecretKey::key_gen(&ikm, &[]).unwrap().sk_to_pk()
                    })
                    .collect::<Vec<_>>();
                let pk_refs = pks.iter().collect::<Vec<_>>();

                let mut aggregate =
                    AggregatePublicKey::aggregate(&pk_refs, false).unwrap();
                let removed =
                    AggregatePublicKey::aggregate(&pk_refs[1..], false)
                        .unwrap();
                let removed_before = removed.to_public_key();

                aggregate.sub_aggregate(&removed);
                assert_eq!(aggregate.to_public_key(), pks[0]);
                assert_eq!(removed.to_public_key(), removed_before);

                let mut identity = removed;
                identity.sub_aggregate(&removed);
                let identity = identity.to_public_key();
                assert!(unsafe { $pk_is_inf(&identity.point) });
            }

            #[test]
            fn test_multiple_agg_sigs() {
                let dst = b"BLS_SIG_BLS12381G2_XMD:SHA-256_SSWU_RO_POP_";
                let num_pks_per_sig = 10;
                let num_sigs = 10;

                let seed = [0u8; 32];
                let mut rng = ChaCha20Rng::from_seed(seed);

                let mut msgs: Vec<Vec<u8>> = vec![vec![]; num_sigs];
                let mut sigs: Vec<Signature> = Vec::with_capacity(num_sigs);
                let mut pks: Vec<PublicKey> = Vec::with_capacity(num_sigs);
                let mut rands: Vec<blst_scalar> = Vec::with_capacity(num_sigs);
                for i in 0..num_sigs {
                    // Create public keys
                    let sks_i: Vec<_> = (0..num_pks_per_sig)
                        .map(|_| gen_random_key(&mut rng))
                        .collect();

                    let pks_i = sks_i
                        .iter()
                        .map(|sk| sk.sk_to_pk())
                        .collect::<Vec<_>>();
                    let pks_refs_i: Vec<&PublicKey> =
                        pks_i.iter().map(|pk| pk).collect();

                    // Create random message for pks to all sign
                    let msg_len = (rng.next_u64() & 0x3F) + 1;
                    msgs[i] = vec![0u8; msg_len as usize];
                    rng.fill_bytes(&mut msgs[i]);

                    // Generate signature for each key pair
                    let sigs_i = sks_i
                        .iter()
                        .map(|sk| sk.sign(&msgs[i], dst, &[]))
                        .collect::<Vec<Signature>>();

                    // Test each current single signature
                    let errs = sigs_i
                        .iter()
                        .zip(pks_i.iter())
                        .map(|(s, pk)| {
                            (s.verify(true, &msgs[i], dst, &[], pk, true))
                        })
                        .collect::<Vec<BLST_ERROR>>();
                    assert_eq!(
                        errs,
                        vec![BLST_ERROR::BLST_SUCCESS; num_pks_per_sig]
                    );

                    let sig_refs_i =
                        sigs_i.iter().map(|s| s).collect::<Vec<&Signature>>();
                    let agg_i =
                        match AggregateSignature::aggregate(&sig_refs_i, false)
                        {
                            Ok(agg_i) => agg_i,
                            Err(err) => panic!("aggregate failure: {:?}", err),
                        };

                    // Test current aggregate signature
                    sigs.push(agg_i.to_signature());
                    let mut result = sigs[i].fast_aggregate_verify(
                        false,
                        &msgs[i],
                        dst,
                        &pks_refs_i,
                    );
                    assert_eq!(result, BLST_ERROR::BLST_SUCCESS);

                    // negative test
                    if i != 0 {
                        result = sigs[i - 1].fast_aggregate_verify(
                            false,
                            &msgs[i],
                            dst,
                            &pks_refs_i,
                        );
                        assert_ne!(result, BLST_ERROR::BLST_SUCCESS);
                    }

                    // aggregate public keys and push into vec
                    let agg_pk_i =
                        match AggregatePublicKey::aggregate(&pks_refs_i, false)
                        {
                            Ok(agg_pk_i) => agg_pk_i,
                            Err(err) => panic!("aggregate failure: {:?}", err),
                        };
                    pks.push(agg_pk_i.to_public_key());

                    // Test current aggregate signature with aggregated pks
                    result = sigs[i].fast_aggregate_verify_pre_aggregated(
                        false, &msgs[i], dst, &pks[i],
                    );
                    assert_eq!(result, BLST_ERROR::BLST_SUCCESS);

                    // negative test
                    if i != 0 {
                        result = sigs[i - 1]
                            .fast_aggregate_verify_pre_aggregated(
                                false, &msgs[i], dst, &pks[i],
                            );
                        assert_ne!(result, BLST_ERROR::BLST_SUCCESS);
                    }

                    // create random values
                    let mut vals = [0u64; 4];
                    vals[0] = rng.next_u64();
                    while vals[0] == 0 {
                        // Reject zero as it is used for multiplication.
                        vals[0] = rng.next_u64();
                    }
                    let mut rand_i = MaybeUninit::<blst_scalar>::uninit();
                    unsafe {
                        blst_scalar_from_uint64(
                            rand_i.as_mut_ptr(),
                            vals.as_ptr(),
                        );
                        rands.push(rand_i.assume_init());
                    }
                }

                let msgs_refs: Vec<&[u8]> =
                    msgs.iter().map(|m| m.as_slice()).collect();
                let sig_refs =
                    sigs.iter().map(|s| s).collect::<Vec<&Signature>>();
                let pks_refs: Vec<&PublicKey> =
                    pks.iter().map(|pk| pk).collect();

                let msgs_rev: Vec<&[u8]> =
                    msgs.iter().rev().map(|m| m.as_slice()).collect();
                let sig_rev =
                    sigs.iter().rev().map(|s| s).collect::<Vec<&Signature>>();
                let pks_rev: Vec<&PublicKey> =
                    pks.iter().rev().map(|pk| pk).collect();

                let mut result =
                    Signature::verify_multiple_aggregate_signatures(
                        &msgs_refs, dst, &pks_refs, false, &sig_refs, true,
                        &rands, 64,
                    );
                assert_eq!(result, BLST_ERROR::BLST_SUCCESS);

                // negative tests (use reverse msgs, pks, and sigs)
                result = Signature::verify_multiple_aggregate_signatures(
                    &msgs_rev, dst, &pks_refs, false, &sig_refs, true, &rands,
                    64,
                );
                assert_ne!(result, BLST_ERROR::BLST_SUCCESS);

                result = Signature::verify_multiple_aggregate_signatures(
                    &msgs_refs, dst, &pks_rev, false, &sig_refs, true, &rands,
                    64,
                );
                assert_ne!(result, BLST_ERROR::BLST_SUCCESS);

                result = Signature::verify_multiple_aggregate_signatures(
                    &msgs_refs, dst, &pks_refs, false, &sig_rev, true, &rands,
                    64,
                );
                assert_ne!(result, BLST_ERROR::BLST_SUCCESS);
            }

            #[test]
            fn test_serialization() {
                let seed = [0u8; 32];
                let mut rng = ChaCha20Rng::from_seed(seed);

                let sk = gen_random_key(&mut rng);
                let sk2 = gen_random_key(&mut rng);

                let pk = sk.sk_to_pk();
                let pk_comp = pk.compress();
                let pk_ser = pk.serialize();

                let pk_uncomp = PublicKey::uncompress(&pk_comp);
                assert_eq!(pk_uncomp.is_ok(), true);
                assert_eq!(pk_uncomp.unwrap(), pk);

                let pk_deser = PublicKey::deserialize(&pk_ser);
                assert_eq!(pk_deser.is_ok(), true);
                assert_eq!(pk_deser.unwrap(), pk);

                let pk2 = sk2.sk_to_pk();
                let pk_comp2 = pk2.compress();
                let pk_ser2 = pk2.serialize();

                let pk_uncomp2 = PublicKey::uncompress(&pk_comp2);
                assert_eq!(pk_uncomp2.is_ok(), true);
                assert_eq!(pk_uncomp2.unwrap(), pk2);

                let pk_deser2 = PublicKey::deserialize(&pk_ser2);
                assert_eq!(pk_deser2.is_ok(), true);
                assert_eq!(pk_deser2.unwrap(), pk2);

                assert_ne!(pk, pk2);
                assert_ne!(pk_uncomp.unwrap(), pk2);
                assert_ne!(pk_deser.unwrap(), pk2);
                assert_ne!(pk_uncomp2.unwrap(), pk);
                assert_ne!(pk_deser2.unwrap(), pk);
            }

            #[test]
            fn test_point_decode_requires_separate_subgroup_validation() {
                fn decode_hex<const N: usize>(encoded: &str) -> [u8; N] {
                    assert_eq!(encoded.len(), 2 * N);
                    let mut out = [0u8; N];
                    for (index, value) in out.iter_mut().enumerate() {
                        *value = u8::from_str_radix(
                            &encoded[2 * index..2 * index + 2],
                            16,
                        )
                        .unwrap();
                    }
                    out
                }

                let g1_order_11 = decode_hex::<48>(
                    "80803f0d09fec09a95f2ee7495323c15c162270c7cceaffa\
                     8566e941c66bcf206e72955d58b3b32e564de3209d672ca5",
                );
                let g2_order_13 = decode_hex::<96>(
                    "808514e6f6c41ddb685cd8e351d7cb2ee891d54d5850d42c\
                     440365b15dd8ab06b53da9ebb22bc31cea2234d30499b09f\
                     01b290cae398ecc1fda3d4beb248af58bd351c3d99a5b6c7\
                     70fdb54bf46b123261b55cf0762aeef2341e35e90608fc31",
                );

                let pk_bytes: &[u8] = if $pk_comp_size == g1_order_11.len() {
                    &g1_order_11
                } else {
                    &g2_order_13
                };
                let pk = PublicKey::uncompress(pk_bytes)
                    .expect("on-curve non-subgroup public-key decode");
                assert_eq!(
                    pk.validate(),
                    Err(BLST_ERROR::BLST_POINT_NOT_IN_GROUP)
                );
                assert_eq!(
                    PublicKey::key_validate(pk_bytes).unwrap_err(),
                    BLST_ERROR::BLST_POINT_NOT_IN_GROUP
                );
                let pk_serialized = pk.serialize();
                let pk_deserialized = PublicKey::deserialize(&pk_serialized)
                    .expect("uncompressed on-curve public-key decode");
                assert_eq!(
                    pk_deserialized.validate(),
                    Err(BLST_ERROR::BLST_POINT_NOT_IN_GROUP)
                );
                let pk_from_bytes = PublicKey::from_bytes(pk_bytes)
                    .expect("from_bytes preserves the decoder contract");
                assert_eq!(
                    pk_from_bytes.validate(),
                    Err(BLST_ERROR::BLST_POINT_NOT_IN_GROUP)
                );

                let sig_bytes: &[u8] = if $sig_comp_size == g1_order_11.len() {
                    &g1_order_11
                } else {
                    &g2_order_13
                };
                let sig = Signature::uncompress(sig_bytes)
                    .expect("on-curve non-subgroup signature decode");
                assert_eq!(
                    sig.validate(false),
                    Err(BLST_ERROR::BLST_POINT_NOT_IN_GROUP)
                );
                assert_eq!(
                    Signature::sig_validate(sig_bytes, false).unwrap_err(),
                    BLST_ERROR::BLST_POINT_NOT_IN_GROUP
                );
                let sig_serialized = sig.serialize();
                let sig_deserialized = Signature::deserialize(&sig_serialized)
                    .expect("uncompressed on-curve signature decode");
                assert_eq!(
                    sig_deserialized.validate(false),
                    Err(BLST_ERROR::BLST_POINT_NOT_IN_GROUP)
                );
                let sig_from_bytes = Signature::from_bytes(sig_bytes)
                    .expect("from_bytes preserves the decoder contract");
                assert_eq!(
                    sig_from_bytes.validate(false),
                    Err(BLST_ERROR::BLST_POINT_NOT_IN_GROUP)
                );

                #[cfg(feature = "serde")]
                {
                    let pk_wire = rmp_serde::encode::to_vec_named(&pk)
                        .expect("serialize non-subgroup public key");
                    let pk_serde: PublicKey =
                        rmp_serde::decode::from_slice(&pk_wire)
                            .expect("serde performs on-curve decode");
                    assert_eq!(
                        pk_serde.validate(),
                        Err(BLST_ERROR::BLST_POINT_NOT_IN_GROUP)
                    );

                    let sig_wire = rmp_serde::encode::to_vec_named(&sig)
                        .expect("serialize non-subgroup signature");
                    let sig_serde: Signature =
                        rmp_serde::decode::from_slice(&sig_wire)
                            .expect("serde performs on-curve decode");
                    assert_eq!(
                        sig_serde.validate(false),
                        Err(BLST_ERROR::BLST_POINT_NOT_IN_GROUP)
                    );
                }
            }

            #[cfg(feature = "serde")]
            #[test]
            fn test_serde() {
                let seed = [0u8; 32];
                let mut rng = ChaCha20Rng::from_seed(seed);

                // generate a sk, pk, and sig, and make sure it signs
                let sk = gen_random_key(&mut rng);
                let pk = sk.sk_to_pk();
                let sig = sk.sign(b"asdf", b"qwer", b"zxcv");
                assert_eq!(
                    sig.verify(true, b"asdf", b"qwer", b"zxcv", &pk, true),
                    BLST_ERROR::BLST_SUCCESS
                );

                // roundtrip through serde
                let pk_ser =
                    rmp_serde::encode::to_vec_named(&pk).expect("ser pk");
                let sig_ser =
                    rmp_serde::encode::to_vec_named(&sig).expect("ser sig");
                let pk_des: PublicKey =
                    rmp_serde::decode::from_slice(&pk_ser).expect("des pk");
                let sig_des: Signature =
                    rmp_serde::decode::from_slice(&sig_ser).expect("des sig");

                // check that we got back the right things
                assert_eq!(pk, pk_des);
                assert_eq!(sig, sig_des);
                assert_eq!(
                    sig.verify(true, b"asdf", b"qwer", b"zxcv", &pk_des, true),
                    BLST_ERROR::BLST_SUCCESS
                );
                assert_eq!(
                    sig_des.verify(true, b"asdf", b"qwer", b"zxcv", &pk, true),
                    BLST_ERROR::BLST_SUCCESS
                );
                assert_eq!(sk.sign(b"asdf", b"qwer", b"zxcv"), sig_des);

                #[cfg(feature = "serde-secret")]
                if true {
                    let sk_ser =
                        rmp_serde::encode::to_vec_named(&sk).expect("ser sk");
                    let sk_des: SecretKey =
                        rmp_serde::decode::from_slice(&sk_ser).expect("des sk");
                    // BLS signatures are deterministic, so this establishes
                    // that sk == sk_des
                    assert_eq!(sk_des.sign(b"asdf", b"qwer", b"zxcv"), sig);
                }
            }

            #[test]
            fn test_multi_point() {
                let dst = b"BLS_SIG_BLS12381G2_XMD:SHA-256_SSWU_RO_POP_";
                let num_pks = 13;

                let seed = [0u8; 32];
                let mut rng = ChaCha20Rng::from_seed(seed);

                // Create public keys
                let sks: Vec<_> =
                    (0..num_pks).map(|_| gen_random_key(&mut rng)).collect();

                let pks =
                    sks.iter().map(|sk| sk.sk_to_pk()).collect::<Vec<_>>();
                let pks_refs: Vec<&PublicKey> =
                    pks.iter().map(|pk| pk).collect();

                // Create random message for pks to all sign
                let msg_len = (rng.next_u64() & 0x3F) + 1;
                let mut msg = vec![0u8; msg_len as usize];
                rng.fill_bytes(&mut msg);

                // Generate signature for each key pair
                let sigs = sks
                    .iter()
                    .map(|sk| sk.sign(&msg, dst, &[]))
                    .collect::<Vec<Signature>>();
                let sigs_refs: Vec<&Signature> =
                    sigs.iter().map(|s| s).collect();

                // create random values
                let mut rands: Vec<u8> = Vec::with_capacity(8 * num_pks);
                for _ in 0..num_pks {
                    let mut r = rng.next_u64();
                    while r == 0 {
                        // Reject zero as it is used for multiplication.
                        r = rng.next_u64();
                    }
                    rands.extend_from_slice(&r.to_le_bytes());
                }

                // Sanity test each current single signature
                let errs = sigs
                    .iter()
                    .zip(pks.iter())
                    .map(|(s, pk)| (s.verify(true, &msg, dst, &[], pk, true)))
                    .collect::<Vec<BLST_ERROR>>();
                assert_eq!(errs, vec![BLST_ERROR::BLST_SUCCESS; num_pks]);

                // sanity test aggregated signature
                let agg_pk = AggregatePublicKey::aggregate(&pks_refs, false)
                    .unwrap()
                    .to_public_key();
                let agg_sig = AggregateSignature::aggregate(&sigs_refs, false)
                    .unwrap()
                    .to_signature();
                let err = agg_sig.verify(true, &msg, dst, &[], &agg_pk, true);
                assert_eq!(err, BLST_ERROR::BLST_SUCCESS);

                // test multi-point aggregation using add
                let agg_pk = pks.add().to_public_key();
                let agg_sig = sigs.add().to_signature();
                let err = agg_sig.verify(true, &msg, dst, &[], &agg_pk, true);
                assert_eq!(err, BLST_ERROR::BLST_SUCCESS);

                // test multi-point aggregation using mult
                let agg_pk = pks.mult(&rands, 64).to_public_key();
                let agg_sig = sigs.mult(&rands, 64).to_signature();
                let err = agg_sig.verify(true, &msg, dst, &[], &agg_pk, true);
                assert_eq!(err, BLST_ERROR::BLST_SUCCESS);
            }
        }
    };
}

#[cfg(feature = "full-rust-bindings")]
pub mod min_pk {
    use super::*;

    sig_variant_impl!(
        "MinPk",
        blst_p1,
        blst_p1_affine,
        blst_p2,
        blst_p2_affine,
        blst_sk_to_pk2_in_g1,
        true,
        blst_hash_to_g2,
        blst_sign_pk2_in_g1,
        blst_p1_affine_is_equal,
        blst_p2_affine_is_equal,
        blst_core_verify_pk_in_g1,
        blst_p1_affine_in_g1,
        blst_p1_to_affine,
        blst_p1_from_affine,
        blst_p1_affine_serialize,
        blst_p1_affine_compress,
        blst_p1_deserialize,
        blst_p1_uncompress,
        48,
        96,
        blst_p2_affine_in_g2,
        blst_p2_to_affine,
        blst_p2_from_affine,
        blst_p2_affine_serialize,
        blst_p2_affine_compress,
        blst_p2_deserialize,
        blst_p2_uncompress,
        96,
        192,
        blst_p1_add_or_double,
        blst_p1_add_or_double_affine,
        blst_p1_cneg,
        blst_p2_add_or_double,
        blst_p2_add_or_double_affine,
        blst_p1_affine_is_inf,
        blst_p2_affine_is_inf,
        blst_p2_in_g2,
    );
}

#[cfg(feature = "full-rust-bindings")]
pub mod min_sig {
    use super::*;

    sig_variant_impl!(
        "MinSig",
        blst_p2,
        blst_p2_affine,
        blst_p1,
        blst_p1_affine,
        blst_sk_to_pk2_in_g2,
        true,
        blst_hash_to_g1,
        blst_sign_pk2_in_g2,
        blst_p2_affine_is_equal,
        blst_p1_affine_is_equal,
        blst_core_verify_pk_in_g2,
        blst_p2_affine_in_g2,
        blst_p2_to_affine,
        blst_p2_from_affine,
        blst_p2_affine_serialize,
        blst_p2_affine_compress,
        blst_p2_deserialize,
        blst_p2_uncompress,
        96,
        192,
        blst_p1_affine_in_g1,
        blst_p1_to_affine,
        blst_p1_from_affine,
        blst_p1_affine_serialize,
        blst_p1_affine_compress,
        blst_p1_deserialize,
        blst_p1_uncompress,
        48,
        96,
        blst_p2_add_or_double,
        blst_p2_add_or_double_affine,
        blst_p2_cneg,
        blst_p1_add_or_double,
        blst_p1_add_or_double_affine,
        blst_p2_affine_is_inf,
        blst_p1_affine_is_inf,
        blst_p1_in_g1,
    );
}

#[cfg(feature = "full-rust-bindings")]
pub trait MultiPoint {
    type Output;

    /// Variable-time multi-scalar multiplication for public scalars only.
    fn mult(&self, scalars: &[u8], nbits: usize) -> Self::Output;
    fn add(&self) -> Self::Output;
    fn validate(&self) -> Result<(), BLST_ERROR> {
        Err(BLST_ERROR::BLST_POINT_NOT_IN_GROUP)
    }
}

#[cfg(feature = "full-rust-bindings")]
include!("pippenger.rs");

#[cfg(all(test, feature = "full-rust-bindings"))]
mod fp12_test {
    use super::*;
    use rand_chacha::ChaCha20Rng;
    use rand_core::{RngCore, SeedableRng};

    #[test]
    fn miller_loop_n() {
        const npoints: usize = 97;
        const nbits: usize = 64;
        const nbytes: usize = (nbits + 7) / 8;

        let mut scalars = Box::new([0u8; nbytes * npoints]);
        ChaCha20Rng::from_seed([0x42; 32]).fill_bytes(scalars.as_mut());

        let mut p1s: Vec<blst_p1> = Vec::with_capacity(npoints);
        let mut p2s: Vec<blst_p2> = Vec::with_capacity(npoints);

        unsafe {
            p1s.set_len(npoints);
            p2s.set_len(npoints);

            for i in 0..npoints {
                blst_p1_mult(
                    &mut p1s[i],
                    blst_p1_generator(),
                    &scalars[i * nbytes],
                    32,
                );
                blst_p2_mult(
                    &mut p2s[i],
                    blst_p2_generator(),
                    &scalars[i * nbytes + 4],
                    32,
                );
            }
        }

        let ps = p1_affines::from(&p1s);
        let qs = p2_affines::from(&p2s);

        let mut naive = blst_fp12::default();
        for i in 0..npoints {
            naive *= blst_fp12::miller_loop(&qs[i], &ps[i]);
        }

        assert_eq!(
            naive,
            blst_fp12::miller_loop_n(qs.as_slice(), ps.as_slice())
        );
    }
}

#[cfg(all(test, feature = "full-rust-bindings"))]
mod sk_test {
    use super::*;
    use rand_chacha::ChaCha20Rng;
    use rand_core::{RngCore, SeedableRng};

    #[test]
    fn inverse() {
        let mut bytes = [0u8; 64];
        ChaCha20Rng::from_entropy().fill_bytes(bytes.as_mut());

        let mut sk = blst_scalar::default();
        let mut p1 = blst_p1::default();
        let mut p2 = blst_p2::default();

        unsafe {
            blst_scalar_from_be_bytes(&mut sk, bytes.as_ptr(), bytes.len());

            blst_p1_mult(&mut p1, blst_p1_generator(), sk.b.as_ptr(), 255);
            blst_sk_inverse(&mut sk, &sk);
            blst_p1_mult(&mut p1, &p1, sk.b.as_ptr(), 255);

            blst_p2_mult(&mut p2, blst_p2_generator(), sk.b.as_ptr(), 255);
            blst_sk_inverse(&mut sk, &sk);
            blst_p2_mult(&mut p2, &p2, sk.b.as_ptr(), 255);
        }

        assert_eq!(p1, unsafe { *blst_p1_generator() });
        assert_eq!(p2, unsafe { *blst_p2_generator() });
    }
}
