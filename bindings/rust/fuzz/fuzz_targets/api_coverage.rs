// Copyright 2026 Category Labs
// SPDX-License-Identifier: Apache-2.0

//! Executable coverage owner for safe Rust APIs.
//!
//! Inputs are either arbitrary byte slices passed to decoding operations or
//! values constructed through safe APIs. The module contains no raw FFI and no
//! unsafe code.

use blst::{
    avx_p1s_in_g1, avx_p2s_in_g2, combine, fft, gt_msm, middle_product,
    min_pk, min_sig, precompute_middle_product_b, update, blst_fp12, blst_fr,
    blst_p1, blst_p1_affine, blst_p2, blst_p2_affine, blst_scalar, p1_affines,
    p2_affines, Pairing, BLST_ERROR,
};

macro_rules! exercise_variant {
    ($variant:ident, $ikm:expr, $message:expr, $dst:expr, $child:expr) => {{
        let sk = $variant::SecretKey::key_gen($ikm, &[])
            .expect("bounded IKM must produce a key");
        let _ = $variant::SecretKey::key_gen_v3($ikm, &[]);
        let _ = $variant::SecretKey::key_gen_v4_5($ikm, &[], &[])
            .expect("bounded v4.5 IKM");
        let _ = $variant::SecretKey::key_gen_v5($ikm, &[], &[])
            .expect("bounded v5 IKM");
        let master = $variant::SecretKey::derive_master_eip2333($ikm)
            .expect("bounded EIP-2333 IKM");
        let _child = master.derive_child_eip2333($child);
        let secret_wire = sk.serialize();
        let _ = $variant::SecretKey::deserialize(&secret_wire)
            .expect("canonical secret key");
        assert_eq!(sk.to_bytes(), secret_wire);
        let _ = $variant::SecretKey::from_bytes(&secret_wire)
            .expect("canonical secret key alias");

        let pk = sk.sk_to_pk();
        pk.validate().expect("derived public key");
        let pk_compressed = pk.compress();
        let pk_serialized = pk.serialize();
        assert_eq!(pk.to_bytes(), pk_compressed);
        let _ = $variant::PublicKey::uncompress(&pk_compressed)
            .expect("canonical compressed public key");
        let _ = $variant::PublicKey::deserialize(&pk_serialized)
            .expect("canonical serialized public key");
        let _ = $variant::PublicKey::from_bytes(&pk_compressed)
            .expect("canonical public key alias");
        let _ = $variant::PublicKey::key_validate(&pk_compressed)
            .expect("validated public key bytes");

        let mut aggregate_pk =
            $variant::AggregatePublicKey::from_public_key(&pk);
        assert_eq!(aggregate_pk.to_public_key(), pk);
        let aggregate_pk_one =
            $variant::AggregatePublicKey::aggregate(&[&pk], true)
                .expect("one validated public key");
        let aggregate_pk_random =
            $variant::AggregatePublicKey::aggregate_with_randomness(
                &[pk], &[1], 8, true,
            )
            .expect("one weighted public key");
        let _ = aggregate_pk_random.to_public_key();
        let pk_wire: &[u8] = &pk_compressed;
        let _ = $variant::AggregatePublicKey::aggregate_serialized(
            &[pk_wire], true,
        )
        .expect("one serialized public key");
        aggregate_pk.add_aggregate(&aggregate_pk_one);
        aggregate_pk.sub_aggregate(&aggregate_pk_one);
        aggregate_pk
            .add_public_key(&pk, true)
            .expect("validated public-key addition");
        let _ = $variant::PublicKey::from_aggregate(&aggregate_pk_one);

        let sig = sk.sign($message, $dst, &[]);
        sig.validate(true).expect("derived non-infinity signature");
        let sig_compressed = sig.compress();
        let sig_serialized = sig.serialize();
        assert_eq!(sig.to_bytes(), sig_compressed);
        let _ = $variant::Signature::uncompress(&sig_compressed)
            .expect("canonical compressed signature");
        let _ = $variant::Signature::deserialize(&sig_serialized)
            .expect("canonical serialized signature");
        let _ = $variant::Signature::from_bytes(&sig_compressed)
            .expect("canonical signature alias");
        let _ = $variant::Signature::sig_validate(&sig_compressed, true)
            .expect("validated signature bytes");
        assert!(sig.subgroup_check());
        assert_eq!(
            sig.verify(true, $message, $dst, &[], &pk, true),
            BLST_ERROR::BLST_SUCCESS
        );
        assert_eq!(
            sig.aggregate_verify(
                true,
                &[$message],
                $dst,
                &[&pk],
                true,
            ),
            BLST_ERROR::BLST_SUCCESS
        );
        assert_eq!(
            sig.fast_aggregate_verify(true, $message, $dst, &[&pk]),
            BLST_ERROR::BLST_SUCCESS
        );
        assert_eq!(
            sig.fast_aggregate_verify_pre_aggregated(
                true,
                $message,
                $dst,
                &aggregate_pk_one.to_public_key(),
            ),
            BLST_ERROR::BLST_SUCCESS
        );
        let rand_ref: &blst_scalar = (&sk).into();
        let rands = [rand_ref.clone()];
        assert_eq!(
            $variant::Signature::verify_multiple_aggregate_signatures(
                &[$message],
                $dst,
                &[&pk],
                true,
                &[&sig],
                true,
                &rands,
                8,
            ),
            BLST_ERROR::BLST_SUCCESS
        );

        let mut aggregate_sig =
            $variant::AggregateSignature::from_signature(&sig);
        aggregate_sig.validate().expect("derived aggregate signature");
        assert_eq!(aggregate_sig.to_signature(), sig);
        let aggregate_sig_one =
            $variant::AggregateSignature::aggregate(&[&sig], true)
                .expect("one validated signature");
        let aggregate_sig_random =
            $variant::AggregateSignature::aggregate_with_randomness(
                &[sig], &[1], 8, true,
            )
            .expect("one weighted signature");
        let _ = aggregate_sig_random.to_signature();
        let sig_wire: &[u8] = &sig_compressed;
        let _ = $variant::AggregateSignature::aggregate_serialized(
            &[sig_wire], true,
        )
        .expect("one serialized signature");
        aggregate_sig.add_aggregate(&aggregate_sig_one);
        aggregate_sig
            .add_signature(&sig, true)
            .expect("validated signature addition");
        assert!(aggregate_sig.subgroup_check());
        let _ = $variant::Signature::from_aggregate(&aggregate_sig_one);

        (sk, pk, sig, aggregate_pk_one, aggregate_sig_one)
    }};
}

#[allow(dead_code)]
pub fn exercise_decode(data: &[u8]) {
    let _ = min_pk::SecretKey::deserialize(data);
    let _ = min_pk::SecretKey::from_bytes(data);
    let _ = min_pk::PublicKey::uncompress(data);
    let _ = min_pk::PublicKey::deserialize(data);
    let _ = min_pk::PublicKey::from_bytes(data);
    let _ = min_pk::PublicKey::key_validate(data);
    let _ = min_pk::Signature::uncompress(data);
    let _ = min_pk::Signature::deserialize(data);
    let _ = min_pk::Signature::from_bytes(data);
    let _ = min_pk::Signature::sig_validate(data, true);

    let _ = min_sig::SecretKey::deserialize(data);
    let _ = min_sig::SecretKey::from_bytes(data);
    let _ = min_sig::PublicKey::uncompress(data);
    let _ = min_sig::PublicKey::deserialize(data);
    let _ = min_sig::PublicKey::from_bytes(data);
    let _ = min_sig::PublicKey::key_validate(data);
    let _ = min_sig::Signature::uncompress(data);
    let _ = min_sig::Signature::deserialize(data);
    let _ = min_sig::Signature::from_bytes(data);
    let _ = min_sig::Signature::sig_validate(data, true);

    let _ = rmp_serde::from_slice::<min_pk::SecretKey>(data);
    let _ = rmp_serde::from_slice::<min_pk::PublicKey>(data);
    let _ = rmp_serde::from_slice::<min_pk::Signature>(data);
    let _ = rmp_serde::from_slice::<min_sig::SecretKey>(data);
    let _ = rmp_serde::from_slice::<min_sig::PublicKey>(data);
    let _ = rmp_serde::from_slice::<min_sig::Signature>(data);
}

#[allow(dead_code)]
pub fn exercise_state(data: &[u8]) {
    let mut ikm = [0u8; 64];
    let copied = data.len().min(ikm.len());
    ikm[..copied].copy_from_slice(&data[..copied]);
    let message = &data[..data.len().min(1024)];
    let dst = b"CATEGORY-BLST-RUST-API-COVERAGE";
    let child = u32::from_le_bytes([
        data.first().copied().unwrap_or(0),
        data.get(1).copied().unwrap_or(0),
        data.get(2).copied().unwrap_or(0),
        data.get(3).copied().unwrap_or(0),
    ]);

    let (_sk_pk, pk, sig, _, _) =
        exercise_variant!(min_pk, &ikm, message, dst, child);
    let _ = exercise_variant!(min_sig, &ikm, message, dst, child);

    let scalar = [data.first().copied().unwrap_or(1) | 1];
    let g1 = blst_p1::default();
    let _ = g1.mult(&scalar, 8);
    let _ = g1.mult_vartime(&scalar, 8);
    let _ = g1.unchecked_mult(&scalar, 8);
    let g2 = blst_p2::default();
    let _ = g2.mult(&scalar, 8);
    let _ = g2.mult_vartime(&scalar, 8);
    let _ = g2.unchecked_mult(&scalar, 8);
    let _ = blst_scalar::hash_to(message, dst);

    let q = blst_p2_affine::default();
    let p = blst_p1_affine::default();
    let miller = blst_fp12::miller_loop(&q, &p);
    let miller_n = blst_fp12::miller_loop_n(&[q], &[p]);
    assert!(blst_fp12::finalverify(&miller, &miller_n));
    let gt = miller.final_exp();
    assert!(gt.in_group());
    let _ = gt.gt_mult(&scalar, 8);
    let _ = gt.gt_mult_vartime(&scalar, 8);
    let _ = gt.to_bendian();

    let mut values = [blst_fr::default()];
    fft::fr_fft(&mut values, false).expect("size-one FFT");
    let domain_len = fft::fr_fft_domain_len(1).expect("size-one domain");
    let mut domain = vec![blst_fr::default(); domain_len];
    fft::fr_fft_domain_init(&mut domain, 1).expect("size-one domain init");
    fft::fr_fft_with_domain(&mut values, &domain, true)
        .expect("size-one prepared FFT");
    let mut coefficients = [blst_fr::default()];
    fft::fr_interpolate(&values, &values, &mut coefficients)
        .expect("one-point interpolation");
    let scratch_len =
        fft::fr_interpolate_fft_domain_suffix_scratch_len(1, 1)
            .expect("size-one interpolation scratch");
    let mut scratch = vec![blst_fr::default(); scratch_len];
    fft::fr_interpolate_fft_domain_suffix_with_scratch(
        &[0],
        &values,
        0,
        1,
        &mut coefficients,
        &mut scratch,
    )
    .expect("size-one suffix interpolation");

    let b = [blst_p2::default(), blst_p2::default()];
    let precomputed = precompute_middle_product_b(&b)
        .expect("valid middle-product B shape");
    let also_precomputed =
        blst::MiddleProductPrecomputedB::new(&b).expect("valid constructor");
    assert_eq!(precomputed.len(), also_precomputed.len());
    assert!(!precomputed.is_empty());
    let a = [blst_p1::default()];
    let mut full = [blst_fp12::default()];
    middle_product(&a, &precomputed, &mut full).expect("middle product");
    let mut constant = [blst_fp12::default()];
    update(&a[0], &precomputed, &mut constant).expect("middle update");
    let mut combined = [blst_fp12::default()];
    combine(&full, &constant, &mut combined).expect("middle combine");

    let raw_pk: blst_p1_affine = pk.into();
    let raw_sig: blst_p2_affine = sig.into();
    let projective_pk: blst_p1 =
        min_pk::AggregatePublicKey::from_public_key(&pk).into();
    let projective_sig: blst_p2 =
        min_pk::AggregateSignature::from_signature(&sig).into();
    let p1_batch = p1_affines::from(&[projective_pk]);
    assert_eq!(p1_batch.as_slice().len(), 1);
    let _ = p1_batch.mult(&scalar, 8);
    let _ = p1_batch.add();
    let p2_batch = p2_affines::from(&[projective_sig]);
    assert_eq!(p2_batch.as_slice().len(), 1);
    let _ = p2_batch.mult(&scalar, 8);
    let _ = p2_batch.add();

    let mut pairing = Pairing::new(true, dst);
    assert_eq!(
        pairing.aggregate(&raw_pk, true, Some(&raw_sig), true, message, &[]),
        BLST_ERROR::BLST_SUCCESS
    );
    assert_eq!(
        pairing.mul_n_aggregate(
            &raw_pk,
            true,
            Some(&raw_sig),
            true,
            &scalar,
            8,
            message,
            &[],
        ),
        BLST_ERROR::BLST_SUCCESS
    );
    pairing.commit();
    let _ = pairing.miller_result();
    let _ = pairing.finalverify(None);
    let mut external = blst_fp12::default();
    Pairing::aggregated(&mut external, &raw_sig);
    let mut other = Pairing::new(true, dst);
    other.raw_aggregate(&raw_sig, &raw_pk);
    other.commit();
    assert_eq!(pairing.merge(&other), BLST_ERROR::BLST_SUCCESS);
    pairing.init(false, message);

    assert!(blst::uniq(&[message]));
    let seed = [0x42; 32];
    let _ = avx_p1s_in_g1(&[], seed);
    let _ = avx_p2s_in_g2(&[], seed, &mut []);
    let _ = gt_msm(&[], &[]);
}
