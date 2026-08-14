// Copyright 2026 Category Labs
// SPDX-License-Identifier: Apache-2.0

#![no_main]

use blst::{
    combine, fft, middle_product, min_pk, min_sig, precompute_middle_product_b,
    update, blst_fp12, blst_fr, blst_p1, blst_p1_affine, blst_p2,
    blst_p2_affine, Pairing, BLST_ERROR,
};
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let mut ikm = [0u8; 64];
    let copied = data.len().min(ikm.len());
    ikm[..copied].copy_from_slice(&data[..copied]);
    let message = &data[..data.len().min(1024)];
    let dst = b"CATEGORY-BLST-RUST-STATE-FUZZ";

    let sk_pk = min_pk::SecretKey::key_gen(&ikm, &[]).expect("bounded IKM");
    let pk_pk = sk_pk.sk_to_pk();
    let sig_pk = sk_pk.sign(message, dst, &[]);
    let aggregate_pk =
        min_pk::AggregatePublicKey::aggregate(&[&pk_pk], true)
            .expect("validated key aggregate");
    let aggregate_sig =
        min_pk::AggregateSignature::aggregate(&[&sig_pk], true)
            .expect("validated signature aggregate");
    assert_eq!(
        aggregate_sig.to_signature().verify(
            true,
            message,
            dst,
            &[],
            &aggregate_pk.to_public_key(),
            true,
        ),
        BLST_ERROR::BLST_SUCCESS
    );
    let pk_wire = pk_pk.to_bytes();
    let sig_wire = sig_pk.to_bytes();
    let _ = min_pk::AggregatePublicKey::aggregate_serialized(&[&pk_wire], true);
    let _ = min_pk::AggregateSignature::aggregate_serialized(&[&sig_wire], true);

    let sk_sig = min_sig::SecretKey::key_gen(&ikm, &[]).expect("bounded IKM");
    let pk_sig = sk_sig.sk_to_pk();
    let sig_sig = sk_sig.sign(message, dst, &[]);
    assert_eq!(
        sig_sig.verify(true, message, dst, &[], &pk_sig, true),
        BLST_ERROR::BLST_SUCCESS
    );

    // Safe state ownership and failure paths.  Affine zero is the documented
    // infinity sentinel, hence it satisfies the raw aggregate precondition.
    let mut pairing = Pairing::new(true, dst);
    pairing.raw_aggregate(&blst_p2_affine::default(), &blst_p1_affine::default());
    pairing.commit();
    let _ = pairing.miller_result();
    let _ = pairing.finalverify(None);
    pairing.init(false, message);

    let fft_size = 1usize << data.first().copied().unwrap_or(0).min(5);
    let values = vec![blst_fr::default(); fft_size];
    let mut allocating = values.clone();
    fft::fr_fft(&mut allocating, false).expect("valid FFT shape");
    fft::fr_fft(&mut allocating, true).expect("valid inverse FFT shape");
    assert_eq!(allocating, values);
    let domain_len = fft::fr_fft_domain_len(fft_size).expect("valid domain size");
    let mut domain = vec![blst_fr::default(); domain_len];
    fft::fr_fft_domain_init(&mut domain, fft_size).expect("valid domain storage");
    let mut prepared = values.clone();
    fft::fr_fft_with_domain(&mut prepared, &domain, false)
        .expect("prepared FFT");

    // Projective zero is the library's canonical infinity.  It is in each
    // subgroup and gives a small, domain-valid middle-product state case.
    let a = vec![blst_p1::default()];
    let b = vec![blst_p2::default(), blst_p2::default()];
    if let Ok(precomputed) = precompute_middle_product_b(&b) {
        let mut full = vec![blst_fp12::default(); 1];
        middle_product(&a, &precomputed, &mut full).expect("middle product");
        let mut constant = vec![blst_fp12::default(); 1];
        update(&a[0], &precomputed, &mut constant).expect("middle update");
        let mut combined = vec![blst_fp12::default(); 1];
        combine(&full, &constant, &mut combined).expect("GT combine");
    }
});
