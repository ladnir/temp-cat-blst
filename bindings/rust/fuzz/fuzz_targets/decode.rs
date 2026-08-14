// Copyright 2026 Category Labs
// SPDX-License-Identifier: Apache-2.0

#![no_main]

mod api_coverage;

use blst::{min_pk, min_sig, BLST_ERROR};
use libfuzzer_sys::fuzz_target;

fn exercise_min_pk(data: &[u8]) {
    let _ = min_pk::SecretKey::from_bytes(data);
    if let Ok(pk) = min_pk::PublicKey::from_bytes(data) {
        let _ = pk.validate();
        let encoded = pk.to_bytes();
        assert_eq!(
            min_pk::PublicKey::from_bytes(&encoded)
                .expect("canonical public key")
                .to_bytes(),
            encoded
        );
    }
    if let Ok(sig) = min_pk::Signature::from_bytes(data) {
        let _ = sig.validate(true);
        let encoded = sig.to_bytes();
        assert_eq!(
            min_pk::Signature::from_bytes(&encoded)
                .expect("canonical signature")
                .to_bytes(),
            encoded
        );
    }

    if let Ok(pk) = rmp_serde::from_slice::<min_pk::PublicKey>(data) {
        let encoded = rmp_serde::to_vec(&pk).expect("serialize decoded key");
        let _: min_pk::PublicKey =
            rmp_serde::from_slice(&encoded).expect("round-trip public key");
    }
    if let Ok(sig) = rmp_serde::from_slice::<min_pk::Signature>(data) {
        let encoded = rmp_serde::to_vec(&sig).expect("serialize decoded signature");
        let _: min_pk::Signature =
            rmp_serde::from_slice(&encoded).expect("round-trip signature");
    }
}

fn exercise_min_sig(data: &[u8]) {
    let _ = min_sig::SecretKey::from_bytes(data);
    if let Ok(pk) = min_sig::PublicKey::from_bytes(data) {
        let _ = pk.validate();
        let encoded = pk.to_bytes();
        assert_eq!(
            min_sig::PublicKey::from_bytes(&encoded)
                .expect("canonical public key")
                .to_bytes(),
            encoded
        );
    }
    if let Ok(sig) = min_sig::Signature::from_bytes(data) {
        let _ = sig.validate(true);
        let encoded = sig.to_bytes();
        assert_eq!(
            min_sig::Signature::from_bytes(&encoded)
                .expect("canonical signature")
                .to_bytes(),
            encoded
        );
    }
}

fuzz_target!(|data: &[u8]| {
    api_coverage::exercise_decode(data);
    exercise_min_pk(data);
    exercise_min_sig(data);

    let mut ikm = [0u8; 64];
    let copied = data.len().min(ikm.len());
    ikm[..copied].copy_from_slice(&data[..copied]);
    let message = &data[..data.len().min(1024)];
    let dst = b"CATEGORY-BLST-RUST-FUZZ";

    let sk_pk = min_pk::SecretKey::key_gen(&ikm, &[]).expect("bounded IKM");
    let pk_pk = sk_pk.sk_to_pk();
    let sig_pk = sk_pk.sign(message, dst, &[]);
    assert_eq!(
        sig_pk.verify(true, message, dst, &[], &pk_pk, true),
        BLST_ERROR::BLST_SUCCESS
    );

    let sk_sig = min_sig::SecretKey::key_gen(&ikm, &[]).expect("bounded IKM");
    let pk_sig = sk_sig.sk_to_pk();
    let sig_sig = sk_sig.sign(message, dst, &[]);
    assert_eq!(
        sig_sig.verify(true, message, dst, &[], &pk_sig, true),
        BLST_ERROR::BLST_SUCCESS
    );
});
