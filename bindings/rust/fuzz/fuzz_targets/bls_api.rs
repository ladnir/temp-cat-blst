// Copyright 2026 Category Labs
// SPDX-License-Identifier: Apache-2.0

#![no_main]

use blst::{min_pk, min_sig, BLST_ERROR};
use libfuzzer_sys::fuzz_target;

fn derive_ikm(data: &[u8], domain: u8) -> [u8; 32] {
    let mut result = [0u8; 32];
    for (index, byte) in result.iter_mut().enumerate() {
        let source = if data.is_empty() {
            index as u8
        } else {
            data[(index + domain as usize) % data.len()]
        };
        *byte =
            source ^ domain ^ (index as u8).wrapping_mul(29).wrapping_add(1);
    }
    result
}

fn altered(data: &[u8], selector: u8) -> Vec<u8> {
    let mut result = data.to_vec();
    if result.is_empty() {
        result.push(if selector == 0 { 1 } else { selector });
    } else {
        let index = selector as usize % result.len();
        result[index] ^= 1u8 << (selector % 8);
    }
    result
}

macro_rules! exercise_variant {
    ($variant:ident, $ikm_a:expr, $ikm_b:expr, $message_a:expr,
     $message_b:expr, $dst:expr, $aug:expr, $altered_a:expr,
     $altered_b:expr $(,)?) => {{
        let sk_a = $variant::SecretKey::key_gen($ikm_a, &[])
            .expect("fixed-width IKM A");
        let sk_b = $variant::SecretKey::key_gen($ikm_b, &[])
            .expect("fixed-width IKM B");
        let pk_a = sk_a.sk_to_pk();
        let pk_b = sk_b.sk_to_pk();
        let sig_a = sk_a.sign($message_a, $dst, $aug);
        let sig_b = sk_b.sign($message_b, $dst, $aug);

        assert_eq!(
            sig_a.verify(true, $message_a, $dst, $aug, &pk_a, true),
            BLST_ERROR::BLST_SUCCESS
        );
        assert_eq!(
            sig_b.verify(true, $message_b, $dst, $aug, &pk_b, true),
            BLST_ERROR::BLST_SUCCESS
        );
        assert_ne!(
            sig_b.verify(true, $altered_b, $dst, $aug, &pk_b, true),
            BLST_ERROR::BLST_SUCCESS
        );

        // Aggregate verification has no augmentation argument. Keep these
        // signatures separate from the augmented direct-verification case.
        let aggregate_sig_a = sk_a.sign($message_a, $dst, &[]);
        let aggregate_sig_b = sk_b.sign($message_b, $dst, &[]);
        let aggregate_signature = $variant::AggregateSignature::aggregate(
            &[&aggregate_sig_a, &aggregate_sig_b],
            true,
        )
        .expect("two validated signatures");
        let reverse_signature = $variant::AggregateSignature::aggregate(
            &[&aggregate_sig_b, &aggregate_sig_a],
            true,
        )
        .expect("reversed validated signatures");
        assert_eq!(
            aggregate_signature.to_signature(),
            reverse_signature.to_signature()
        );
        assert_eq!(
            aggregate_signature.to_signature().aggregate_verify(
                true,
                &[$message_a, $message_b],
                $dst,
                &[&pk_a, &pk_b],
                true,
            ),
            BLST_ERROR::BLST_SUCCESS
        );
        assert_ne!(
            aggregate_signature.to_signature().aggregate_verify(
                true,
                &[$message_a, $altered_b],
                $dst,
                &[&pk_a, &pk_b],
                true,
            ),
            BLST_ERROR::BLST_SUCCESS
        );

        let mut incremental_signature =
            $variant::AggregateSignature::from_signature(&aggregate_sig_a);
        incremental_signature
            .add_signature(&aggregate_sig_b, true)
            .expect("validated incremental signature");
        assert_eq!(
            incremental_signature.to_signature(),
            aggregate_signature.to_signature()
        );
        let sig_a_wire = aggregate_sig_a.to_bytes();
        let sig_b_wire = aggregate_sig_b.to_bytes();
        let serialized_signature =
            $variant::AggregateSignature::aggregate_serialized(
                &[sig_a_wire.as_slice(), sig_b_wire.as_slice()],
                true,
            )
            .expect("canonical serialized signatures");
        assert_eq!(
            serialized_signature.to_signature(),
            aggregate_signature.to_signature()
        );

        let aggregate_public_key =
            $variant::AggregatePublicKey::aggregate(&[&pk_a, &pk_b], true)
                .expect("two validated public keys");
        let reverse_public_key =
            $variant::AggregatePublicKey::aggregate(&[&pk_b, &pk_a], true)
                .expect("reversed validated public keys");
        assert_eq!(
            aggregate_public_key.to_public_key(),
            reverse_public_key.to_public_key()
        );
        let mut incremental_public_key =
            $variant::AggregatePublicKey::from_public_key(&pk_a);
        incremental_public_key
            .add_public_key(&pk_b, true)
            .expect("validated incremental public key");
        assert_eq!(
            incremental_public_key.to_public_key(),
            aggregate_public_key.to_public_key()
        );
        let pk_a_wire = pk_a.to_bytes();
        let pk_b_wire = pk_b.to_bytes();
        let serialized_public_key =
            $variant::AggregatePublicKey::aggregate_serialized(
                &[pk_a_wire.as_slice(), pk_b_wire.as_slice()],
                true,
            )
            .expect("canonical serialized public keys");
        assert_eq!(
            serialized_public_key.to_public_key(),
            aggregate_public_key.to_public_key()
        );

        let same_a = sk_a.sign($message_a, $dst, &[]);
        let same_b = sk_b.sign($message_a, $dst, &[]);
        let same_aggregate =
            $variant::AggregateSignature::aggregate(&[&same_a, &same_b], true)
                .expect("same-message signatures");
        assert_eq!(
            same_aggregate.to_signature().fast_aggregate_verify(
                true,
                $message_a,
                $dst,
                &[&pk_a, &pk_b],
            ),
            BLST_ERROR::BLST_SUCCESS
        );
        assert_eq!(
            same_aggregate
                .to_signature()
                .fast_aggregate_verify_pre_aggregated(
                    true,
                    $message_a,
                    $dst,
                    &aggregate_public_key.to_public_key(),
                ),
            BLST_ERROR::BLST_SUCCESS
        );
        assert_ne!(
            same_aggregate.to_signature().fast_aggregate_verify(
                true,
                $altered_a,
                $dst,
                &[&pk_a, &pk_b],
            ),
            BLST_ERROR::BLST_SUCCESS
        );

        let weighted_public_key =
            $variant::AggregatePublicKey::aggregate_with_randomness(
                &[pk_a, pk_b],
                &[1, 2],
                8,
                true,
            )
            .expect("weighted public keys");
        let weighted_signature =
            $variant::AggregateSignature::aggregate_with_randomness(
                &[aggregate_sig_a, aggregate_sig_b],
                &[1, 2],
                8,
                true,
            )
            .expect("weighted signatures");
        weighted_public_key
            .to_public_key()
            .validate()
            .expect("weighted public key validation");
        weighted_signature
            .validate()
            .expect("weighted signature validation");
    }};
}

fuzz_target!(|data: &[u8]| {
    let data = &data[..data.len().min(512)];
    let selector = data.first().copied().unwrap_or(0);
    let body = data.get(1..).unwrap_or_default();
    let split = body.len() / 2;
    let message_a = &body[..split];
    let message_b = &body[split..];
    let dst = &data[..data.len().min(64)];
    let aug = &data[data.len().saturating_sub(64)..];
    let altered_message_a = altered(message_a, selector ^ 0x5a);
    let altered_message_b = altered(message_b, selector);
    let ikm_a = derive_ikm(data, 0x31);
    let ikm_b = derive_ikm(data, 0xa7);

    exercise_variant!(
        min_pk,
        &ikm_a,
        &ikm_b,
        message_a,
        message_b,
        dst,
        aug,
        &altered_message_a,
        &altered_message_b,
    );
    exercise_variant!(
        min_sig,
        &ikm_a,
        &ikm_b,
        message_a,
        message_b,
        dst,
        aug,
        &altered_message_a,
        &altered_message_b,
    );
});
