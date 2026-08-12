// Copyright Supranational LLC and Category Labs
// Licensed under the Apache License, Version 2.0, see LICENSE for details.
// SPDX-License-Identifier: Apache-2.0

#![cfg(feature = "full-rust-bindings")]

use blst::{min_pk, BLST_ERROR};

#[test]
fn exact_g1_and_g2_compressed_infinity_encodings_are_accepted() {
    let mut g1 = [0u8; 48];
    g1[0] = 0xc0;
    let public_key =
        min_pk::PublicKey::uncompress(&g1).expect("valid G1 infinity encoding");
    assert_eq!(public_key.validate(), Err(BLST_ERROR::BLST_PK_IS_INFINITY));

    let mut g2 = [0u8; 96];
    g2[0] = 0xc0;
    let signature =
        min_pk::Signature::uncompress(&g2).expect("valid G2 infinity encoding");
    assert_eq!(
        signature.validate(true),
        Err(BLST_ERROR::BLST_PK_IS_INFINITY)
    );
}
