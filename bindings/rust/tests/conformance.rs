use blst::*;

const G1_COMPRESSED: &str = concat!(
    "97f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac5",
    "86c55e83ff97a1aeffb3af00adb22c6bb"
);
const G2_COMPRESSED: &str = concat!(
    "93e02b6052719f607dacd3a088274f65596bd0d09920b61ab5da61bbdc7f504",
    "9334cf11213945d57e5ac7d055d042b7e024aa2b2f08f0a91260805272dc510",
    "51c6e47ad4fa403b02b4510b647ae3d1770bac0326a805bbefd48056c8c121b",
    "db8"
);

// draft-irtf-cfrg-pairing-friendly-curves-13, Appendix B.1. The draft's
// e_0..e_11 tower coefficients are reordered into BLST Fp12 serialization.
const GT_OUTPUT_CHUNKS: [&str; 12] = [
    "11619b45f61edfe3b47a15fac19442526ff489dcda25e59121d9931438907dfd448299a87dde3a649bdba96e84d54558",
    "153ce14a76a53e205ba8f275ef1137c56a566f638b52d34ba3bf3bf22f277d70f76316218c0dfd583a394b8448d2be7f",
    "01ecfcf31c86257ab00b4709c33f1c9c4e007659dd5ffc4a735192167ce197058cfb4c94225e7f1b6c26ad9ba68f63bc",
    "08890726743a1f94a8193a166800b7787744a8ad8e2f9365db76863e894b7a11d83f90d873567e9d645ccf725b32d26f",
    "095668fb4a02fe930ed44767834c915b283b1c6ca98c047bd4c272e9ac3f3ba6ff0b05a93e59c71fba77bce995f04692",
    "16deedaa683124fe7260085184d88f7d036b86f53bb5b7f1fc5e248814782065413e7d958d17960109ea006b2afdeb5f",
    "0e61c752414ca5dfd258e9606bac08daec29b3e2c57062669556954fb227d3f1260eedf25446a086b0844bcd43646c10",
    "0fe63f185f56dd29150fc498bbeea78969e7e783043620db33f75a05a0a2ce5c442beaff9da195ff15164c00ab66bdde",
    "09c92cf02f3cd3d2f9d34bc44eee0dd50314ed44ca5d30ce6a9ec0539be7a86b121edc61839ccc908c4bdde256cd6048",
    "111061f398efc2a97ff825b04d21089e24fd8b93a47e41e60eae7e9b2a38d54fa4dedced0811c34ce528781ab9e929c7",
    "10900338a92ed0b47af211636f7cfdec717b7ee43900eee9b5fc24f0000c5874d4801372db478987691c566a8c474978",
    "1454814f3085f0e6602247671bc408bbce2007201536818c901dbd4d2095dd86c1ec8b888e59611f60a301af7776be3d",
];

// Independent cross-implementation vectors from @noble/curves 2.2.0
// (npm gitHead 043905dfc90e75fe77bbbd60c791baa1e5177801).
#[cfg(feature = "full-rust-bindings")]
const REFERENCE_SECRET_KEY: &str = "411602cb19a6ccc34494d79d98ef1e7ed5af25f7";
#[cfg(feature = "full-rust-bindings")]
const REFERENCE_G1_SIGNATURE: &str = concat!(
    "af6305904771c5883edf9f1142b0df66b74f2fcd5b512e0ef1adb197424f2a1",
    "52e69bd0e3fa13ab4ae1457431717dd00"
);
#[cfg(feature = "full-rust-bindings")]
const REFERENCE_G2_SIGNATURE: &str = concat!(
    "b38c2b4675bbc8118304a58717acf03fcf4f4416be9d3c57756009427d0566ed",
    "1bdd7c90476a0fe8d1ffc3e37594703a121cb7c540f2d02721fd7684536063dfa",
    "d528fc113ce8900ec316eba0f22713afabddacd931de3fa13abeb3cb53b51dc"
);

fn hex_nibble(value: u8) -> u8 {
    match value {
        b'0'..=b'9' => value - b'0',
        b'a'..=b'f' => value - b'a' + 10,
        _ => panic!("invalid lowercase hex digit"),
    }
}

fn decode_hex<const N: usize>(encoded: &str) -> [u8; N] {
    assert_eq!(encoded.len(), N * 2);
    let mut out = [0u8; N];
    for (i, byte) in out.iter_mut().enumerate() {
        *byte = (hex_nibble(encoded.as_bytes()[2 * i]) << 4)
            | hex_nibble(encoded.as_bytes()[2 * i + 1]);
    }
    out
}

fn draft_pairing_vector() -> blst_fp12 {
    let mut raw = blst_fp12 {
        fp6: [blst_fp6::default(); 2],
    };
    let mut chunk = 0;
    for fp6_component in 0..3 {
        for fp12_component in 0..2 {
            for fp2_component in 0..2 {
                let coefficient = decode_hex::<48>(GT_OUTPUT_CHUNKS[chunk]);
                unsafe {
                    blst_fp_from_bendian(
                        &mut raw.fp6[fp12_component].fp2[fp6_component].fp
                            [fp2_component],
                        coefficient.as_ptr(),
                    );
                }
                chunk += 1;
            }
        }
    }
    raw
}

fn expected_production_pairing() -> blst_fp12 {
    let raw = draft_pairing_vector();
    // Appendix B.1 defines e_lib(P,Q) = e_pseudocode(P,Q)^3 for BLS12-381.
    raw * raw * raw
}

#[test]
fn pairing_friendly_curves_draft13_scalar() {
    unsafe {
        let mut g1 = [0u8; 48];
        let mut g2 = [0u8; 96];
        blst_p1_compress(g1.as_mut_ptr(), blst_p1_generator());
        blst_p2_compress(g2.as_mut_ptr(), blst_p2_generator());
        assert_eq!(g1, decode_hex::<48>(G1_COMPRESSED));
        assert_eq!(g2, decode_hex::<96>(G2_COMPRESSED));

        let pairing = blst_fp12::miller_loop(
            &*blst_p2_affine_generator(),
            &*blst_p1_affine_generator(),
        )
        .final_exp();
        assert_eq!(pairing, expected_production_pairing());

        let raw_bytes = draft_pairing_vector().to_bendian();
        for (actual, expected) in
            raw_bytes.chunks_exact(48).zip(GT_OUTPUT_CHUNKS)
        {
            assert_eq!(actual, decode_hex::<48>(expected));
        }
    }
}

#[cfg(feature = "avx")]
#[test]
fn pairing_friendly_curves_draft13_avx() {
    unsafe {
        assert!(blst_avx_enabled(), "AVX feature did not enable its backend");
        let g1_generator = blst_avx_p1_generator();
        let g2_generator = blst_avx_p2_generator();
        assert!(!g1_generator.is_null());
        assert!(!g2_generator.is_null());

        let mut g1 = [0u8; 48];
        let mut g2 = [0u8; 96];
        assert!(blst_avx_p1_compress(g1.as_mut_ptr(), g1_generator));
        assert!(blst_avx_p2_compress(g2.as_mut_ptr(), g2_generator));
        assert_eq!(g1, decode_hex::<48>(G1_COMPRESSED));
        assert_eq!(g2, decode_hex::<96>(G2_COMPRESSED));

        let mut avx_pairing = blst_avx_fp12::default();
        let mut scalar_pairing = blst_fp12::default();
        assert!(blst_avx_pairing(
            &mut avx_pairing,
            g1_generator,
            g2_generator
        ));
        assert!(blst_fp12_from_avx(&mut scalar_pairing, &avx_pairing));
        assert_eq!(scalar_pairing, expected_production_pairing());
    }
}

#[cfg(feature = "full-rust-bindings")]
#[test]
fn bls_draft07_basic_signature_vectors() {
    const MESSAGE: &[u8] = b"sample";
    const G1_DST: &[u8] = b"BLS_SIG_BLS12381G1_XMD:SHA-256_SSWU_RO_NUL_";
    const G2_DST: &[u8] = b"BLS_SIG_BLS12381G2_XMD:SHA-256_SSWU_RO_NUL_";

    let mut secret_key = [0u8; 32];
    secret_key[12..].copy_from_slice(&decode_hex::<20>(REFERENCE_SECRET_KEY));

    let min_sig_sk = min_sig::SecretKey::from_bytes(&secret_key).unwrap();
    let min_sig_pk = min_sig_sk.sk_to_pk();
    let g1_signature = min_sig_sk.sign(MESSAGE, G1_DST, &[]);
    assert_eq!(
        g1_signature.compress(),
        decode_hex::<48>(REFERENCE_G1_SIGNATURE)
    );
    assert_eq!(
        g1_signature.verify(true, MESSAGE, G1_DST, &[], &min_sig_pk, true),
        BLST_ERROR::BLST_SUCCESS
    );

    let min_pk_sk = min_pk::SecretKey::from_bytes(&secret_key).unwrap();
    let min_pk_pk = min_pk_sk.sk_to_pk();
    let g2_signature = min_pk_sk.sign(MESSAGE, G2_DST, &[]);
    assert_eq!(
        g2_signature.compress(),
        decode_hex::<96>(REFERENCE_G2_SIGNATURE)
    );
    assert_eq!(
        g2_signature.verify(true, MESSAGE, G2_DST, &[], &min_pk_pk, true),
        BLST_ERROR::BLST_SUCCESS
    );
}
