use blst::{
    blst_fp12, blst_p1, blst_p1_affine, blst_p1_generator, blst_p1_in_g1, blst_p1_mult,
    blst_p1_to_affine, blst_p2, blst_p2_affine, blst_p2_generator, blst_p2_in_g2, blst_p2_mult,
    blst_p2_to_affine, combine, middle_product, precompute_middle_product_b, update,
};

fn scalar_bytes(value: u64) -> [u8; 32] {
    let mut out = [0u8; 32];
    out[..8].copy_from_slice(&value.to_le_bytes());
    out
}

fn raw_pairing_smoke() {
    unsafe {
        let p_scalar = scalar_bytes(3);
        let q_scalar = scalar_bytes(11);

        let mut p = blst_p1::default();
        let mut q = blst_p2::default();
        blst_p1_mult(&mut p, blst_p1_generator(), p_scalar.as_ptr(), 256);
        blst_p2_mult(&mut q, blst_p2_generator(), q_scalar.as_ptr(), 256);

        assert!(blst_p1_in_g1(&p));
        assert!(blst_p2_in_g2(&q));

        let mut p_affine = blst_p1_affine::default();
        let mut q_affine = blst_p2_affine::default();
        blst_p1_to_affine(&mut p_affine, &p);
        blst_p2_to_affine(&mut q_affine, &q);

        let gt = blst_fp12::miller_loop(&q_affine, &p_affine).final_exp();
        assert!(gt.in_group());
        assert!(blst_fp12::finalverify(&gt, &gt));
    }
}

fn middle_product_smoke() {
    unsafe {
        let mut a_full = [blst_p1::default()];
        let a_without_constant = [blst_p1::default()];
        let mut b = [blst_p2::default(), blst_p2::default()];

        let a_scalar = scalar_bytes(5);
        let b0_scalar = scalar_bytes(7);
        let b1_scalar = scalar_bytes(13);

        blst_p1_mult(&mut a_full[0], blst_p1_generator(), a_scalar.as_ptr(), 256);
        blst_p2_mult(&mut b[0], blst_p2_generator(), b0_scalar.as_ptr(), 256);
        blst_p2_mult(&mut b[1], blst_p2_generator(), b1_scalar.as_ptr(), 256);

        let precomputed = precompute_middle_product_b(&b).expect("precompute(b) failed");
        let mut expected = [blst_fp12::default()];
        let mut mid = [blst_fp12::default()];
        let mut cnst = [blst_fp12::default()];
        let mut combined = [blst_fp12::default()];

        middle_product(&a_full, &precomputed, &mut expected).expect("middle_product failed");
        middle_product(&a_without_constant, &precomputed, &mut mid)
            .expect("middle_product without constant failed");
        update(&a_full[0], &precomputed, &mut cnst).expect("update failed");
        combine(&mid, &cnst, &mut combined).expect("combine failed");

        assert!(expected[0].in_group());
        assert_eq!(expected[0], combined[0]);
    }
}

#[cfg(feature = "full-rust-bindings")]
fn full_bindings_smoke() {
    use blst::min_pk::SecretKey;

    let ikm = [7u8; 32];
    let dst = b"BLS_SIG_BLS12381G2_XMD:SHA-256_SSWU_RO_NUL_";
    let msg = b"blst smoke example";

    let sk = SecretKey::key_gen(&ikm, &[]).expect("key_gen failed");
    let pk = sk.sk_to_pk();
    let sig = sk.sign(msg, dst, &[]);
    let err = sig.verify(true, msg, dst, &[], &pk, true);
    assert_eq!(err, blst::BLST_ERROR::BLST_SUCCESS);
}

#[cfg(not(feature = "full-rust-bindings"))]
fn full_bindings_smoke() {
    eprintln!("full-rust-bindings not enabled; skipping sign/verify smoke");
}

fn main() {
    raw_pairing_smoke();
    middle_product_smoke();
    full_bindings_smoke();
    println!("blst smoke example passed");
}
