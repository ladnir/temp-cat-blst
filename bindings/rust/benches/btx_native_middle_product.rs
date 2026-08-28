// Copyright (C) 2025 Category Labs, Inc.

use std::time::Duration;

use blst::{
    blst_avx_gt, blst_avx_p1, blst_avx_p1_generator, blst_avx_p1_mult,
    blst_avx_p2, blst_avx_p2_generator, blst_avx_p2_mult,
    middle_product_avx_native_with_workspace,
    NativeAvxMiddleProductPrecomputedB, NativeAvxMiddleProductWorkspace,
};
use criterion::{
    black_box, criterion_group, criterion_main, BenchmarkId, Criterion,
};

fn scalar_bytes(value: u64) -> [u8; 32] {
    let mut scalar = [0u8; 32];
    scalar[..8].copy_from_slice(&value.to_le_bytes());
    scalar
}

fn g1_multiple(value: u64) -> blst_avx_p1 {
    let generator = unsafe { blst_avx_p1_generator() };
    assert!(!generator.is_null());
    let mut output = blst_avx_p1::default();
    assert!(unsafe {
        blst_avx_p1_mult(
            &mut output,
            generator,
            scalar_bytes(value).as_ptr(),
            64,
        )
    });
    output
}

fn g2_multiple(value: u64) -> blst_avx_p2 {
    let generator = unsafe { blst_avx_p2_generator() };
    assert!(!generator.is_null());
    let mut output = blst_avx_p2::default();
    assert!(unsafe {
        blst_avx_p2_mult(
            &mut output,
            generator,
            scalar_bytes(value).as_ptr(),
            64,
        )
    });
    output
}

fn bench_native_middle_product(c: &mut Criterion) {
    let mut group = c.benchmark_group("btx_native_middle_product");
    group.sample_size(10);
    group.warm_up_time(Duration::from_millis(500));
    group.measurement_time(Duration::from_secs(2));

    for count in [32usize, 64, 128, 256] {
        let input = (0..count)
            .map(|index| g1_multiple(index as u64 + 1))
            .collect::<Vec<_>>();
        let precompute_input = (0..2 * count)
            .map(|index| g2_multiple(2 * index as u64 + 3))
            .collect::<Vec<_>>();
        let precomputed =
            NativeAvxMiddleProductPrecomputedB::new(&precompute_input)
                .expect("native middle-product precomputation");
        let mut workspace = NativeAvxMiddleProductWorkspace::new(&precomputed)
            .expect("native middle-product workspace");
        let mut output = vec![blst_avx_gt::default(); count];

        group.bench_function(BenchmarkId::from_parameter(count), |b| {
            b.iter(|| {
                middle_product_avx_native_with_workspace(
                    black_box(&input),
                    black_box(&precomputed),
                    black_box(&mut workspace),
                    black_box(&mut output),
                )
                .expect("native middle product");
                black_box(&output);
            });
        });
    }

    group.finish();
}

criterion_group!(benches, bench_native_middle_product);
criterion_main!(benches);
