// Copyright (C) 2026 Category Labs, Inc.

use std::time::Duration;

use blst::avx::{
    G1Projective, G2Projective, MillerLoopResult, PairingWorkspace,
    PreparedG2Batch,
};
use criterion::{
    black_box, criterion_group, criterion_main, BenchmarkId, Criterion,
};
use rayon::{prelude::*, ThreadPool};

fn evaluate_prepared_ranges(
    batch: &PreparedG2Batch,
    g1: &[G1Projective],
    indices: &[usize],
    grain: usize,
) {
    let mut product: Option<MillerLoopResult> = None;
    for (points, prepared_indices) in
        g1.chunks(grain).zip(indices.chunks(grain))
    {
        let range_product = batch
            .miller_loop_product_indices(points, prepared_indices)
            .expect("prepared Miller range");
        product = Some(match product {
            Some(accumulator) => accumulator
                .multiply(&range_product)
                .expect("Miller-product multiplication"),
            None => range_product,
        });
    }
    black_box(
        product
            .expect("the benchmark batch is nonempty")
            .final_exponentiation()
            .expect("one final exponentiation"),
    );
}

fn evaluate_prepared_ranges_parallel(
    pool: &ThreadPool,
    batch: &PreparedG2Batch,
    g1: &[G1Projective],
    indices: &[usize],
    grain: usize,
) {
    let product = pool
        .install(|| {
            g1.par_chunks(grain)
                .zip(indices.par_chunks(grain))
                .map(|(points, prepared_indices)| {
                    batch
                        .miller_loop_product_indices(points, prepared_indices)
                        .expect("prepared Miller range")
                })
                .reduce_with(|left, right| {
                    left.multiply(&right)
                        .expect("Miller-product multiplication")
                })
        })
        .expect("the benchmark batch is nonempty");
    black_box(
        product
            .final_exponentiation()
            .expect("one final exponentiation"),
    );
}

fn bench_prepared_pairing_product(c: &mut Criterion) {
    let g1_generator =
        G1Projective::generator().expect("native AVX G1 generator");
    let g2_generator =
        G2Projective::generator().expect("native AVX G2 generator");
    let pool2 = rayon::ThreadPoolBuilder::new()
        .num_threads(2)
        .build()
        .expect("two-worker benchmark pool");
    let pool4 = rayon::ThreadPoolBuilder::new()
        .num_threads(4)
        .build()
        .expect("four-worker benchmark pool");
    let mut group = c.benchmark_group("btx_prepared_pairing_product");
    group.sample_size(10);
    group.warm_up_time(Duration::from_millis(500));
    group.measurement_time(Duration::from_secs(2));

    for count in [32usize, 64, 128, 256] {
        let g1 = vec![g1_generator; count];
        let g2 = vec![g2_generator; count];
        let indices: Vec<_> = (0..count).rev().collect();
        let prepared =
            PreparedG2Batch::new(&g2).expect("prepared native G2 batch");
        let mut legacy =
            PairingWorkspace::for_pairs(count).expect("pairing workspace");

        group.bench_with_input(
            BenchmarkId::new("unprepared", count),
            &count,
            |b, _| {
                b.iter(|| {
                    legacy.clear();
                    for (g1, g2) in g1.iter().zip(&g2) {
                        legacy
                            .push(black_box(g1), black_box(g2))
                            .expect("stage pairing");
                    }
                    black_box(legacy.evaluate().expect("pairing product"));
                });
            },
        );
        group.bench_with_input(
            BenchmarkId::new("prepared_full", count),
            &count,
            |b, _| {
                b.iter(|| {
                    black_box(
                        prepared
                            .miller_loop_product_indices(
                                black_box(&g1),
                                black_box(&indices),
                            )
                            .expect("prepared Miller product")
                            .final_exponentiation()
                            .expect("final exponentiation"),
                    );
                });
            },
        );

        for grain in [1usize, 2, 4, 8, 9, 10, 11, 12, 16, 32] {
            group.bench_with_input(
                BenchmarkId::new(format!("prepared_grain_{grain}"), count),
                &count,
                |b, _| {
                    b.iter(|| {
                        evaluate_prepared_ranges(
                            &prepared,
                            black_box(&g1),
                            black_box(&indices),
                            grain,
                        );
                    });
                },
            );
        }

        for (workers, pool) in [(2usize, &pool2), (4, &pool4)] {
            group.bench_with_input(
                BenchmarkId::new(format!("parallel_{workers}_grain_10"), count),
                &count,
                |b, _| {
                    b.iter(|| {
                        evaluate_prepared_ranges_parallel(
                            pool,
                            &prepared,
                            black_box(&g1),
                            black_box(&indices),
                            10,
                        );
                    });
                },
            );
        }
    }

    group.finish();
}

fn bench_prepared_pairing_range(c: &mut Criterion) {
    let g1_generator =
        G1Projective::generator().expect("native AVX G1 generator");
    let g2_generator =
        G2Projective::generator().expect("native AVX G2 generator");
    let g1 = vec![g1_generator; 16];
    let g2 = vec![g2_generator; 16];
    let indices: Vec<_> = (0..16).rev().collect();
    let prepared = PreparedG2Batch::new(&g2).expect("prepared native G2 batch");
    let mut group = c.benchmark_group("btx_prepared_pairing_range");
    group.sample_size(20);
    group.warm_up_time(Duration::from_millis(500));
    group.measurement_time(Duration::from_secs(2));

    for grain in [8usize, 9, 10, 11, 12, 16] {
        group.bench_with_input(
            BenchmarkId::from_parameter(grain),
            &grain,
            |b, &grain| {
                b.iter(|| {
                    black_box(
                        prepared
                            .miller_loop_product_indices(
                                black_box(&g1[..grain]),
                                black_box(&indices[..grain]),
                            )
                            .expect("prepared Miller range"),
                    );
                });
            },
        );
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_prepared_pairing_product,
    bench_prepared_pairing_range
);
criterion_main!(benches);
