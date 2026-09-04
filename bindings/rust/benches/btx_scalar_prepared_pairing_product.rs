// Copyright (C) 2026 Category Labs, Inc.

use std::time::Duration;

use blst::scalar::{
    G1Affine, G1Projective, G2Affine, G2Projective, MillerLoopResult,
    PreparedG2Batch, Scalar,
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
    let product = g1
        .chunks(grain)
        .zip(indices.chunks(grain))
        .map(|(points, prepared_indices)| {
            batch
                .miller_loop_product_indices(points, prepared_indices)
                .expect("prepared scalar Miller range")
        })
        .fold(MillerLoopResult::identity(), |left, right| left * right);
    black_box(product.final_exponentiation());
}

fn evaluate_prepared_ranges_parallel(
    pool: &ThreadPool,
    batch: &PreparedG2Batch,
    g1: &[G1Projective],
    indices: &[usize],
    grain: usize,
) {
    let product = pool.install(|| {
        g1.par_chunks(grain)
            .zip(indices.par_chunks(grain))
            .map(|(points, prepared_indices)| {
                batch
                    .miller_loop_product_indices(points, prepared_indices)
                    .expect("prepared scalar Miller range")
            })
            .reduce(MillerLoopResult::identity, |left, right| left * right)
    });
    black_box(product.final_exponentiation());
}

fn bench_scalar_prepared_pairing_product(c: &mut Criterion) {
    let g1_value = G1Projective::generator() * Scalar::from(17);
    let g2_value = G2Projective::generator() * Scalar::from(29);
    let pool2 = rayon::ThreadPoolBuilder::new()
        .num_threads(2)
        .build()
        .expect("two-worker benchmark pool");
    let pool4 = rayon::ThreadPoolBuilder::new()
        .num_threads(4)
        .build()
        .expect("four-worker benchmark pool");
    let mut group = c.benchmark_group("btx_scalar_prepared_pairing_product");
    group.sample_size(10);
    group.warm_up_time(Duration::from_millis(500));
    group.measurement_time(Duration::from_secs(2));

    for count in [32usize, 64, 128, 256] {
        let g1 = vec![g1_value; count];
        let g2 = vec![g2_value; count];
        let indices: Vec<_> = (0..count).rev().collect();
        let prepared =
            PreparedG2Batch::new(&g2).expect("prepared scalar G2 batch");
        let mut g1_affine = vec![G1Affine::identity(); count];
        let mut g2_affine = vec![G2Affine::identity(); count];

        group.bench_with_input(
            BenchmarkId::new("unprepared", count),
            &count,
            |b, _| {
                b.iter(|| {
                    for (output, point) in g1_affine.iter_mut().zip(&g1) {
                        *output = black_box(point).to_affine();
                    }
                    for (output, point) in g2_affine.iter_mut().zip(&g2) {
                        *output = black_box(point).to_affine();
                    }
                    black_box(
                        MillerLoopResult::from_pairs(
                            black_box(&g1_affine),
                            black_box(&g2_affine),
                        )
                        .expect("unprepared scalar Miller product")
                        .final_exponentiation(),
                    );
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
                            .expect("prepared scalar Miller product")
                            .final_exponentiation(),
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

        for grain in [2usize, 4, 8, 10, 16] {
            for (workers, pool) in [(2usize, &pool2), (4, &pool4)] {
                group.bench_with_input(
                    BenchmarkId::new(
                        format!("parallel_{workers}_grain_{grain}"),
                        count,
                    ),
                    &count,
                    |b, _| {
                        b.iter(|| {
                            evaluate_prepared_ranges_parallel(
                                pool,
                                &prepared,
                                black_box(&g1),
                                black_box(&indices),
                                grain,
                            );
                        });
                    },
                );
            }
        }
    }

    group.finish();
}

fn bench_scalar_prepared_pairing_range(c: &mut Criterion) {
    let g1 = vec![G1Projective::generator() * Scalar::from(17); 16];
    let g2 = vec![G2Projective::generator() * Scalar::from(29); 16];
    let indices: Vec<_> = (0..16).rev().collect();
    let prepared = PreparedG2Batch::new(&g2).expect("prepared scalar G2 batch");
    let mut group = c.benchmark_group("btx_scalar_prepared_pairing_range");
    group.sample_size(20);
    group.warm_up_time(Duration::from_millis(500));
    group.measurement_time(Duration::from_secs(2));

    for grain in [1usize, 2, 3, 4, 8, 9, 10, 11, 12, 16] {
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
                            .expect("prepared scalar Miller range"),
                    );
                });
            },
        );
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_scalar_prepared_pairing_product,
    bench_scalar_prepared_pairing_range
);
criterion_main!(benches);
