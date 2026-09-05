// Copyright 2026 Category Labs
// SPDX-License-Identifier: Apache-2.0

use std::time::Duration;

use blst::{
    avx::{G1MsmWorkspace, G1Projective},
    Scalar,
};
use criterion::{
    black_box, criterion_group, criterion_main, BenchmarkId, Criterion,
};

fn benchmark_scalars(count: usize) -> Vec<Scalar> {
    (0..count)
        .map(|index| {
            let mut wide = [0u8; 64];
            for (byte_index, byte) in wide.iter_mut().enumerate() {
                *byte = (index as u8)
                    .wrapping_mul(131)
                    .wrapping_add((byte_index as u8).wrapping_mul(29))
                    .wrapping_add(17);
            }
            Scalar::from_bytes_wide_le(&wide)
        })
        .collect()
}

fn bench_avx_g1_msm(c: &mut Criterion) {
    let generator = G1Projective::generator().expect("native AVX G1 generator");
    for count in [2usize, 4] {
        let workspace =
            G1MsmWorkspace::for_points(count).expect("small MSM workspace");
        println!(
            "AVX_G1_MSM_SCRATCH points={count} bytes={}",
            workspace.capacity_bytes()
        );
    }
    let mut group = c.benchmark_group("avx_g1_msm");
    group.sample_size(10);
    group.warm_up_time(Duration::from_millis(500));
    group.measurement_time(Duration::from_secs(2));

    for count in [
        2usize, 4, 8, 16, 24, 32, 40, 48, 64, 80, 96, 128, 160, 192, 224, 256,
    ] {
        let point_scalars = benchmark_scalars(count);
        let points = point_scalars
            .iter()
            .map(|scalar| {
                generator
                    .multiply(scalar)
                    .expect("native AVX G1 point construction")
            })
            .collect::<Vec<_>>();
        let scalars = benchmark_scalars(count + 1)
            .into_iter()
            .skip(1)
            .collect::<Vec<_>>();
        let scalars_le =
            scalars.iter().map(Scalar::to_bytes_le).collect::<Vec<_>>();
        let mut constant_time_workspace =
            G1MsmWorkspace::for_points(count).expect("constant-time workspace");
        let mut variable_time_workspace =
            G1MsmWorkspace::for_points(count).expect("variable-time workspace");
        let mut pippenger_workspace =
            G1MsmWorkspace::for_points(count).expect("Pippenger workspace");
        let mut constant_time_output =
            G1Projective::identity().expect("constant-time output");
        let mut variable_time_output =
            G1Projective::identity().expect("variable-time output");
        constant_time_workspace
            .secret_msm_into(&mut constant_time_output, &points, &scalars)
            .expect("constant-time reference");
        variable_time_workspace
            .msm_into(&mut variable_time_output, &points, &scalars)
            .expect("variable-time reference");
        assert_eq!(constant_time_output, variable_time_output);
        let mut pippenger_output =
            G1Projective::identity().expect("Pippenger output");
        pippenger_workspace
            .pippenger_msm_into(&mut pippenger_output, &points, &scalars_le)
            .expect("Pippenger reference");
        assert_eq!(constant_time_output, pippenger_output);

        group.bench_with_input(
            BenchmarkId::new("constant_time", count),
            &count,
            |b, _| {
                b.iter(|| {
                    constant_time_workspace
                        .secret_msm_into(
                            black_box(&mut constant_time_output),
                            black_box(&points),
                            black_box(&scalars),
                        )
                        .expect("constant-time G1 MSM");
                });
            },
        );
        group.bench_with_input(
            BenchmarkId::new("public_vartime", count),
            &count,
            |b, _| {
                b.iter(|| {
                    variable_time_workspace
                        .msm_into(
                            black_box(&mut variable_time_output),
                            black_box(&points),
                            black_box(&scalars),
                        )
                        .expect("public variable-time G1 MSM");
                });
            },
        );
        group.bench_with_input(
            BenchmarkId::new("pippenger", count),
            &count,
            |b, _| {
                b.iter(|| {
                    pippenger_workspace
                        .pippenger_msm_into(
                            black_box(&mut pippenger_output),
                            black_box(&points),
                            black_box(&scalars_le),
                        )
                        .expect("Pippenger G1 MSM");
                });
            },
        );
    }

    group.finish();
}

criterion_group!(benches, bench_avx_g1_msm);
criterion_main!(benches);
