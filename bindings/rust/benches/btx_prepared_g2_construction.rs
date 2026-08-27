// Copyright (C) 2025 Category Labs, Inc.

use std::time::Duration;

use blst::avx::G2Projective;
use blst::PreparedG2Batch;
use criterion::{
    black_box, criterion_group, criterion_main, BenchmarkId, Criterion,
};

fn bench_prepared_g2_construction(c: &mut Criterion) {
    let generator = G2Projective::generator().expect("native AVX G2 generator");
    let mut group = c.benchmark_group("btx_prepared_g2_construction");
    group.sample_size(10);
    group.warm_up_time(Duration::from_millis(500));
    group.measurement_time(Duration::from_secs(3));

    for count in [256usize, 1024] {
        // The input is resident native state and remains outside the timed
        // region. Construction measures only allocation, affine conversion,
        // Miller-line preparation, and release of the completed batch.
        let points = vec![generator; count];
        group.bench_function(BenchmarkId::from_parameter(count), |b| {
            b.iter(|| {
                let prepared = PreparedG2Batch::new(black_box(&points))
                    .expect("native prepared-G2 batch construction");
                assert_eq!(prepared.len(), count);
                black_box(prepared);
            });
        });
    }

    group.finish();
}

criterion_group!(benches, bench_prepared_g2_construction);
criterion_main!(benches);
