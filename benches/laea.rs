use criterion::{black_box, criterion_group, criterion_main, Criterion};
use geomatic::laea::{backward, backward_ref, forward, forward_ref};
use geomatic::{Point3035, Point4326};

pub fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("forward_fast", |b| b.iter(|| forward(black_box(Point4326::new(50.0, 5.0)))));
    c.bench_function("forward_ref", |b| b.iter(|| forward_ref(black_box(Point4326::new(50.0, 5.0)))));
    c.bench_function("backward", |b| b.iter(|| backward(black_box(Point3035::new(3962799.45, 2999718.85)))));
    c.bench_function("backward_ref", |b| b.iter(|| backward_ref(black_box(Point3035::new(3962799.45, 2999718.85)))));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
