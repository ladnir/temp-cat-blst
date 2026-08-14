#!/usr/bin/env bash
# Copyright 2026 Category Labs
# SPDX-License-Identifier: Apache-2.0

set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
out_dir="$repo_root/out/fuzz-inventory"
scalar_config="$out_dir/config-scalar"
portable_config="$out_dir/config-portable-avx"
rust_target="$out_dir/rust-target"

mkdir -p "$out_dir"

cmake -S "$repo_root" -B "$scalar_config" -G Ninja \
    -DBLST_TEST=OFF \
    -DBLST_BENCHMARK=OFF \
    -DBLST_ASM=OFF \
    -DBLST_DEFAULT_ASM=OFF \
    -DBLST_ADX=OFF \
    -DBLST_AVX512=OFF \
    -DBLST_AVX512_PORTABLE=OFF \
    -DBLST_CUDA=OFF

cmake -S "$repo_root" -B "$portable_config" -G Ninja \
    -DCMAKE_BUILD_TYPE=Debug \
    -DBLST_TEST=OFF \
    -DBLST_BENCHMARK=OFF \
    -DBLST_ASM=OFF \
    -DBLST_DEFAULT_ASM=OFF \
    -DBLST_ADX=OFF \
    -DBLST_AVX512=OFF \
    -DBLST_AVX512_PORTABLE=ON \
    -DBLST_DEFAULT_AVX=ON \
    -DBLST_CUDA=OFF

ast_kinds='(CXXMethodDecl|CXXConstructorDecl|CXXConversionDecl|CXXDestructorDecl|FunctionDecl|UsingDecl|UsingShadowDecl|CXXRecordDecl)'

clang++ -std=c++20 \
    -I"$repo_root" -I"$scalar_config" -I"$scalar_config/blst" \
    -Xclang -ast-dump -Xclang -ast-dump-filter=blst \
    -fsyntax-only "$repo_root/fuzz/inventory/cpp_surface.cpp" \
    | grep -E "$ast_kinds" > "$out_dir/cpp-scalar-compact.txt"

clang++ -std=c++20 \
    -I"$repo_root" -I"$portable_config" -I"$portable_config/blst" \
    -Xclang -ast-dump -Xclang -ast-dump-filter=blst::avx \
    -fsyntax-only "$repo_root/fuzz/inventory/cpp_surface.cpp" \
    | grep -E "$ast_kinds" > "$out_dir/cpp-avx-compact.txt"

{
    printf '# profile: scalar\n'
    cat "$out_dir/cpp-scalar-compact.txt"
    printf '# profile: avx\n'
    cat "$out_dir/cpp-avx-compact.txt"
} > "$out_dir/cpp-compact.txt"

export CARGO_BUILD_JOBS=2
cargo +nightly-2026-01-15 rustdoc \
    --manifest-path "$repo_root/bindings/rust/Cargo.toml" \
    --target-dir "$rust_target" \
    --lib \
    --features full-rust-bindings,portable,serde-secret \
    -- -Z unstable-options --output-format json

python3 "$repo_root/scripts/fuzz/api_inventory.py" \
    --root "$repo_root" --require-generated --write-snapshots --write
python3 "$repo_root/scripts/fuzz/api_inventory.py" \
    --root "$repo_root" --require-generated
