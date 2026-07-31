# Phase 5 Uniqueness-Tree Hardening Evidence

## Record

- Date: 2026-07-30
- Baseline commit: `0b99644`
- Implementation commit: `1e43e0e`
- State: validated and committed
- Scope: `AUD-011` and conceptual upstream row `06a5711`
- Benchmarks: disabled and not run
- CUDA: disabled

## Root cause

The retained uniqueness tree accepted an allocation sized for a caller-known
node count but did not store or enforce that count. A caller mistake or extra
insertion therefore indexed beyond the trailing node allocation. Its fixed
path arrays were also written without a depth check, allocation arithmetic
could overflow, message length was shifted into a colour bit without a range
check, and unequal lengths narrowed `size_t` through `long`.

The boolean insertion result further collapsed duplicate and internal failure,
making fail-closed handling ambiguous.

## Breaking API

The unsafe legacy interface was removed:

- `blst_uniq_sizeof(capacity)` returns zero for zero capacity or exact-size
  overflow;
- `blst_uniq_alignof()` reports caller-storage alignment;
- `blst_uniq_init(tree, capacity)` validates storage and stores capacity;
- `blst_uniq_insert()` returns inserted, duplicate, capacity exceeded, depth
  exceeded, or invalid input; and
- `blst_uniq_test` no longer exists.

The public header documents initialization, external synchronization,
null/length rules, and the requirement that borrowed message bytes remain
readable and unchanged for the tree lifetime. Rust bindgen declarations and
the retained safe `uniq` helper were migrated in the same commit.

## Storage and traversal proof

The exact layout is `offsetof(rb_tree, nodes) + capacity * sizeof(node)`, with
addition and multiplication checked before returning a size. Insertion checks
capacity before indexing the trailing array. Duplicate comparison occurs first
so a duplicate remains distinguishable when full.

Red-black height is at most `2 * log2(n + 1)`. The guarded path has 64 entries:

- on 32-bit ABIs, a 12-byte header and 16-byte node limit any representable
  allocation to fewer than `2^28` nodes, hence height below 56;
- on 64-bit ABIs, reaching depth greater than 64 requires more than `2^32`
  nodes in the worst case and is rejected before either path-array write.

This preserves the original 64-bit stack footprint while fixing the feasible
32-bit overrun. The insertion algorithm remains iterative and allocation-free.

## Validation

The C regressions cover:

- inserted, duplicate, capacity-exceeded, and invalid-input outcomes;
- duplicate lookup after reaching capacity;
- zero and `SIZE_MAX` capacities plus a binary-searched exact overflow
  boundary;
- null and misaligned tree storage;
- null nonempty messages, null empty messages, and lengths above
  `SIZE_MAX >> 1`;
- 32,768 ascending eight-byte messages, exercising adversarial insertion order
  and subsequent duplicate probes; and
- exact storage alignment and exported-symbol behavior.

| Profile | Configuration | Result |
| --- | --- | --- |
| Windows scalar direct ADX | MSVC Debug, `BLST_PORTABLE=OFF`, `BLST_ADX=ON` | 273 passed, nine expected AVX-only skips, zero failed (282 registered) |
| Windows scalar portable dispatch | MSVC Debug, `BLST_PORTABLE=ON`, `BLST_ADX=OFF` | 315 passed, nine expected AVX-only skips, zero failed (324 registered) |
| Windows AVX-512 portable | MSVC Debug, portable AVX-512 default | Three focused uniqueness tests passed; Phase 4 full baseline remains 144/144 |
| Linux scalar portable dispatch | GCC 13 Debug, ASan and leak detection | Three focused uniqueness tests passed with no sanitizer or leak report |
| Rust native | Full bindings | 71 unit tests and one integration test passed |
| Rust portable | Full bindings plus `portable` | Two focused uniqueness tests passed |
| Rust no-std compile | Full bindings, no default features | `cargo check` passed |
| Public headers | GCC 13 C11 and C++23, `-Wall -Wextra -Werror` | Standalone umbrella-header consumers compiled |
| 32-bit model | GCC 13 `-m32 -std=c++23` | `rb_tree.cpp` compiled to an ELF i386 object |

MSVC archive inspection found the four new C symbols and no legacy
`blst_uniq_test` symbol. `cargo fmt --check`, `git diff --check`, and the audit
integrity checker passed. No benchmark was run.

## Residual validation

- ARM64 and WebAssembly were not executed.
- The checked-in Rust FFI artifact compiled and passed cross-language tests,
  but bindgen regeneration was not run because the bindgen CLI is unavailable
  on this host. `bindings/rust/src/bindings.rs` remains pending for its full
  generator/artifact file review.
- A tree deep enough to return `BLST_UNIQ_DEPTH_EXCEEDED` is not practically
  allocatable on the validation host; the guard is directly before both array
  writes and the architecture bounds are reviewed above.
- The tree intentionally borrows message data. C cannot enforce that lifetime,
  so misuse remains a documented caller error; Rust's slice borrow lasts
  through the complete helper call.
- This phase fully reviews `blst/common/rb_tree.cpp` and its new internal
  header. It reviews only the uniqueness-related slices of the larger C and
  Rust API files, which remain pending in the file-by-file manifest.
