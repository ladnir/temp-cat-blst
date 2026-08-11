# File Review: `blst/fft/ctFft.h`

## Record

| Field | Value |
| --- | --- |
| Path | `blst/fft/ctFft.h` |
| Roadmap phase / queue | Phase 17, core FFT owner after middle-product audit |
| Reviewed commit | `9d1e2bd60e610e88b44c1be36daea7be793898c0` |
| Reviewer | Codex audit team |
| Date | 2026-08-11 |
| Disposition | `in_scope_production` |
| Subsystem | Generic radix-two FFT, truncated transform plans, and polynomial FFT composition |
| Profiles | scalar and AVX-selected native C++; C/Rust Fr adapters use the same engine |
| Generation/provenance | Locally composed handwritten C++ templates over BLST/CatBLST arithmetic contexts |
| Upstream references | `not_applicable: upstream BLST has no corresponding FFT engine` |
| Related files | `blst/fft/context.h`, `blst/fft/domain.h`, `blst/fft/middle_product.h`, `blst/fft/fr_interpolation.cpp`, `blst/fft_c_api.cpp`, `bindings/blst_fft.h`, `bindings/rust/src/fft.rs`, `test/Fft_api_tests.cpp`, `test/Fft_tests.cpp`, frozen BTx/BTE callers |

## Purpose and trust boundary

The file owns the one retained CPU radix-two transform, reusable public-shape
plans for zero-padded forward transforms and selected inverse output ranges,
and an internal cyclic polynomial multiplication composition. It operates on
generic coefficient/group types through a compile-time context and performs no
runtime backend selection.

Transform values may be public or secret-derived. The transform schedule,
indices, roots, and plan contents depend only on public lengths/ranges. The
arithmetic context remains responsible for the side-channel contract of its
`plus`, `minus`, `mult`, `inverse`, and `from_u64` operations.

## Contracts and invariants

- Transform length is a nonzero power of two no greater than the largest power
  of two representable by the signed-index engine (`2^30` on supported hosts).
- A directional root table contains at least `n` ordered powers starting at
  one. `fft(..., invert=true)` receives inverse roots; it only adds inverse
  scaling. The prepared-forward-root operation derives reciprocal indices.
- Mutable data, root storage, and a prepared inverse-size object do not overlap.
- Forward truncation requires `a[xn..n)` to be additive identity and publishes
  all n frequency outputs.
- Inverse range truncation consumes all n frequency inputs, publishes only the
  normalized requested range, and leaves every other element unspecified.
- Plan builders validate public shapes before allocation and may throw
  `std::invalid_argument` or an allocation exception. Executors require an
  exact immutable matching plan; changing schedule fields violates the narrow
  native contract.
- Polynomial FFT multiplication uses equal, disjoint buffers. A skipped B FFT
  means B already contains the complete evaluation-domain transform.
- Hot transforms allocate no memory. Arithmetic exceptions, if supplied by a
  generic context, propagate and may leave the in-place buffer partially
  transformed. Supported native contexts do not throw in their hot arithmetic.
- Optional timing is caller-owned, unsynchronized, and diagnostic. A null sink
  performs no key construction; enabled allocation failure is contained.

The C and safe Rust Fr surfaces are documented in their own reviewed adapters.
They validate untrusted lengths and translate exceptions before crossing the C
ABI. These generic native templates are internal expert interfaces rather than
separately exported ABI symbols.

## Dependency and call-path map

- `blst/fft/domain.h` prepares reusable forward roots and inverse-size storage.
- `blst/fft_c_api.cpp` uses ordinary and prepared Fr transforms after C-boundary
  validation.
- `blst/fft/fr_interpolation.cpp` uses ordinary Fr transforms and polynomial
  multiplication.
- `blst/fft/middle_product.h` stores and reuses forward/inverse truncated plans
  for scalar and AVX-selected G1/G2/GT work.
- Frozen BTx/BTE source calls only ordinary, polynomial, and truncated APIs. It
  never called the removed radix-4/radix-8 or reversal helpers and was not
  modified or audited in this packet.

## Upstream/specification comparison

Upstream BLST contains no corresponding general FFT module. The ordinary
transform is the standard in-place radix-two Cooley-Tukey DIT construction
with binary bit reversal. The deleted radix-4/radix-8 functions were local
experiments introduced in `9525694`; later pruning removed their dispatcher but
left the unused kernels. A prime-field oracle showed they were not equivalent
to the retained transform, so no compatibility surface justified repair.

## Checklist evidence

### 1. Identity and purpose

The entire 897-line pre-remediation file, history, production callers, tests,
benchmark references, and read-only frozen BTx/BTE uses were inventoried.

### 2. Contracts and invariants

Size, roots, direction, zero-padding, output publication, aliasing, plan
immutability, allocation, exception, and destructive-update contracts are
attached to the retained declarations and summarized above.

### 3. Upstream and specification lineage

The module is local with no BLST counterpart. Experimental history and Apache
2.0/Category Labs provenance are recorded in the source.

### 4. Correctness

Active radix-two indexing and terminal-shift handling were reviewed. Reusable
forward/inverse plans match full transforms over prefix, interior, suffix, and
full shapes. Dead incorrect high-radix paths were deleted.

### 5. Memory and lifetime safety

Hot transforms allocate nothing. Plan allocation occurs before mutation.
Overflow-safe range validation precedes plan indexing. Plans own their vectors
and executors retain no caller pointers.

### 6. Undefined and implementation-defined behavior

The shared size predicate protects the signed `int` transform indices.
Subtraction-form range validation eliminates unsigned wrap. Root and plan
shape checks are Debug preconditions in hot expert templates; supported C
boundaries validate in Release.

### 7. Side channels and cryptographic hygiene

All control flow and memory indices depend on public transform shapes and root
tables. Generic arithmetic timing is delegated explicitly to `CoeffCtx`; no
data inspection, randomization, or secret-dependent pruning occurs.

### 8. Randomness and secret lifecycle

N/A: no entropy or persistent secret storage. The transform is in place and
does not make hidden copies of coefficient arrays.

### 9. API and ABI

Templates have no exported binary ABI. Supported C and Rust Fr adapters route
to this engine but retain their separate validation and error contracts.

### 10. Concurrency

There is no global state. Completed plans and roots may be shared by concurrent
readers; data buffers, mutable plans, and one timing sink require caller
synchronization.

### 11. Resource exhaustion and denial of service

Plan builders reject invalid/unsupported public sizes before vector
allocation. Valid maximum plans can still be very large; public C APIs own
their documented allocation/failure translation.

### 12. Portability and build behavior

Supported hosts are 64-bit. The signed-index maximum is explicit. MSVC 19.50
Debug rebuilt scalar-default, portable AVX headers, and unchanged frozen BTx/BTE
consumers. CUDA is outside scope.

### 13. Performance-sensitive structure

The radix-two inner loop, skipped `w == 1` multiplication, compile-time context
dispatch, and precomputed plan schedules are retained. Added hot checks are
`BLST_ASSERT` only. Metric key construction occurs only for a nonnull sink and
is contained by its no-throw barrier.

### 14. Tests and evidence

Focused tests cover invalid and overflowed shapes, Debug hot preconditions,
plan accounting, reusable plan/full-transform equivalence, prepared-domain
equivalence, and optional timing. Existing middle-product compilation covers
scalar and portable-selected context instantiation.

### Specialized additions

The active transform is radix two only. Introducing another radix or
permutation requires an independent oracle, production caller, focused tests,
and performance evidence before inclusion.

## Findings

`AUD-032`, `AUD-046`, `AUD-285`, and `AUD-292` through `AUD-296`.

## Validation

- `cmake --build out/build/x64-Debug-Scalar --target blst_tests --config Debug --parallel 2`: passed.
- `blst_tests.exe --gtest_filter=FftCore.*:FftApi.*:Timer.*`: 8/8 passed.
- MSVC standalone composed-key allocation-failure probe: passed.
- Structural scan: no `fft_radix4`, `fft_radix8`, `derive_ntt_of_reversed`,
  `fft_stockham`, or `DefaultCoeffCtx` remains outside generated Cargo output.
- No benchmark, CUDA, full suite, sanitizer, native AVX-512, remote, or Apple
  run was performed.

## Automated evidence

| Tool/profile | Version and configuration | Result or finding IDs | Blind spots / N/A reason |
| --- | --- | --- | --- |
| CodeQL | Earlier Phase 12 database | `AUD-032`, `AUD-046` | Does not prove FFT mathematics or constant-time scheduling |
| Clang Static Analyzer / clang-tidy | Earlier retained profiles | No additional mapped diagnostic | Current packet validated with MSVC only |
| Clippy | N/A | N/A | Native template engine; Rust adapter reviewed separately |
| ASan/UBSan | Earlier focused Phase 17 ASan | No mapped active radix-two failure | Not rerun after this cleanup |
| Fuzzing | Not yet run | Phase 17 residual | Differential FFT fuzzing remains planned |

## Residual risk and decision

The mathematical validity/order of a caller-supplied native root table is a
narrow precondition rather than an O(n) runtime check. Hot Release templates
likewise do not validate plan internals or overlap. Public C/Rust adapters
validate their exposed shapes. The generic names remain in the global namespace
for compatibility with frozen BTx/BTE callers; namespace cleanup is deferred
until that integration can move without editing frozen code. Final manifest
closure awaits later sanitizer/native-AVX milestone evidence.
