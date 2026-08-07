# File Review: `bindings/blst_aux.h`

## Record

| Field | Value |
| --- | --- |
| Path | `bindings/blst_aux.h` |
| Roadmap phase | Phase 18 public C/C++/Rust surfaces |
| Reviewed implementation commit | `98faac3fb001f59f72fd1b131c278c101953f1a0` |
| Reviewer | Codex audit team |
| Date | 2026-08-06 |
| Disposition | in-scope production |
| Profiles | C99-or-later, native C++, scalar and AVX declarations |
| Included surfaces | scalar, AVX, middle product, FFT |

This record reviews the exact implementation commit above.

## Purpose and include contract

The historical BLST auxiliary header is now CatBLST's compatibility aggregator
for the split public C headers. It must expose each supported split operation
surface without duplicating declarations, imposing include order on consumers,
or changing the configuration and representation policy owned by those files.

It intentionally includes `blst_scalar.h` before the expert AVX and algorithm
headers. Every included header independently obtains `blst_types.h`, so the
order does not provide a hidden type dependency. Repeated inclusion is safe.

## Complete review

1. The entire file and every include edge were inspected.
2. `BLST_AUX_H` is project-prefixed and not a reserved identifier.
3. The file contains no declaration, definition, executable behavior, mutable
   state, allocation, exception path, or secret-dependent operation.
4. Scalar, AVX, middle-product, and FFT headers have independent guards and no
   edge back to `blst_aux.h`; the include graph is acyclic.
5. Adding FFT only exposes declarations already shipped and linked by CatBLST;
   it neither creates symbols nor changes ABI or runtime selection.
6. The local split predates the FFT header. `AUD-185` resulted when commit
   `9daaf27` added `blst_fft.h` without extending the earlier `c71496b`
   compatibility list.
7. Direct split-header inclusion remains recommended for consumers that want
   the smallest parse surface.

## Validation and residual risk

The strict standalone C11/C++20 compile, link, and execution evidence is shared
with the `bindings/blst.h` review because that file reaches this aggregator
directly. The binding-header regression proves the transitive FFT edge. Full
packaging and non-MSVC supported-compiler consumers remain Phase 18 evidence;
the semantic contracts of each included surface remain with their own reviews.
