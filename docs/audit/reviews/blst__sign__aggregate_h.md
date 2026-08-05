# File Review: Native pairing context storage

## Record

| Field | Value |
| --- | --- |
| Path | `blst/sign/aggregate.h` |
| Roadmap phase / queue | Phase 15 curves and protocols, tenth packet |
| Reviewed commit | `091dfb2cae096ee8ed801036e98827b3253b2e84` |
| Reviewer | Codex audit team |
| Date | 2026-08-05 |
| Disposition | in-scope production |
| Subsystem | native opaque pairing-context representation |
| Profiles | scalar and selected AVX implementation |
| Upstream reference | pinned BLST `src/aggregate.c` and `src/server.c` |

## Review

The header contains only internal native storage and the two declarations needed
across CatBLST translation units. `PAIRING` is never a public field-access API:
the C surface exposes queried size/alignment and opaque caller storage. Its
selected-backend `MillerResult`, independent live G1/G2 signature objects, and
fixed arrays have valid C++ lifetimes. The typed batch capacity is exactly eight
and cannot be altered by an external preprocessor definition. Scalar and selected
portable-vector builds agree on the context queried by their matching library;
no cross-configuration binary layout promise is made.

The only exported native accessor commits and writes a pre-final Miller product
to caller storage. It neither allocates nor enlarges the context. The full state
machine and field meanings are documented adjacent to the structure.

## Evidence and disposition

Scalar MSVC Debug/Release and selected portable-AVX Debug build and pass the six
focused aggregate regressions. Static assertions pin capacity eight and the C
getter signature. Structural checks find no `N_MAX` or old accessor. `AUD-152`,
`AUD-157`, `AUD-158`, and `AUD-159` are validated. Native AVX-512 hardware,
sanitizer/fuzz, Apple, and full milestone execution remain release evidence.

Complete Phase 15 file review at implementation commit
`091dfb2cae096ee8ed801036e98827b3253b2e84`.
