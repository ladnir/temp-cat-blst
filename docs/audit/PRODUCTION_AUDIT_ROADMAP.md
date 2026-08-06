# Production Audit Roadmap

## Program boundary

Phases 0-11 were preparation: upstream inventory and translation, urgent
correctness/security hardening, generator repair, unsupported-surface cleanup,
platform scope control, and conformance evidence. The main production audit is
Phases 12-20. It contains most of the remaining work.

The program has four continuous workstreams:

1. Manual file-by-file review and contract reconstruction.
2. Automated static analysis, sanitizers, fuzzing, and differential/oracle
   coverage.
3. Finding triage, remediation, regression tests, and performance review.
4. Durable evidence, CI, platform qualification, packaging, and release gates.

No workstream substitutes for another. A static-analysis baseline does not
complete a file; a clean fuzz run does not prove arithmetic or constant-time
behavior; manual review does not replace executed UB and boundary evidence.

## Subsystem review cycle

Every production subsystem follows the same cycle:

1. Freeze the phase queue at an exact commit and map callers, callees,
   generators, public consumers, and scalar/AVX relationships.
2. Run the applicable automated baseline and triage results into findings,
   duplicates, or documented false positives.
3. Review files individually using `CHECKLIST.md` and
   `FILE_REVIEW_TEMPLATE.md`.
4. Implement findings in separable commits with deterministic regressions and
   performance/code-generation review where hot paths are affected.
5. Expand sanitizer, fuzz, property, independent-oracle, cross-language, and
   scalar/AVX coverage around the now-understood contracts.
6. Run one batched subsystem gate and record exact evidence in `PHASE_*.md`.
7. Mark files complete only when their manifest rows and review records name
   the same reviewed commit and pass audit integrity.

If a later change touches a completed file's audited behavior, the review must
be reopened or receive a linked delta review before the release claim can use
it. Shared analyzer/test results may support several files, but every record
must explain the file-specific coverage and blind spots.

## Main phases

| Phase | Scope | Primary automation/evidence | Closure result |
| --- | --- | --- | --- |
| 12. Audit platform | CodeQL, CodeChecker/Clang Static Analyzer, clang-tidy, Clippy, ASan/UBSan, fuzz infrastructure, corpora/artifacts, triage policy, immutable CI provenance | Initial repository baseline and public byte-boundary fuzzers | Reproducible platform; baseline fully triaged; tools remain active |
| 13. Fundamental scalar trust chain | Common types, limbs, bytes/endian, configuration, portable backend, assembly generators/outputs, CPU selection | Integer/conversion oracles, UB/bounds tests, generation proof, portable/assembly differential checks | Fundamental representations and scalar backend files complete |
| 14. Scalar and extension fields | Fr, Fp, Fp2, Fp6/Fp12, reductions, inversion, square roots, hash-to-field support | Algebraic properties, modulus-boundary fuzzing, independent field oracles | Field files complete across applicable scalar backends |
| 15. Curves and protocols | G1/G2, serialization, subgroup checks, map/hash to curve, pairing/GT, keygen/signatures/aggregation, scalar MSM, uniqueness tree | Malformed-input and state fuzzing, group/pairing properties, resource/bounds sanitizers | Scalar curve/protocol files complete |
| 16. AVX-512 | Primitives, RNS/reduction, fields, curves, pairing/GT, MSM, conversions, dispatch/gating | Exhaustive exposed-operation scalar differential work, edge/lane/alignment fuzzing, native AVX-512 execution | AVX files complete with native hardware evidence |
| 17. FFT and middle product | Scalar/AVX domains, roots, permutations, normalization, aliasing, scratch/sizing, APIs | Simple bounded reference implementations, identity/differential fuzzing, sanitizers | FFT/middle-product files complete |
| 18. C/C++/Rust surfaces | Public headers and symbols, native C++, C contracts, Rust FFI/safe wrappers, unsafe/ownership/features/errors, package inputs | Cross-language differential tests, standalone consumers, FFI sanitizer/fuzz coverage, Clippy feature matrices | Supported language/API files complete |
| 19. Support and release evidence | Tests, vectors, build/config, CI, policy/audit scripts, docs, examples, packaging, residual gaps | Full static/sanitizer/fuzz corpus and supported scalar/AVX/Apple/C/C++/Rust/package matrix | All remaining in-scope files and release evidence complete |
| 20. CatBLST identity and ABI | Hard rename plus build-time `blst_` or `catblst_` C prefix generation | Separate symbol, consumer, package, and supported test matrices for both prefix modes | Rename and selectable ABI validated without runtime overhead |

Phase 12 is the first platform milestone, not the end of automated work.
`AUTOMATED_ANALYSIS_PLAN.md` remains controlling through Phase 19. Phases
13-18 own both their hand review and the subsystem-specific automation needed
to close them.

## Current work volume

The generated manifest currently contains 498 in-scope production, generated,
and support files. Fifty are complete and 448 remain pending or in review. This rollup exposes
the scale of the main program; `FILE_MANIFEST.csv` remains authoritative.

| Main work group | Roadmap phases | Manifest subsystems rolled up | Total | Complete | Not complete |
| --- | --- | --- | ---: | ---: | ---: |
| Scalar foundations, fields, curves, and protocols | 13-15 | `common`, `core`, `scalar_backend`, `fields`, `curve`, `assembly_generator`, `generated_assembly` | 192 | 30 | 162 |
| AVX-512 | 16 | `avx512`, `avx512_generated` | 42 | 1 | 41 |
| FFT and middle product | 17 | `fft_middle_product`, `fft_middle_product_avx512` | 10 | 0 | 10 |
| C/C++/Rust APIs and bindings | 18 | `c_api`, `rust_binding`, `rust_binding_support`, `rust_ffi` | 28 | 2 | 26 |
| Support, tests, builds, CI, audit process, and release controls | 12 and 19 | `tests`, `test_vectors`, `build_configuration`, `ci`, `repository`, `audit_process`, `api_contracts`, `command_demo` | 226 | 17 | 209 |
| **Total** | **12-19** | | **498** | **50** | **448** |

The table is a planning rollup, not a shortcut around one-row-per-file review.
Before each phase begins, its exact path queue and dependency order must be
recorded. Files that span phases are assigned to the earliest phase that can
review their complete production contract; later phases cite that record and
add a delta review if new behavior is exercised.

## Parallel review policy

Independent leaf files or clearly separated file families may be reviewed in
parallel after their shared dependencies and automated baseline are stable.
Each review has one named owner and exact commit. Findings that alter shared
types, representations, generators, or public contracts stop dependent review
queues until the change is integrated and affected records are rebased or
reopened. Subsystem closure is always an integrated gate, never the union of
unreconciled parallel notes.

## Program completion gate

The production audit is complete only when:

- all in-scope production, generated, and support files have complete or
  explicitly accepted-risk records;
- all findings are resolved or accepted with durable rationale;
- generators and generated outputs have exact provenance;
- scalar and AVX-512 implementations have independent and differential
  evidence, including native AVX-512 hardware;
- C, native C++, and Rust APIs agree and pass their supported matrices;
- every retained C ABI entry point is proven nonthrowing across its transitive
  call graph; legacy `void` entries are allocation-free, and fallible entries
  translate errors without partially publishing output;
- Apple Silicon and all other claimed platforms have native release evidence;
- the full configured static-analysis, sanitizer, and fuzz-regression corpus is
  clean with every diagnostic dispositioned;
- benchmarks remain outside the audit except separately approved performance
  gates, with no concurrent benchmark execution;
- frozen BTx/BTE and experimental CUDA/non-Apple ARM64/CHERI remain outside the
  production claim; and
- audit integrity, packaging, security policy, and release controls pass at the
  exact release candidate commit.
