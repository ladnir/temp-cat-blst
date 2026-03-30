# VROOM Backend Reset Plan

## Purpose

This document defines a phased plan to:

1. remove backend code that no longer has a clear role,
2. reduce the packed AVX path to the pieces we still believe in,
3. delete the current VROOM integration work,
4. reintroduce VROOM as a literal in-tree copy,
5. harden that copied code until it is fit to live inside `blst`.

`Vendored` means the source code is copied into this repository and built here, rather than referenced as an external checkout.

## Why We Are Doing This

We currently have three overlapping SIMD or AVX directions:

- single-operation `avxbls` integration
- packed-limb batched AVX integration
- VROOM CRNS integration

That is too much overlap, and the performance data does not support keeping all three alive.

The current conclusions are:

- packed AVX limb arithmetic is very strong at base field level
- that advantage does not survive the climb to packed `Fp12`
- VROOM is materially better for `Fp12` and higher extension-field work
- the current in-tree VROOM work should not be treated as a permanent foundation because it mixes scaffolding, experiments, and integration glue

So this reset exists to replace a mixed experimental state with a deliberate backend split and a clean rebuild path.

## Current Direction

The current evidence supports this backend split:

- packed AVX limb backend for `Fp`, `Fp2`, likely `G1`, and maybe `G2`
- VROOM CRNS backend for `Fp12`, `GT`, and pairing-extension work

The existing single-operation `avxbls` path does not currently justify continued investment unless a specific benchmark later proves otherwise.

## Goals

- end with a clean backend story instead of three overlapping experiments
- preserve the strong packed-AVX base field work
- rebuild VROOM in-tree from a literal copy, not a hand port
- add robust tests as each layer is introduced
- make the resulting code maintainable enough for long-term work inside `blst`

## Final Intended State

At the end of this plan, the intended state is:

- one packed AVX backend that is intentionally limited in scope
  - owns `Fp`
  - owns `Fp2`
  - may later own `G1`
  - may later own `G2` if benchmarks justify it
- one VROOM-derived backend that is intentionally focused on higher-tower work
  - owns `Fp12`
  - owns `GT`
  - owns pairing-extension arithmetic
- no active single-operation `avxbls` tower backend
- no active hand-ported VROOM arithmetic core

## Definitions

- `parity`: same algorithmic path, same benchmark shape, and performance that matches upstream within normal benchmark noise
- `benchmark noise`: small run-to-run variation that does not indicate a real implementation difference
- `gate`: a required condition that must be met before the next phase can begin
- `literal copy`: source copied with only the minimal edits required for include paths, build integration, and symbol collision avoidance
- `robust tests`: deterministic tests, edge-case tests, randomized parity tests, and regression tests sufficient to make later refactors safe

## Execution Rules

- no phase may begin until the previous phase gate is satisfied
- no “cleanup” edits inside copied VROOM hot code before parity is established
- if a copied VROOM file must be changed for build reasons, the change must be minimal and called out explicitly
- every imported arithmetic layer must be benchmarked before the next layer is added
- every imported arithmetic layer must gain tests before the next layer is added

## Parity Target

For VROOM, `parity` is not informal. It means all of the following:

- same upstream VROOM revision is identified and recorded
- same machine class is used for the comparison
- same compiler family and version is used for the comparison
- same benchmark binary shape is used, or a direct benchmark with equivalent timed work
- same preloaded or boundary condition is used
- measured runtime is within a small tolerance band across repeated runs

### Required Parity Comparison Set

- upstream VROOM direct base-field benchmark
- in-tree copied VROOM direct base-field benchmark
- upstream VROOM direct `Fp12 mul` benchmark
- in-tree copied VROOM direct `Fp12 mul` benchmark

### Required Parity Reporting

Every parity checkpoint must record:

- machine type
- compiler version
- benchmark command
- benchmark output
- whether the result is considered pass or fail

### Default Pass Rule

Unless replaced by a stricter rule later, parity should mean:

- median runtime within `5%` of upstream over repeated pinned-core runs, and
- no persistent directional drift that suggests code-generation differences

If the median is outside `5%`, parity is not achieved and the phase is blocked.

## Test Bar

`Robust tests` in this document means:

- deterministic known-answer tests
- deterministic edge-case tests
- randomized parity tests against a trusted reference
- regression tests for any bug fixed during the migration
- explicit coverage of both supported parameter sets where applicable
- tests that fail with a nonzero exit code and are easy to run from one documented command

`Robust tests` does not mean:

- only a few random smoke tests
- benchmark-only validation
- tests that depend on manual inspection of printed output

## Phase Gate Template

Each phase must explicitly answer these questions before the next phase starts:

1. What code was added or deleted?
2. What tests prove the current layer is correct?
3. What benchmark proves the current layer matches expectations?
4. What known gaps remain?
5. Why is it safe to move to the next phase?

## Non-Goals

- no API integration into the main library yet
- no refactor for elegance before parity
- no algorithm changes to VROOM during the literal-port phase
- no attempt to preserve obsolete intermediate experiments

## Guiding Rules

- copy VROOM literally first, optimize or clean later
- delete dead or unclear code aggressively
- each phase must end in a passing test and benchmark state
- every added layer must have parity tests against a trusted reference
- new tests must be deterministic where practical

## Phase 0: Freeze The Target

### Deliverables

- record the exact upstream VROOM revision that will be copied
- record the exact benchmark configuration used for parity
- record the intended backend ownership:
  - packed AVX owns `Fp` and `Fp2`
  - VROOM owns `Fp12` and above

### Exit Criteria

- one chosen upstream commit hash
- one chosen compiler/toolchain baseline
- one benchmark command set for parity
- one written record of the intended backend ownership split

## Phase 1: Remove Obsolete Backend Code

### Scope

Delete code that no longer matches the target architecture.

### Deletion Decision Rules

Delete a file, target, or test if any of the following is true:

- it exists only to support the single-operation `avxbls` tower path
- it exists only to support packed AVX above `Fp2`
- it exists only to support the hand-ported CRNS or VROOM path
- it duplicates functionality that will be reintroduced by the literal VROOM copy
- it is benchmark-only scaffolding for a backend we no longer intend to keep
- it is dead code, abandoned experiment code, or disconnected build glue

Keep a file, target, or test if all of the following are true:

- it supports packed AVX `Fp` or packed AVX `Fp2`
- it is needed by benchmarks we still care about
- it is not tied to a discarded tower design
- it has a clear role in the final intended state

Replace rather than delete when:

- the file has the right responsibility, but the current implementation is the wrong backend
- the benchmark harness is still useful, but one of the benchmarked paths is obsolete
- the test concept is still useful, but its current reference or target implementation is being removed

### Survival Test

For any file or target under review, ask these questions:

1. Is this part of the final intended backend split?
2. If this disappeared today, would we still want to recreate it later?
3. Does it serve `Fp` or `Fp2` packed AVX, or high-tower VROOM?
4. Is it a thin adapter, a test, or a benchmark we still need?
5. Is it carrying only migration history rather than future value?

If the answers do not clearly justify the file, delete it.

### Work Items

- inventory all current AVX backend files and classify them as:
  - keep
  - replace
  - delete
- remove dead CMake targets and stale benchmark targets
- remove tests that only validate code we are deleting
- keep only the benchmark paths needed for future side-by-side comparison

### Exit Criteria

- the repo builds without the deleted backends
- packed AVX `Fp` and `Fp2` still compile and test cleanly
- there is no active hand-ported VROOM code left in-tree
- the remaining backend inventory is simple enough to describe in one page

## Phase 2: Create The New In-Tree VROOM Landing Zone

### Scope

Add a clean destination for the literal VROOM copy before importing code.

### Proposed Layout

- `blst/vroom/`
- `blst/vroom/cpu/vector/`
- `blst/vroom/cpu/reduction/`
- `blst/vroom/cpu/precompute/`
- `blst/vroom/src/`
- `blst/vroom/tests/`
- `blst/vroom/bench/`
- `blst/vroom/adapters/`

### Rules

- preserve upstream folder structure as much as possible
- keep all `blst`-specific adaptation code outside the copied core
- do not mix copied VROOM code with `blst` field types in the same files

### Exit Criteria

- the directory structure exists
- build wiring exists but may still be stubbed
- copied code and adapters have a clear boundary
- there is a written rule for which files are copied core versus `blst`-owned adapters

## Phase 3: Literal Import Of VROOM Core

### Scope

Delete the old VROOM integration and copy the upstream VROOM implementation directly into the new in-tree location.

### Rules

- copy-paste, do not reinterpret
- preserve filenames where practical
- preserve call graph and data layout
- preserve accumulation order and reduction order
- keep comments from upstream

### Allowed Changes

- include-path adjustments
- namespace wrapping if needed
- build-system glue
- symbol renaming only when required to avoid collisions

### Forbidden Changes In This Phase

- no “cleanup” refactors
- no formula rewrites
- no representation changes
- no style-driven rewrites

### Exit Criteria

- the copied VROOM code compiles in-tree
- the copied base-field direct benchmark matches upstream within benchmark noise on the same machine/config
- no handwritten arithmetic replacement remains in the imported base-field path

## Phase 4: Harden The Base Field Layer

### Scope

Before climbing the tower, make `Fp`-level VROOM code production-grade enough to trust.

### Required Tests

- deterministic known-answer tests for conversion to and from RNS/Montgomery
- deterministic known-answer tests for `modmul`
- batch-size coverage for supported batch sizes
- randomized parity tests against trusted scalar reference
- fallback-path correctness tests if fallback build remains supported

### Required Benchmarking

- direct VROOM-vs-upstream parity benchmark
- benchmark reproducibility script or documented command set

### Cleanup Targets

- remove or quarantine known-broken code paths
- document invariants for:
  - representation form
  - valid bounds
  - expected reduction state

### Exit Criteria

- base `Fp` tests are deterministic and broad
- parity to upstream is confirmed
- no unexplained benchmark drift remains at base field
- all known-broken or unused base-field paths are either removed or clearly quarantined

## Phase 5: Add `Fp2`

### Scope

Introduce `Fp2` on top of the literal VROOM `Fp` layer.

### Required Tests

- deterministic algebraic identity tests
- parity tests against trusted reference implementation
- edge-value tests:
  - zero
  - one
  - negative-one
  - maximal reduced inputs

### Required Benchmarking

- direct `Fp2 mul` and `Fp2 sqr`
- compare in-tree results to upstream VROOM where applicable

### Exit Criteria

- `Fp2` is correct and benchmarked
- no design drift from upstream has been introduced
- `Fp2` test coverage is good enough that later cleanup can be done safely

## Phase 6: Add `Fp12 Mul` Literally

### Scope

Reintroduce only the `Fp12 mul` stack, still as a literal VROOM import plus thin adapters.

### Required Tests

- deterministic fixed-vector `Fp12 mul` tests
- randomized parity tests against trusted scalar reference
- parity tests against upstream VROOM
- round-trip boundary tests through the adapter layer

### Required Benchmarking

- direct VROOM-in-tree `Fp12 mul`
- upstream VROOM `Fp12 mul`
- scalar baseline
- current packed AVX baseline if still relevant

### Success Condition

The literal in-tree VROOM `Fp12 mul` must match upstream VROOM performance on the same machine, compiler, and benchmark shape.

### Exit Criteria

- `Fp12 mul` parity established
- benchmark parity established
- adapter layer remains thin and isolated
- there is a one-command benchmark path for:
  - upstream VROOM
  - in-tree copied VROOM
  - scalar baseline
  - packed AVX baseline if still relevant

## Phase 7: Expand Only What We Need

### Scope

After `Fp12 mul` is stable, add the next required operations one by one.

### Candidate Order

- `Fp12 sqr`
- cyclotomic square
- sparse line multiplication
- Miller-step helpers
- `GT` operations needed by application paths

### Rule

Do not import additional VROOM functionality until the previous layer has:

- deterministic tests
- parity tests
- benchmark parity
- documented invariants

## Phase 8: Raise Codebase Quality To `blst` Standard

### Scope

Once the literal import is stable, improve maintainability without changing semantics.

### Quality Targets

- clearer file-level documentation
- stronger test organization
- predictable build targets
- removal of dead experimental branches
- explicit naming of hot-path assumptions

### Recommended Improvements

- add dedicated test entrypoints for:
  - base field
  - extension field
  - adapters
- separate benchmark-only utilities from arithmetic code
- isolate BLST-dependent conversion helpers from pure VROOM arithmetic

### Exit Criteria

- a new engineer can identify:
  - copied upstream code
  - `blst`-owned adapter code
  - tests validating each layer

## Phase 9: Decide External Versus Permanent Vendor

### Decision Point

After the in-tree literal port is stable, decide whether to:

- keep VROOM vendored permanently, or
- treat the vendored copy as a temporary staging area before switching back to an external dependency

### Decision Criteria

- how often we need to patch internals
- how much local cleanup was required
- whether upstream remains stable and useful
- whether our tests now materially exceed upstream’s test coverage

## Immediate Next Actions

1. produce a deletion inventory for the existing AVX and CRNS code
2. delete the current VROOM integration and hand-port remnants
3. create the new `blst/vroom/` landing zone
4. copy the chosen VROOM revision literally into that landing zone
5. wire only the base-field direct benchmark first
6. add deterministic base-field tests before any tower code
7. climb to `Fp2`
8. climb to `Fp12 mul`

## Acceptance Bar

The reset is successful only if all of the following are true:

- the backend story is simple
- dead code is gone
- VROOM is copied literally, not reinterpreted
- every imported layer has strong tests
- in-tree VROOM performance matches upstream for the same benchmark
- packed AVX remains focused on the parts where it still has a clear advantage

## What This Document Should Enable

A new engineer should be able to pick up this document and answer:

- what we are deleting
- what we are keeping
- why we are doing the reset
- what “literal VROOM copy” means
- what parity means
- what test quality is required
- what blocks the next phase
- what the end state is supposed to look like
