# Phase 6B Go Support Removal Evidence

## Record

- Date: 2026-07-31
- Implementation commit: `aeea65edf5d3777b853baf52925ffd3caea13186`
- State: validated and committed
- Scope: Go language support, Go-only upstream rows, and `AUD-019`
- Benchmarks: excluded and not run
- CUDA: disabled and not exercised

## Decision

Go is not a supported CatBLST language surface. The project has no external
users or compatibility requirement for the binding, so Phase 6B removes it
rather than auditing, repairing, or publishing it. Supported language surfaces
remain C, native C++, and Rust.

The Go-only upstream rows remain `not_applicable`. Cross-language defect
signals were not discarded: compressed-infinity, empty-salt, uniqueness-tree,
Pippenger-planner, and scalar-length issues were already translated to and
validated on retained C, native C++, or Rust paths in earlier phases.

## Removed surface

The implementation commit deletes:

- all 14 files under `bindings/go/`, including generated and template Go,
  tests, the generator, CGo assembly/server shims, and documentation;
- root `go.mod` and `.golangci.yml`;
- the Go-only lint workflow;
- the unused Go source-root helper `build/srcroot.go`; and
- the orphaned `blst/mains/server_.cpp` amalgamation.

It also removes Go classification rules from the audit manifest generator and
updates public project guidance to list only C, native C++, and Rust as
supported surfaces.

## Orphaned amalgamation finding

The inventory found that `blst/mains/server_.cpp` retained the repository's
only `__BLST_CGO__` conditional. No build, package, or source file referenced
the amalgamation. A direct compile also showed it was already stale after the
fork's directory refactor: it includes the nonexistent
`blst/curve/GT.cpp`. The deleted CGo shim independently included a nonexistent
`server.c`, so the local Go path was not a viable build surface before this
phase.

Deleting the orphan is safer than retaining a nominal production file that is
outside every supported target and cannot compile. This is recorded as
`AUD-019` rather than silently disappearing from the production inventory.

## Validation

| Surface | Configuration | Result |
| --- | --- | --- |
| C and native C++ | Windows MSVC Debug, direct ADX scalar | 282 registered: 273 passed, nine expected AVX-only skips, zero failed |
| AVX-512 | Windows MSVC Debug, portable AVX-512 default | `blst_tests` target remained buildable |
| Rust | `cargo test --features full-rust-bindings` | 71 unit tests and one compressed-infinity integration test passed |
| Rust packaging | `cargo package --list --allow-dirty` | No Go, CGo, module, or source-root helper path is packaged |
| Repository | Filename and content reference searches | No retained non-audit Go/CGo support reference or Go source file remains |
| Formatting/integrity | Cargo formatting and Git whitespace | Passed |

## Residual validation

- The removed Go suite was not run because Go is outside the supported surface
  and the local CGo server path was already incomplete.
- No ARM64 or WebAssembly build was run for this removal-only batch.
- Unsupported SWIG, C#, Emscripten, Java, Node.js, and Python surfaces remain
  scheduled for separate reviewed removal.
