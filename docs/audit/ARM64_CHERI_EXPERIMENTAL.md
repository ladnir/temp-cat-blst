# Non-Apple ARM64 and CHERI Experimental Target Policy

## State

Apple Silicon macOS is a supported scalar target and does not require an
experimental override. Other ARM64 and Morello/CHERI sources and generated
assembly are retained because they may work and remain useful for future
qualification. Those targets are not fully supported, production-ready, or
release-qualified CatBLST targets.

The x86-64 production profiles cover scalar and AVX-512. Apple Silicon covers
the supported scalar profile and requires its own release evidence. Other
ARM64 and CHERI evidence does not inherit from x86-64 or Apple testing,
deterministic assembly generation, or source review.

## Default build policy

- CMake accepts Apple Silicon macOS without an experimental override.
- CMake rejects every other detected ARM64 target unless
  `BLST_ALLOW_UNTESTED_ARM64=ON` is explicitly provided.
- A CHERI pure-capability compiler additionally requires
  `BLST_ALLOW_UNTESTED_CHERI=ON`; both acknowledgements are mandatory.
- An opted-in CMake configuration emits a warning that the target may work but
  is unsupported and unqualified.
- Generated `config.h` repeats the compile-time guards so the restriction is
  enforced when sources or assembly are compiled.
- Cargo accepts Apple Silicon macOS normally. It rejects other ARM64 targets
  before using either a source build or a prebuilt native archive. Rust callers
  opt in with the corresponding environment variables set to a true value such
  as `1` or `ON`.
- The opt-ins acknowledge risk only. They do not grant ABI stability,
  packaging support, security support, or production-readiness claims.

The target sources must remain reproducible and receive relevant source-level
security fixes while retained. Experimental status is not permission to let
generated files drift.

The audit manifest therefore keeps the 9 Apple ARM64 Mach-O artifacts in the
production set and classifies the 36 retained ELF, COFF, Win64, and CHERI ARM
artifacts as experimental. They continue to regenerate from the same owners;
this classification does not delete them or imply a separate implementation.

Commit `61d3b66` also corrected the opted-in MSVC ARM64 source selection to use
the validated target architecture rather than pointer width. The correction is
structurally reviewed, but native Windows ARM64 configuration and execution
remain promotion evidence because the local MSVC installation has no ARM64
toolset.

## Promotion gate

Non-Apple ARM64 remains experimental until a dedicated phase records, at
minimum:

1. native AArch64 assembly and portable builds with documented compilers and
   operating systems;
2. focused inverse tests plus the complete scalar C/native C++ and Rust suites;
3. sanitizer coverage where supported;
4. package and public-header consumers; and
5. repeatable CI on maintained AArch64 infrastructure.

CHERI additionally requires:

1. assembly acceptance by a current Morello-capable LLVM toolchain;
2. pure-capability execution under CheriBSD on QEMU-CHERI, the Morello FVP, or
   Morello hardware;
3. capability-fault-free 256- and 384-bit inverse boundary tests with
   multiplicative-identity checks; and
4. an explicit review of pointer provenance, bounds, ABI, and dependency
   compatibility across the retained C and native C++ surfaces.

Removing either default block requires a separately reviewed policy change
linked to that evidence. Merely obtaining a successful compilation is not
sufficient.

Apple Silicon remains in supported scope, but release readiness still requires
a real macOS ARM64 CMake build, the scalar C/native C++ and Rust suites, public
header/package consumers, and repeatable Apple CI. Its exemption is a support
decision, not fabricated test evidence.

Current native evidence is GitHub Actions run `30754724252` at commit
`0b7ec47`: the `macos-15` ARM64 runner used AppleClang 17 with
`BLST_ALLOW_UNTESTED_ARM64=OFF`, passed all 336 C/C++ tests, and passed the Rust
default and full-binding test matrix. Public header and package-consumer gates
remain part of the Phase 19 release qualification rather than a condition for
the existing Apple Silicon support designation.
