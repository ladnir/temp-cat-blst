# ARM64 and CHERI Experimental Target Policy

## State

ARM64 and Morello/CHERI sources and generated assembly are retained because
they may work and remain useful for future qualification. They are not fully
supported, production-ready, or release-qualified CatBLST targets.

The production-qualified architecture is currently x86-64, covering the
primary scalar and AVX-512 profiles. ARM64 and CHERI evidence does not inherit
from x86-64 testing, deterministic assembly generation, or source review.

## Default build policy

- CMake rejects every detected ARM64 target unless
  `BLST_ALLOW_UNTESTED_ARM64=ON` is explicitly provided.
- A CHERI pure-capability compiler additionally requires
  `BLST_ALLOW_UNTESTED_CHERI=ON`; both acknowledgements are mandatory.
- An opted-in CMake configuration emits a warning that the target may work but
  is unsupported and unqualified.
- Generated `config.h` repeats the compile-time guards so the restriction is
  enforced when sources or assembly are compiled.
- Cargo rejects ARM64 before using either a source build or a prebuilt native
  archive. Rust callers opt in with the corresponding environment variables
  set to a true value such as `1` or `ON`.
- The opt-ins acknowledge risk only. They do not grant ABI stability,
  packaging support, security support, or production-readiness claims.

The target sources must remain reproducible and receive relevant source-level
security fixes while retained. Experimental status is not permission to let
generated files drift.

## Promotion gate

ARM64 remains experimental until a dedicated phase records, at minimum:

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
