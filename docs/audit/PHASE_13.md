# Phase 13 - Foundational Scalar Trust Chain

## Scope and reviewed state

Phase 13 begins the file-by-file hand audit at the configuration and common
representation boundary. This packet completes `blst/config.h.in` and
`blst/common/defines.h` at exact reviewed commit
`7719303ad809dd7a738421c721be8ebe91bf6ce5`. It remediates dependent defects in
the C ABI and FFT boundary without claiming complete reviews of those larger
files.

The opening review baseline was
`33914f23da23fced86ca51ce9917cf52b59bd96e`. Six planned findings were recorded
as `AUD-042` through `AUD-047`. Required boundary tests exposed two additional
correctness defects, `AUD-048` and `AUD-049`; both were fixed in the same
responsible subsystem packet. The retained-C-ABI inventory also found the
separate allocation-policy gap `AUD-050`, which remains open for a dedicated
keygen review.

## Remediation results

### C ABI exception safety

`blst_miller_loop_n` no longer constructs `std::vector` at the C boundary.
Scalar-default builds statically require identical public and scalar affine
representations and call the existing compact-pointer implementation directly.
AVX-default builds convert at most 16 pairs at a time into fixed stack arrays,
run the existing contiguous AVX implementation, and accumulate each batch in
`PublicGt`. Output is published only after all batches complete.

The required infinity regression found `AUD-048`: both scalar multi-pair
implementations accumulated infinity terms for `n > 1`. They now skip those
pairs while preserving the fixed batch structure. Explicit and compact arrays
at sizes 1, 15, 16, 17, and 32, malformed heads, and infinity cases passed four
focused C API tests. The exact Rust `miller_loop_n` wrapper test also passed.
Function-level scalar and AVX wrapper disassembly contains no allocation call,
conversion heap, runtime dispatch, or threading.

A read-only inventory covered retained `extern "C"` definitions in `c_api.cpp`,
`avx_c_api.cpp`, `fft_c_api.cpp`, `middle_product_c_api.cpp`, `exports.cpp`, and
`aggregate.cpp`, plus their retained C/Rust consumers. BTx/BTE, CUDA,
benchmarks, and unsupported bindings were excluded. No second C++ exception
escape was found: FFT catches and translates invalid-domain/allocation failure
to `false`, interpolation and middle-product boundaries return explicit status,
and the remaining retained wrappers expose no source-level throwing allocation,
`.at`, string construction, or callback.

The inventory did identify `AUD-050`: legacy `void` keygen/EIP-2333 APIs still
allocate an HKDF `info` buffer and clear the output key if allocation fails.
That path does not throw and therefore is not part of `AUD-042`, but it violates
the allocation-free legacy-void policy. It is intentionally left open for a
keygen-specific change and vector packet.

### Common helpers and FFT boundaries

`log2ceil(0)` now asserts in Debug and returns zero in Release. Empty native
root preparation throws its existing invalid-domain exception before sizing;
the public C FFT boundary continues translating it to `false`. `ceil_div` now
uses quotient/remainder arithmetic, asserts for a zero divisor in Debug, and
returns the documented fail-closed zero result in Release. Supported callers
were inventoried and either provide nonzero divisors or explicitly handle
failure.

The size-one root regression exposed `AUD-049`: generic and Fr-specialized
primitive-root checks skipped their first unique prime factor and underflowed
the loop bound for a factor list of size one. Both checks now traverse the full
factor list with an ordinary zero-based bound. Focused Release and Debug helper
and root tests passed; the Debug build exercises the two assertion contracts,
and the C prepared-domain API rejects an empty domain.

### Constant-time value-barrier evidence

`test/audit_launder_codegen.cpp` provides noinline probes for conditional
select, swap, conditional zero, secure zero, and a representative portable
arithmetic mask. Optimized output was inspected for secret-dependent branches,
addresses, calls, and variable spills. Secret operations were branchless and
used fixed addresses on all three compilers. The only branch is the public
length loop in secure zero; GCC's portable-mask alias checks depend only on
public pointer relationships.

| Toolchain and options | Evidence SHA-256 |
| --- | --- |
| MSVC 19.50.35728, `/O2 /Ob3 /DNDEBUG`, object | `71C049CF92E85E394F665A6C1E573791861799D1EA6ED5BE48BB3369D3003292` |
| MSVC 19.50.35728, `/O2 /Ob3 /DNDEBUG`, assembly | `1BB80AB1865DDC7F7E7633A1FCE50F1C5FEC8B9A892C214D78C615738A22AB19` |
| Clang 18.1.3, optimized assembly | `A998A1B52BE63722C5C40E747CE84F0FE884F23E093E5C39458CD7995B578D10` |
| GCC 13.3, optimized assembly | `BB50BACD6F29AD1311F2F6F435329BB0A3C0C1ADFFF618F7892B0FC72127B5AE` |

Artifacts remain local under ignored `out/audit/barrier/`. `BLST_launder` was
not changed, so no benchmark was required.

### Warning, provenance, and 64-bit configuration

`defines.h` now carries the repository's Apache-2.0/SPDX header and explicitly
states that it is locally composed with lineage from upstream BLST `src/vect.h`.
Translation-unit-wide Intel, GCC `-Wpedantic`, and MSVC 4127/4189 suppressions
were removed. The exposed diagnostics were fixed at their source: unsigned
mask construction, a GNU statement-expression most-significant-bit helper, an
Fp12 narrowing conversion, compiler-resistant byte clearing, FFT index-range
validation, and an unused FFT parameter. No replacement suppression was added.

The most-significant-bit helper remains forced-inline and the optimized strict
objects contain no helper call. FFT retains `int` loop indices after one cold
entry range check, so no hot-loop widening was introduced.

`BLST_LIMB_SIZE` is now a CMake cache string whose only supported value is
`64`; every other value fails configuration. `defines.h` independently rejects
non-64-bit generated configuration and defines `limb_t` unconditionally as
`uint64_t`. Two dormant 32-bit branches remain in `backend_portable.h` at the
reviewed commit and are deferred to that file's Phase 13 review rather than
mechanically deleting unaudited backend code.

## C ABI exception policy

No C++ exception may cross a CatBLST C ABI boundary. Legacy `void` entry points
must be allocation-free and nonthrowing for inputs satisfying their documented
preconditions. Fallible entry points must translate failure to an explicit
status, `bool`, or null result and must not partially update caller output.
Adding `noexcept` is not a substitute for proving this contract because an
escaping exception would terminate the process.

## Focused validation

All compiler work was run sequentially with one build job. No full suite,
sanitizer sweep, portable-AVX Release build, benchmark, CUDA build, or BTx/BTE
target was run.

| Area | Focused result |
| --- | --- |
| C ABI Miller | Scalar `CApi.MillerLoopN*` 4/4; exact Rust wrapper 1/1; native-AVX Debug C API object compiled; scalar/native wrappers have no allocation calls |
| Helpers and FFT | Release and portable-AVX Debug focused helper/root tests 2/2 each; empty prepared-domain C test 1/1 |
| Warning policy | Standalone header and exact backend/C API/FFT C API objects passed MSVC `/W4 /WX`, Clang and GCC `-Wall -Wextra -Wpedantic -Werror` |
| Limb configuration | Explicit `64` configuration succeeded; explicit `32` configuration failed with the CatBLST 64-bit-only diagnostic |
| Constant time | MSVC, Clang, and GCC optimized probes inspected and hashed as above |
| Audit process | `git diff --check` and `scripts/audit/check_audit_manifest.ps1` passed at closure |
| Public zero shifts | Fr/Fp/Fp2 C ABI identity passed out-of-place and in-place; affected native C++/C test objects compile; constant-count Release wrappers add no branch |
| Key derivation | Four published EIP-2333 master/child vectors passed, including the 255-block Lamport expansion; HKDF test mode compiles; no allocation reference remains |
| Portable arithmetic | Exact lost-borrow regression and zero/one/top-bit scan passed; portable inverse matched assembly for zero, one, and modulus-minus-one |
| ADX typing | 10,072 portable-versus-ADX cases passed under MSVC, GCC, and Clang; MSVC Release emitted direct arithmetic with no conversion temporary |

## Closure and residual work

`AUD-042` through `AUD-049` are validated. `blst/config.h.in` and
`blst/common/defines.h` have complete per-file review records at the same exact
commit. The larger C API, scalar GT, FFT, and primitive-root files remain
pending for their complete Phase 15, 17, or 18 reviews even though their linked
findings are resolved.

The next Phase 13 dependency is `blst/backend/backend_portable.h` followed by
`backend_portable.cpp`. That review should remove or deliberately retain its
two dormant 32-bit branches, map each constant-time primitive to the barrier
probe evidence, and preserve current inlining and fixed-width arithmetic.
`AUD-050` should be handled as a separate keygen packet before the public C ABI
phase can close.

The subsequent portable-backend review upgraded `AUD-050` to high severity:
Lamport expansion failure is not merely silent zeroization because child
derivation continues into key generation and can publish an incorrect key.
That review also opened `AUD-051` through `AUD-054` for zero-count public
shifts, ADX intrinsic output aliasing, lost-borrow portable division, and the
MSVC zero-input bit scan. These are remediation packets, not completed file
reviews; exact validation evidence will be recorded before closure.

The remediation packets are complete. `AUD-050` through `AUD-054` are
validated. Public zero shifts are identity before backend selection; keygen and
EIP-2333 are allocation-free; portable division propagates borrow without limb
overflow; the MSVC zero bit scan is defined and branchless; and ADX intrinsic
outputs use an exact matching private word type without changing `limb_t` or a
public representation. The remaining 32-bit backend branches, broad Clang
warning suppression, and unused secret printers were removed. Complete review
records now cover both portable-backend files at exact commit
`5606e84262401225c4e02b3ec786b849c85adfcc`. Field and keygen files remain
pending for their broader Phase 14/15/18 reviews despite closure of these
specific findings.

The subsequent `blst/backend/backend.h` review reopened `AUD-051`: its first
remediation covered the public field and C boundaries but not every supported
native C++ template operator. The same review opened `AUD-055` because the
selected SHA-256 compression wrapper has unconditionally used the portable
implementation since CUDA refactor `4e8775d1`, even in default-assembly host
builds. Both findings require focused remediation before the selector file can
be closed.

The selector remediation is complete. Zero-count identity now lives once at
all four `_v256`/`_v384` shift boundaries, and duplicate field checks were
removed. Native/default-assembly focused tests passed 6/6, including the
deterministic SHA comparison; no-assembly shift tests passed 5/5. MSVC
`/O2 /Ob3` probes show direct tail calls for all constant nonzero shifts and
the intended native versus portable SHA symbols. `AUD-051` and `AUD-055` are
validated, and `blst/backend/backend.h` is reviewed at exact commit
`dc07be802c39a3c04a17347a5cec64f74d95ecd7`.

The subsequent assembly-interface review reopened `backend.h` and opened
`AUD-056` and `AUD-057`. Zero-length vector predicates and selected SHA
operations can reach raw assembly routines whose count contract is strictly
positive; the generated x86-64 and AArch64 implementations load before proving
the count and then underflow their loops. Current production callers use
nonzero fixed sizes or guard the count, but the supported native C++ boundary
is unsafe. The same review found that persistent provider-renaming macros make
the apparent q/x symbol surface include-order-dependent, conceal which Apple
provider is tested, and prevent combined builds from directly validating all
built providers. Remediation must remain forced-inline and compile-time, with
code-generation proof that constant nonzero hot calls gain no branch.
