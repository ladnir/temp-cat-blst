# Phase 14 - Scalar and Extension Fields

## Scope

Phase 14 reviews the scalar `Fr` and base `Fp` implementations before moving
upward into extension fields. The first wave reviewed every line of
`blst/scalar/fields/Fr.h`, `Fr.cpp`, `Fp.h`, and `Fp.cpp` at
`c3096482cdab88edc4f05e0c6583c7896fd75bab`, together with only the directly
corresponding C and Rust declarations, wrappers, and tests. CUDA code was
excluded. No build, benchmark, or source edit occurred during the review.

## First-wave findings

The review opened `AUD-067` through `AUD-074`. The approved straightforward
packet covers zero exponentiation, repeated-square boundary safety, public
square-root failure output, secret-boundary cleanup, restoration of retained
scalar declarations, bulk-size hardening, and file documentation/provenance.

`AUD-069` and `AUD-070` were deliberately excluded from that implementation
packet. The then-current C boundary relied on ABI-compatible but nonidentical C
and C++ function/object types, while several inherited C storage idioms did not
establish valid C++ object lifetime. Their repair required an explicit
zero-copy representation design and optimized scalar/AVX code-generation
evidence; routine boundary copies are not acceptable.

## Validation policy

The straightforward packet uses only focused Fr/Fp/C-boundary tests,
standalone C/Rust declaration checks, assembly regeneration when a generator
changes, and small optimized code-generation probes. Full suites, sanitizer
sweeps, portable-AVX Release, benchmarks, CUDA, and CI pushes remain separate
permission-gated evidence.

## Straightforward remediation result

Commit `7c0e399` changes the x86 q/x repeated-square generators and all eight
affected ELF, COFF, Mach-O, and Win64 products to preserve full `size_t`
counters. Commit `b6a51bb` implements the remaining approved straightforward
field packet. `AUD-067`, `AUD-068`, `AUD-072`, `AUD-073`, and `AUD-074` are
validated. `AUD-071` is fixed with focused C evidence and remains pending a
retained Rust execution check. `AUD-069` is fixed pending focused runtime,
Rust, and Apple evidence; `AUD-070` is validated at `61c8cd4`.

Validation remained focused. One existing Debug target was rebuilt with one
Ninja worker after the field headers invalidated dependent objects; only three
public C/header tests were executed and passed. A standalone scalar field probe
passed zero-exponent, zero-count, exact-alias, scalar conversion, square-root
failure, and oversized-count cases. The scalar conversion header compiled as
strict C11 under MSVC `/W4 /WX`; GCC 13 and Clang 18 strict syntax checks passed
for the probes and modified C++ boundaries. Assembly regeneration verified all
100 products. MSVC 19.50 Release code generation shows constant-count
`sqr_n_mul` remains one direct native call without the alias slow path, while
constant-zero exponent calls reduce to the required identity copy. No full
test sweep, sanitizer run, benchmark, CUDA job, portable-AVX Release build, CI
push, or Apple run was performed.

## AUD-069 representation remediation

The approved representation is a single conditional public definition rather
than a copied boundary object. C sees the existing data-only `blst_fr`,
`blst_fp`, extension-field, and point structures. C++ sees standard-layout,
trivially-copyable classes with the same public storage names and the existing
member functions; `blst::scalar::Fr`, `Fp`, `Fp2`, `Fp6`, and `Fp12` are exact
aliases of the corresponding `blst_*` classes. The generic scalar point types
are likewise the exact specializations named by `blst_p1`, `blst_p2`, and their
affine forms. The public C layout is unchanged.

The same public-header cleanup establishes C99 as the minimum supported C
language mode. Standard headers now own `size_t`, the fixed-width integer
types, and `bool`; CatBLST no longer redeclares those names or maps legacy C
`bool` to ABI-incompatible `int`. Library-owned global C typedefs use the
`blst_` namespace (`blst_byte` and `blst_limb_t`), while internal C++ retains
only the namespaced conveniences `blst::byte` and `blst::limb_t`. Public
binding-header guards use the conventional non-reserved `BLST_*_H` spelling,
and Rust bindgen is invoked explicitly in C99 mode.

## AUD-070 lifetime remediation

CatBLST now rejects non-little-endian targets at compile time with
`std::endian`; this matches the supported x86-64 and Apple Silicon platform
matrix. Canonical `blst_scalar` remains a 32-byte little-endian value distinct
from Montgomery-form `blst_fr`, and their fixed-width conversions document an
exact-or-disjoint overlap contract with no partial overlap.

The inherited `radix384` union is replaced by one real typed `Fp12` identity.
Fp and Fp2 identities and zero values are references to its live nested
subobjects. Reciprocal retains one 96-byte `UFp` scratch object and addresses
its two halves only as limbs, preserving the native inverse, reduction, and
Montgomery-multiply call structure without changing stack size. Focused MSVC,
GCC, and Clang compilation, the scalar library build, exact-alias and identity
boundary checks, reciprocal correctness, and optimized MSVC code inspection
pass at `61c8cd4`.

The follow-up review found that `blst_t.hpp` was an isolated generic-modulus
shim rather than the supported fixed BLS12-381 native C++ core. It duplicated
low-level arithmetic and retained separate lifetime violations, while only
three dedicated tests consumed it. Commit `a720184` removes that header and
those tests; canonical native C++ field, curve, and protocol types remain.

## Fp2/Fp12 hand-audit wave

The complete scalar Fp2 and Fp12 source review at `f0a9d9e` opened
`AUD-079` through `AUD-087`. The inherited tower arithmetic and selected
widened multiplication structure remain suitable, but the supported surface
is not ready to close. High-severity defects cover zero/alias exponent
boundaries, unchecked public Frobenius degrees, conflation of the cyclotomic
and order-r GT domains, undeclared variable-time GT scalar processing, and
untrusted precomputed-schedule lengths. Medium findings cover Fp2 C square-root
failure publication, noncanonical unchecked GT decoding, and an advertised
but undefined native GLV operation. The cleanup finding covers reserved/dead
implementation artifacts, warnings, provenance, and missing negative tests.

No behavior was changed in this review packet. The GT timing policy requires
an explicit design decision before remediation: safe constant-time default
with a named public-scalar fast path, or an explicitly variable-time-only API
plus a separate constant-time operation. Focused compiler probes only were
run; no product build, broad test, sanitizer, benchmark, CUDA, CI, or Apple job
was started.

That strict all-header check also opened and fixed `AUD-075`: MSVC C rejected
the AVX opaque-shell alignment when its `__declspec(align(...))` argument was
spelled as an enum identifier. The fixed expert-ABI sizes and 64-byte alignment
are now prefixed preprocessing constants, leaving their values and layouts
unchanged.

The field and point API documentation moved with those declarations. The
canonical public header retains the Montgomery representation, constant-time,
failure-output, aliasing, subgroup, encoding, GLV, and scalar contracts;
implementation-only commentary remains in the implementation headers. The
real packed AVX API remains separate in `blst::avx` and `blst_avx.h`; four
scalar `Fp12::*_avx` placeholders that asserted and fell back to scalar code
were removed during final review rather than being presented as working AVX.

`exports.cpp` and the retained C boundaries now include their public headers
and define the exact pointer, Boolean, fixed-width, and opaque-object function
types. Redundant Fr/Fp/Fp12/G1/G2 reinterpretation was removed from the scalar,
FFT, middle-product, and AVX conversion boundaries. Intentional opaque pairing,
uniqueness-tree, AVX shell, and caller scratch-buffer casts remain separate
review subjects; none aliases a duplicate public field or point type.

Focused evidence currently marks `AUD-069` fixed, not yet validated. MSVC parses
the complete retained scalar library source list, native and portable AVX C
boundaries, and the affected scalar/public tests. GCC 13 and Clang 18 compile
the standalone header as C11 and C++20, and GCC parses the scalar exports, C
API, FFT, and middle-product boundaries. Compile-time checks require exact type
identity, standard layout, trivial copying, public sizes, alignments, and first
member offsets. MSVC optimized Fr/Fp/Fp2 arithmetic wrappers contain direct
backend calls or tail calls and no `memcpy`, `memmove`, conversion loop, or
temporary lane object. A focused linked/runtime check, retained Rust compile,
and Apple execution remain before validation. No full build, sanitizer sweep,
benchmark, CUDA build, or CI run was started.

Final review also completed the `blst_byte`/`blst_limb_t` migration in retained
BTx/TRX boundary sources and public-header tests. BTx/TRX behavior was not
changed or audited. `C_api_tests.cpp` now includes the public umbrella header;
its retained declarations are checked against that header by the compiler, and
the keygen, EIP-2333, pairing, uniqueness-tree, and fixed-width call sites use
the exact public types.
