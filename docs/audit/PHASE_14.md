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

The approved representation uses one core C++ object type and no copied
boundary object. The native field and point headers own the standard-layout,
trivially-copyable `Fr`, `Fp`, `Fp2`, `Fp6`, `Fp12`, `Point`, and
`AffinePoint` definitions. In C++ mode `bindings/blst_types.h` includes those
owners and makes every `blst_*` name an exact alias. In C mode it exposes the
corresponding data-only structures. The public C layout is unchanged, and the
core no longer depends on a binding header for its mathematical types.

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

The inherited `radix384` union was first replaced by typed storage. The later
core-ownership cleanup gives Fp, Fp2, and Fp12 independent correctly typed
zero/one constants, removing the dependency that made the base-field header
require a complete Fp12 type. Reciprocal retains one 96-byte `UFp` scratch
object and addresses its two halves only as limbs, preserving the native
inverse, reduction, and Montgomery-multiply call structure without changing
stack size. Focused compilation and value-based identity checks preserve the
same mathematical constants without relying on coincident subobject addresses.

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

The first remediation packet subsequently fixed `AUD-079`, accepted and
documented the trusted-schedule and failed-Fp2-sqrt-output contracts in
`AUD-083`/`AUD-084`, removed the obsolete `glv_split_4` declaration for
`AUD-086`, and completed the `AUD-087` cleanup. Caller inventory narrowed the
historical `cyclotomic_exp` contract to GT. Follow-up work fixed `AUD-080` by
making the existing C/Rust Frobenius boundary reject powers outside 1-3 without
touching output; the direct internal primitive remains limited to the only
powers actually used. `AUD-081` is fixed: focused comparison at r, r+1, 2r-1,
2r, and 2^256-1 proves the GLV decomposition accepts the full 256-bit scalar
domain modulo the GT order, so the false Debug-only reduced-scalar assertion
was removed. The exact final quotient bound is two and is now asserted in both
scalar and native AVX splitters; the same five cases pass on Peach under native
AVX-512/IFMA. `AUD-085` now explicitly requires both checked `deserialize` and
`deserialize_unchecked` for protocols that defer one group check across an
aggregate. The next section records the subsequent `AUD-082`/`AUD-088`
constant-time remediation and performance experiment. Validation of the
earlier field packet was limited to directly affected objects and focused
tests, strict MSVC syntax, and optimized code inspection; no benchmark,
sanitizer, CUDA, portable-AVX Release, or broad suite ran in that packet.

## Primary 256-bit group multiplication

Commits `0ee24de` and `5bebdc0` validate `AUD-088` and `AUD-089`. The public
contracts are now distinct across C, native C++, and Rust:

- `mult` is subgroup-optimized and constant-time for canonical scalars at a
  fixed public `nbits`; noncanonical inputs are accepted modulo r, while timing
  may reveal that reduction was required;
- `mult_vartime` has the same point precondition and is restricted to public
  scalars;
- `unchecked_mult` is a generic constant-time operation for any on-curve point
  and never invokes subgroup endomorphisms.

Release builds deliberately do not validate those point preconditions. Debug
builds assert the bit-count, scalar-pointer, curve, subgroup, and GT-membership
invariants. Ordinary Pippenger is subgroup-only and its one-point fallback uses
`mult_vartime`; FFT, middle product, signing, and validated aggregation likewise
use the appropriate subgroup operations. Secret scalar preparation masks at
most 256 little-endian bits, skips classification through 254 bits, and
otherwise uses one branchless subtraction followed by a single
canonical/noncanonical branch. The noncanonical path reuses that difference
and performs one additional branchless conditional subtraction because the
complete input range is below 3r. Secret-derived scalar copies are explicitly
wiped. Public scalar preparation uses the same cheap canonical fast path
without secure wiping. The
post-correctness `AUD-090` pass extended that rule through the native AVX G1/G2
decomposition and executor: division scratch is erased after extraction,
superseded split objects are erased after their final copy, and the fixed Booth
schedule plus executor scalars are erased after the final group operation.
Public-scalar variable-time and MSM paths receive no additional wipe overhead.

Scalar checked multiplication restores the established constant-time GLV
crossover at 176 bits for G1 and GLS crossover at 144 bits for G2. Smaller
inputs, including the secondary 128-bit path, retain generic fixed-window
multiplication. AVX checked multiplication canonicalizes before its inherited
fixed-schedule decomposition. Unchecked multiplication always uses the generic
fixed-window path.

Sequential Peach candidate measurements used several deterministic dense
canonical 256-bit scalars and one benchmark process at a time. Scalar G1 and
G2 selected odd-table wNAF-5 over Booth width five by 6.62% and 12.20%
geometric-mean improvements. Native AVX G1 selected wNAF-5 by 6.08%; native
AVX G2 selected the simpler wNAF-4 path because it improved over Booth width
five by 11.71%, while wNAF-5 beat wNAF-4 by only 2.78% and did not satisfy the
3% complexity threshold. Constant-time scalar and AVX G1/G2 retained their
established width-five kernels because no candidate met the 3%/2% selection
rule. All losing forced candidates were removed.

GT constant-time multiplication retains four joint components, a final
16-entry owning table, 64 cyclotomic squarings, and a full constant-time table
scan. Its splitter now performs fixed-time reduction modulo r followed by one
`div_by_zz` and two `div_by_z` operations. The four unsigned radix-z digits use
conjugated odd Frobenius bases because the GT eigenvalue is -z. Debug builds
assert the canonical bound, r < z^4, each digit below z, and exact
recomposition. The unsigned splitter improved scalar variable-time GT by only
about 1.1%, so scalar retains the simpler balanced splitter; native AVX gained
about 7% and uses the unsigned splitter for its variable-time path.

Final sequential 256-bit production means, in microseconds, were:

| Backend | G1 CT / VT | G2 CT / VT | GT CT / VT |
| --- | ---: | ---: | ---: |
| scalar | 68.01 / 60.17 | 133.79 / 111.81 | 387.71 / 156.79 |
| native AVX-512/IFMA | 32.22 / 29.53 | 63.91 / 50.60 | 71.39 / 55.41 |

The final seven-test target passes under local scalar Release, local
portable-AVX Debug, and Peach native AVX-512/IFMA Release. It covers 0, 1,
r-1, r, r+1, 2r-1, 2r, and 2^256-1; exact output aliasing; infinity and GT
identity; ignored bits above non-byte-aligned lengths; scalar/AVX agreement;
and on-curve non-subgroup fixtures generated without cofactor clearing.
Checked Debug calls reject the latter while unchecked multiplication matches
an independent bitwise reference. The focused Rust CT/VT method test passes.

GCC 15.2 optimized native inspection finds no AVX gathers, allocation
references, or out-of-line `_vv256`/`_vv512` wrappers in `blst_gt_mult`.
The two stack-probe pages plus final 0x700-byte subtraction record a 0x2700
(9,984-byte) native GT frame, well below the previous 54,144-byte acceptance
cap. Source and optimized inspection show fixed decomposition, table-scan, and
executor loop bounds; statistical leakage testing remains a future evidence
upgrade. No full suite, sanitizer, CUDA, BTx/BTE, Apple CI, or portable-AVX
Release job ran for this packet.

The `AUD-090` follow-up reused the `vec_zero` primitive already validated by
`AUD-045`; plain `memset` is not used for secret erasure. Portable-AVX Debug and
Peach native AVX-512/IFMA Release focused group tests pass. GCC 15.2 optimized
inspection retains the exact expected fixed-size volatile erase regions: G1
32/32/104/32 bytes for division scratch, split output, Booth schedule, and
executor scalars; G2 retains 32/16/16/32/16/16/104/32 bytes for its three
division scratches, three split objects, schedule, and executor scalars. The
erase-loop branches depend only on compile-time sizes. This evidence covers
explicit object storage; neither C++ nor the platform ABI promises erasure of
all transient compiler-created register copies.

The subsequent full `blst/avx/curve/scalar_mult.h` review made the independent
correctness oracle permanent. Twenty boundary, crossover, and deterministic
full-width scalar cases now compare G1/G2 CT, VT, and unchecked multiplication
against simple double-and-add, and compare GT CT/VT against basic cyclotomic
square/multiply exponentiation. The focused target passes 8/8 in scalar Release,
portable-AVX Debug, and Peach native AVX-512/IFMA Release. Native G1/G2 CT symbol
sizes remain exactly `0x13e9` and `0x3033`; no benchmark was needed.

Caller inventory opened `AUD-091`: the general AVX Pippenger entry silently
selects GLV/GLS for every 256-bit input even though its public contract does not
require subgroup points or canonical scalars. The BLST split accepts a 255-bit
value and has only a 128-bit quotient; `2^256-1` produces a 129-bit quotient.
Existing MSM tests use subgroup points and mask the top scalar nibble to `0x0f`,
so they miss both violations. This requires a dedicated MSM contract/design
packet; canonical reduction alone cannot repair the non-subgroup case.

The final group-domain recalibration validates `AUD-091`, `AUD-092`,
`AUD-093`, and `AUD-095`. Unsuffixed scalar and AVX MSM is consistently the
subgroup operation for every bit width and point count. Full-width GLV/GLS
preparation reduces public scalars before division and consumes caller-owned
split scratch; plain Pippenger does not preprocess scalars because exact and
reduced multiplication agree on subgroup points. Projective one-point MSM calls
`mult_vartime` and retains its no-scratch path. The experimental unchecked MSM
family was removed; only upstream's two legacy single-point `unchecked_mult`
exceptions remain.

The GT-only GLV4 operation is hard-named `blst_avx_gt_mult_pippenger` with a
semantic `blst_avx_gt` storage alias. The CatBLST-only Fp12-named spellings were
removed rather than retained as misleading compatibility forwards. Rust
`gt_msm` calls the GT-named entry under an explicit previously-validated-GT
precondition. The dedicated target excludes frozen BTx/BTE and covers
one/two-point subgroup modulo-r consistency, 257-bit rejection, complete
membership predicates, and GT-only naming. Local MSVC portable-AVX Debug passes
5/5. That run exposed a generic affine G2 membership defect: an Fp2 projective
`z` was initialized from G1-sized `Fp::one()` storage; the implementation now
uses `Field::one()`. The structural contract check covers 242 public symbols in
51 disjoint families and rejects every removed unchecked-MSM and Fp12-MSM
spelling. Native AVX/Rust/codegen/benchmark evidence for the recalibrated fast
path is complete. Peach GCC 15.2 Release passes the 10-test scalar packet and
the 5-test native AVX-512/IFMA contract packet. Focused Rust multiplication
passes and the full-Rust-binding AVX library type-check succeeds. Optimized
scalar and native AVX G1 code has the approved branchless classification, one
canonical/reduction branch, branchless slow-path selection, and retained
zeroization; selected-vector conversion wrappers are absent.

Sequential CPU-pinned benchmarks compare the exact pre-packet commit
`2744efc` with commits `b72f33c`, `036313a`, and `a383828` under identical
flags. Canonical G1/G2/GT multiplication and MSM sizes 1/8/64/256 meet the 2%
gate. The largest regression is GT MSM at one point (`+1.193%`); every other
primary regression is at most `+0.458%`. G1/G2 one-point MSM improves by
39%/53%, and GT CT multiplication improves by about 2.9%. An initial redundant
GT reduction caused a 2.1-3.8% regression and was rejected; the balanced GLV4
splitter's existing complete-256-bit modulo-r preparation is now used directly.
No CUDA, portable-AVX Release, broad suite, sanitizer, Apple job, or push ran.

That strict all-header check also opened and fixed `AUD-075`: MSVC C rejected
the AVX opaque-shell alignment when its `__declspec(align(...))` argument was
spelled as an enum identifier. The fixed expert-ABI sizes and 64-byte alignment
are now prefixed preprocessing constants, leaving their values and layouts
unchanged.

The field and point API documentation lives with the core declarations in
their owning field and curve headers. The C binding header retains the C ABI
storage and aliasing contract, while the core headers retain the Montgomery
representation, constant-time, failure-output, aliasing, subgroup, encoding,
GLV, and scalar contracts. The real packed AVX API remains separate in
`blst::avx` and `blst_avx.h`; four
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

## Incremental API documentation pass

The audited field and focused group families now carry declaration-level
contracts on the C, selected native C++, and safe Rust surfaces. The pass adds
language-specific and audit-closure rules to `docs/API_CONTRACTS.md`, tracks 14
families in `PUBLIC_API_DOCUMENTATION.csv`, and adds a build-free structural
checker. Generated raw Rust FFI inherits the corresponding C contract and stays
deferred to the Phase 18 generator review rather than receiving
hand-maintained comments.

The pass fixes the documentation portion of `AUD-100`, states the current
incomplete-addition preconditions without closing the `AUD-096` design choice,
and records `AUD-103`. Scalar GT conditional negation now matches the AVX
backend's constant-time contract by selecting an unconditionally computed
conjugate through the audited conditional-copy primitive. The scalar
variable-time MSM helper retains its public signed-digit branch locally, so the
API correction does not add a conditional-copy cost to that hot path.

No ABI, generated binding, benchmark, CUDA, BTx, or runtime dispatch changed in
this packet.

## Generated Fp addition chains

The complete p-2 reciprocal and (p-3)/4 reciprocal-square-root schedules now
have retained symbolic verification under `AUD-181`. Every scratch read,
hexadecimal exponent trace, final exponent, and operation count is checked.
Focused scalar Debug/Release tests compare both schedules with independent
generic exponentiation and cover exact aliasing and the actual zero semantics.

`AUD-182` documents that native fixed inversion maps zero to zero because
branchless infinity signing relies on it, while the public C inverse retains a
nonzero precondition. `AUD-183` moves the Fp-specific helpers from the global
generic namespace into `blst::scalar::detail`; `AUD-184` replaces byte copying
with typed assignment and pins original BLST commits. The arithmetic schedules,
880/832-byte MSVC frames, and 86/85 arithmetic call counts are unchanged; the
typed copy removes one useless entry pointer comparison from each function.

Proposed broad scratch erasure (`AUD-180`) was explicitly declined. Generic
constant-time arithmetic does not promise comprehensive transient erasure, and
no generic hot-path zeroing was added. No benchmark or broad test ran.

## Selected field-facade closure

Wave 1 completed the small selected native C++ field facades `Fp.h`, `Fp2.h`,
`Fp12.h`, and `Fr.h` at commit
`aa4c7ed187f5a63496f0e94f9284108cbd46c7f1`. They are exact aliases and
declaration re-exports over the owning scalar implementations; none adds a
wrapper, conversion, dispatch layer, allocation, or executable arithmetic.
`AUD-196` replaces the reserved Fp2 include guard, and the locally composed Fr
facade now carries its license/provenance header. The `AUD-070` exact-alias
follow-up is implemented at `02f5e5ff078b251ca3b5f86f4f83ea64af664bc7`.

Strict C++20 parsing and the standalone Release field-boundary executable pass.
The owning scalar field implementations retain responsibility for arithmetic,
constant-time, sanitizer, and cross-platform evidence.

## Wave 2 scalar and extension-field audit

Wave 2 completed the line-by-line rereview of scalar `Fr`, `Fp`, `Fp2`, and
`Fp12` at implementation commit `314c186`. Findings `AUD-198` through
`AUD-207` record exponent/timing contracts, REDC bounds, missing Fp2 transforms,
inverse timing evidence, square-root publication, extreme shift work, Fp12
domain wording, a bad cyclotomic test, and sparse packed-operand semantics.
`AUD-060` was extended because portable modular conditional negation disagreed
with assembly for non-bit flags.

The approved repairs preserve hot structure. Portable `cneg` uses branchless
truth normalization. Shift counts through 256 use the original selected
backend; only larger public counts use O(log count) exponentiation. Optimized
MSVC shows constant count three still calls the original 256/384-bit backend
directly. Focused scalar Release tests pass 7/7 across portable/assembly cneg,
extreme counts, Fp2 transforms/sqrt, and valid-domain cyclotomic squaring.
Scalar Debug and Release libraries compile serially. No benchmark, CUDA, BTx,
portable-AVX Release, sanitizer, Apple job, or broad runtime sweep ran.

`AUD-201` remains open: Fp/Fp2 reciprocal includes a result-verification branch
and rare fixed-chain fallback, so the conservative variable-time contract is
retained until selected-platform code generation is reviewed. `AUD-199` remains
fixed rather than validated until its Debug assertion has a retained death
regression. Broader independent field oracles, fuzzing, sanitizers, Rust, and
supported-platform evidence remain Phase 14/19 work.
