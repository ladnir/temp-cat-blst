# Group-Domain API Contract Audit

Date: 2026-08-04

State: subgroup API packet validated locally and on native AVX-512/IFMA

Scope: supported C, native C++, and Rust surfaces; scalar and AVX backends

Excluded: CUDA, BTx/BTE, benchmark file review, and removed/unsupported bindings

## Purpose

This audit began by recording the behavior of the API and implementation before
semantic changes. The approved recalibration gives a deliberately small basic
primitive layer an explicit on-curve contract and keeps advanced G1/G2
algorithms subgroup-only. It preserves upstream's two legacy
`unchecked_mult` exceptions, removes the experimental unchecked-MSM family,
and hard-names the GT-only AVX MSM without Fp12 compatibility forwards. It also
makes public membership predicates complete curve-plus-subgroup checks. The
independent memory-safety defect `AUD-101` remains remediated.

The repeated source of ambiguity is that these three statements are different:

1. an operation validates subgroup membership;
2. an operation assumes subgroup membership was validated earlier;
3. an operation remains mathematically correct outside the subgroup.

A Debug assertion records a precondition. It is not Release validation. A
subgroup-relation test is also not necessarily a combined curve-and-subgroup
decoder.

## Inventory control

`GROUP_API_CONTRACTS.csv` classifies the public C symbols from
`bindings/blst_scalar.h` and `bindings/blst_avx.h` into 53 disjoint contract
families. Each row records:

- mathematical input and output domains;
- Release validation and Debug assertions;
- scalar interpretation and timing contract;
- invalid-input behavior;
- the trust boundary at which the operation is safe to use.

`scripts/audit/check_group_api_contracts.ps1` extracts the public declarations
and fails if a relevant symbol is missing, matches multiple families, or a
contract family becomes empty. It is intentionally separate from the full file
manifest checker so it can run after every API edit without a build.

## Terminology used by this audit

Normative contract language is defined in
[`../API_CONTRACTS.md`](../API_CONTRACTS.md). The table below is the subset of
mathematical-domain shorthand used frequently in this audit; it does not
override the public vocabulary.

| Term | Meaning |
| --- | --- |
| general Fp12 | Any ordinary Fp12 field element. |
| cyclotomic | An Fp12 element satisfying the cyclotomic relation; not necessarily in the order-r target group. |
| GT | A nonzero element of the prime-order pairing target subgroup. |
| valid representation | Coordinates/limbs obey the in-memory representation rules; this says nothing by itself about curve or subgroup membership. |
| on-curve | A point satisfies the curve equation, including the library's infinity representation. |
| subgroup | An on-curve point in the prime-order G1/G2 subgroup. |
| relation-only subgroup test | An internal endomorphism relation that is meaningful only after on-curve status is established. Public membership predicates compose the curve equation with this relation. |
| validating operation | Rejects invalid input in Release and reports status. |
| assumed-valid operation | Performs arithmetic under a documented precondition; it may assert only in Debug. |

Public point-domain preconditions use only **on-curve** and **subgroup**.
Subgroup implies on-curve. Validators need no prior point-domain precondition.
Infinity belongs to both domains; protocol-specific nonidentity rules remain
separate. `POINT_DOMAIN_MATRIX.md` records the exact approved C/C++/Rust
families.

## Cross-surface results

### Point decoding and membership

Scalar and AVX `deserialize`/`uncompress` enforce canonical field encodings,
flags, infinity encoding, and the curve equation. They do not perform a full
prime-order subgroup check. The scalar decoder additionally rejects the known
`x == 0` exceptional point, but that is not a general subgroup test.

`blst_p1_in_g1`, `blst_p2_in_g2`, their affine forms, and the AVX equivalents
now evaluate the curve equation before the fast subgroup endomorphism relation.
Those individual predicates are complete for the mathematical property named
by the API, while protocol-specific infinity policy remains separate.

The AVX expert batch predicates validate every point and shift as on-curve
before evaluating caller-specified randomized relations. They are sound only
when the shift is in the subgroup and coefficient rows satisfy the documented
challenge model. Zero rounds and all-zero rows remain valid expert inputs, so
these entries are not general validation boundaries.

`blst_avx_p1s_in_g1` and `blst_avx_p2s_in_g2` provide the friendly boundary.
They validate every curve equation, obtain one fresh unpredictable 32-byte seed
after the immutable batch is fixed, and derive 81 independent uniform ternary
rows with domain-separated SHA-256. CatBLST fixes the subgroup shift. Because
the selected curves have no nontrivial 2-torsion, conditioning on all other
coefficients leaves at most one coefficient in {-1,0,+1} that can cancel a
fixed nonzero quotient component. Under the SHA-256 PRF model, false acceptance
is therefore at most `(1/3)^81 < 2^-128`. `AUD-161` is validated. Native Peach
AVX-512/IFMA Release passes the friendly acceptance, rejection, entropy, and
identity-MSM cases; the focused G1/G2 contract test averages 6.53 milliseconds.

The normal safe sequence for adversarial point bytes is:

1. deserialize/uncompress, which establishes canonical encoding and on-curve;
2. call the relevant complete `in_g*` predicate, or a higher-level validator;
3. use subgroup-only arithmetic.

This sequence is now stated for the focused low-level membership family. The
high-level Rust decoding policy remains `AUD-094` and will be documented when
those containing files receive their complete review.

### Point arithmetic and multiplication

The scalar `Point::add` overloads use incomplete addition formulas requiring
distinct, nonopposite finite inputs (with infinity handled). Opposite finite
points have the same affine x coordinate and negated affine y coordinates, so
their mathematical sum is infinity. The C functions
`blst_p1_add`, `blst_p2_add`, and their affine-add variants select those
formulas. The C and native comments now state the current incomplete-formula
preconditions and direct callers to `add_or_double` when equality is possible.
`AUD-096` is accepted: CatBLST retains the inherited names and fast narrow
formula, while `add_or_double`, subtraction, and their corresponding complete
paths handle equality/opposites. Safe high-level Rust already selects the
complete operation; generated raw Rust inherits the C declaration contract.

The approved wider contract is intentionally limited to basic single-point
primitives: add/subtract/double, negation, representation conversion,
observation, encoding, and the two inherited exact-scalar multiplications.
Advanced batch conversion/addition, fixed-table multiplication, MSM, pairing,
and protocol algorithms remain subgroup-only even if their current formulas
happen to evaluate correctly on some wider inputs. This keeps the supported
surface small and lets users construct specialized wider algorithms from the
basic primitives when needed.

Single-point `mult` and `mult_vartime` now use subgroup endomorphisms for their
primary widths. They require previously validated on-curve subgroup points and
reduce the complete 256-bit scalar domain modulo r. `unchecked_mult` is
constant-time generic curve multiplication: it requires on-curve input but is
correct outside the prime-order subgroup and applies the exact unsigned scalar.
The name "unchecked" therefore currently describes a wider algebraic domain,
not merely omitted validation.

Scalar batch conversion/addition, Pippenger/wbits, and native-Fr MSM are
subgroup operations. Some current
low-level kernels happen to use formulas that also work on wider curve inputs,
but CatBLST does not promise or test that incidental domain. The one-point
Pippenger fallback calls subgroup `mult_vartime`.

The AVX MSM contracts are now point-count and bit-width independent:

- unsuffixed affine and projective Pippenger entries require previously
  validated subgroup points and interpret at most 256 scalar bits modulo r;
  full-width GLV/GLS preparation canonicalizes every scalar before splitting;
- every fallible boundary rejects `nbits > 256`; Debug builds assert the stated
  point-domain preconditions;
- projective one-point fallbacks call subgroup `mult_vartime` and require no
  scratch. There is no CatBLST unchecked MSM family.

G1/G2/GT `mult` accepts noncanonical values modulo r. Canonical inputs take a
cheap classification/copy path and retain the fixed-`nbits` constant-time
contract. Noncanonical inputs take a reduction path whose timing may reveal
only that reduction was required. Public `mult_vartime`, weighted pairing, and
G1/G2 decomposed MSM use the corresponding non-wiping public-scalar
preparation. The balanced GT GLV4 splitter already consumes the complete
256-bit domain modulo r, so GT MSM does not perform a redundant standalone
reduction. Plain Pippenger does not preprocess scalars when no decomposition
needs it.

### Fp12, cyclotomic, and GT

The storage types do not enforce these domains. Both scalar and AVX GT wrappers
can contain arbitrary Fp12 data.

- General `fp12_mul`, `fp12_sqr`, equality, conversion, and Frobenius operations
  are valid for general Fp12.
- `fp12_cyclotomic_sqr` requires a cyclotomic input. The C and native C++
  declarations now state that precondition, state that Release does not
  validate it, and direct general Fp12 callers to `fp12_sqr`.
- GT subtraction/negation rely on conjugation as inversion, doubling may use a
  cyclotomic shortcut, and optimized exponentiation relies on order-r
  identities. These are GT-only operations with no Release membership check.
- Final exponentiation accepts a nonzero general Fp12/Miller-loop value and
  produces GT. Zero is not rejected through a fallible API.
- `blst_avx_gt_mult_pippenger` and the semantic `blst_avx_gt` storage alias
  state the GT-only contract. The CatBLST-only Fp12 spelling was removed rather
  than retained as a misleading compatibility forward. Rust `gt_msm` calls the
  GT-named symbol and explicitly requires previously
  validated GT elements; it performs no repeated membership check.

There is no identified use case for generic-Fp12 GLV4 MSM. General Fp12 remains
necessary for Miller outputs, multiplication/squaring, layout conversion, and
final exponentiation.

### Pairing and verification

Raw Miller-loop and raw pairing aggregation routines assume subgroup points and
perform no subgroup validation. Their output is a general Fp12 Miller product,
not GT until final exponentiation. The underlying formulas may evaluate on a
wider curve domain, but CatBLST deliberately does not promise that advanced
contract.

The default aggregate-verification functions are intentionally asymmetric:

| Entry | PK infinity | PK subgroup | Signature subgroup |
| --- | --- | --- | --- |
| `blst_pairing_aggregate_pk_in_g*` | rejected | not checked | checked |
| `blst_pairing_mul_n_aggregate_pk_in_g*` | rejected | not checked | checked |
| `blst_pairing_chk_n_*` | rejected | complete membership when requested | complete membership when requested |
| `blst_core_verify_pk_in_g*` | rejected | checked | checked |
| `blst_pairing_raw_aggregate` | not checked | not checked | not checked |
| `blst_avx_pairing` | not checked | not checked | not checked |

The default assumes public keys were validated once before aggregation, while
checking the possibly aggregated signature contribution. CatBLST retains this
upstream-compatible policy because repeated PK checks add avoidable work. A
false explicit-check flag transfers the subgroup-point precondition to the
caller; it does not opt the operand into an on-curve-only domain. Weighted
aggregation additionally reaches subgroup-only point multiplication, making
prior subgroup validation a correctness precondition and a protocol-security
precondition.

`blst_avx_pairing` and raw Miller helpers remain expert subgroup-only
operations. Their result flags report boundary or backend availability, not
mathematical validation. CatBLST adds no duplicate validated/raw pairing
family: callers use the explicit check flags or the one-shot core verifier when
validation is required at this boundary.

The weighted C entry points use subgroup `mult_vartime`, reject `scalar ==
nullptr` with nonzero `nbits`,
and reject `nbits > 256` with `BLST_BAD_SCALAR` before changing the pairing
context. Safe Rust performs the same bound check before native entry. The
shared fixed-size canonicalizer also fails closed to scalar zero in Release if
an internal or legacy void caller violates those shape preconditions, closing
the memory-safety defect independently of higher-level validation (`AUD-101`).

### Hashing, mapping, signing, and public-key derivation

`encode_to_g*`, `hash_to_g*`, and `map_to_g*` all clear the cofactor and return
prime-order subgroup points. The C comments do not currently state this useful
postcondition.

The low-level signing and public-key functions consume a 32-byte scalar through
the GLV/GLS splitter without Release validation. Public declarations now
require a canonical nonzero key produced by keygen or accepted by
`blst_sk_check`; signing additionally requires an already validated subgroup
hash point. Native scalar and AVX implementations encode those invariants with
Debug assertions only. The inherited `*2` byte outputs are explicitly
documented as full-width pre-decorated buffers, not directly usable encodings.

## Native C++ surface

| Family | Current behavior |
| --- | --- |
| `Point::add` | Basic on-curve primitive: incomplete distinct/nonopposite addition. |
| `Point::dadd` | Basic on-curve primitive: complete add-or-double path. |
| `Point::mult` / `mult_vartime` | Subgroup-only, scalar reduced modulo r, Debug membership assertions, no Release validation. |
| `Point::unchecked_mult` | Generic on-curve constant-time multiplication; no Release validation. |
| `Point::in_subgroup` / affine equivalent | Complete curve-equation plus subgroup-relation predicate. |
| `AffinePoint::deserialize` / `uncompress` | Canonical and on-curve decoding, but not a complete subgroup check. |
| `Point::sign` / `to_pk` | Assumes a valid canonical secret scalar; signing also assumes a subgroup hash point. |
| `Fp12::cyclotomic_*` | Cyclotomic or GT preconditions vary by operation and are documented in the native header. |
| `GT` construction | `GT(Fp12)` accepts arbitrary Fp12 and establishes no invariant. |
| `GT` arithmetic | Assumes GT where conjugation, cyclotomic squaring, or scalar decomposition is used. |
| `GT::deserialize` | Canonical, failure-atomic decoder that establishes GT membership; returns `false` for any coefficient at least p or a canonical non-GT value. |
| `GT::deserialize_unchecked` | Canonical, failure-atomic Fp12 decoder that intentionally defers only GT membership; candidates must be combined as general Fp12 values and validated before ordinary GT operations. |

The native types are high-performance semantic wrappers, not validated
capability types. That can remain a valid design, but every operation must say
whether it accepts the representation, on-curve, cyclotomic, or subgroup
domain.

Selected scalar and AVX `GT::cneg` implementations are constant-time for a
zero-or-one condition bit (`AUD-103`). Scalar `cneg` constructs the conjugate
unconditionally and uses the audited `ccopy`/`vec_select` primitive. The scalar
variable-time MSM helper keeps its public signed-digit branch locally, avoiding
a full conditional copy in that hot path without weakening the selected API.

## Rust surface

The generated/raw Rust FFI mirrors C and is `unsafe`; its mathematical
contracts are exactly the C contracts.

The low-level conveniences in `primitives.rs` enforce scalar slice length and
`nbits <= 256`, but they do not validate the receiver. Their multiplication
documentation states the subgroup/on-curve preconditions and exact
little-endian bit interpretation. Miller-loop helpers require previously
validated subgroup points, matching the C/native advanced-operation policy;
they do not promise the wider on-curve domain merely because the underlying
formula may evaluate there. Returned cryptographic values and decisions are
marked must-use.

The higher-level BLS types behave as follows:

| Rust operation | Current behavior |
| --- | --- |
| `PublicKey::deserialize` / `from_bytes` / Serde `Deserialize` | Canonical and on-curve only; does not establish subgroup or reject infinity. |
| `PublicKey::validate` / `key_validate` | Rejects infinity and evaluates complete curve-plus-subgroup membership. |
| `Signature::deserialize` / `from_bytes` / Serde `Deserialize` | Canonical and on-curve only; does not establish subgroup. |
| `Signature::validate` / `sig_validate` | Evaluates complete curve-plus-subgroup membership and optionally rejects infinity. |
| `AggregatePublicKey` / `AggregateSignature` constructors and `From` conversions | Preserve whatever validation state the source value had; types do not encode it. |
| aggregate/add methods | Boolean parameters decide whether each input receives subgroup validation. |
| verify methods | Boolean parameters are forwarded to the checked C pairing entries; false transfers responsibility to the caller. |
| fast aggregate verify | Deliberately assumes proof-of-possession-validated public keys. |
| `Pairing::aggregate` | Exposes independent PK/signature validation flags. |
| `Pairing::mul_n_aggregate` | Exposes the same flags; rejects scalar slices shorter than `nbits` and rejects `nbits > 256` before FFI. |
| AVX `gt_msm` | Safe Rust wrapper over GT-only GLV4; accepts storage-compatible `blst_fp12` under an explicit previously-validated-GT precondition. |

Rust memory safety does not imply cryptographic-domain validation. In
particular, the safe weighted-pairing bit-length defect is a memory-safety issue
and is not merely a mathematical precondition.

## Recorded findings

| Finding | Summary | State |
| --- | --- | --- |
| `AUD-091` | Automatic 256-bit AVX affine MSM silently changed point/scalar domain. | validated |
| `AUD-092` | AVX/Rust Fp12 MSM was actually GT-only. | validated |
| `AUD-093` | Public subgroup predicates were relation-only and did not establish on-curve. | validated |
| `AUD-094` | Point decoders and Rust deserializers establish on-curve, not subgroup; contracts are not explicit. | validated |
| `AUD-095` | AVX projective MSM changed point domain at one point and lacked an nbits bound. | validated |
| `AUD-096` | Inherited unsuffixed add remains the fast incomplete operation under explicit finite distinct/nonopposite preconditions; complete add-or-double remains available. | accepted |
| `AUD-097` | Pairing retains the validate-once-PK policy under explicit complete-membership and raw subgroup contracts. | validated |
| `AUD-098` | Signing, weighted pairing, and AVX multiplication rely on explicit narrow subgroup/SK preconditions without Release validation. | validated |
| `AUD-099` | GT is represented by Fp12 storage under explicit GT-operation preconditions and a separate MillerResult native type. | accepted risk |
| `AUD-100` | C cyclotomic squaring states its cyclotomic precondition and lack of Release validation. | validated |
| `AUD-101` | Weighted pairing permitted oversized nbits into a fixed 32-byte scalar buffer; safe Rust could reach it. | validated |
| `AUD-103` | Selected scalar/AVX GT conditional negation is constant-time across both backends. | validated |
| `AUD-160` | Basic on-curve primitives and subgroup-only advanced operations lacked one bounded cross-language policy. | validated |
| `AUD-161` | Expert AVX batch soundness is caller-controlled; friendly G1/G2 APIs now own fresh seeded challenge construction and target a false-accept bound below 2^-128. | validated |

## Incremental public documentation control

`PUBLIC_API_DOCUMENTATION.csv` tracks the declaration coverage of API families
whose implementations have already received a complete or focused hand audit.
It currently accounts for 23 field, scalar, group, and expert-AVX families
across C, selected native C++, safe Rust, and generated raw Rust FFI.

`scripts/audit/check_public_api_documentation.ps1` verifies every ledger source,
anchor, scope, and disposition. Generated raw Rust declarations explicitly
inherit the C contract; complete Rust generator/comment propagation remains a
Phase 18 review item. Deferred entries do not claim that their containing file
is audited. The per-file checklist now makes cross-language declaration
documentation an audit completion gate.

## Remaining decision queue

The selected point policy is recorded in `POINT_DOMAIN_MATRIX.md`: deserialize
and uncompress establish on-curve only; the basic single-point primitive layer
explicitly accepts on-curve inputs; advanced algorithms require subgroup
inputs; and no new validation-state types or generic-MSM families are added.
The remaining questions require separate review:

1. Should the default pairing aggregate names preserve the upstream
   "signature checked, PK prevalidated" policy, or should validated and raw
   variants be named explicitly?
2. Should signing/public-key C entry points remain precondition-only void APIs,
   or gain fallible validated variants?
3. What supported-native performance does the new friendly AVX batch
   membership boundary deliver, and should its 81-row construction later be
   supplemented by a separately reviewed wide-coefficient MSM design?
`AUD-101` was independent of those semantic choices and is remediated. Legacy
void arithmetic APIs retain documented 256-bit preconditions, while the shared
canonicalizer now prevents a violated precondition from becoming memory
corruption in Release.

## Validation performed

- Complete static trace from both public C headers through `c_api.cpp`,
  `avx_c_api.cpp`, scalar/AVX point and GT implementations, pairing aggregation,
  native headers, and Rust safe/raw wrappers.
- Focused MSVC scalar Debug tests use an order-11 G1 point and order-13 G2 point
  to exercise incomplete/complete addition, subtraction, doubling, conditional
  negation, representation conversion, observation, compression, serialization,
  and separate subgroup validation. The three selected group-contract cases pass.
- Focused portable-AVX Debug tests exercise the same low-order conversion,
  complete addition, negation, observation, compressed/uncompressed decoding,
  and a nonzero challenged batch rejection; both selected cases pass.
- Both safe Rust variants pass compressed, uncompressed, `from_bytes`, and Serde
  decoding of those on-curve low-order fixtures and reject them through the
  complete individual validators.
- The contract inventory validator passes with 243 public C symbols in 53
  nonempty disjoint families. Structural checks find no unchecked MSM symbol
  and no legacy AVX Fp12-MSM symbol in retained source.
- Focused scalar group multiplication passes 10/10 under local MSVC Debug and
  Release and under Peach GCC 15.2 native AVX-default Release. The cases cover
  `0`, `1`, `r-1`, `r`, `r+1`, `2r-1`, `2r`, `2^256-1`, exact output aliasing,
  ignored high bits, crossover widths, subgroup MSM, complete membership, and
  independent double-and-add/cyclotomic exponentiation oracles.
- The dedicated no-BTx/BTE AVX contract target passes 5/5 under local portable
  AVX Debug and Peach native AVX-512/IFMA Release. It covers one/two-point G1/G2
  subgroup MSM, noncanonical public scalars, 257-bit rejection, off-curve batch
  rejection, and the GT-only name against the scalar backend. The maximal
  noncanonical GT scalar also validates the balanced splitter's integrated
  modulo-r contract.
- The focused Rust multiplication test passes, and
  `cargo check --features "full-rust-bindings,avx"` passes with two build jobs.
- GCC 15.2 optimized G1 and native AVX G1 inspection shows a branchless
  four-limb `sub`/`sbb` classification, one canonical-versus-reduction branch,
  a branchless second-subtraction selection, and retained fixed-size
  zeroization loops. The native AVX object has no unresolved `vv512`, `cc512`,
  or selected-vector conversion wrapper. G2 has the same scalar-preparation
  structure.
- Sequential, CPU-pinned Peach benchmarks compared commit `2744efc` with the
  candidate using identical GCC 15.2 Release AVX-512/IFMA flags and five
  repetitions. Canonical G1/G2/GT multiplication and MSM sizes 1, 8, 64, and
  256 all satisfy the 2% gate: the largest regression is GT MSM at one point
  (`+1.193%`), every other primary regression is at most `+0.458%`, G1/G2
  one-point MSM improves by 39%/53%, and GT CT multiplication improves by
  approximately 2.9%.
- AUD-101 focused validation: scalar Release C/C++ boundary and fail-closed
  helper tests pass; three full-Rust-binding scalar-shape tests pass.
- Native compilation repeats the already recorded low-severity `AUD-102`
  `-Wignored-attributes` warnings for `std::array<__m512i>`; this packet does not
  broaden warning suppressions or change vector layout.
- Code commits: `b72f33c` (scalar preparation and focused benchmark fixture),
  `036313a` (reduced unchecked API and GT hard rename), and `a383828` (complete
  membership predicates and focused tests).
- No sanitizer sweep, portable-AVX Release build, full suite, CUDA, BTx/BTE,
  Apple CI, or push ran.
