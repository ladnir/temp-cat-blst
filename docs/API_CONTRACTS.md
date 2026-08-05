# CatBLST API Contract Vocabulary

This document defines the vocabulary used by CatBLST API contracts. Public C,
C++, and Rust documentation should use these terms consistently. A function's
individual contract remains authoritative for its inputs and outputs; a name
such as `fast`, `checked`, or `unchecked` is never a substitute for that
contract.

CatBLST is a low-level cryptographic library. Fast operations may deliberately
have narrow preconditions. The objective is precise, predictable semantics,
not implicit validation at every arithmetic boundary.

## Contract language

| Term | Meaning |
| --- | --- |
| precondition | A condition the caller must establish before calling the operation. Unless the contract explicitly says the operation validates it, violating a precondition is outside the operation's supported domain. |
| postcondition | A condition the operation establishes on successful return when all preconditions hold. |
| invariant | A condition that remains true across the stated lifetime or sequence of operations. |
| validates | Checks the stated property in Release builds and reports failure through the documented result. |
| assumes | Requires the stated property as a precondition and does not necessarily check it in Release. |
| preserves | Produces an output with the stated property when the corresponding input preconditions hold. It does not imply that the property was checked. |
| rejects | Detects the stated condition in Release and returns the documented failure result. |
| Debug assertion | A development-time check of an invariant or precondition. It is not Release validation and is not a recoverable error mechanism. |
| fallible | Reports whether the operation succeeded, normally through `BLST_ERROR`, `bool`, or a null result. |
| infallible | Has no failure result for inputs satisfying its preconditions. This does not mean every bit pattern is accepted. |
| fail closed | On a detected invalid condition, returns a safe documented failure result and does not accidentally continue with broader authority or corrupted state. |

Normative words have their usual meanings:

- **must** or **requires** states a precondition or mandatory behavior;
- **must not** states a forbidden input, overlap, or use;
- **may** states permitted behavior;
- **undefined by CatBLST** means the library makes no semantic guarantee for
  that precondition violation, even when the C or C++ language would not itself
  classify the behavior as undefined behavior.

## Representation and encoding

| Term | Meaning |
| --- | --- |
| valid in-memory representation | The object satisfies its C/C++ layout and internal representation requirements. This alone establishes no curve, subgroup, or protocol validity. |
| canonical field element | A field element represented in the unique range required by the stated interface, normally `[0, p)` or `[0, r)`. Internal Montgomery or lazy representations are canonical only when their specific contract says so. |
| canonical encoding | A unique accepted byte encoding with required range, flag, and padding rules. Canonical decoding does not by itself imply subgroup membership. |
| deserialize | Parse the documented byte representation and establish only the postconditions explicitly listed by that function. |
| serialize | Produce the documented byte representation from an input satisfying the serializer's preconditions. Serialization does not validate an arbitrary in-memory input unless explicitly stated. |
| exact aliasing | An output begins at the same address as an input and covers the same object. |
| disjoint | The relevant input and output storage ranges do not overlap. |
| partial overlap | Storage ranges overlap without exact aliasing. Partial overlap is forbidden unless explicitly supported. |

## Curve and group domains

| Term | Meaning |
| --- | --- |
| infinity | The identity point of the elliptic-curve group, represented using CatBLST's documented projective or affine convention. |
| on-curve | The point satisfies the curve equation or is the accepted infinity representation. This does not imply membership in the prime-order subgroup. |
| prime-order subgroup | An on-curve point belonging to the order-`r` subgroup used by BLS12-381 protocols. |
| G1 / G2 element | A point in the corresponding prime-order subgroup, including its identity, unless a narrower protocol contract excludes infinity. |
| subgroup relation | The fast endomorphism relation used to test subgroup membership after on-curve status is known. A relation-only test does not independently establish on-curve. |
| complete group validation | Establishes every property named by the contract, normally valid representation, on-curve, subgroup membership, and any required infinity policy. |
| generic curve operation | Mathematically correct for every on-curve input in its stated curve domain, including non-subgroup points unless explicitly excluded. |
| subgroup-optimized operation | May use endomorphisms or reduction modulo `r` and therefore requires prime-order subgroup inputs. |
| incomplete addition | An addition formula with documented exceptional input pairs, such as equal or opposite finite points. |
| complete addition | Handles every pair of valid on-curve inputs covered by its infinity policy. |

Protocol validity is separate from mathematical group validity. For example, a
valid G1 element can still be forbidden as a public key because it is infinity.

Ordinary G1 and G2 arithmetic uses the subgroup domain by default. A function
need not repeat `subgroup` in its name to require a G1/G2 element, and Release
arithmetic normally assumes rather than validates that precondition. CatBLST
does not create parallel `unchecked` families for MSM, addition, or other
ordinary operations. The two inherited `blst_p1_unchecked_mult` and
`blst_p2_unchecked_mult` symbols are legacy upstream exceptions: they accept
on-curve non-subgroup points and apply exact scalars without endomorphisms.

## Fp12, cyclotomic, and GT domains

| Term | Meaning |
| --- | --- |
| general Fp12 | Any ordinary element of the extension field Fp12. |
| nonzero Fp12 | A general Fp12 element other than field zero; required where inversion is used. |
| Miller-loop value | A general Fp12 value produced by Miller-loop accumulation before final exponentiation. It is not assumed to be in GT. |
| cyclotomic element | An Fp12 element in the cyclotomic subgroup. This is a broader domain than GT. |
| GT element | A nonzero element of the prime-order order-`r` target subgroup. GT is contained in the cyclotomic subgroup. |
| final exponentiation | The operation mapping a valid nonzero Miller-loop/general-Fp12 input into GT. |
| cyclotomic operation | An operation whose shortcut is valid only for cyclotomic inputs; it must not be documented as general Fp12 arithmetic. |
| GT operation | An operation whose algebra, scalar reduction, Frobenius identities, or conjugation-as-inversion requires a GT input. |

The current C and raw-Rust storage representation uses `blst_fp12` for general
Fp12, Miller-loop, and GT values. Storage compatibility does not widen an
operation's mathematical domain. Native C++ uses distinct `MillerResult` and
`GT` types; only `GT::final_exp` maps the former into the latter.

## Scalar domains

| Term | Meaning |
| --- | --- |
| scalar bitstring | An unsigned little-endian byte string with exactly the stated `nbits` low bits meaningful. Bits above `nbits` are ignored only when the function explicitly says so. |
| canonical scalar | An integer in `[0, r)`. A canonical secret key normally has the narrower protocol range `[1, r)`. |
| reduced scalar | The input integer reduced modulo `r`. |
| exact integer multiplication | Multiplies by the unsigned integer represented by the input bits without reducing modulo `r`. This distinction matters for non-subgroup points. |
| subgroup scalar multiplication | Interprets the scalar modulo `r`, because multiplication is defined on an order-`r` input group. |
| fixed public `nbits` | The bit count may affect control flow or loop bounds and therefore is public even when scalar bits are secret. |
| secret scalar | Scalar contents must not influence branches, memory addresses, or variable-time operations under the stated constant-time contract. |
| public scalar | Scalar contents may be used by a variable-time implementation. |

Unless an individual API says otherwise, supported scalar bitstrings contain at
most 256 meaningful bits, and a scalar pointer may be null only when `nbits` is
zero.

## Validation and trust-boundary names

| Name | Required meaning |
| --- | --- |
| checked | Validates the properties explicitly listed by the function in Release and reports failure. It does not imply that every conceivable property is checked. |
| unchecked | Omits one or more explicitly named validations. All remaining preconditions still apply. Documentation must say what is unchecked and what remains required. |
| raw | Expert representation-level operation. It performs only the validation explicitly listed and normally assumes its mathematical domain. |
| validate | Establishes the complete set of properties listed in its postcondition. A bare `validate` name should not be used for a relation-only test. |
| fast / optimized | Performance description only. It has no validation or mathematical-domain meaning. |
| safe Rust | Memory-safe to call under Rust's type and slice rules. It does not imply cryptographic validation or subgroup membership. |

An unchecked decoder and an unchecked arithmetic operation need not omit the
same checks. Each must name its actual retained preconditions and
postconditions.

CatBLST reserves new `checked` names for functions that perform the documented
validation in Release. The inherited `unchecked_mult` spelling does not imply
that every ordinary operation needs a checked/unchecked pair and will not be
extended to new API families.

## Timing language

| Term | Meaning |
| --- | --- |
| constant-time for fixed public sizes | Control flow and memory access do not depend on secret values when documented sizes, including `nbits`, are held fixed. |
| constant-time for canonical scalars | Constant-time for fixed public sizes when the scalar lies in `[0,r)`. If the contract accepts noncanonical values, it must state what timing distinction reduction introduces. |
| variable-time | Control flow or memory access may depend on the stated input values. Secret inputs must not be used. |
| public-data operation | Intended only for data that is already public, such as verification scalars or encoded points. |

Constant-time is an implementation property that requires code-generation and
platform evidence. It does not mean the operation takes identical wall-clock
time under every system condition.

## Failure and output publication

Every fallible contract must state one of these output behaviors:

- **unchanged on failure**: no caller-visible output or context mutation occurs;
- **defined failure value**: the documented identity, zero, null, or other value
  is written;
- **unspecified on failure**: the caller must discard the output.

No C ABI function may allow a C++ exception to cross the boundary. Legacy
`void` APIs must be allocation-free and nonthrowing for supported inputs.
Fallible APIs must translate failure to their documented result without
publishing a partially initialized output unless partial publication is itself
part of the contract.

## Standard function-contract form

New and revised public functions should document, in this order:

1. **Operation**: the mathematical result.
2. **Preconditions**: representation and mathematical domain, scalar bounds,
   nullability, aliasing, and required scratch.
3. **Validation**: properties actually checked in Release.
4. **Postconditions**: properties established or preserved on success.
5. **Timing**: secret/public inputs and constant-time or variable-time status.
6. **Failure**: result code and output/context publication behavior.

Example:

> Computes variable-time subgroup scalar multiplication. `p` must be an
> on-curve G1 element. `scalar` is a public unsigned little-endian bitstring
> with `nbits <= 256`; it may be null only when `nbits == 0`. The operation
> does not validate `p` in Release, reduces the scalar modulo `r`, and produces
> a G1 element. Exact output aliasing with `p` is supported; partial overlap is
> forbidden.

The contract should say **assumes subgroup membership**, not merely "does no
group check," because the former states the mathematical domain required for a
correct result.

## Language surfaces

The declaration on each supported language surface is the user-facing source
of truth. Audit tables track coverage and cross-language agreement; they do not
replace declaration documentation.

### C

Public C declarations use Doxygen comments. A family-level contract may state
shared representation, aliasing, timing, and nullability rules, but a function
must restate every exception to those rules. In particular:

- a `void` result does not imply that arbitrary inputs are accepted;
- validation performed only by a Debug assertion must be described as a
  precondition, not validation;
- a fallible function must say whether output is unchanged, defined, or
  unspecified on failure; and
- every supported call is allocation-free and nonthrowing unless its contract
  explicitly documents a translated failure result.

### Native C++

The selected `blst::` names are the main native C++ API. When a selected name
is an alias of a scalar or AVX implementation type, the implementation type's
Doxygen contract is authoritative and the selected header must identify that
relationship. Public methods use `@pre`, `@post`, `@return`, `@warning`, and
aliasing notes where applicable. A method is not documented as `noexcept`
unless its complete call chain has been reviewed as nonthrowing.

### Rust

Safe Rust methods document their mathematical domain even when violating that
domain cannot violate Rust memory safety. Use ordinary Rustdoc sections such as
`# Panics`, `# Errors`, and `# Timing`; reserve `# Safety` for obligations of an
`unsafe` API. Safe slice checks do not establish field, curve, subgroup, or
protocol validity.

The checked-in raw FFI declarations are generated artifacts. Their semantic
contract is the matching C declaration plus Rust's pointer/alignment rules.
They must not accumulate hand-maintained per-function promises that can diverge
from the C header; generation and inherited Rustdoc are reviewed together in
the Rust binding phase.

## Audit closure rule

For an externally reachable family, a file review cannot complete until:

- every supported surface implemented by the reviewed files has declaration
  documentation covering all applicable items in the standard form;
- equivalent C, C++, and Rust operations agree on their mathematical domain,
  validation, scalar interpretation, timing, and failure publication;
- intentionally absent or generated-only surfaces are recorded explicitly;
  and
- unresolved behavior is linked to a finding instead of being documented as a
  stable promise.

Documentation may be completed for an audited family without marking its much
larger containing public header complete. The audit manifest must distinguish
focused contract coverage from a complete file review.
