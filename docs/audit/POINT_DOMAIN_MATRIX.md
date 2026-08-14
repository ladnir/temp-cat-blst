# Point-Domain Contract Matrix

Date: 2026-08-10

State: policy applied and focused cross-language and native AVX validation
passing; `AUD-161` validated

Scope: supported scalar and AVX C APIs, selected native C++, safe Rust helpers,
and generated raw Rust FFI. CUDA and frozen BTx/BTE are excluded.

## Policy

Point-taking APIs use only two mathematical preconditions:

1. **on-curve**: the point satisfies the curve equation; subgroup membership
   is not required;
2. **subgroup**: the point is in G1 or G2 and is therefore also on-curve.

Validators have no prior point-domain precondition because their purpose is to
establish one of these properties. Infinity belongs to both domains; a protocol
may separately reject it where required. Pointer, length, aliasing, scalar, and
scratch requirements remain ordinary mechanical preconditions.

The wider on-curve contract is opt-in. It is promised only for the basic
primitive families below. Advanced algorithms remain subgroup-only even when a
current implementation happens to work for some non-subgroup inputs.

In the tables, `pN` expands independently to `p1` and `p2`, `gN` to `g1` and
`g2`, `GNPoint` to `G1Point` and `G2Point`, and `GNAffinePoint` to
`G1AffinePoint` and `G2AffinePoint`. Generated raw Rust declarations inherit
the matching C contract exactly.

`GROUP_API_CONTRACTS.csv` is the machine-checked exact C inventory: all 243
relevant scalar and AVX symbols match exactly one of its 53 disjoint families.
This matrix supplies the selected C++ and Rust mapping and the approved policy
used to migrate declarations.

## On-curve primitives and validators

| Family | C symbols | Selected/native C++ | Rust surface | Contract |
| --- | --- | --- | --- | --- |
| Decode | `blst_pN_deserialize`; `blst_pN_uncompress`; `blst_avx_pN_deserialize`; `blst_avx_pN_uncompress` | `GNAffinePoint::deserialize`; `GNAffinePoint::uncompress`; AVX equivalents | Raw FFI equivalents; high-level `PublicKey::deserialize`/`from_bytes`; `Signature::deserialize`/`from_bytes`; Serde decode | No point input. Success establishes on-curve only; it does not establish subgroup membership. Failure publication remains function-specific and must be documented. |
| Curve predicate | `blst_pN_on_curve`; `blst_pN_affine_on_curve` | `GNPoint::on_curve`; `GNAffinePoint::on_curve` | Raw FFI | Validator: no prior on-curve or subgroup precondition. Returns whether the input is on-curve. |
| Subgroup predicate | `blst_pN_in_gN`; `blst_pN_affine_in_gN`; `blst_avx_pN_in_gN` | `GNPoint::in_subgroup`; `GNAffinePoint::in_subgroup`; AVX equivalents | Raw FFI; `PublicKey::validate`; `Signature::validate`; `subgroup_check` helpers | Validator: no prior on-curve or subgroup precondition. Performs the curve check before the subgroup relation. Protocol validators may separately reject infinity. |
| Incomplete addition | `blst_pN_add`; `blst_pN_add_affine` | `GNPoint::add` projective/projective and projective/affine overloads | Raw FFI | Inputs are on-curve; subgroup membership is not required. Equal and opposite finite inputs remain outside the formula contract. |
| Complete basic arithmetic | `blst_pN_sub`; `blst_pN_add_or_double`; `blst_pN_sub_affine`; `blst_pN_add_or_double_affine`; `blst_pN_double` | `GNPoint::dadd`; `GNPoint::doubling`; negation plus complete addition for subtraction | Raw FFI | Inputs are on-curve; subgroup membership is not required. Output remains on-curve. If every input is in the subgroup, the output remains in it. |
| Conditional negation | `blst_pN_cneg`; `blst_pN_affine_cneg`; `blst_avx_pN_cneg` | `GNPoint::cneg`/`neg`; `GNAffinePoint::cneg`; AVX equivalents | Raw FFI | Input is on-curve; subgroup membership is not required. Output is the same point or its inverse. |
| Single-point representation conversion | `blst_pN_to_affine`; `blst_pN_from_affine`; `blst_pN_from_jacobian`; scalar/AVX `blst_avx_pN_from_scalar`, `blst_pN_from_avx`, and affine peers; AVX `to_affine`/`from_affine` | Corresponding `GNPoint`/`GNAffinePoint` conversion operations | Raw FFI | Input is on-curve; subgroup membership is not required. Conversion preserves the represented mathematical point and therefore preserves, but does not establish, subgroup membership. |
| Observation | `blst_pN_is_equal`; `blst_pN_affine_is_equal`; `blst_pN_is_inf`; `blst_pN_affine_is_inf`; AVX equivalents | `equal`; `is_inf`; AVX equivalents | Raw FFI; raw-type `PartialEq` | Inputs are on-curve; subgroup membership is not required. These are not validators. |
| Encoding | `blst_pN_serialize`; `blst_pN_compress`; `blst_pN_affine_serialize`; `blst_pN_affine_compress`; AVX equivalents | `serialize`; `compress`; AVX equivalents | Raw FFI; high-level `serialize`, `compress`, and `to_bytes` | Input is on-curve; subgroup membership is not required. Encoding does not validate or upgrade the input domain. |
| Exact scalar multiplication | `blst_p1_unchecked_mult`; `blst_p2_unchecked_mult` | `G1Point::unchecked_mult`; `G2Point::unchecked_mult` | Raw FFI; `blst_p1::unchecked_mult`; `blst_p2::unchecked_mult` | Input is on-curve; subgroup membership is not required. Applies the exact unsigned scalar without modulo-r subgroup semantics or GLV/GLS. These are the only inherited generic-multiplication exceptions. |
| Friendly batch subgroup validation | `blst_avx_pNs_in_gN` | Internal AVX batch relation and owned challenge stream | Safe and raw AVX Rust APIs | Validator: candidates have no prior point-domain precondition. The implementation checks every curve equation, obtains one fresh unpredictable seed after the batch is fixed, and owns 81 uniform ternary rows and the subgroup shift. Under the documented SHA-256 PRF and no-2-torsion assumptions, false acceptance is below 2^-128. |
| Expert batch subgroup relation | `blst_avx_pN_batch_subgroup_check` | Internal AVX batch relation implementation | Raw expert AVX FFI | Candidate points are on-curve and the shift is a subgroup point. Caller-provided challenge width and rows determine soundness. Zero width and all-zero rows deliberately remain possible; this is not the general validation boundary. |

## Subgroup-only families

| Family | C symbols | Selected/native C++ | Rust surface | Contract |
| --- | --- | --- | --- | --- |
| Single-point optimized multiplication | `blst_pN_mult`; `blst_pN_mult_vartime`; `blst_avx_pN_mult` | `GNPoint::mult`; `GNPoint::mult_vartime`; AVX equivalents | Raw FFI; `blst_pN::mult`; `blst_pN::mult_vartime` | Input is in G1/G2. Scalar reduction and GLV/GLS may be used. |
| Batch conversion and addition | `blst_pNs_to_affine`; `blst_pNs_add` | `points_to_affine`; `points_add` | Raw FFI; high-level aggregate helpers where applicable | Inputs are subgroup points. CatBLST deliberately does not promise the wider domain even when the implementation uses basic curve formulas. Users needing wider batching can compose the single-point primitives. |
| Fixed-table and MSM | Every `blst_pNs_mult_wbits*`, `blst_pNs_mult_pippenger*`, `blst_pNs_tile_pippenger`, `blst_pNs_mult_fr*`, and `blst_pNs_mult_batch_affine_fr*` entry | Scalar/selected wbits, Pippenger, and batch-affine implementations | Raw FFI and `MultiPoint` helpers | Inputs are subgroup points. Precomputation helpers share the domain of the evaluator that consumes their output. No generic MSM promise exists. |
| Expert AVX MSM | Every `blst_avx_pN_mult_pippenger*` and `blst_avx_pNs_mult_pippenger*` entry | AVX GLV/GLS and Pippenger implementations | Raw/expert AVX FFI | Inputs are subgroup points. No incidental full-curve behavior is promised. |
| Raw Miller and pairing | `blst_miller_loop`; `blst_miller_loop_n`; `blst_miller_loop_lines`; `blst_precompute_lines`; `blst_pairing_raw_aggregate`; `blst_avx_pairing` | Selected Miller-loop, precomputed-line, and raw pairing implementations | Raw FFI and expert pairing wrappers | Point inputs are subgroup points. Although underlying formulas may evaluate on a wider curve domain, CatBLST does not promise that advanced contract. |
| Protocol arithmetic | Pairing aggregation, signing, public-key derivation from points, aggregate operations, and verification entries | Selected sign/pairing implementations | `PublicKey`, `Signature`, aggregate, `Pairing`, and verification methods | Point operands are subgroup points unless the function explicitly validates them. A false group-check flag transfers that subgroup precondition to the caller. Protocol-specific infinity rules remain explicit per function. |

## Constructors and postconditions

Generators, hash-to-curve, encode-to-curve, and direct map-to-curve have no
point-input domain. Their successful outputs are subgroup points. A successful
subgroup validator likewise establishes the subgroup condition; ordinary
arithmetic only preserves a condition already supplied by the caller.

## Migration rule

The reviewed scalar/AVX C declarations, selected native C++ point declarations,
and safe Rust point/protocol declarations now state this contract. Generated raw
Rust inherits the C declaration. Future containing-file reviews must apply the
same vocabulary to any additional point-taking surface. A wider promise requires
focused non-subgroup evidence. A narrower subgroup promise needs no artificial
runtime check and remains valid even if the current implementation is
incidentally full-curve-correct.

The inherited AVX batch entries remain expert relation evaluators. Zero rounds
accept without inspecting candidates, and all-zero caller coefficient rows
ignore every candidate and test only the shift. General callers use
`blst_avx_p1s_in_g1` or `blst_avx_p2s_in_g2`; these entries own challenge
generation and cannot be configured into either bypass. Their probabilistic
bound requires a fresh unpredictable seed after the immutable batch is fixed.
