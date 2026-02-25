# AVXBLS Porting Survey V2 (Logical Order, Full Detail)

## 1) Scope and Objective

This document is a complete, ordered survey of the in-tree `avxbls` implementation and how to port it into main `blst` in phases.

Repository locations surveyed:

- `C:\Users\peter\repo\blst\avxbls\src`
- `C:\Users\peter\repo\blst\avxbls\baseline`

Required porting phases:

1. `Fp`
2. `Fp2`
3. `Fp4`
4. `Fp6`
5. `Fp12`
6. `G1`
7. `G2`

This V2 includes:

- full file-level AVX deltas,
- complete vector-type contracts/invariants/layouts,
- full function behavior reference,
- pre/post conditions for all function families,
- operation-type inventories,
- phased integration plan with Fp-first execution strategy.

---

## 2) Source Delta (What AVXBLS Adds/Changes)

## 2.1 New files in `avxbls/src`

- `fp12_avx.h`
- `fp12_avx.c`
- `intrin.h`
- `intrin_compat.h`

## 2.2 Modified files vs `baseline/`

- `fields.h`
- `fp12_tower.c`
- `pairing.h`
- `pairing.c`
- `vect.h`
- `main.c` (test/benchmark expansion)

## 2.3 Integration pattern introduced by AVXBLS

- Scalar and vector implementations coexist in several modules.
- Header macros route generic symbols to vector versions (notably in `fields.h`, `pairing.h`).
- Wrappers in `fp12_tower.c` and `pairing.c` do scalar<->packed conversion around AVX kernels.

---

## 3) Core Data Model and Notation

## 3.1 SIMD register model

- `__m512i` is treated as 8 lanes x 64 bits.
- Emulation definition (`intrin_compat.h`): `struct { uint64_t lane[8]; }`.
- `VSET(X7,...,X0)` maps `X0 -> lane0`, `X7 -> lane7`.

## 3.2 Limb/radix model

From `fp12_avx.h`:

- `BRADIX = 48`
- `NWORDS = 8`
- `VWORDS = 4`
- `SWORDS = 6`
- `BMASK = 0xFFFFFFFFFFFFULL`
- `BALIGN = 4`

Meaning:

- AVX internal arithmetic is radix-`2^48`.
- Single-length Fp = 8 x 48-bit limbs.
- Scalar wrapper boundary still uses 6 x 64-bit limb form.

## 3.3 Layout naming convention

Examples:

- `FpPack8x1`: 8-way channels, 1 word per channel per vector slot.
- `FpPack4x2`: 4-way channels with lane-pair packing.
- `fp2_2x2x2w`: 2-way, each Fp2 component split across lane-pairs.

Dimension product maps to 8 logical lanes.

---

## 4) Complete Vector Type Contracts (Representation, Invariants, Layout)

## 4.1 Global invariants

- Field values are Montgomery-domain unless explicitly documented otherwise.
- Scalar limb order is little-endian by limb index.
- Packed radix-48 normalized limbs satisfy `< 2^48` after carry/normalization.
- "x2" types are widened accumulators (not reduced).

## 4.2 Scalar/base types

- `vec384`: one Fp value, 6x64 limbs.
- `vec768`: widened 12x64-limb container.
- `vec384x`: Fp2 pair `[c0, c1]`.
- `vec384fp6`: `[c0, c1, c2]` of Fp2.
- `vec384fp12`: `[d0, d1]` of Fp6.

## 4.3 Fp packed types

- `FpPack8x1` (`__m512i[8]`)
  - Represents 8 parallel Fp values.
  - For channel `ch`, limb `i`: `elem[ch].limb[i] = a[i].lane[ch]`.
- `UFpPack8x1` (`__m512i[16]`)
  - Widened product form of `FpPack8x1`.
- `FpPack4x2` (`__m512i[4]`)
  - 4 parallel Fp values with lane-pair packing.
  - Channel `ch` uses lanes `(2*ch, 2*ch+1)` for low/high limb halves.
- `UFpPack4x2` (`__m512i[12]`)
  - Widened product form of `FpPack4x2`.

## 4.4 Fp2 packed types

- `fp2_8x1x1w` (`FpPack8x1[2]`): two Fp planes (`c0`, `c1`) of 8-way packed values.
- `fp2x2_8x1x1w` (`UFpPack8x1[2]`): widened equivalent.
- `fp2_4x2x1w` (`__m512i[8]`): 4-way Fp2 in lane-pair format.
- `fp2x2_4x2x1w` (`__m512i[16]`): widened equivalent.
- `fp2_2x4x1w` (`__m512i[8]`): 2-way/4-lane-group intermediate format.
- `fp2x2_2x4x1w` (`__m512i[16]`): widened equivalent.
- `fp2_2x2x2w` (`__m512i[4]`): 2-way Fp2 line-kernel layout.
- `fp2x2_2x2x2w` (`__m512i[12]`): widened equivalent.
- `fp2_1x4x2w` (`__m512i[4]`): 1-way / 4 subgroup representation.
- `fp2x2_1x4x2w` (`__m512i[12]`): widened equivalent.

## 4.5 Fp4 packed types

- `fp4_2x2x2x1w` (`__m512i[8]`)
- `fp4x2_2x2x2x1w` (`__m512i[16]`)
- `fp4_1x2x2x2w` (`__m512i[4]`)
- `fp4x2_1x2x2x2w` (`__m512i[12]`)

Used by cyclotomic and higher-tower kernels.

## 4.6 Fp6/Fp12 packing strategy

- No dedicated `fp6_*` typedef family.
- Fp6/Fp12 vector kernels are composed from packed Fp2/Fp4 blocks plus lane transforms.
- Wrapper boundaries in `fp12_tower.c` and `pairing.c` explicitly build/tear down packings using `VSET` + conversion helpers.
- Important: these are schedule kernels, not always lane-independent kernels.
  - internal `perm_*`, `blend_*`, and `shuf_*` operations intentionally mix slots to reuse partial products.
  - therefore a naive "lane i equals scalar multiply i" test is often invalid unless the exact schedule mapping is documented.

---

## 5) Function Catalog in Dependency Order (Behavior + Preconditions + Postconditions)

Legend:

- `Pre-L`: input layout is valid for the function's packed/scalar type.
- `Pre-M`: inputs are Montgomery-domain field values.
- `Pre-N`: single-length radix-48 normalized inputs.
- `Pre-W`: widened (unreduced) inputs.
- `Post-L`: output respects declared layout.
- `Post-R`: reduced field output (Montgomery).
- `Post-N`: normalized radix-48 limbs (or after stated carry helper).
- `Post-W`: widened unreduced output.
- `Post-P`: permutation/blend only.
- `Post-G2`: updates G2 state and returns corresponding line coefficients.

## 5.1 Header utility and conversion functions (`fp12_avx.h`, static)

- `blend_0x0F`: blends high lanes from `a` and low lanes from `b`.
  - Pre: `Pre-L`; Post: `Post-L`, `Post-P`.
- `blend_0x33`: fixed blend mask `0x33`.
  - Pre: `Pre-L`; Post: `Post-L`, `Post-P`.
- `blend_0xC0`: fixed blend mask `0xC0`.
  - Pre: `Pre-L`; Post: `Post-L`, `Post-P`.
- `blend_0xC3`: fixed blend mask `0xC3`.
  - Pre: `Pre-L`; Post: `Post-L`, `Post-P`.
- `perm_var`: variable lane permutation over 8-lane vectors.
  - Pre: `Pre-L`; Post: `Post-L`, `Post-P`.
- `blend_0x0F_hl`: half-length blend variant.
  - Pre: `Pre-L`; Post: `Post-L`, `Post-P`.
- `perm_var_hl`: half-length variable permutation variant.
  - Pre: `Pre-L`; Post: `Post-L`, `Post-P`.
- `conv_64to48_fp`: one scalar Fp conversion 64-bit-limb form -> radix-48.
  - Pre: scalar little-endian 64-bit limbs; Post: scalar radix-48 limbs.
- `FpPack8x1::conv_64to48`: packed conversion (8x1w) to radix-48.
  - Pre: `Pre-L`; Post: `Post-L`, `Post-N`.
- `trace_dump__m512i`: debug register dump.
  - Pre: any `__m512i`; Post: no state change.
- `conv_64to48_FpPack4x2`: packed conversion (4x2w) to radix-48.
  - Pre: `Pre-L`; Post: `Post-L`, `Post-N`.
- `conv_64to48_mpi`: generic multiprecision 64->48 conversion.
  - Pre: valid buffer sizes; Post: converted limb stream.
- `conv_48to64_fp`: one scalar Fp conversion radix-48 -> 64-bit-limb form.
  - Pre: scalar radix-48 limbs; Post: scalar 64-bit-limb form.
- `FpPack8x1::conv_48to64`: packed conversion (8x1w) back to 64-bit-limb packing.
  - Pre: `Pre-L`; Post: `Post-L`.
- `conv_48to64_FpPack4x2`: packed conversion (4x2w) back to 64-bit-limb packing.
  - Pre: `Pre-L`; Post: `Post-L`.
- `conv_48to64_mpi`: generic multiprecision 48->64 conversion.
  - Pre: valid buffers; Post: converted limb stream.
- `carryp_mpi48`: scalar radix-48 forward carry propagation.
  - Pre: radix-48 style array; Post: normalized carry.
- `carryp_FpPack4x2`: packed radix-48 carry propagation for 4x2w.
  - Pre: `Pre-L`; Post: `Post-L`, `Post-N`.
- `get_channel_8x1w`: extracts one lane-channel from `FpPack8x1`.
  - Pre: channel index in range; Post: scalar channel limbs.
- `get_channel_4x2w`: extracts one packed channel from `FpPack4x2`.
  - Pre: channel index in range; Post: scalar channel limbs.

## 5.2 SIMD transform helpers (`fp12_avx.c`, static)

### 5.2.1 Single-length shuffle/permute/blend

- `shuf_00`, `shuf_01`, `shuf_11`, `shuf_z1`
- `perm_0123`, `perm_0300`, `perm_1010`, `perm_1032`, `perm_1200`, `perm_10zz`, `perm_2000`, `perm_2332`, `perm_3322`, `perm_zz32`
- `blend_0x03`, `blend_0x11`, `blend_0x30`, `blend_0x55`, `blend_0x77`
- `shuf_01_fp2_8x1x1w`, `perm_0123_fp2_8x1x1w`, `perm_0300_fp2_8x1x1w`, `perm_1032_fp2_8x1x1w`, `perm_1200_fp2_8x1x1w`, `perm_2000_fp2_8x1x1w`, `blend_0x11_fp2_8x1x1w`, `blend_0x33_fp2_8x1x1w`, `blend_0x55_fp2_8x1x1w`, `blend_0x77_fp2_8x1x1w`

Behavior: deterministic lane rearrangement, duplication, zero-masked selection, and mixing used to set up arithmetic schedules.

Pre/Post for all above: `Pre-L` -> `Post-L`, `Post-P`.

### 5.2.2 Double-length shuffle/permute/blend

- `shuf_00_dl`, `shuf_01_dl`, `shuf_11_dl`
- `perm_1010_dl`, `perm_1032_dl`, `perm_1100_dl`, `perm_10zz_dl`, `perm_2323_dl`, `perm_3000_dl`, `perm_3232_dl`, `perm_zz32_dl`, `perm_var_dl`
- `blend_0x01_dl`, `blend_0x03_dl`, `blend_0x0F_dl`, `blend_0x33_dl`, `blend_0x55_dl`, `blend_0xC0_dl`
- `shuf_00_fp2x2_8x1x1w`, `shuf_01_fp2x2_8x1x1w`, `shuf_11_fp2x2_8x1x1w`
- `perm_1032_fp2x2_8x1x1w`, `perm_2323_fp2x2_8x1x1w`, `perm_3000_fp2x2_8x1x1w`, `perm_var_fp2x2_8x1x1w`
- `blend_0x01_fp2x2_8x1x1w`, `blend_0x03_fp2x2_8x1x1w`, `blend_0x33_fp2x2_8x1x1w`, `blend_0x55_fp2x2_8x1x1w`

Behavior: same as single-length family but for widened (`x2`) vector stacks and widened Fp2 stacks.

Pre/Post for all above: `Pre-L` -> `Post-L`, `Post-P`.

### 5.2.3 Half-lane and layout-bridge transforms

- `perm_1010_hl`, `perm_1032_hl`, `perm_3232_hl`, `perm_zz32_hl`, `blend_0x55_hl`
- `conv_sltohl`, `perm_var_vl`, `blend_0x0F_vl`, `conv_dltovl`

Behavior: repacks between internal layout families used by line kernels, cyclotomic kernels, and multiply schedules.

Pre/Post: `Pre-L` -> `Post-L`, `Post-P`.

## 5.3 Fp kernels and helpers

### 5.3.1 Static Fp/Fpx2 helpers

- `FpPack8x1::add`, `FpPack8x1::sub`, `FpPack8x1::asx4`, `FpPack8x1::asax6`, `FpPack8x1::assx6`
- `add_FpPack4x2`, `sub_FpPack4x2`, `asx2_FpPack4x2`
- `add_UFpPack8x1`, `sub_UFpPack8x1`, `aaasx2_UFpPack8x1`, `aassx2_UFpPack8x1`, `asx4_UFpPack8x1`, `ssaax2_UFpPack8x1`, `ssasx2_UFpPack8x1`
- `add_UFpPack4x2`, `sub_UFpPack4x2`, `asx2_UFpPack4x2`, `aass_UFpPack4x2`

Behavior:

- `add/sub` families perform modular add/sub in packed layouts.
- `as*`/`aa*`/`ss*` families are schedule-specific linear-combination helpers used by multiply/reduction kernels.

Pre/Post:

- Single-length helpers: `Pre-L`, `Pre-M` (typically `Pre-N`) -> `Post-L`, `Post-N`.
- Widened helpers: `Pre-L`, `Pre-W` -> `Post-L`, `Post-W`.

### 5.3.2 Public Fp kernels

- `mul_UFpPack8x1` (blst alias; implementation is AVXBLS `mul_UFpPack8x1_v2`)
- `mul_UFpPack8x1_v2`
- `mul_UFpPack4x2` (blst alias; implementation is AVXBLS `mul_UFpPack4x2_v1`)
- `mul_UFpPack4x2_v1`
- `redc_UFpPack8x1`
- `redc_UFpPack4x2`

Behavior:

- `mul_fpx2_*`: widened Montgomery-compatible multiplication in packed radix-48 domain.
- `redc_fpx2_*`: Montgomery reduction from widened product back to single-length packed Fp.

Pre/Post:

- `mul_fpx2_*`: `Pre-L`, `Pre-M`, `Pre-N` -> `Post-L`, `Post-W`.
- `redc_fpx2_*`: `Pre-L`, `Pre-W` -> `Post-L`, `Post-R` (caller may carry-normalize where required).

## 5.4 Fp2 kernels and helpers

### 5.4.1 Static Fp2/Fp2x2 helpers

- `add_fp2_8x1x1w`
- `add_fp2_4x2x1w`, `sub_fp2_4x2x1w`, `aass_fp2_4x2x1w`, `assa_fp2_4x2x1w`, `sass_fp2_4x2x1w`, `ssaa_fp2_4x2x1w`, `mul_by_u_plus_1_fp2_4x2x1w`
- `mul_fp2_2x4x1w`
- `add_fp2_2x2x2w`, `sub_fp2_2x2x2w`, `as_fp2_2x2x2w`, `sa_fp2_2x2x2w`, `mul_by_u_plus_1_fp2_2x2x2w`
- `mul_fp2_1x4x2w`
- `add_fp2x2_8x1x1w`, `sub_fp2x2_8x1x1w`, `asx4_fp2x2_8x1x1w`, `mul_by_u_plus_1_fp2x2_8x1x1w`
- `add_fp2x2_4x2x1w`, `sub_fp2x2_4x2x1w`, `sax2_fp2x2_4x2x1w`, `mul_by_u_plus_1_fp2x2_4x2x1w`
- `redc_fp2x2_4x2x1w`
- `redc_fp2x2_2x2x2w`
- `add_fp2x2_2x2x2w`, `sub_fp2x2_2x2x2w`, `mul_by_u_plus_1_fp2x2_2x2x2w`

Behavior:

- Arithmetic helpers perform add/sub and schedule combinations in each packing.
- `mul_by_u_plus_1_*` multiplies by Fp2 non-residue `(u+1)`.
- `redc_fp2x2_*` reduces widened Fp2x2 to single-length Fp2.

Pre/Post:

- single-length helpers: `Pre-L`, `Pre-M` -> `Post-L`, `Post-R`/`Post-N`.
- widened helpers: `Pre-L`, `Pre-W` -> `Post-L`, `Post-W`.
- reduction helpers: `Pre-L`, `Pre-W` -> `Post-L`, `Post-R`.

### 5.4.2 Public Fp2/Fp2x2 kernels

- `mul_fp2_4x2x1w`, `sqr_fp2_4x2x1w`
- `mul_fp2_2x2x2w`, `sqr_fp2_2x2x2w`
- `mul_fp2x2_8x1x1w`
- `mul_fp2x2_4x2x1w`, `sqr_fp2x2_4x2x1w`
- `mul_fp2x2_2x4x1w`
- `mul_fp2x2_2x2x2w`, `sqr_fp2x2_2x2x2w`
- `mul_fp2x2_1x4x2w`

Behavior:

- `mul_fp2_*`/`sqr_fp2_*`: reduced Fp2 multiply/square in declared layout.
- `mul_fp2x2_*`/`sqr_fp2x2_*`: widened Fp2 multiply/square in declared layout.

Pre/Post:

- reduced kernels: `Pre-L`, `Pre-M`, `Pre-N` -> `Post-L`, `Post-R`.
- widened kernels: `Pre-L`, `Pre-M` -> `Post-L`, `Post-W`.

## 5.5 Fp4 kernels

- `sqr_fp4_2x2x2x1w_v1`
- `sqr_fp4_1x2x2x2w_v1`
- `sqr_fp4_2x2x2x1w_v2` (static)
- `sqr_fp4_2x2x2x1w_v3` (static)
- `sqr_fp4_1x2x2x2w_v2` (static)

Behavior:

- Fp4 square in different packed schedules/layouts; `v1` variants are selected primary APIs.

Pre/Post:

- `Pre-L`, `Pre-M` -> `Post-L`, `Post-R`.

## 5.6 Fp6 kernels

- `mul_fp6x2_4x2x1x1w` (static)
- `mul_fp6x2_2x2x2x1w`
- `mul_fp6x2_1x2x2x2w`
- `mul_by_xy00z0_fp6x2_2x2x2x1w` (static)
- `mul_by_xy00z0_fp6x2_2x4x1x1w` (static)

Behavior:

- widened Fp6 multiplication in multiple packings,
- sparse `xy00z0` specialization for line-function multiplies.

Pre/Post:

- `Pre-L`, `Pre-M` -> widened packed outputs (`Post-W`).
- additional contract note:
  - output slot positions are schedule-dependent for `4x2x1x1w`, `2x2x2x1w`, and sparse `xy00z0` kernels.
  - treat outputs as schedule blocks unless a slot map is explicitly defined.
  - recommended validation during bring-up is scalar differential testing with broadcast/degenerate packed inputs; add strict slot-map tests once mapping is derived.

## 5.7 Fp12 and line kernels

- `cyclotomic_sqr_fp12_vec_v1`
- `cyclotomic_sqr_fp12_vec_v2`
- `compressed_cyclotomic_sqr_fp12_vec_v1`
- `mul_fp12_vec_v1`
- `mul_fp12_vec_v2`
- `mul_fp12_vec_v3`
- `mul_fp12_vec_v4`
- `mul_by_xy00z0_fp12_vec_v1`
- `mul_by_xy00z0_fp12_vec_v2`
- `sqr_fp12_vec_v1`
- `line_dbl_vec_v1`
- `line_dbl_vec_v2`
- `line_add_vec_v1`
- `line_by_Px2_4x2x1w`
- `line_by_Px2_2x2x2w`
- `trace_dump_fp2_2x2x2w_state`

Behavior:

- Fp12 multiply/square/cyclotomic operations in packed layouts.
- Sparse `xy00z0` multiply for Miller line multiplication.
- G2 line add/double kernels and line scaling by `Px2`.

Pre/Post:

- Fp12 kernels: `Pre-L`, `Pre-M` -> `Post-L`, reduced/split outputs as declared.
- compressed cyclotomic: updates compressed BC-part representation.
- line kernels: packed G2 input state -> updated state and line outputs (`Post-G2`).
- `line_add_vec_v1` specific preconditions:
  - packed tuples must encode valid on-curve mixed-add operands `R=(X1,Y1,Z1)` and affine `Q=(X2,Y2)`,
  - kernel is the non-exceptional `madd-2007-bl` schedule and requires `H = X2*Z1^2 - X1 != 0`,
  - synthetic random-limb tuples are outside this contract and are not valid scalar-differential oracle inputs.
- packing contract note:
  - this kernel expects `X1Y1 = X1 | Y1` (the opposite order from some other schedule helpers),
  - with that packing, strict scalar differential checks are valid for `line_add_vec_v1`.

## 5.8 Scalar/vector wrappers in `fp12_tower.c`

- `mul_fp2x2`, `sqr_fp2x2`, `sqr_fp4`, `mul_fp6x2`
- `mul_fp12_scalar`, `mul_fp12_vector`
- `mul_by_xy00z0_fp12_scalar`, `mul_by_xy00z0_fp12_vector`
- `sqr_fp12_scalar`, `sqr_fp12_vector`
- `cyclotomic_sqr_fp12_scalar`, `cyclotomic_sqr_fp12_vector`
- `compressed_cyclotomic_sqr_fp12_scalar`
- `conjugate_fp12`, `inverse_fp12`, `frobenius_map_fp12`
- `back_cyclotomic_sim_fp12`

Behavior:

- Scalar functions implement tower arithmetic directly over `vec384*` structures.
- Vector versions pack scalar inputs -> call AVX kernels -> unpack outputs.
- Final-exponent support functions (`conjugate/inverse/frobenius`) remain scalar flow control around chosen multiply/square backend.

Pre/Post:

- Scalar wrappers: scalar Montgomery-domain structs in/out.
- Vector wrappers: scalar boundary inputs/outputs, with internal pack/unpack invariants.

## 5.9 Pairing flow functions in `pairing.c`

- `line_add_scalar`, `line_add_vector`
- `line_dbl_scalar`, `line_dbl_vector`
- `line_by_Px2`
- `start_dbl_n`
- `add_n_dbl_n_scalar`, `add_n_dbl_n_vector`
- `miller_loop_n`
- `mul_n_sqr_scalar`, `mul_n_sqr_vector`
- `raise_to_z_div_by_2`, `raise_to_z`
- `final_exp`
- `optimal_ate_pairing`

Behavior:

- Implements Miller loop and final exponentiation.
- Scalar and vector line-eval/mul-sqr blocks coexist; macros in `pairing.h` select vector path by default.

Pre/Post:

- Valid pairing inputs (`Q`, `P`) and intermediate projective state.
- Produces pairing value in `vec384fp12`.

---

## 6) Operation-Type Inventory (Cross-Check View)

## 6.1 Fp

Types:

- `vec384`, `vec768`
- `FpPack8x1`, `UFpPack8x1`, `FpPack4x2`, `UFpPack4x2`

Primary kernels:

- `mul_UFpPack8x1` / `mul_UFpPack8x1_v2` (AVXBLS v2 kernel)
- `mul_UFpPack4x2_v1`
- `redc_UFpPack8x1`, `redc_UFpPack4x2`

## 6.2 Fp2

Types:

- `vec384x`, `vec384fp2`
- `fp2_8x1x1w`, `fp2x2_8x1x1w`
- `fp2_4x2x1w`, `fp2x2_4x2x1w`
- `fp2_2x4x1w`, `fp2x2_2x4x1w`
- `fp2_2x2x2w`, `fp2x2_2x2x2w`
- `fp2_1x4x2w`, `fp2x2_1x4x2w`

Primary kernels:

- `mul_fp2_*`, `sqr_fp2_*`
- `mul_fp2x2_*`, `sqr_fp2x2_*`

## 6.3 Fp4

Types:

- `vec384fp4`
- `fp4_2x2x2x1w`, `fp4x2_2x2x2x1w`
- `fp4_1x2x2x2w`, `fp4x2_1x2x2x2w`

Primary kernels:

- `sqr_fp4_2x2x2x1w_v1`
- `sqr_fp4_1x2x2x2w_v1`

## 6.4 Fp6

Types:

- `vec384fp6`, `vec768fp6`
- composed packed Fp2/Fp4 blocks

Primary kernels:

- `mul_fp6x2_2x2x2x1w`
- `mul_fp6x2_1x2x2x2w`

## 6.5 Fp12

Types:

- `vec384fp12`
- packed split forms using `fp2_*` + `fp4_*`

Primary kernels:

- `mul_fp12_vec_v1..v4`
- `sqr_fp12_vec_v1`
- `cyclotomic_sqr_fp12_vec_v1/v2`
- `compressed_cyclotomic_sqr_fp12_vec_v1`
- `mul_by_xy00z0_fp12_vec_v1/v2`

Wrappers:

- `mul_fp12_vector`, `sqr_fp12_vector`, `cyclotomic_sqr_fp12_vector`, `mul_by_xy00z0_fp12_vector`

## 6.6 G1

Status:

- no dedicated AVX kernel additions found.
- G1 remains scalar/backend path.

## 6.7 G2

Packed working types:

- `fp2_2x2x2w`, `fp2_4x2x1w`

Primary kernels:

- `line_dbl_vec_v1`, `line_dbl_vec_v2`
- `line_add_vec_v1`
- `line_by_Px2_4x2x1w`, `line_by_Px2_2x2x2w`

Pairing wrappers:

- `line_dbl_vector`, `line_add_vector`, `add_n_dbl_n_vector`, `mul_n_sqr_vector`

---

## 7) File-Level Integration Notes

## 7.1 `fields.h`

- Adds scalar/vector split declarations for Fp12 operations.
- Macro-routing:
  - `sqr_fp12 -> sqr_fp12_vector`
  - `cyclotomic_sqr_fp12 -> cyclotomic_sqr_fp12_vector`
  - `mul_fp12 -> mul_fp12_vector`
  - `mul_by_xy00z0_fp12 -> mul_by_xy00z0_fp12_vector`
- Exposes helper types for benchmarking (`vec768x`, `vec768fp6`, `vec384fp4`) and helper APIs under `BENCHMARK`.

## 7.2 `pairing.h`

- Includes `fp12_avx.h`.
- Macro-routing:
  - `line_dbl -> line_dbl_vector`
  - `line_add -> line_add_vector`
  - `add_n_dbl_n -> add_n_dbl_n_vector`
  - `mul_n_sqr -> mul_n_sqr_vector`

## 7.3 `fp12_tower.c`

- Keeps scalar baseline functions and introduces vector wrappers.
- Uses explicit pack/unpack lane maps around AVX kernels.

## 7.4 `pairing.c`

- Keeps scalar line formulas and introduces vectorized wrappers and block functions.
- Miller loop uses macro-selected path (vector by default with current header mapping).

---

## 8) Porting Plan (Phase Order and Dependencies)

## 8.1 Phase 0: Infrastructure

- Port intrinsics abstraction layer (`intrin.h` plus compatibility strategy).
- Port constants/macros and base conversion helpers from `fp12_avx.h`.
- Validate 64<->48 pack/unpack and carry normalization in isolation.

## 8.2 Phase 1: Fp (first concrete target)

Port now:

- Types: `FpPack8x1`, `UFpPack8x1`, `FpPack4x2`, `UFpPack4x2`
- Kernels: `mul_UFpPack8x1`/`mul_UFpPack8x1_v2` (AVXBLS v2), `mul_UFpPack4x2`/`mul_UFpPack4x2_v1`
- Reductions: `redc_UFpPack8x1`, `redc_UFpPack4x2`
- Required helpers: Fp/Fpx2 add-sub/combination helpers + needed lane transforms + conversions.

Test gates:

- `64->48->64` round-trip tests for both `8x1w` and `4x2w`.
- widened mul consistency:
  - AVX widened -> scalar reference,
  - `8x1w` and `4x2w` path equivalence after canonical unpack.
- reduction correctness vs scalar Montgomery reduction.
- randomized and edge vectors near `p`, `0`, carry-heavy patterns.

## 8.3 Phase 2: Fp2

- Port reduced and widened Fp2 mul/sqr families.
- Validate across all active layouts (`4x2x1w`, `2x2x2w`, `8x1x1w`).

## 8.4 Phase 3: Fp4

- Port selected Fp4 square kernels.
- Validate against scalar Fp4 square.

## 8.5 Phase 4: Fp6

- Port `mul_fp6x2_*` and sparse helper dependencies.
- Validate vs scalar Fp6 multiplication.

## 8.6 Phase 5: Fp12

- Port `mul_fp12_vec_*`, `sqr_fp12_vec_v1`, cyclotomic/compressed kernels.
- Enable wrapper-level vector path only after kernel parity is proven.

## 8.7 Phase 6: G1

- No AVX kernels expected; verify no regressions from integration changes.

## 8.8 Phase 7: G2/Pairing

- Port line kernels and pairing block wrappers.
- Validate full Miller + final exp parity for single and batched cases.

---

## 9) Practical Integration Policy

- Keep scalar and AVX paths side-by-side during bring-up.
- Do not switch generic macro-routing to AVX globally before per-phase tests pass.
- Preserve bisectability: one phase, one test gate, then enable next phase.
- Add phase-specific tests in `test/` aligned with existing BLST style.
