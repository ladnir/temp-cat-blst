# `redc_core` K=2 ASM Schedule (Current)

## Goal
Provide a spill-free, explicit Linux AVX512 schedule for two independent `UFpPack8x1::redc_core` streams.

## High-Level Structure

- Function: `blst_redc_core_k2_asm`
- Location: `benchmark/redc_core_k2_asm.S`
- Inputs:
  - `out0`, `out1` (`LFpPack8x1*`)
  - `in0`, `in1` (`UFpPack8x1*`)
- Arithmetic model:
  - Montgomery reduction over radix-48 limbs
  - lazy high queue per stream
  - final carry-only normalization to 8 limbs

## Register Layout

- Stream 0 active limbs: `zmm0..zmm7`
- Stream 1 active limbs: `zmm8..zmm15`
- Stream 0 lazy queue (7 slots): `zmm16..zmm22`
- Stream 1 lazy queue (7 slots): `zmm23..zmm29`
- Temporaries: `zmm30`, `zmm31`

No stack spills are used in the kernel body.

## Round Schedule

For each round:

1. Compute `u = low(pivot * montw) & bmask`
2. Apply `p0..p6` low/high updates to active window
3. Carry from pivot and fold `y0 << 4` into next limb
4. Load next wide limb and apply `p7` low/high
5. Eagerly fold `p7` high into the loaded next limb
6. Advance queue mapping by caller-side register rotation

Both streams run each round in sequence with fixed register mapping.

## Constant-Feed Revision

Latest change in this schedule:

- Each `p[i]` is loaded once into a temp register and reused for both `vpmadd52lo`/`vpmadd52hi` of that stream.
- This replaces repeated memory-form IFMAs for each lo/hi pair.

Observed effect (Linux host):

- Lower load-port pressure (`uops_dispatched.port_2_3_10` dropped materially)
- Lower memory-bound slots (`topdown.memory_bound_slots` dropped materially)
- Higher IPC
- Throughput improved to about `23.2M .. 23.3M redc_core/s` for `BM_Fp8x1RedcCore_AsmK2`

## Remaining Limiter

The kernel is still backend-bound after spill removal and constant-feed cleanup:

- carry/shift/add dependency chains
- IFMA scheduling overlap limits

So further gains need schedule-level arithmetic improvements, not just spill removal.
