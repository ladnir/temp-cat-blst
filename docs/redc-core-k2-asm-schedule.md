# `redc_core` K=2 ASM Schedule (Promoted)

## Goal
Provide a spill-free, explicit Linux AVX512 schedule for two independent `UFpPack8x1::redc_core` streams.

## High-Level Structure

- Function: `blst_fp8x1_redc_core_k2`
- Location: `blst/fields/Fp_avx_redc_core_k2_asm.S`
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

Current register layout:
- stream limbs: `zmm0..zmm15`
- stream0 queue (6): `zmm16..zmm21`
- stream1 queue (6): `zmm22..zmm27`
- temps: `zmm28..zmm31`

This is the former benchmark `v3` schedule promoted to the library.

## Round Schedule

For each round:

1. Compute `u = low(pivot * montw) & bmask`
2. Apply `p0..p6` low/high updates to active window
3. Carry from pivot and fold `y0 << 4` into next limb
4. Load next wide limb and apply `p7` low/high
5. Eagerly fold `p7` high into the loaded next limb
6. Advance queue mapping by caller-side register rotation

Both streams run each round in sequence with fixed register mapping.

Current round strategy:
- compute `u0`/`u1` first,
- apply `p0` low first, then add:
  - carry from pivot,
  - old queue-head,
  - `p0` high directly into `a[r+1]` / `b[r+1]`,
- shift 6-slot queues,
- apply `p1..p6` with high terms accumulated into queue slots,
- apply `p7` with eager high-fold into loaded next limb.

Key change from old variants:
- for `p0` and `p7` high folds, stream 0 uses temp `T` and stream 1 uses register `C` after constant use in that sub-step.
- this removes the single-temp reuse chain without adding spill pressure.

## Constant-Feed Revision

Latest change in this schedule:

- Each `p[i]` is loaded once into a temp register and reused for both `vpmadd52lo`/`vpmadd52hi` of that stream.
- This replaces repeated memory-form IFMAs for each lo/hi pair.
- `montw`/`bmask` in the `u` derivation are back to memory-form operations (avoids extra `vmovdqa64` instructions).
- Final carry phase loads `bmask` once and reuses it for both 8-limb carry chains.
- Historical note: earlier benchmark variants (`v2`, `v3`) were A/B-only names; `v3` is now the promoted library symbol.

Observed effect (Linux host):

- Lower load-port pressure (`uops_dispatched.port_2_3_10` dropped materially)
- Lower memory-bound slots (`topdown.memory_bound_slots` dropped materially)
- Higher IPC
- Throughput improved to about `25.1M redc_core/s` for `BM_Fp8x1RedcCore_AsmK2`
- promoted kernel reaches about `26.3M..26.4M redc_core/s` in aggregate runs.

Recent A/B (Linux host, `--benchmark_min_time=1s`, 10 reps):

- `BM_Fp8x1RedcCore_AsmK2` (promoted kernel): `~26.34M/s`
- gain: about `+1.8%`

Pinned perf snapshot (`taskset -c 0`, 1 rep, 2s min time):

- old `v2`: `~25.70M/s`, `uops_executed.stall_cycles ~418.8M`, `topdown.memory_bound_slots ~1.55B`
- promoted kernel (old `v3`): `~25.98M/s`, `uops_executed.stall_cycles ~328.2M`, `topdown.memory_bound_slots ~1.20B`

Static instruction footprint (`objdump`, symbol body):

- old `v2`: total `771`, `vpmadd52` `272`, shift `90`, add `90`
- promoted kernel (old `v3`): total `773`, `vpmadd52` `272`, shift `90`, add `90`

## Remaining Limiter

The kernel is still backend-bound after spill removal and constant-feed cleanup:

- carry/shift/add dependency chains
- IFMA scheduling overlap limits

The promoted schedule confirms there is still headroom in dependency scheduling without changing IFMA count.
Further gains likely require deeper carry/queue scheduling changes or algorithmic changes above `redc_core`.
