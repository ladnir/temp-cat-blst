

# Plan 2 — High-ROI VROOM/CRNS competitiveness experiment (macro-kernel first)

### Goal

Answer: **Does CRNS/RNS beat your best limb+IFMA throughput approach at the operations you actually care about (Fp12 mul / cyclotomic sqr / pairing inner step)?**

The paper itself indicates Fp mul may be near parity while Fp12 and pairing see large wins, because CRNS amortizes reduction across tower arithmetic.
(We already saw this in the uploaded MIT paper: Fq mul ~1.02×, F12 mul ~5.27×, pairing ~3.29× vs blst in their setup.) 

So: don’t waste time trying to beat your Fp×8 kernel first. Beat it where it matters.

## 1. Define the minimal CRNS prototype scope

You need only enough CRNS machinery to implement:

* `Fp_mul` (as a building block, can be slower)
* `Fp2_mul` (optional if you jump directly)
* **Fp12_mul**
* **Fp12_cyclo_sqr**
* (optional) sparse mul used in Miller loop

But you can prototype faster if you do:

* **Fp2 and Fp12 directly in RNS** with batched base extension calls, and treat Fp mul as “internal”.

## 2. Pick a single CPU target and match the playing field

VROOM is very sensitive to:

* AVX-512IFMA availability,
* number of zmm regs,
* microarchitecture (ICL vs SPR).

Fix a target (e.g., SPR). Use identical:

* compiler flags,
* frequency governor,
* bench harness.

**Deliverable**

* A benchmark harness that runs:

  * your current SoA-limb macro-kernels
  * CRNS macro-kernels
    under identical batching, pinned cores, warmed caches.

## 3. Implement CRNS in the most “throughput-honest” way

CRNS is basically a structured matrix-vector multiply plus correction; it’s IFMA-friendly if:

* residues are stored SoA (one modulus per vector stream),
* constant tables are laid out contiguous and aligned,
* you minimize gathers and avoid shuffles.

**Do not** implement CRNS in an AoS “array of residues” way — you’ll lose to memory.

**Representation sketch**

* Choose `t` moduli (e.g., ~52-bit primes) per base.
* Store an element as:

  * `resM[j]` for j in 0..t-1 (SoA over j)
  * each `resM[j]` is a `__m512i` across 8 lanes
* Same for base N if using dual-base.

**Deliverables**

* `struct Rns8 { __m512i r[t]; };`
* `crns_extend_k` that base-extends K elements in flight (like Plan 1’s k-way ILP)

## 4. The decisive experiment: Fp12 mul and cyclotomic sqr throughput

Implement *one* Fp12 multiplication and one cyclotomic squaring in CRNS form, with:

* aggressive lazy reduction: delay CRNS correction as long as algebra allows
* fused subcomputations: compute all needed residue products, then do the minimum number of CRNS steps

Then compare against your best limb+IFMA macro-kernel:

* cycles per Fp12 mul
* cycles per cyclo sqr
* bytes loaded/stored per op (approx)

**Why this is high-ROI**
If CRNS doesn’t beat you here, it probably won’t win end-to-end.
If it does, it likely wins big in pairing and GT exponentiation.

## 5. If CRNS wins at Fp12: one “pairing inner step” macro-kernel

Implement the smallest pairing-relevant unit:

* “line evaluation * accumulator update”
  which includes:
* sparse multiplication into Fp12 accumulator
* a squaring (depending on step)
* minimal normalizations

Benchmark:

* throughput for B pairings batched (B large)
* compare to your existing batched pipeline

## 6. Interpretation guide: what outcomes mean

* **CRNS loses at Fp12**: stay with limb+IFMA, focus on Plan 1 fusion/lazy reduction and memory layout.
* **CRNS wins at Fp12 but loses at pairing inner step**: likely your sparse mul / data layout / correction frequency is wrong; fix fusion.
* **CRNS wins at pairing inner step**: it’s worth considering a serious CRNS backend for GT/pairing paths, possibly leaving G1 scalar ops on limb backend.
