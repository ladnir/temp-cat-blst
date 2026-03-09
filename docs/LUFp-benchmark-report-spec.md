# LUFp Benchmark Report Spec

## Purpose

Each benchmark run must produce one Markdown report with a fixed structure.

The report must:

- present the same information in the same order every time;
- make it obvious when the run is not trustworthy;
- provide both raw timings and accounting-derived constituent timings;
- make the top-level conclusions easy to compare across runs.

This spec defines the required report structure and the required red-flag rules.

## Output Format

- Markdown only
- one report per run

Recommended filename:

- `out/benchmarks/<timestamp>-<commit>-report.md`

## Report Order

The report must contain these sections in this exact order:

1. `Red Flags`
2. `Run Metadata`
3. `Benchmark Coverage`
4. `Raw Results`
5. `Speed-of-Light Summary`
6. `Accounting Breakdown`
7. `Per-Operation Summary`
8. `Full Comparisons`
9. `Interpretation`
10. `Notes`

## 1. Red Flags

This section is mandatory and must appear first.

It must contain either:

- `None`

or

- a flat list of red flags

Each red flag line must include:

- severity: `fatal` or `warning`
- benchmark family
- short explanation
- suggested action

### Meaning

- `fatal`
  - the run must not be treated as trustworthy for conclusions
- `warning`
  - the run may still be useful, but interpretation must be cautious

### Mandatory red flag conditions

A run must emit a red flag if any of the following occur.

#### Accounting integrity

- accounting stage deltas do not sum to the full timing within tolerance
- accounting family is missing a required stage
- per-op averages are derived from missing or inconsistent op counts

#### Variance / stability

- coefficient of variation exceeds threshold
- repetition count is lower than required
- benchmark run time per case is below the configured minimum

#### Environment quality

- machine was not sufficiently idle
- build and benchmark overlapped
- benchmark command or environment metadata is missing

#### Correctness / validity

- correctness gate was disabled
- a known-incorrect path was benchmarked
- mocked or fake behavior was benchmarked
- compared variants were not on equal footing
  - different sink shape
  - different output storage behavior
  - different visibility/inlining footing

#### Coverage / drift

- expected benchmark family missing from the run
- benchmark registration drift
  - expected benchmark not present in binary
- code/report version mismatch

#### Speed-of-light integrity

- lower-bound estimate missing for a required family
- lower-bound estimate presented without explicit assumptions
- report states or implies a causal explanation for a gap without evidence

### Required thresholds

These thresholds must be applied consistently.

#### Accounting sum tolerance

For any accounting family:

- if `abs(sum(deltas) - full) / full > 3%`
  - emit `fatal`

#### Coefficient of variation

- if `cv > 2%`
  - emit `warning`
- if `cv > 5%`
  - emit `fatal`

#### Environment idleness

The report must record the load average before the run.

If the load average is above the agreed threshold for the benchmark host, emit:

- `warning` for mild excess
- `fatal` for significant excess

The exact host-specific idle threshold should be documented alongside the runner configuration.

#### Missing family

- any missing required benchmark family is `fatal`

#### Correctness disabled

- any run with correctness disabled is `fatal`
  - unless the report explicitly marks the run as exploratory

#### Mocked or fake benchmark

- any mocked/fake benchmark is `fatal`
  - unless the report explicitly marks the run as exploratory

#### Speed-of-light explanation discipline

- if a report presents speculation as fact when discussing the gap between lower bound and measurement
  - emit `fatal`
- if a report presents an unexplained gap without clearly marking it unexplained
  - emit `warning`

## 2. Run Metadata

This section is mandatory.

It must include:

- timestamp
- git commit
- dirty/clean status
- hostname
- server IP if remote
- CPU model
- load average before run
- benchmark binary path
- benchmark command line
- repetition count
- minimum benchmark time
- any CPU pinning / affinity settings

This section should be concise and stable in formatting.

## 3. Benchmark Coverage

This section lists every expected benchmark family and whether it ran.

Expected families should include, at minimum:

- `UFp add`
- `UFp add accounting`
- `UFp Fp6 mul`
- `UFp Fp6 mul accounting`
- `LUFp add`
- `LUFp add accounting`
- `LUFp renorm`
- `LUFp Fp6 mul`
- `LUFp Fp6 mul accounting`

For each family, report:

- `present`
- `missing`
- `skipped intentionally`

Any `missing` family must also appear in `Red Flags`.

## 4. Raw Results

This section contains one table of raw benchmark numbers.

Required columns:

- benchmark name
- mean
- median
- stddev
- cv

No interpretation belongs here.

This is the raw measurement record.

## 5. Speed-of-Light Summary

This section is mandatory for every important benchmark family.

For each family, report:

- lower-bound estimate
- assumptions used to derive that estimate
- measured timing
- ratio of measured to lower bound
- known causes of the gap
- unknown causes of the gap

The section must explicitly separate:

- `Facts`
- `Inferences`
- `Speculation`

Requirements:

- facts must be supported by measured data or directly known code structure
- inferences must be labeled as inferences
- speculation must be labeled as speculation
- if the gap is not explained, the report must say so explicitly

Recommended format per family:

- `Lower bound: ...`
- `Measured: ...`
- `Gap: ...`
- `Facts: ...`
- `Inferences: ...`
- `Speculation: ...`

This section is not a place to hand-wave away a discrepancy.

## 6. Accounting Breakdown

This is the authoritative constituent-timing section.

For each accounting family, report:

- cumulative stage time
- delta from previous stage
- percent of full

Example:

| Stage | Cumulative | Delta | % Full |
|---|---:|---:|---:|
| RawAdd | ... | ... | ... |
| RawAddBias | ... | Bias delta | ... |
| RawAddBiasCarry | ... | Carry delta | ... |
| Full | ... | Final repair delta | ... |

Requirements:

- stages must appear in declared order
- deltas must be computed from consecutive cumulative stages
- the final line must be the full benchmark
- the sum of deltas must match the full benchmark within tolerance

## 7. Per-Operation Summary

This section provides derived average costs for larger functions.

For each relevant family, report:

- op counts by kind
- average time per op kind
- formulas used

Example for `Fp6 mul`:

- narrow pre-add count
- `Fp2 mul` count
- wide add/sub count
- renorm count
- average time per `Fp2 mul`
- average time per wide add/sub
- average time per renorm

The report must explicitly state the op counts used.

If an op count is assumed rather than measured from topology, that must be stated.

## 8. Full Comparisons

This section reports direct top-level comparisons.

Examples:

- `UFp add` vs `LUFp add`
- `UFp Fp6 mul` vs `LUFp Fp6 mul`

For each comparison, report:

- baseline time
- comparison time
- absolute delta
- percent delta

This section should remain concise.

## 9. Interpretation

This section contains short conclusions only.

It should answer:

- what clearly improved
- what clearly regressed
- what remains inconclusive

This section must not hide red flags or missing coverage.

If the run has any `fatal` red flag, this section must begin with:

- `Run not trustworthy for final conclusions`

This section must also avoid converting speed-of-light gaps into unsupported explanations.

## 10. Notes

This section is for non-fatal contextual details, such as:

- unusual machine conditions
- intentionally skipped families
- follow-up experiments needed

This section should not duplicate `Red Flags`.

## Consistency Rules

Every report must:

- use the same section order
- use the same benchmark family names
- use the same op-count formulas
- explicitly mark missing families
- explicitly mark exploratory runs

Silent omission is not allowed.

## Decision Rule

When raw microbenchmarks and accounting benchmarks disagree:

- accounting benchmarks are the source of truth for constituent attribution

Raw microbenchmarks are supporting evidence only.

When a speed-of-light estimate and the measured result disagree:

- the measured result is the fact;
- the lower-bound estimate remains a model;
- the gap must be reported honestly;
- unsupported explanation is not allowed.
