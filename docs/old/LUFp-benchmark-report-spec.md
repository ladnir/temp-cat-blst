# LUFp Benchmark Report Spec

## Purpose

Each benchmark run must produce one Markdown report with a fixed structure.

The report must:

- make it obvious what run profile was executed;
- make it obvious when the run is not trustworthy;
- preserve raw timing data;
- report accounting-derived constituent timings for the families that were actually requested;
- keep canonical add-family performance guardrails visible.

## Output Format

- Markdown only
- one report per run

Recommended filename:

- `docs/report/<timestamp>-<commit>-report.md`

## Run Profiles

Every report must declare one profile:

- `canonical-add-only`
- `canonical-core`
- `canonical-full`
- `exploratory`

Required family coverage is profile-dependent.

Missing required families are `fatal`.

Families outside the chosen profile may be:

- absent;
- present as optional;
- skipped intentionally.

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

### Automated red flag conditions

The current report script is expected to flag the following conditions automatically.

#### Coverage / drift

- missing required family for the selected profile
- expected benchmark not present in the binary for the selected profile

#### Stability

- coefficient of variation exceeds threshold
- repetition count is lower than required
- benchmark time per case is below the configured minimum

#### Environment

- machine was not sufficiently idle
- benchmark command or metadata is missing from the report

#### Validity

- mocked or fake behavior was benchmarked in a canonical run

#### Accounting integrity

- required stage missing from an accounting family
- accounting deltas derived from inconsistent stage ordering

#### Canonical add guardrail

- canonical add-family benchmark exceeds the accepted `3%` baseline threshold
- canonical 1-step or 4-step repair canary delta exceeds the accepted `3%` threshold

### Manual review expectations

These checks still matter, but they are not fully enforced automatically by the current script and must be verified by the person reviewing the run:

- code/report version mismatch
- build and benchmark overlapped
- correctness gate was disabled for a canonical run
- compared variants were not on equal footing
- accounting full does not match a separate corresponding full benchmark within tolerance when both are intentionally included in the same profile

### Thresholds

- if `cv > 2%`: emit `warning`
- if `cv > 5%`: emit `fatal`
- if accounting/full mismatch exceeds `3%`: emit `fatal`
- if canonical add-family or repair canary exceeds `3%` from accepted baseline: emit `fatal`

## 2. Run Metadata

This section is mandatory.

It must include:

- timestamp
- run profile
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

## 3. Benchmark Coverage

This section lists the families relevant to the selected profile.

For each family, report one of:

- `present`
- `missing`
- `optional`
- `skipped intentionally`

Any `missing` required family must also appear in `Red Flags`.

### Required families by profile

#### `canonical-add-only`

- `UFp add accounting`
- `LUFp v2 add accounting`

#### `canonical-core`

- `UFp add accounting`
- `LUFp v2 add accounting`
- `LUFp v2 renorm accounting`

#### `canonical-full`

- `UFp add accounting`
- `LUFp v2 add accounting`
- `LUFp v2 renorm accounting`
- `UFp Fp6 mul accounting`
- canonical `LUFp v2 Fp6 mul accounting`

#### `exploratory`

- no implicit required families
- the report must explicitly list the families that were intentionally targeted

## 4. Raw Results

This section contains one table of raw benchmark numbers.

Required columns:

- benchmark name
- mean
- median
- stddev
- cv

No interpretation belongs here.

## 5. Speed-of-Light Summary

This section is mandatory for each important family present in the run.

For each reported family, include:

- lower-bound estimate
- assumptions
- measured timing
- ratio of measured to lower bound
- known causes of the gap
- unknown causes of the gap

The section must explicitly separate:

- `Facts`
- `Inferences`
- `Speculation`

Unsupported explanation is not allowed.

## 6. Accounting Breakdown

This is the authoritative constituent-timing section.

For each accounting family present in the report, include:

- cumulative stage time
- derived delta
- percent of full

### Add-family rule for alternative terminal stages

For canonical `LUFp v2` add and renorm families:

- `Full` is the canonical exact-add terminal stage
- `FullSign4` is an alternative terminal stage

Therefore:

- bias delta is measured from `RawAdd` or `Bias`
- carry delta is measured from the immediately preceding carry predecessor
- 1-step repair delta is measured from `RawAddBiasCarry` or `BiasCarry`
- 4-step repair delta is measured from the same shared predecessor, not from `Full`

The report must not present `FullSign4 - Full` as the repair cost.

## 7. Per-Operation Summary

This section provides derived average costs where the topology is known.

For each relevant family, report:

- op counts by kind
- average time per op kind
- formulas used

If an op count is assumed rather than directly known, say so explicitly.

## 8. Full Comparisons

This section reports direct top-level comparisons that are meaningful for the selected profile.

Examples:

- `UFp add` vs `LUFp v2 add`
- `UFp Fp6 mul` vs canonical `LUFp v2 Fp6 mul`

For each comparison, report:

- baseline time
- comparison time
- absolute delta
- percent delta

## 9. Interpretation

This section contains short conclusions only.

It should answer:

- what clearly improved
- what clearly regressed
- what remains inconclusive

If the run has any `fatal` red flag, this section must begin with:

- `Run not trustworthy for final conclusions`

## 10. Notes

This section is for non-fatal context, such as:

- unusual machine conditions
- intentionally skipped optional families
- follow-up experiments needed

## Consistency Rules

Every report must:

- use the same section order
- declare the run profile
- use stable family names
- mark missing required families explicitly
- mark exploratory runs explicitly

Silent omission is not allowed.

## Decision Rule

When raw microbenchmarks and accounting benchmarks disagree:

- accounting benchmarks are the source of truth for constituent attribution.

When a speed-of-light estimate and the measured result disagree:

- the measured result is the fact;
- the lower-bound estimate is a model;
- the gap must be reported honestly;
- unsupported explanation is not allowed.
