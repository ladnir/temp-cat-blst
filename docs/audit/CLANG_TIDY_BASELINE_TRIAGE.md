# clang-tidy Baseline Triage

## Scope and configuration

The first report-only clang-tidy baseline ran at commit `3474eff1af99` with
Ubuntu LLVM clang-tidy 18.1.3. It used independent scalar and portable-AVX
Debug compilation databases for the retained `blst_static` production target.
Tests, benchmarks, frozen BTx/BTE, experimental CUDA, assembly, and removed
bindings were excluded. Portable-AVX Release was not configured or compiled.

The curated checker expression was:

```text
-*,bugprone-*,performance-*,portability-*,-bugprone-easily-swappable-parameters,-bugprone-reserved-identifier
```

The existing Clang Static Analyzer baseline already covers path-sensitive
analyzer checks. Identifier policy and easily-swappable-parameter checks were
excluded because they generate broad style/API churn without establishing a
defect. No fix-it was applied.

Ignored raw logs and fix-it YAML are in local run
`3474eff1af99-20260802T230337Z`. The reproducible runner is
`scripts/audit/run_clang_tidy_local.sh`.

## Result summary

Repeated diagnostics from shared headers were deduplicated by canonical file,
line, column, and rule. Scalar produced 47 unique diagnostics. Portable AVX
produced the same 47 plus 26 AVX-specific diagnostics, for 73 unique locations
overall.

| Rule | Unique | Disposition |
| --- | ---: | --- |
| `bugprone-implicit-widening-of-multiplication-result` | 34 | Non-defects. Operations use fixed small loop bounds, compile-time dimensions, checked sizes, or constants whose products fit before conversion. The scalar MSM branch containing a caller-sized product can execute only when `npoints < 32`; unsigned wrap before that later conjunct cannot select the branch. |
| `bugprone-macro-parentheses` | 14 | Current-use non-defects. Backend wide-arithmetic macros receive simple object addresses, `ROTR` receives fixed shift literals, and `TO_BYTES` receives integer literals. Parenthesization is optional hardening but not evidence of a current incorrect expansion. |
| `bugprone-narrowing-conversions` | 9 | Bounded conversions. Results are booleans, JSF digits in `{-1,0,1}`, a mask below 2^52, debug-print widths capped at 256, or iterator differences bounded by an existing object. |
| `performance-enum-size` | 8 | Declined. Public C enums retain their ABI representation; changing internal enums offers no demonstrated hot-path or storage benefit and can create layout churn. |
| `bugprone-assignment-in-if-condition` | 2 | Intentional assignments in pairing commit and SHA buffering. Both assigned values are consumed by the selected branch; the SHA expression intentionally initializes the cached offset while testing both conditions. |
| `bugprone-signed-char-misuse` | 2 | False positives. The GLV-JSF representation deliberately stores signed digits and promotes `-1`, `0`, or `1` to `int`. Converting through unsigned char would change the algorithm. |
| `performance-avoid-endl` | 2 | Debug-only bound-check diagnostics. They are outside arithmetic hot paths and flushing is useful immediately before a failed validation returns. |
| `performance-no-int-to-ptr` | 1 | Required alignment boundary. The AVX C adapter converts a scratch pointer to `uintptr_t`, rounds it up to the statically required alignment, and converts it back inside caller-provided storage. |
| `bugprone-too-small-loop-variable` | 1 | `AUD-039`. The template has no declared `table_log_size` bound, while an 8-bit one-based loop and booth index cannot safely represent a table size at or above 256. |

## Resolved finding

### AUD-039 - unconstrained AVX fixed-window table width

`PrecomputedTable<Curve, Ring, table_log_size>` stores its booth index and
gather loop variable in `uint8_t`. All current production instantiations use
the default `table_log_size == 4`, but the template accepts arbitrary explicit
values. At `table_log_size == 8`, the gather loop's inclusive upper bound is
256, the counter wraps after 255, and the booth result cannot represent 256.
Larger values also make the integer shift and table sizing unsafe.

The approved remediation:

- add a compile-time contract requiring `table_log_size > 0` and
  `table_log_size < 8`;
- keep the 8-bit booth representation and fixed loop structure unchanged;
- add compile-time acceptance/rejection coverage plus focused G1/G2 scalar
  multiplication differential tests; and
- verify optimized native AVX code generation is unchanged for the production
  width of four.

The contract is preferable to widening the hot loop because the compact booth
representation is intentional and no supported caller needs a 256-entry
fixed-window table.

The template now rejects widths outside 1 through 7 at compile time. Boundary
assertions cover accepted widths 1, 4, and 7 plus rejected widths 0 and 8.
Focused portable-AVX Debug tests exercised the production width-four G1 GLV
and G2 GLS entry points against the scalar backend and passed 2/2. An optimized
native AVX-512 Release probe compiled the same width-four G1/G2 entry points at
the pre-change and remediated revisions; after normalizing source-path and
line-marker metadata, all 102,762 assembly lines were identical. No benchmark
was needed because code generation did not change. `AUD-039` is validated.

## Residual limits

clang-tidy does not prove arithmetic bounds, constant-time behavior,
cryptographic correctness, or native AVX-512 instruction execution. The
bounded-conversion dispositions must be revisited during each owning file's
manual audit. This packet establishes native code generation but does not
replace Phase 16 native AVX-512 execution and broader scalar differential
evidence.
