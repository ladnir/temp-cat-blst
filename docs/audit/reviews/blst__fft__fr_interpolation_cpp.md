# File Review: `blst/fft/fr_interpolation.cpp`

## Record

| Field | Value |
| --- | --- |
| Primary path | `blst/fft/fr_interpolation.cpp` |
| Interface path | `blst/fft/fr_interpolation.h` |
| Roadmap phase | Phase 17 Packet 17F |
| Implementation commit | `1821c58b80904f07cabaae9b9482d495f344256c` |
| Final documented interface commit | `baf3e456cd695eae65cc66cd85a0a17e91b4d203` |
| Reviewer | Codex audit team |
| Date | 2026-08-11 |
| Disposition | `in_scope_production` |
| Profiles | scalar Fr; scalar-default builds may also compile portable AVX support |
| Upstream reference | `not_applicable: local CatBLST interpolation implementation` |

## Purpose and trust boundary

The module owns arbitrary-point Fr interpolation and a specialized
allocation-free coefficient-suffix interpolation for samples selected from a
root-of-unity domain. It is reached by native C++, C, and safe Rust entry points.
Points, indices, lengths, and requested ranges are public scheduling inputs.
Values and recovered coefficients may be secret-derived.

## Mathematical contracts

- General interpolation accepts n pairwise-distinct Fr points and n values and
  returns the unique degree-<n polynomial in monomial order.
- Its product tree represents `M(X) = product_i (X - x_i)`. Denominators are
  `M'(x_i)`, and fixed-width reconstruction sums the Lagrange numerators.
- FFT-domain suffix interpolation accepts unique domain indices and evaluates
  the high coefficients of `M(X)Q(X) = (X^L - 1)P(X)` to recover
  `P[ell..n)`.
- General output is published only after the complete result exists and may
  overlap either Fr input. Suffix output may overlap values. Scratch cannot
  overlap another argument.
- Failure leaves output unchanged. Suffix scratch becomes unspecified after
  initial validation and contains value-derived intermediates.
- General interpolation may allocate and is variable-time in public points and
  length. For fixed points, value-derived reconstruction uses exact public
  widths and fixed dispatch shapes. Suffix interpolation allocates nothing and
  schedules only from public lengths and indices.

## Construction and dependency map

- `ctFft.h` supplies the single retained radix-two transform.
- `domain.h` constructs one ordered forward-root table.
- The product tree stores `2n-1` topology records in one vector and all node
  polynomials in one offset-addressed Fr arena.
- Remainder descent is point-derived and computes no discarded quotient vector.
- Reconstruction reuses the tree offsets in one separate arena plus one n-element
  temporary. It never trims value-derived coefficients.
- `fft_c_api.cpp` translates typed status to the legacy C `bool` boundary.
- `bindings/blst_fft.h` and `bindings/rust/src/fft.rs` state the same shape,
  alias, timing, allocation, publication, and secret-storage contracts.

## Checklist evidence

### 1. Identity, provenance, and scope

The complete 629-line pre-remediation implementation, header, history, supported
callers, tests, and read-only frozen-source references were inventoried. The
implementation is local and now carries Apache-2.0/Category Labs provenance.

### 2. Preconditions and postconditions

Nonempty/equal spans, supported maximum size, distinct points, power-of-two
domain length, bounded unique indices, suffix length, exact output size, scratch
size, overlap, failure publication, and scratch lifecycle are explicit.

### 3. Correctness

The product-tree derivative/Lagrange construction and the suffix identity were
checked algebraically. Independent Horner-oracle tests recover deterministic
dense-with-zero polynomials at 63, 64, and 65 points, crossing both fast-path
thresholds and an unbalanced tree shape. All-zero values recover all-zero
coefficients without changing reconstruction shapes.

### 4. Memory and lifetime safety

Checked size arithmetic precedes tree reserve and scratch partitioning. Tree
spans are recovered from offsets after arena growth; no retained pointer crosses
a vector operation. The C boundary requires valid arrays of the declared sizes.

### 5. Undefined and implementation-defined behavior

The transform maximum stays below the signed-index engine limit. Power-of-two
growth and polynomial sizes are checked. Fast division rejects an impossible
high-degree remainder instead of truncating it silently; Debug asserts expose
the violated invariant immediately.

### 6. Side channels

General point processing is intentionally variable-time. The value-dependent
half has exact node widths, allocation sizes, dispatch decisions, loop bounds,
and addresses for a fixed point set. Suffix value handling likewise follows a
public index/shape schedule. Fr arithmetic supplies the underlying constant-time
field operations.

### 7. Secret lifecycle

General interpolation uses ordinary heap storage for value-derived weights and
reconstruction data and does not securely erase allocator storage. Suffix
scratch is caller-owned and explicitly documented as value-derived so callers
with an erasure policy can clear it. No broad zeroization overhead was added.

### 8. Exceptions and C ABI

Native entry points are `noexcept` status APIs and translate duplicate points,
allocation failure, invalid sizes, and internal exceptions. No exception crosses
the C ABI. Output publication occurs after every throwing allocation and
calculation.

### 9. API and ABI

Existing C signatures and Rust wrapper signatures remain unchanged. Result-
bearing C declarations are nodiscard. Scratch units remain Fr elements; the
required count is intentionally reduced and callers must query it.

### 10. Concurrency

There is no global state, cache, thread pool, or background work. Every allocation
and scratch region is call-local or caller-owned.

### 11. Resource exhaustion

General interpolation can allocate O(n log n) field storage and rejects sizes
beyond the largest supported internal transform. Suffix interpolation performs
no allocation but retains quadratic work in n; `AUD-301` accepts that explicit
performance tradeoff until deployment evidence justifies a crossover.

### 12. Portability and build behavior

The interface uses C99-compatible public types. Strict MSVC C11 parsing passed.
MSVC Debug compiled with portable AVX support enabled; scalar Release compiled
the implementation and focused oracle object without portable AVX Release.

### 13. Performance-sensitive structure

Persistent per-node ownership allocations, quotient-only temporaries, a third
FFT product buffer, and inverse-root tables were removed. Large multiplication
uses two operand buffers and one forward-root table. No runtime polymorphism,
type erasure, thread scheduling, or hidden hot-path validation was introduced.
No benchmark was run. The simple quadratic suffix is an accepted performance
risk, not an unresolved correctness or memory-safety finding.

### 14. Tests and automated evidence

- Scalar-default plus portable-AVX Debug: eight focused interpolation/C tests.
- Scalar Release: three focused native tests linked from the independently
  compiled FFT test object.
- Rust Release FFT module: 9/9.
- Strict MSVC C11 header parse: `/W4 /WX /Zs` passed.
- Structural scan: no owning product-tree node, discarded quotient path,
  inverse-root suffix table, or value-derived reconstruction trimming remains.

## Findings and decision

`AUD-025`, `AUD-034`, `AUD-297` through `AUD-300`, and `AUD-302` are validated.
`AUD-301` is accepted risk; a replacement still requires explicitly approved
sequential deployment-CPU benchmark evidence.

Milestone sanitizer, fuzzing, native AVX-512, and Apple evidence remains pending.
BTx/BTE was neither audited nor modified; public-header rebuilding exposed only
its already-known frozen ignored-nodiscard warning.
