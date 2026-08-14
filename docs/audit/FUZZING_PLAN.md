# API-Complete Fuzzing Plan

## Objective

The fuzzing program covers every supported C, native C++, and safe Rust API.
Clang libFuzzer runs the primary scalar profile under ASan and UBSan.

Fuzzing searches for memory-safety defects, undefined behavior, contract
violations, and arithmetic inconsistencies. It does not establish
constant-time behavior. It also does not replace proof or hand-audit evidence.

The harness must distinguish adversarial input from contract violations:

- validation and decoding functions receive arbitrary bytes;
- arithmetic functions receive values constructed for their documented domain;
- fuzz-only types enforce subgroup, nonzero, reduced-field, aliasing, and state
  preconditions;
- raw library calls occur only inside one reviewed adapter; and
- a structural gate rejects direct calls and unclassified APIs elsewhere.

Fuzz-only code must not change a production API, representation, dispatch
decision, allocation policy, or hot path.

## Scope

The primary profile is the pure scalar C/C++ backend on x86-64 Linux. Later
profiles cover safe Rust and native AVX-512 differential execution.

The following areas remain excluded:

- frozen BTx/BTE;
- experimental CUDA;
- experimental non-Apple ARM64 and CHERI;
- benchmarks; and
- portable-AVX Release builds.

Portable AVX may replay a bounded corpus in Debug mode. It is developer
evidence, not native AVX or production evidence.

## Phase 0: inventory and enforce preconditions

### Symbol ledger

`FUZZ_API_LEDGER.csv` assigns every supported symbol to one campaign and one
classification. Each row contains:

`surface`, `symbol`, `owner`, `campaign`, `classification`,
`required_domain`, `factory`, `allowed_invalid_cases`, `oracles`, and `notes`.

The classifications have these meanings:

- `direct`: the API accepts arbitrary fuzzer-controlled input;
- `constructed`: a factory establishes every mathematical precondition first;
- `oracle`: the API observes or validates another result;
- `indirect`: another public surface reaches the same implementation;
- `regression_only`: deterministic testing is more appropriate; and
- `excluded`: the row gives an explicit scope reason.

The inventory gate extracts all C declarations from the scalar, AVX, FFT, and
middle-product headers. It extracts native C++ methods from the reviewed public
owners through Clang AST JSON. It extracts safe Rust methods through rustdoc
JSON from the pinned fuzzing toolchain. Generated raw Rust FFI functions are
`indirect` through their corresponding C rows.

The gate fails if a supported API is missing, duplicated, or assigned to an
unknown campaign. Zero-argument accessors, constants, and ABI layout queries
are normally `oracle` or `regression_only`. Scratch queries are `constructed`
because shape campaigns consume their results.

### Domain values

The fuzz support library defines `DomainValue<Storage, Tag>` with private
constructors. Only reviewed factories can construct or refine a domain value.

The required domains are:

- `ScalarBits`: at most 256 meaningful little-endian bits;
- `CanonicalScalar`: a scalar in `[0,r)`;
- `SecretScalar`: a scalar in `[1,r)`;
- `ReducedFr`, `ReducedFp`, `ReducedFp2`, and `ReducedFp12`;
- `Nonzero<T>`;
- `FieldValidG1Candidate` and `FieldValidG2Candidate`;
- affine and projective `OnCurveG1` and `OnCurveG2`;
- affine and projective `SubgroupG1` and `SubgroupG2`;
- `IncompleteAddPair<T>` for finite, distinct, non-opposite inputs;
- `MillerValue` and `GtElement`;
- `PreparedFftDomain` and `LiveMiddleProductHandle`; and
- `PairingState<Orientation, PendingState>`.

Factories implement the following refinements:

- arbitrary bytes to a successfully decoded on-curve point;
- an on-curve point to a subgroup point after membership succeeds;
- arbitrary bytes to a reduced scalar or field value;
- a canonical scalar to a nonzero secret scalar;
- reduced coordinates to an on-curve point after the curve predicate succeeds;
- a Miller value to GT through final exponentiation.

A failed refinement skips only the narrower operation. It is not a library
failure.

### Structural boundary

LibFuzzer entrypoints include only the byte-span driver interface. Only the
internal adapter includes public library declarations or native C++ arithmetic
headers. Adapter functions accept domain values rather than raw mathematical
objects.

The structural gate rejects:

- direct `blst_*` or native arithmetic calls outside the adapter;
- public library headers outside the adapter;
- `reinterpret_cast` in fuzz code;
- raw-byte fabrication of a domain value;
- forged opaque handles; and
- adapter calls absent from the symbol ledger.

Compile-time tests verify that a broader domain cannot invoke a narrower
adapter.

### Pointer and shape domains

Typed call-frame builders represent:

- valid disjoint storage;
- documented exact aliasing;
- null pointers paired with zero length;
- documented invalid sizes;
- documented short scratch; and
- documented overlap rejection.

Harnesses never generate arbitrary addresses, unreadable fixed-size buffers,
unsupported partial overlap, forged handles, use-after-free, or unsupported
misalignment.

C point decoding always receives the largest documented readable buffer. Rust
slice APIs may receive arbitrary lengths because their contracts check length.

## Phase 1: byte boundaries

### `blst_fuzz_decode`

The target directly covers G1/G2 `uncompress`, `deserialize`, and aggregate
wire decoding. Predicates, representation conversions, equality, serialization,
and compression serve as oracles.

The target checks:

- rejected calls preserve output;
- successful decoding produces an on-curve point;
- re-encoding is canonical;
- canonical decoding is idempotent;
- compressed `deserialize` agrees with `uncompress`;
- conversions preserve subgroup status; and
- aggregate-wire functions agree with decoding, membership, and complete
  addition.

The input limit is 192 bytes. The corpus contains generators, infinity,
retained non-subgroup fixtures, modulus boundaries, malformed flags, and prior
decoding regressions.

### `blst_fuzz_bytes`

The target covers scalar and Fp fixed-width and variable-width imports in both
byte orders. It also covers 32-bit and 64-bit word forms, hexadecimal parsing,
scalar validity, secret-key validity, scalar/Fr conversions, and batch
conversion.

Hexadecimal input is copied into bounded NUL-terminated storage. The byte
reduction limit is 1,024 bytes.

An independent bit-serial modulus oracle uses `uint32_t` limbs and no library
arithmetic. Boundary tests establish the oracle before fuzz comparison.

The target checks endian symmetry, canonical output, exact scalar round trips,
Fp residues, zero-length behavior, and Boolean/nonzero agreement.

## Phase 2: constructed arithmetic

### `blst_fuzz_fields`

The target covers every public Fr, Fp, Fp2, general Fp12, and GT arithmetic
function. Ordinary arithmetic receives reduced Montgomery inputs. Inversion
receives nonzero inputs. Cyclotomic and GT functions receive GT inputs.

The target checks independent reduction, algebraic identities, inverse
identities, square-root postconditions, exponentiation, Frobenius bounds,
canonical exports, and documented exact aliasing.

General Fp12 and GT remain distinct domains. A general Fp12 value cannot reach
a GT or cyclotomic adapter.

### `blst_fuzz_groups`

The target covers every public single-point G1/G2 operation and validator.
Factories establish reduced coordinates, on-curve status, subgroup status, or
incomplete-addition conditions as required.

The target checks complete and incomplete addition, doubling, negation,
representation conversion, equality, serialization, checked and unchecked
multiplication on subgroup inputs, and documented non-subgroup primitives.

Native C++ methods that share the C implementation are `indirect`. Unique
native methods execute through the same typed campaign.

## Phase 3: protocols, state, and batches

### `blst_fuzz_protocol`

The target covers SHA-256, XMD expansion, all key-generation variants,
EIP-2333 derivation, hash/encode/map to G1/G2, public-key derivation, signing,
decorated-output helpers, and one-shot verification.

The target uses these bounds:

- IKM: 256 bytes;
- salt and info: 256 bytes each;
- message: 1,024 bytes;
- DST: 512 bytes;
- augmentation: 512 bytes; and
- XMD output: through 8,161 bytes.

The target checks determinism, invalid-input clearing, subgroup membership,
canonical signatures and public keys, sign/verify agreement, and failure
atomicity.

### `blst_fuzz_pairing_state`

The target covers Miller loops, prepared lines, pairing contexts, aggregate
variants, weighted aggregate, merge, commit, final verification, and external
Miller contributions.

The state model records orientation, pending or committed state, borrowed DST
lifetime, subgroup trust, and external Miller contributions. Rejected
transitions preserve both context bytes and model state.

The target compares one-shot verification, incremental aggregation, raw Miller
products, prepared-line products, and merged committed contexts.

### `blst_fuzz_uniq`

The target covers uniqueness-tree size, alignment, initialization, insertion,
duplicate detection, capacity, and empty messages. The harness owns all
borrowed messages until tree destruction.

Capacity is at most 64 nodes. Total retained message storage is at most 4,096
bytes.

### `blst_fuzz_fft_middle`

The target covers FFT, prepared domains, interpolation, and middle products.
It uses reduced Fr vectors and subgroup points.

FFT size is at most 64. General interpolation size is at most 32.
Middle-product length is at most 16.

The target checks FFT round trips, prepared and allocating agreement, bounded
DFT results, interpolation reconstruction, invalid size tags, overlap rules,
failure atomicity, direct pairing formulas, update/combine agreement, and
handle lifetime.

### `blst_fuzz_msm`

The target covers scalar G1/G2 Wbits, Pippenger, Fr-scalar, batch-affine, tile,
scratch-query, and batch-conversion APIs. It uses at most 16 subgroup points.

Every result is compared with a direct sum of single-point multiplications.
Only documented invalid shapes are generated.

## Phase 4: AVX and Rust

### `blst_fuzz_avx`

The target replays decoder corpora through AVX G1/G2 decoders. It compares
every nontrivial public AVX conversion, point, MSM, Fp12, GT, and pairing API
with the scalar implementation.

Native AVX-512/IFMA runs only on Peach with the supported compiler. Portable
AVX runs only as bounded Debug replay. No native AVX campaign runs without
explicit approval.

Batch subgroup checks must accept subgroup inputs and reject off-curve inputs.
Randomized non-subgroup false acceptance is not an absolute oracle because the
API states a probability bound.

### Safe Rust

The pinned `cargo-fuzz` package provides two targets:

- `blst_rust_fuzz_decode` covers byte constructors, validation, Serde decoding,
  and both `min_pk` and `min_sig`;
- `blst_rust_fuzz_state` covers aggregation, signing, verification, pairing
  ownership and state, FFT, and middle-product wrappers.

Raw Rust FFI remains `indirect` through C campaigns. Safe Rust targets slice
bounds, ownership, panic freedom, practical allocation-failure translation,
and agreement with C. Native C/C++ code reached through Rust uses the same
sanitizer profile.

## Phase 5: execution and closure

### Build profile

`BLST_FUZZ=ON` is an advanced Clang-only option. It requires `BLST_ASAN=ON`
and `BLST_UBSAN=ON`. It adds libFuzzer only to fuzz targets and excludes them
from ordinary builds.

The primary scalar profile disables assembly, ADX, AVX, CUDA, tests, and
benchmarks. Discovery uses `RelWithDebInfo`. Debug builds replay retained
corpora with assertions.

Builds use at most two workers. Fuzzing uses exactly one worker and one job.
No large run starts without explicit approval. An approved default discovery
packet runs each target sequentially for 15 minutes.

Per-input limits are:

- two seconds for decode, byte, field, and group targets;
- ten seconds for protocol, pairing, uniqueness, FFT/middle-product, MSM, and
  AVX targets; and
- 2 GiB RSS for every target.

### Corpus and findings

Each target owns a reviewed corpus and optional dictionary. Every confirmed
defect receives:

1. the next available `AUD-*` entry;
2. a minimized reproducer;
3. a deterministic non-fuzz regression test; and
4. a retained seed when the seed preserves material coverage.

Crash artifacts remain under `out/` until triage. CI never modifies the
checked-in corpus.

### Coverage and automation

Runs record libFuzzer edge, feature, corpus, and execution counts. Counts are
comparable only for identical compiler and build profiles.

LLVM source coverage is descriptive. It must show each named success,
rejection, boundary, aliasing, and round-trip path. The initial program has no
percentage gate.

Pull requests run structural gates and deterministic corpus replay only. Long
discovery is manual. Peach AVX is manual. Apple Silicon performs
release-candidate replay, not long discovery. Workflow concurrency prevents
simultaneous fuzz jobs. Fuzz workflows never run benchmarks.

## Acceptance criteria

Phase 0 is complete when:

- the ledger exactly matches every supported C symbol and reviewed C++/Rust
  surface;
- structural self-tests reject missing rows, unknown calls, forbidden includes,
  and direct library calls;
- compile-time tests prevent broad domains from invoking narrow adapters;
- factory tests cover every refinement result; and
- call-frame tests cover valid disjoint storage and every permitted invalid
  shape.

Each campaign must pass deterministic corpus replay before discovery. The two
initial targets must replay nonempty corpora under Release ASan/UBSan and Debug
assertions. Independent modulus oracles must pass retained boundary vectors
before fuzz use.

Every mismatch is classified as a harness defect, library finding, duplicate,
or documented contract rejection.
