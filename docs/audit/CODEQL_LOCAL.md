# Local CodeQL Procedure

## Purpose and license boundary

The repository is private and does not currently have GitHub Code Security.
The checked-in GitHub workflow therefore remains disabled unless the repository
variable `CODEQL_ENABLED=true` is explicitly set after Code Security is
enabled. Phase 12 runs CodeQL manually and locally instead of bypassing the
private-repository service boundary through CI artifacts.

The CodeQL bundle is distributed under the GitHub CodeQL Terms and Conditions.
Downloading it signifies acceptance of those terms. Confirm that the intended
use is eligible before downloading or running the CLI. CatBLST carries the
Apache-2.0 license, but this procedure does not assert a separate legal
conclusion and must not be converted into automated CI without confirming the
applicable CodeQL license.

## Pinned bundle

- Release: `codeql-bundle-v2.26.2`
- Platform: Linux x64
- Archive: `codeql-bundle-linux64.tar.gz`
- Size: 822,920,884 bytes
- SHA-256:
  `cb361567fa1bdb9d322da4240f621b36f245e4d7bb97db3c3a2ad7f743c8e8e7`
- URL:
  `https://github.com/github/codeql-action/releases/download/codeql-bundle-v2.26.2/codeql-bundle-linux64.tar.gz`

The gzip archive avoids adding a local `zstd` prerequisite. Allow several GiB
for the extracted bundle, databases, query cache, and results.

## Installation in WSL

Choose an explicit WSL-native installation directory. The following example
uses `/opt/catblst-codeql/2.26.2`:

```bash
archive=/var/tmp/codeql-bundle-linux64-2.26.2.tar.gz
install_root=/opt/catblst-codeql/2.26.2

curl --fail --location \
  https://github.com/github/codeql-action/releases/download/codeql-bundle-v2.26.2/codeql-bundle-linux64.tar.gz \
  --output "$archive"

printf '%s  %s\n' \
  cb361567fa1bdb9d322da4240f621b36f245e4d7bb97db3c3a2ad7f743c8e8e7 \
  "$archive" |
  sha256sum --check -

mkdir -p "$install_root"
tar -xzf "$archive" -C "$install_root"

"$install_root/codeql/codeql" version
"$install_root/codeql/codeql" resolve packs
```

Do not execute an unverified archive. Do not replace an existing installation
directory; choose a new versioned directory for a future CodeQL release.

## Running an analysis

Invoke the runner from WSL at the repository root:

```bash
CODEQL_EXE=/opt/catblst-codeql/2.26.2/codeql/codeql \
  bash scripts/audit/run_codeql_local.sh all
```

The optional mode is `cpp-scalar`, `cpp-avx`, `cpp`, `rust`, or `all`.
`cpp` and `all` create independent scalar and portable-AVX databases rather
than combining both compiler profiles. Worker counts default to two and can be
controlled without editing the scripts:

```bash
BLST_CODEQL_THREADS=4 \
BLST_CODEQL_BUILD_JOBS=2 \
CODEQL_EXE=/opt/catblst-codeql/2.26.2/codeql/codeql \
  bash scripts/audit/run_codeql_local.sh cpp
```

The runner refuses to overwrite prior databases. It writes databases under
the WSL-native `/var/tmp/catblst-codeql/<run-id>/` by default and writes the
following ignored evidence under `out/audit/codeql/<run-id>/`:

- `cpp-scalar.sarif`, `cpp-avx-portable.sarif`, and/or `rust.sarif`;
- the complete run log;
- the analyzed Git commit and exact query suites; and
- CodeQL, CMake, compiler, and Rust versions.

Before extraction, the runner snapshots tracked and non-ignored untracked
files from the current working tree. The snapshot includes local source edits
but excludes ignored build products such as `out/` and
`bindings/rust/target/`. Its file list, archive hash, and working-tree status
are captured with the evidence so the analyzed input is identifiable.

C/C++ database extraction builds only `blst_static` in Debug mode. Portable
scalar and portable AVX are extracted into separate databases so a partial AVX
extraction cannot obscure scalar coverage or be mistaken for a complete
combined result. Assembly, tests, benchmarks, frozen BTx/BTE, experimental
CUDA, and removed bindings are not added to the scan. Rust uses CodeQL's
build-free extraction mode.

## Query suites

- C/C++:
  `codeql/cpp-queries@1.8.0:codeql-suites/cpp-security-and-quality.qls`
- Rust:
  `codeql/rust-queries@0.1.39:codeql-suites/rust-security-and-quality.qls`
- SARIF format: fixed `sarifv2.1.0`

These query-pack versions are bundled with CodeQL 2.26.2 and match the first
hosted Phase 12 run.

## Triage and cleanup

Preserve databases until extraction warnings and errors have been understood.
Every applicable diagnostic or alert is mapped to an `AUD-*` finding,
duplicate, or documented false positive according to
`AUTOMATED_ANALYSIS_PLAN.md`.

After triage, remove only the exact versioned run directory printed by the
runner. Raw databases and SARIF remain ignored and are not committed. Durable
findings, tool versions, commands, exclusions, and residual gaps belong in the
audit records.

## Compensating Clang extraction

CodeQL extractor diagnostics are not treated as successful analysis. Run the
pinned compiler-native Clang Static Analyzer profiles documented in
`AUTOMATED_ANALYSIS_PLAN.md` to cover buildable scalar and AVX template
instantiations that CodeQL does not fully model:

```bash
bash scripts/audit/run_clang_analyzer_local.sh all
```

The runner requires Ubuntu Clang 18.1.3 plus `scan-build`, snapshots the same
tracked and non-ignored source set, and writes ignored SARIF evidence beneath
`out/audit/clang-analyzer/<run-id>/`. Scalar and portable AVX are separate
analyzer runs. Tests, benchmarks, assembly, BTx/BTE, and CUDA remain excluded.
