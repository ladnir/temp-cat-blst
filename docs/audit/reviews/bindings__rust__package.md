# File Review: Rust package metadata and guides

## Record

| Field | Value |
| --- | --- |
| Paths | `bindings/rust/Cargo.toml`, `Cargo.lock`, `README.md`, `docs/overview.md`, `docs/full_bindings.md` |
| Roadmap phase / queue | Phase 18 Rust packaging and supported-user documentation |
| Reviewed commits | Package/MSRV/lock `7956fb53166729089305317f4f0e2e81b2621814`; dependency automation `a71068a47c361c9a1edae4a706b75435023756f6`; guides `778f2f60d1fd4ba19569ba0dd1d4e3f5f1a50b88` |
| Reviewer | Codex audit team |
| Date | 2026-08-12 |
| Disposition | in-scope production metadata and support documentation |
| Profiles | Rust scalar, native AVX-512, no-std alloc; Windows portable AVX development only |
| Upstream references | Common base `f48500c1fdbefa7c0bf9800bccd65d28236799c1`; local support policy intentionally narrowed |

## Reviewed behavior

Cargo declares Rust 1.78 as the enforced MSRV and tracks the version-4 lockfile
used by CI. Serde disables default features and opts into `alloc`/`derive`;
`serde-secret` remains explicit. Feature docs
distinguish portable scalar from native AVX and Windows's developer emulation.
Optimized emulation requires explicit selection and is not part of the ordinary
release matrix. Cross builds require a toolchain, prebuilt use is explicit, Apple
Silicon is supported, other ARM64/CHERI is experimental, and CUDA is outside
the retained release claim.

The overview maps each module to its owner. The full-binding guide defines the
validation-state model, decode-versus-subgroup postconditions, aggregate flow,
proof-of-possession assumption, batch randomness, caller-thread execution, and
the nonzero `SecretKey` invariant.

## Checklist and validation

- API/portability: feature/MSRV/no-std claims agree with Cargo and build policy.
- Security: guides do not imply that safe Rust performs hidden subgroup checks.
- Concurrency/resources: no internal scheduler; allocations are described at
  their safe owners.
- Performance: portable AVX is never presented as deployment or native proof.
- Rust/Cargo 1.78 locked metadata and the combined supported-feature check
  pass; no-std Serde check passes; rustdoc passes
  with warnings denied; `check_rust_api_contracts.ps1` enforces the critical
  guide sections.

## Findings and residual risk

`AUD-328`, `AUD-330`, `AUD-331`, `AUD-333`, and `AUD-351` are validated. Exact generated-
binding regeneration (`AUD-324`) remains separately validated. Publishing/package
installation on every supported platform remains Phase 19 release evidence;
the reviewed metadata and guides are complete.
