# File Review: Rust package metadata and guides

## Record

| Field | Value |
| --- | --- |
| Paths | `bindings/rust/Cargo.toml`, `README.md`, `docs/overview.md`, `docs/full_bindings.md` |
| Roadmap phase / queue | Phase 18 Rust packaging and supported-user documentation |
| Reviewed commits | Package/README `b2c93e9d3ec4b4478edca51c9cd67187f859eba4`; guides `778f2f60d1fd4ba19569ba0dd1d4e3f5f1a50b88` |
| Reviewer | Codex audit team |
| Date | 2026-08-12 |
| Disposition | in-scope production metadata and support documentation |
| Profiles | Rust scalar, native AVX-512, no-std alloc; Windows portable AVX Debug development only |
| Upstream references | Common base `f48500c1fdbefa7c0bf9800bccd65d28236799c1`; local support policy intentionally narrowed |

## Reviewed behavior

Cargo declares Rust 1.65 as the enforced MSRV. Serde disables default features
and opts into `alloc`/`derive`; `serde-secret` remains explicit. Feature docs
distinguish portable scalar from native AVX and Windows's developer-only Debug
emulation. Cross builds require a toolchain, prebuilt use is explicit, Apple
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
- `cargo metadata` reports Rust 1.65; no-std Serde check passes; rustdoc passes
  with warnings denied; `check_rust_api_contracts.ps1` enforces the critical
  guide sections.

## Findings and residual risk

`AUD-328`, `AUD-330`, `AUD-331`, `AUD-333` are validated. Exact generated-
binding regeneration (`AUD-324`) remains separate and open. Publishing/package
installation on every supported platform remains Phase 19 release evidence;
the reviewed metadata and guides are complete.
