# Copyright 2026 Category Labs
# SPDX-License-Identifier: Apache-2.0

[CmdletBinding()]
param()

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '..\..')).Path
$workflowRoot = Join-Path $repoRoot '.github\workflows'
$workflowFiles = @(Get-ChildItem -LiteralPath $workflowRoot -Filter '*.yml' |
    Sort-Object Name)

if ($workflowFiles.Count -eq 0) {
    throw 'No GitHub Actions workflows were found.'
}

foreach ($workflow in $workflowFiles) {
    $text = Get-Content -LiteralPath $workflow.FullName -Raw

    if ($text -notmatch '(?m)^concurrency:\s*$') {
        throw "$($workflow.Name): missing workflow-level concurrency policy."
    }
    if ($text -notmatch '(?m)^permissions:\s*$') {
        throw "$($workflow.Name): missing workflow-level permissions policy."
    }

    $jobsOffset = $text.IndexOf("`njobs:")
    if ($jobsOffset -lt 0) {
        throw "$($workflow.Name): missing jobs mapping."
    }
    $jobsText = $text.Substring($jobsOffset)
    $jobCount = [regex]::Matches($jobsText, '(?m)^  [A-Za-z0-9_-]+:\s*$').Count
    $timeoutCount = [regex]::Matches($jobsText, '(?m)^    timeout-minutes:\s*[0-9]+\s*$').Count
    if ($jobCount -eq 0 -or $timeoutCount -ne $jobCount) {
        throw "$($workflow.Name): every job must have one timeout-minutes value."
    }

    foreach ($match in [regex]::Matches($text, '(?m)^\s*uses:\s*([^\s#]+)')) {
        $reference = $match.Groups[1].Value
        if ($reference.StartsWith('./')) {
            continue
        }
        if ($reference -notmatch '^[^@\s]+@[0-9a-f]{40}$') {
            throw "$($workflow.Name): external action is not pinned to a full commit SHA: $reference"
        }
    }

    if ($text -match '(?m)^\s*runs-on:\s*(ubuntu-latest|windows-latest|macos-latest)\s*$') {
        throw "$($workflow.Name): mutable latest runner alias is not allowed."
    }
}

$ciPath = Join-Path $workflowRoot 'ci.yml'
$ci = Get-Content -LiteralPath $ciPath -Raw
$expectedPresets = @(
    'linux-scalar',
    'linux',
    'macos-apple-silicon',
    'x64-Release'
)
$actualPresets = @([regex]::Matches($ci, '(?m)^\s+preset:\s*([^\s#]+)\s*$') |
    ForEach-Object { $_.Groups[1].Value })
if (($actualPresets -join '|') -ne ($expectedPresets -join '|')) {
    throw "ci.yml: supported CMake preset matrix drifted: $($actualPresets -join ', ')"
}
if ($ci -match '(?m)^\s*\*\*/Cargo\.lock\s*$') {
    throw 'ci.yml: Cargo.lock must never be restored from an Actions cache.'
}
if ($ci -notmatch "hashFiles\('bindings/rust/Cargo.lock'\)") {
    throw 'ci.yml: the Cargo cache key must include the committed lockfile hash.'
}
$cargoLockPath = Join-Path $repoRoot 'bindings\rust\Cargo.lock'
if (-not (Test-Path -LiteralPath $cargoLockPath -PathType Leaf)) {
    throw 'ci.yml: locked Rust CI requires the tracked bindings/rust/Cargo.lock.'
}
$cargoManifestPath = Join-Path $repoRoot 'bindings\rust\Cargo.toml'
$cargoManifest = Get-Content -LiteralPath $cargoManifestPath -Raw
if ($cargoManifest -notmatch '(?m)^rust-version\s*=\s*"1\.78"\s*$') {
    throw 'Cargo.toml: the tested MSRV must remain Rust 1.78.'
}
$cargoLock = Get-Content -LiteralPath $cargoLockPath -Raw
if ($cargoLock -notmatch '(?m)^version\s*=\s*4\s*$') {
    throw 'Cargo.lock: the reviewed Rust 1.78 CI graph must use lockfile version 4.'
}

$cargoCommands = [regex]::Matches(
    $ci,
    '(?m)^\s+cargo(?: \+[^\s]+)? (test|run|check) .*$')
if ($cargoCommands.Count -eq 0) {
    throw 'ci.yml: no retained Cargo test, smoke, or MSRV check commands were found.'
}
foreach ($command in $cargoCommands) {
    if ($command.Value -notmatch ' --locked(?:\s|$)') {
        throw "ci.yml: Cargo command is not locked: $($command.Value.Trim())"
    }
}
foreach ($required in @(
    'rustup toolchain install 1.78.0 --profile minimal',
    'cargo +1.78.0 check --locked --lib',
    '--features full-rust-bindings,portable,serde-secret'
)) {
    if (-not $ci.Contains($required)) {
        throw "ci.yml: missing Rust 1.78 MSRV contract '$required'."
    }
}

$portablePath = Join-Path $workflowRoot 'portable-avx-debug.yml'
$portable = Get-Content -LiteralPath $portablePath -Raw
foreach ($required in @(
    '-DCMAKE_BUILD_TYPE=Debug',
    '-DBLST_ASM=OFF',
    '-DBLST_ADX=OFF',
    '-DBLST_AVX512=OFF',
    '-DBLST_AVX512_PORTABLE=ON',
    '-DBLST_DEFAULT_AVX=ON'
)) {
    if (-not $portable.Contains($required)) {
        throw "portable-avx-debug.yml: missing required developer-emulation setting $required"
    }
}

$fuzzPath = Join-Path $workflowRoot 'fuzz.yml'
$fuzz = Get-Content -LiteralPath $fuzzPath -Raw
foreach ($required in @(
    'group: category-blst-fuzz',
    'python3 scripts/fuzz/check_structure.py --boundary-only --self-test',
    'python3 scripts/fuzz/check_coverage_paths.py --check-config --self-test',
    'cargo +nightly-2026-01-15 check',
    'cargo +nightly-2026-01-15 install cargo-fuzz --version 0.13.2 --locked',
    '--locked',
    'bash scripts/fuzz/replay.sh',
    'bash scripts/fuzz/replay_rust.sh',
    'BLST_FUZZ_DISCOVERY_APPROVED: YES',
    'bash scripts/fuzz/discovery.sh',
    'bash scripts/fuzz/discovery_rust.sh',
    'bash scripts/fuzz/coverage_replay.sh',
    'runs-on: macos-15',
    'bash scripts/fuzz/configure_avx.sh out/fuzz-avx-portable portable'
)) {
    if (-not $fuzz.Contains($required)) {
        throw "fuzz.yml: missing required bounded-fuzz policy '$required'."
    }
}
if ($fuzz -match '(?m)^\s*run:.*benchmark') {
    throw 'fuzz.yml: fuzz workflows must not invoke benchmarks.'
}

$dependabotPath = Join-Path $repoRoot '.github\dependabot.yml'
$dependabot = Get-Content -LiteralPath $dependabotPath -Raw
foreach ($required in @(
    'package-ecosystem: github-actions',
    'directory: /',
    'interval: weekly',
    'groups:',
    '      github-actions:',
    'package-ecosystem: cargo',
    'directory: /bindings/rust',
    '      rust-dependencies:'
)) {
    if (-not $dependabot.Contains($required)) {
        throw "dependabot.yml: missing retained action-update policy '$required'."
    }
}

$allWorkflowText = ($workflowFiles | ForEach-Object {
    Get-Content -LiteralPath $_.FullName -Raw
}) -join "`n"
$pinnedActionCount = [regex]::Matches(
    $allWorkflowText,
    '(?m)^\s*uses:\s*[^\s#]+@[0-9a-f]{40}').Count

$summary = (
    'CI policy validation passed: {0} workflows, {1} pinned external action uses, ' +
    '{2} supported CMake presets, and {3} locked Cargo commands.') -f
    $workflowFiles.Count,
    $pinnedActionCount,
    $actualPresets.Count,
    $cargoCommands.Count
Write-Output $summary
