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

$cargoCommands = [regex]::Matches($ci, '(?m)^\s+cargo (test|run) .*$')
if ($cargoCommands.Count -eq 0) {
    throw 'ci.yml: no retained Cargo test or smoke commands were found.'
}
foreach ($command in $cargoCommands) {
    if ($command.Value -notmatch ' --locked(?:\s|$)') {
        throw "ci.yml: Cargo command is not locked: $($command.Value.Trim())"
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

$dependabotPath = Join-Path $repoRoot '.github\dependabot.yml'
$dependabot = Get-Content -LiteralPath $dependabotPath -Raw
foreach ($required in @(
    'package-ecosystem: github-actions',
    'directory: /',
    'interval: weekly',
    'groups:',
    '      github-actions:'
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
