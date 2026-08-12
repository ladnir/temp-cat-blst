# Copyright 2026 Category Labs, Inc.
# Licensed under the Apache License, Version 2.0, see LICENSE for details.
# SPDX-License-Identifier: Apache-2.0

[CmdletBinding()]
param()

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '..\..')).Path
$libPath = Join-Path $repoRoot 'bindings\rust\src\lib.rs'
$guidePath = Join-Path $repoRoot 'bindings\rust\docs\full_bindings.md'
$lines = [System.IO.File]::ReadAllLines($libPath)
$source = [System.IO.File]::ReadAllText($libPath)
$guide = [System.IO.File]::ReadAllText($guidePath)
$errors = [System.Collections.Generic.List[string]]::new()

# Keep the reviewed safe surface documented even though the generated raw FFI
# declarations intentionally inherit their contracts from the C header.
for ($i = 0; $i -lt $lines.Length; $i++)
{
    if ($lines[$i] -notmatch '^\s*pub\s+(?:fn|struct|trait|mod)\s+([A-Za-z0-9_]+)')
    {
        continue
    }

    $name = $Matches[1]
    if ($name -eq 'gen_random_key')
    {
        # Test-only helper inside the macro's private test module.
        continue
    }

    $previous = $i - 1
    while ($previous -ge 0 -and
           ($lines[$previous] -match '^\s*$' -or
            $lines[$previous] -match '^\s*#'))
    {
        $previous--
    }
    if ($previous -lt 0 -or $lines[$previous] -notmatch '^\s*///')
    {
        $errors.Add(
            "Missing contract documentation before public '$name' at " +
            "bindings/rust/src/lib.rs:$($i + 1)."
        )
    }
}

if ($source -match 'SecretKey::default\s*\(' -or
    $source -match 'impl\s+Default\s+for\s+SecretKey' -or
    $source -match '(?s)#\[derive\([^\)]*\bDefault\b[^\)]*\)\]\s*' +
                  '#\[zeroize\(drop\)\]\s*pub struct SecretKey')
{
    $errors.Add('SecretKey must not expose or construct the invalid zero Default value.')
}

$uninitializedConstructors = [regex]::Matches(
    $source,
    'MaybeUninit::<blst_scalar>::uninit\(\)'
).Count
if ($uninitializedConstructors -lt 7)
{
    $errors.Add(
        "Expected all seven SecretKey construction/decoding paths to avoid an " +
        "invalid initialized value; found $uninitializedConstructors markers."
    )
}

for ($i = 0; $i -lt $lines.Length; $i++)
{
    if (-not $lines[$i].Contains('transmute::<&[_], &['))
    {
        continue
    }

    $start = [Math]::Max(0, $i - 4)
    $context = ($lines[$start..$i] -join "`n")
    if (-not $context.Contains('SAFETY:'))
    {
        $errors.Add(
            "Slice representation conversion lacks a local SAFETY argument at " +
            "bindings/rust/src/lib.rs:$($i + 1)."
        )
    }
}

foreach ($requiredGuideSection in @(
    '## Validation model',
    '## Aggregate verification flow',
    '## Secret-key invariant',
    'proof of possession',
    'unpredictable, nonzero random weights'
))
{
    if (-not $guide.Contains($requiredGuideSection))
    {
        $errors.Add("Full Rust binding guide lacks required contract text: $requiredGuideSection")
    }
}

if ($errors.Count -ne 0)
{
    foreach ($errorMessage in $errors)
    {
        Write-Error $errorMessage
    }
    throw "Rust API contract check failed with $($errors.Count) error(s)."
}

Write-Host (
    'Rust API contracts passed: documented safe surface, non-default SecretKey, ' +
    'and justified slice representation conversions.'
)
