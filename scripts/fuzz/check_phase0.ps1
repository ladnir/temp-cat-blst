# Copyright 2026 Category Labs
# SPDX-License-Identifier: Apache-2.0

[CmdletBinding()]
param(
    [switch]$RequireGeneratedInventory
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

$root = (Resolve-Path (Join-Path $PSScriptRoot '..\..')).Path
$python = Get-Command python -ErrorAction Stop

if ($RequireGeneratedInventory) {
    & $python.Source (Join-Path $root 'scripts\fuzz\api_inventory.py') `
        '--root' $root '--self-test' '--require-generated'
    if ($LASTEXITCODE -ne 0) {
        throw 'Fuzz API inventory gate failed.'
    }
}

$structureArgs = @('--root', $root, '--self-test')
if (-not $RequireGeneratedInventory) {
    $structureArgs += '--boundary-only'
}
& $python.Source (Join-Path $root 'scripts\fuzz\check_structure.py') @structureArgs
if ($LASTEXITCODE -ne 0) {
    throw 'Fuzz structural gate failed.'
}

Write-Host 'Fuzz Phase 0 structural checks passed.'
