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

$inventoryArgs = @(
    (Join-Path $root 'scripts\fuzz\api_inventory.py'),
    '--root', $root,
    '--self-test'
)
if ($RequireGeneratedInventory) {
    $inventoryArgs += '--require-generated'
}
& $python.Source @inventoryArgs
if ($LASTEXITCODE -ne 0) {
    throw 'Fuzz API inventory gate failed.'
}

& $python.Source (Join-Path $root 'scripts\fuzz\check_structure.py') `
    '--root' $root '--self-test'
if ($LASTEXITCODE -ne 0) {
    throw 'Fuzz structural gate failed.'
}

Write-Host 'Fuzz Phase 0 structural checks passed.'
