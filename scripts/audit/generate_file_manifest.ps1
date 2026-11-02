[CmdletBinding()]
param()

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '..\..')).Path
$manifestPath = Join-Path $repoRoot 'docs\audit\FILE_MANIFEST.csv'
$preservedColumns = @(
    'disposition',
    'subsystem',
    'audit_status',
    'profiles',
    'reviewer',
    'review_record',
    'upstream_reference',
    'findings',
    'tests',
    'validation',
    'residual_risk',
    'reviewed_commit',
    'notes'
)

function Get-Classification
{
    param([Parameter(Mandatory)][string]$Path)

    $lower = $Path.ToLowerInvariant()
    $disposition = 'in_scope_support'
    $subsystem = 'repository'
    $profiles = 'process'
    $notes = ''

    $isBenchmark = $lower.StartsWith('benchmark/')

    $activeBtx =
        $lower.StartsWith('bte/btx/') -or
        $lower -eq 'bindings/blst_btx.h' -or
        $lower -like 'test/btx*'

    $frozenBte =
        $lower.StartsWith('bte/') -or
        $lower -eq 'bindings/blst_trx.h'

    $experimentalCuda =
        $lower.StartsWith('blst/cuda/') -or
        $lower.StartsWith('blst/backend/cuda/') -or
        $lower -match '(^|/).*cuda[^/]*\.(c|cc|cpp|cuh|h|hpp)$' -or
        $lower -like 'test/*cuda*' -or
        $lower -in @(
            'blst/exports.h',
            'blst/fields/fptile.h',
            'blst/fields/parallel_prefix.h',
            'blst/common/counting_iterator.h',
            'blst/common/mdspan_compat.h',
            'blst/common/transform_iterator.h'
        )

    if ($isBenchmark)
    {
        $disposition = 'out_of_scope_no_audit'
        $subsystem = 'benchmarks'
        $profiles = 'non_shipping'
        if ($lower -like 'benchmark/btx*')
        {
            $notes = 'Excluded from file-by-file audit; run only for an explicit BTX performance question and never concurrently with another benchmark.'
        }
        elseif ($lower -match 'cuda')
        {
            $notes = 'Excluded benchmark for the retained experimental CUDA subsystem; do not run without an explicit request.'
        }
        else
        {
            $notes = 'Benchmarks are inventoried but explicitly excluded from file-by-file audit.'
        }
    }
    elseif ($activeBtx)
    {
        $disposition = if ($lower.StartsWith('test/')) { 'in_scope_support' } else { 'in_scope_production' }
        $subsystem = 'btx_bte'
        $profiles = 'C|C++|scalar|avx512'
        $notes = 'BTX freeze lifted 2026-08-17; active for performance comparison, specification audit, modification, and extraction.'
    }
    elseif ($frozenBte)
    {
        $disposition = 'out_of_scope_frozen'
        $subsystem = 'btx_bte'
        $profiles = 'frozen_pending_merge'
        $notes = 'PFE, TRX, or unrelated shared BTE code remains frozen; the 2026-08-17 decision unfroze BTX only.'
    }
    elseif ($experimentalCuda)
    {
        $disposition = 'out_of_scope_experimental'
        $subsystem = 'cuda'
        $profiles = 'cuda_experimental'
        $notes = 'Retained experimental, opt-in CUDA surface; outside scalar/AVX-512 production readiness.'
    }
    elseif ($lower.StartsWith('build/') -or
            $lower -eq 'bindings/rust/src/bindings.rs' -or
            $lower -match '^blst/avx/generated.*\.(h|hpp|c|cc|cpp)$')
    {
        $disposition = 'generated_production'
        $subsystem = if ($lower.StartsWith('bindings/rust/')) {
            'rust_ffi'
        } elseif ($lower.StartsWith('blst/avx/')) {
            'avx512_generated'
        } else {
            'generated_assembly'
        }
        $profiles = if ($lower.StartsWith('blst/avx/')) { 'avx512' } else { 'scalar|avx512' }
        $notes = 'Generated production input/output; audit generator and generated artifact together.'
    }
    elseif ($lower.StartsWith('docs/old/') -or
            $lower -match '\.(png|jpg|jpeg|gif|pdf)$')
    {
        $disposition = 'documentation_or_archive'
        $subsystem = 'documentation_archive'
        $profiles = 'non_shipping'
        $notes = 'Catalogued non-shipping documentation or archive.'
    }
    elseif ($lower.StartsWith('blst/'))
    {
        $disposition = 'in_scope_production'
        $profiles = 'scalar'

        if ($lower -match '(^|/)(fft|middle_product)(/|\.|_)')
        {
            $subsystem = if ($lower.StartsWith('blst/avx/')) {
                'fft_middle_product_avx512'
            } else {
                'fft_middle_product'
            }
            $profiles = if ($lower.StartsWith('blst/avx/')) { 'avx512' } else { 'scalar|avx512' }
        }
        elseif ($lower.StartsWith('blst/avx/') -or
            $lower.StartsWith('blst/backend/avx512/') -or
            $lower -eq 'blst/avx_c_api.cpp')
        {
            $subsystem = 'avx512'
            $profiles = 'avx512'
        }
        elseif ($lower.StartsWith('blst/asm/'))
        {
            $subsystem = 'assembly_generator'
            $profiles = 'scalar'
        }
        elseif ($lower.StartsWith('blst/backend/'))
        {
            $subsystem = 'scalar_backend'
        }
        elseif ($lower.StartsWith('blst/fields/'))
        {
            $subsystem = 'fields'
        }
        elseif ($lower.StartsWith('blst/curve/'))
        {
            $subsystem = 'curve'
        }
        elseif ($lower.StartsWith('blst/common/'))
        {
            $subsystem = 'common'
        }
        elseif ($lower.StartsWith('blst/mains/'))
        {
            $subsystem = 'command_demo'
        }
        else
        {
            $subsystem = 'core'
        }
    }
    elseif ($lower -match '^bindings/blst(_.*)?\.(h|hpp)$')
    {
        $disposition = 'in_scope_production'
        $subsystem = if ($lower.EndsWith('.hpp')) { 'cpp_api' } else { 'c_api' }
        $profiles = if ($lower.Contains('_avx')) { 'C|C++|avx512' } else { 'C|C++|scalar|avx512' }
    }
    elseif ($lower.StartsWith('bindings/rust/'))
    {
        if ($lower.StartsWith('bindings/rust/src/') -or
            $lower -eq 'bindings/rust/build.rs' -or
            $lower -eq 'bindings/rust/cargo.toml')
        {
            $disposition = 'in_scope_production'
            $subsystem = 'rust_binding'
            $profiles = 'Rust|scalar|avx512'
        }
        else
        {
            $disposition = 'in_scope_support'
            $subsystem = 'rust_binding_support'
            $profiles = 'Rust|scalar|avx512'
        }
    }
    elseif ($lower.StartsWith('bindings/vectors/'))
    {
        $disposition = 'in_scope_support'
        $subsystem = 'test_vectors'
        $profiles = 'scalar|avx512|C|C++|Rust'
    }
    elseif ($lower.StartsWith('test/'))
    {
        $disposition = 'in_scope_support'
        $subsystem = 'tests'
        $profiles = 'scalar|avx512'
    }
    elseif ($lower.StartsWith('benchmark/'))
    {
        $disposition = 'in_scope_support'
        $subsystem = 'benchmarks'
        $profiles = 'scalar|avx512'
    }
    elseif ($lower.StartsWith('.github/'))
    {
        $disposition = 'in_scope_support'
        $subsystem = 'ci'
        $profiles = 'process'
    }
    elseif ($lower.StartsWith('scripts/audit/') -or $lower.StartsWith('docs/audit/'))
    {
        $disposition = 'in_scope_support'
        $subsystem = 'audit_process'
        $profiles = 'process'
    }
    elseif ($lower -match '(^|/)cmakelists\.txt$' -or
            $lower -eq 'cmakepresets.json')
    {
        $disposition = 'in_scope_support'
        $subsystem = 'build_configuration'
        $profiles = 'scalar|avx512|C|C++|Rust'
    }

    $status = switch ($disposition)
    {
        'out_of_scope_frozen' { 'frozen' }
        'out_of_scope_experimental' { 'experimental' }
        'out_of_scope_planned_removal' { 'scope_review_required' }
        'out_of_scope_no_audit' { 'excluded' }
        'documentation_or_archive' { 'catalogued' }
        default { 'pending' }
    }

    return [pscustomobject]@{
        disposition = $disposition
        subsystem = $subsystem
        audit_status = $status
        profiles = $profiles
        notes = $notes
    }
}

$oldRows = @{}
if (Test-Path -LiteralPath $manifestPath)
{
    foreach ($row in @(Import-Csv -LiteralPath $manifestPath))
    {
        if (-not [string]::IsNullOrWhiteSpace($row.path))
        {
            $oldRows[$row.path] = $row
        }
    }
}

$paths = [System.Collections.Generic.HashSet[string]]::new(
    [System.StringComparer]::Ordinal
)
foreach ($path in @(& git -C $repoRoot ls-files --cached --others --exclude-standard))
{
    if (-not [string]::IsNullOrWhiteSpace($path))
    {
        [void]$paths.Add($path.Replace('\', '/'))
    }
}

[void]$paths.Add('docs/audit/FILE_MANIFEST.csv')

$rows = foreach ($path in @($paths) | Sort-Object)
{
    $classification = Get-Classification -Path $path
    $values = [ordered]@{
        path = $path
        disposition = $classification.disposition
        subsystem = $classification.subsystem
        audit_status = $classification.audit_status
        profiles = $classification.profiles
        reviewer = ''
        review_record = ''
        upstream_reference = ''
        findings = ''
        tests = ''
        validation = ''
        residual_risk = ''
        reviewed_commit = ''
        notes = $classification.notes
    }

    if ($oldRows.ContainsKey($path))
    {
        $old = $oldRows[$path]
        $liftedFreeze =
            [string]$old.disposition -eq 'out_of_scope_frozen' -and
            [string]$classification.disposition -ne 'out_of_scope_frozen'
        foreach ($column in $preservedColumns)
        {
            if ($liftedFreeze -and $column -in @(
                    'disposition',
                    'subsystem',
                    'audit_status',
                    'profiles',
                    'notes'
                ))
            {
                continue
            }
            if ($old.PSObject.Properties.Name -contains $column)
            {
                $oldValue = [string]$old.$column
                if (-not [string]::IsNullOrWhiteSpace($oldValue))
                {
                    $values[$column] = $oldValue
                }
            }
        }
    }

    [pscustomobject]$values
}

$manifestDirectory = Split-Path -Parent $manifestPath
if (-not (Test-Path -LiteralPath $manifestDirectory))
{
    New-Item -ItemType Directory -Path $manifestDirectory | Out-Null
}

$csvLines = @($rows | ConvertTo-Csv -NoTypeInformation)
[System.IO.File]::WriteAllLines(
    $manifestPath,
    $csvLines,
    [System.Text.UTF8Encoding]::new($false)
)

$counts = $rows | Group-Object disposition | Sort-Object Name
Write-Output "Wrote $($rows.Count) rows to $manifestPath"
foreach ($count in $counts)
{
    Write-Output ("  {0}: {1}" -f $count.Name, $count.Count)
}
