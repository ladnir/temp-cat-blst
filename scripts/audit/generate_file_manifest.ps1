[CmdletBinding()]
param()

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '..\..')).Path
$manifestPath = Join-Path $repoRoot 'docs\audit\FILE_MANIFEST.csv'
$preservedColumns = @(
    'audit_status',
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

    $plannedRemoval =
        $lower.StartsWith('bte/') -or
        $lower -eq 'bindings/blst_btx.h' -or
        $lower -eq 'bindings/blst_trx.h' -or
        $lower -like 'test/btx*' -or
        $lower.StartsWith('bindings/go/') -or
        $lower.StartsWith('bindings/c#/') -or
        $lower.StartsWith('bindings/emscripten/') -or
        $lower.StartsWith('bindings/java/') -or
        $lower.StartsWith('bindings/node.js/') -or
        $lower.StartsWith('bindings/python/') -or
        $lower -eq 'bindings/blst.swg' -or
        $lower -eq 'go.mod' -or
        $lower -eq '.golangci.yml'

    $cudaRemoval =
        $lower.StartsWith('blst/cuda/') -or
        $lower.StartsWith('blst/backend/cuda/') -or
        $lower -match '(^|/).*cuda[^/]*\.(c|cc|cpp|cuh|h|hpp)$' -or
        $lower -like 'test/*cuda*'

    if ($isBenchmark)
    {
        $disposition = 'out_of_scope_no_audit'
        $subsystem = 'benchmarks'
        $profiles = 'non_shipping'
        if ($lower -like 'benchmark/btx*')
        {
            $notes = 'Excluded from file-by-file audit; this BTx benchmark is also removed when BTx/BTE moves.'
        }
        elseif ($lower -match 'cuda')
        {
            $notes = 'Excluded from file-by-file audit; this CUDA benchmark is also removed with CUDA.'
        }
        else
        {
            $notes = 'Benchmarks are inventoried but explicitly excluded from file-by-file audit.'
        }
    }
    elseif ($plannedRemoval -or $cudaRemoval)
    {
        $disposition = 'out_of_scope_planned_removal'
        $profiles = 'planned_removal'

        if ($cudaRemoval)
        {
            $subsystem = 'cuda'
            $notes = 'Planned removal: CUDA is not a production target.'
        }
        elseif ($lower.StartsWith('bte/') -or $lower -match '(^|/)btx' -or
                $lower -eq 'bindings/blst_trx.h')
        {
            $subsystem = 'btx_bte'
            $notes = 'Planned removal: BTx/BTE will move to a separate repository.'
        }
        elseif ($lower.StartsWith('bindings/go/') -or $lower -eq 'go.mod' -or
                $lower -eq '.golangci.yml')
        {
            $subsystem = 'go_binding'
            $notes = 'Planned removal: Go is not a supported language surface.'
        }
        else
        {
            $subsystem = 'unsupported_binding'
            $notes = 'Planned removal: only C, C++, and Rust bindings are supported.'
        }
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
        foreach ($column in $preservedColumns)
        {
            if ($column -in @('audit_status', 'notes') -and
                [string]$old.disposition -ne [string]$classification.disposition)
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
