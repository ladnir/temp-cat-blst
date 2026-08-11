[CmdletBinding()]
param()

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '..\..')).Path
$manifestPath = Join-Path $repoRoot 'docs\audit\FILE_MANIFEST.csv'
$upstreamManifestPath = Join-Path $repoRoot 'docs\audit\UPSTREAM_MANIFEST.md'
$allowedDispositions = @(
    'in_scope_production',
    'in_scope_support',
    'generated_production',
    'out_of_scope_frozen',
    'out_of_scope_experimental',
    'out_of_scope_planned_removal',
    'out_of_scope_no_audit',
    'documentation_or_archive'
)
$allowedStatuses = @(
    'pending',
    'in_review',
    'blocked',
    'complete',
    'accepted_risk',
    'frozen',
    'experimental',
    'scope_review_required',
    'excluded',
    'catalogued'
)
$requiredColumns = @(
    'path',
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
$errors = [System.Collections.Generic.List[string]]::new()

if (-not (Test-Path -LiteralPath $manifestPath))
{
    throw "Missing manifest: $manifestPath"
}

$rows = @(Import-Csv -LiteralPath $manifestPath)
if ($rows.Count -eq 0)
{
    throw "Manifest contains no rows: $manifestPath"
}

foreach ($column in $requiredColumns)
{
    if ($rows[0].PSObject.Properties.Name -notcontains $column)
    {
        $errors.Add("Missing required column '$column'.")
    }
}

$expectedPaths = [System.Collections.Generic.HashSet[string]]::new(
    [System.StringComparer]::Ordinal
)
foreach ($path in @(& git -C $repoRoot ls-files --cached --others --exclude-standard))
{
    if (-not [string]::IsNullOrWhiteSpace($path))
    {
        [void]$expectedPaths.Add($path.Replace('\', '/'))
    }
}
foreach ($path in @(& git -C $repoRoot ls-files --deleted))
{
    if (-not [string]::IsNullOrWhiteSpace($path))
    {
        [void]$expectedPaths.Remove($path.Replace('\', '/'))
    }
}

[void]$expectedPaths.Add('docs/audit/FILE_MANIFEST.csv')

$seenPaths = [System.Collections.Generic.HashSet[string]]::new(
    [System.StringComparer]::Ordinal
)
foreach ($row in $rows)
{
    $path = [string]$row.path
    if ([string]::IsNullOrWhiteSpace($path))
    {
        $errors.Add('Manifest contains an empty path.')
        continue
    }
    if (-not $seenPaths.Add($path))
    {
        $errors.Add("Duplicate manifest path: $path")
    }
    if (-not $allowedDispositions.Contains([string]$row.disposition))
    {
        $errors.Add("Invalid disposition '$($row.disposition)' for $path")
    }
    if (-not $allowedStatuses.Contains([string]$row.audit_status))
    {
        $errors.Add("Invalid audit status '$($row.audit_status)' for $path")
    }
    if ([string]::IsNullOrWhiteSpace([string]$row.subsystem))
    {
        $errors.Add("Missing subsystem for $path")
    }
    if ([string]::IsNullOrWhiteSpace([string]$row.profiles))
    {
        $errors.Add("Missing profiles for $path")
    }

    if ($row.disposition -in @(
            'out_of_scope_frozen',
            'out_of_scope_experimental',
            'out_of_scope_planned_removal',
            'out_of_scope_no_audit'
        ) -and
        [string]::IsNullOrWhiteSpace([string]$row.notes))
    {
        $errors.Add("Out-of-scope row lacks a reason: $path")
    }
    if ($row.disposition -eq 'out_of_scope_no_audit' -and
        $row.audit_status -ne 'excluded')
    {
        $errors.Add("No-audit row must have status 'excluded': $path")
    }
    if ($row.disposition -eq 'out_of_scope_frozen' -and
        $row.audit_status -ne 'frozen')
    {
        $errors.Add("Frozen row must have status 'frozen': $path")
    }
    if ($row.audit_status -eq 'frozen' -and
        $row.disposition -ne 'out_of_scope_frozen')
    {
        $errors.Add("Status 'frozen' requires disposition 'out_of_scope_frozen': $path")
    }
    if ($row.disposition -eq 'out_of_scope_experimental' -and
        $row.audit_status -ne 'experimental')
    {
        $errors.Add("Experimental row must have status 'experimental': $path")
    }
    if ($row.audit_status -eq 'experimental' -and
        $row.disposition -ne 'out_of_scope_experimental')
    {
        $errors.Add("Status 'experimental' requires disposition 'out_of_scope_experimental': $path")
    }

    if ($row.audit_status -in @('complete', 'accepted_risk'))
    {
        foreach ($column in @(
            'reviewer',
            'review_record',
            'upstream_reference',
            'findings',
            'tests',
            'validation',
            'residual_risk',
            'reviewed_commit'
        ))
        {
            if ([string]::IsNullOrWhiteSpace([string]$row.$column))
            {
                $errors.Add("$path is $($row.audit_status) but '$column' is empty.")
            }
        }

        if (-not ([string]$row.reviewed_commit -match '^[0-9a-fA-F]{40}$'))
        {
            $errors.Add("$path has a non-full reviewed commit: $($row.reviewed_commit)")
        }

        if (-not [string]::IsNullOrWhiteSpace([string]$row.review_record))
        {
            foreach ($reviewRecord in ([string]$row.review_record -split '\|'))
            {
                $reviewPath = Join-Path $repoRoot $reviewRecord.Replace('/', '\')
                if (-not (Test-Path -LiteralPath $reviewPath -PathType Leaf))
                {
                    $errors.Add("$path references missing review record: $reviewRecord")
                }
            }
        }

        if ($row.audit_status -eq 'accepted_risk' -and
            -not ([string]$row.findings -match 'AUD-\d{3}'))
        {
            $errors.Add("$path accepts risk without a finding ID.")
        }
    }
}

foreach ($path in $expectedPaths)
{
    if (-not $seenPaths.Contains($path))
    {
        $errors.Add("Repository file missing from manifest: $path")
    }
}
foreach ($path in $seenPaths)
{
    if (-not $expectedPaths.Contains($path))
    {
        $errors.Add("Manifest path does not exist in the repository inventory: $path")
    }
}

if (-not (Test-Path -LiteralPath $upstreamManifestPath))
{
    $errors.Add("Missing upstream manifest: $upstreamManifestPath")
}
else
{
    $upstreamText = Get-Content -LiteralPath $upstreamManifestPath -Raw
    foreach ($boundary in @(
        'f48500c1fdbefa7c0bf9800bccd65d28236799c1',
        '54e6e55674722fc2797ebb4bbb71b26d881eb4b8'
    ))
    {
        if (-not $upstreamText.Contains($boundary))
        {
            $errors.Add("Upstream manifest is missing boundary commit $boundary")
        }
    }

    $upstreamRows = @(
        Get-Content -LiteralPath $upstreamManifestPath |
            Where-Object { $_ -match '^\|\s*\d+\s*\|' }
    )
    if ($upstreamRows.Count -ne 63)
    {
        $errors.Add("Upstream manifest must contain 63 commit rows; found $($upstreamRows.Count).")
    }

    $seenUpstreamCommits = [System.Collections.Generic.HashSet[string]]::new(
        [System.StringComparer]::OrdinalIgnoreCase
    )
    for ($index = 0; $index -lt $upstreamRows.Count; $index++)
    {
        $cells = @($upstreamRows[$index].Split('|') | ForEach-Object { $_.Trim() })
        if ($cells.Count -lt 10)
        {
            $errors.Add("Malformed upstream manifest row: $($upstreamRows[$index])")
            continue
        }

        $number = 0
        if (-not [int]::TryParse($cells[1], [ref]$number) -or
            $number -ne ($index + 1))
        {
            $errors.Add("Upstream row sequence error near '$($cells[1])'.")
        }

        $commit = $cells[2].Trim('`')
        if ($commit -notmatch '^[0-9a-fA-F]{7,40}$')
        {
            $errors.Add("Invalid upstream commit identifier '$commit'.")
        }
        elseif (-not $seenUpstreamCommits.Add($commit))
        {
            $errors.Add("Duplicate upstream commit identifier '$commit'.")
        }

        if ($cells[4] -notin @('direct', 'conceptual', 'removal_candidate'))
        {
            $errors.Add("Invalid upstream relevance '$($cells[4])' for $commit.")
        }
        if ($cells[6] -notin @(
            'pending',
            'port',
            'translate',
            'already_equivalent',
            'not_applicable',
            'reject',
            'defer'
        ))
        {
            $errors.Add("Invalid upstream decision '$($cells[6])' for $commit.")
        }
        if ($cells[7] -notin @(
            'not_started',
            'reviewing',
            'implemented',
            'validated',
            'closed'
        ))
        {
            $errors.Add("Invalid upstream status '$($cells[7])' for $commit.")
        }
        if ($cells[7] -eq 'closed' -and $cells[6] -eq 'pending')
        {
            $errors.Add("Closed upstream row still has a pending decision: $commit")
        }
        if ($cells[7] -eq 'closed' -and [string]::IsNullOrWhiteSpace($cells[8]))
        {
            $errors.Add("Closed upstream row lacks evidence: $commit")
        }
    }
}

if ($errors.Count -gt 0)
{
    Write-Error ("Audit manifest validation failed with {0} error(s):`n - {1}" -f
        $errors.Count, ($errors -join "`n - "))
    exit 1
}

Write-Output "Audit manifest validation passed: $($rows.Count) files accounted for."
Write-Output "Upstream manifest validation passed: 63 commits accounted for."
foreach ($group in @($rows | Group-Object disposition | Sort-Object Name))
{
    Write-Output ("  {0}: {1}" -f $group.Name, $group.Count)
}
foreach ($group in @($rows | Group-Object audit_status | Sort-Object Name))
{
    Write-Output ("  status/{0}: {1}" -f $group.Name, $group.Count)
}
