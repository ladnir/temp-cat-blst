param(
    [string]$Root = (Resolve-Path (Join-Path $PSScriptRoot "..\.."))
)

$ErrorActionPreference = "Stop"

$ledgerPath = Join-Path $Root "docs\audit\PUBLIC_API_DOCUMENTATION.csv"
if (-not (Test-Path -LiteralPath $ledgerPath)) {
    throw "missing public API documentation ledger: $ledgerPath"
}

$rows = @(Import-Csv -LiteralPath $ledgerPath)
if ($rows.Count -eq 0) {
    throw "public API documentation ledger is empty"
}

$allowedScopes = @("audited", "focused")
$allowedStatuses = @(
    "documented",
    "documented_current_contract",
    "documented_focused",
    "documented_safe",
    "documented_safe_and_generated",
    "documented_safe_and_raw",
    "generated_c_contract",
    "deferred_containing_file",
    "deferred_raw_ffi"
)
$documentedStatuses = @(
    "documented",
    "documented_current_contract",
    "documented_focused",
    "documented_safe",
    "documented_safe_and_generated",
    "documented_safe_and_raw"
)

$ids = @{}
$errors = [System.Collections.Generic.List[string]]::new()

function Test-SourceAnchors {
    param(
        [pscustomobject]$Row,
        [string]$Surface
    )

    $sourceField = "${Surface}_source"
    $anchorField = "${Surface}_anchors"
    $statusField = "${Surface}_status"
    $status = $Row.$statusField

    if ($allowedStatuses -notcontains $status) {
        $errors.Add("$($Row.id): invalid $Surface status '$status'")
        return
    }

    $sources = @($Row.$sourceField -split ";" | Where-Object { $_ })
    $anchors = @($Row.$anchorField -split "\|" | Where-Object { $_ })
    if ($sources.Count -eq 0 -or $anchors.Count -eq 0) {
        $errors.Add("$($Row.id): $Surface source and anchors are required")
        return
    }

    $combined = ""
    foreach ($source in $sources) {
        $path = Join-Path $Root $source
        if (-not (Test-Path -LiteralPath $path)) {
            $errors.Add("$($Row.id): missing $Surface source '$source'")
            continue
        }
        $combined += "`n" + (Get-Content -LiteralPath $path -Raw)
    }

    foreach ($anchor in $anchors) {
        $index = $combined.IndexOf($anchor, [StringComparison]::Ordinal)
        if ($index -lt 0) {
            $errors.Add("$($Row.id): $Surface anchor '$anchor' was not found")
            continue
        }

        if ($documentedStatuses -contains $status) {
            $start = [Math]::Max(0, $index - 1600)
            $prefix = $combined.Substring($start, $index - $start)
            if ($prefix -notmatch '(?s)(/\*\*|///|//!)') {
                $errors.Add("$($Row.id): $Surface anchor '$anchor' lacks nearby declaration documentation")
            }
        }
    }
}

foreach ($row in $rows) {
    if (-not $row.id -or $ids.ContainsKey($row.id)) {
        $errors.Add("missing or duplicate ledger id '$($row.id)'")
        continue
    }
    $ids[$row.id] = $true

    if ($allowedScopes -notcontains $row.audit_scope) {
        $errors.Add("$($row.id): invalid audit scope '$($row.audit_scope)'")
    }

    Test-SourceAnchors -Row $row -Surface "c"
    Test-SourceAnchors -Row $row -Surface "cpp"
    Test-SourceAnchors -Row $row -Surface "rust"

    if ($row.rust_status -eq "generated_c_contract" -and
        $documentedStatuses -notcontains $row.c_status) {
        $errors.Add("$($row.id): generated Rust contract requires a documented C source")
    }
}

if ($errors.Count -ne 0) {
    $errors | ForEach-Object { Write-Error $_ }
    exit 1
}

Write-Host "Public API documentation ledger: $($rows.Count) audited/focused families accounted for."
Write-Host "Public API documentation check passed."
