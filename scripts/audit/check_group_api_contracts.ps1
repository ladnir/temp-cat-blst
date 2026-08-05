param(
    [string]$RepositoryRoot = ""
)

$ErrorActionPreference = "Stop"

if ([string]::IsNullOrWhiteSpace($RepositoryRoot)) {
    $RepositoryRoot = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
} else {
    $RepositoryRoot = (Resolve-Path $RepositoryRoot).Path
}

$contractPath = Join-Path $RepositoryRoot "docs\audit\GROUP_API_CONTRACTS.csv"
if (-not (Test-Path -LiteralPath $contractPath)) {
    throw "Missing group API contract inventory: $contractPath"
}

$rows = @(Import-Csv -LiteralPath $contractPath)
if ($rows.Count -eq 0) {
    throw "Group API contract inventory is empty"
}

$requiredColumns = @(
    "family",
    "surface",
    "source",
    "match_kind",
    "symbol_pattern",
    "input_domain",
    "release_validation",
    "output_domain",
    "trust_boundary",
    "status"
)

$families = @{}
foreach ($row in $rows) {
    foreach ($column in $requiredColumns) {
        if ([string]::IsNullOrWhiteSpace($row.$column)) {
            throw "Contract row '$($row.family)' has an empty '$column' field"
        }
    }
    if ($families.ContainsKey($row.family)) {
        throw "Duplicate contract family '$($row.family)'"
    }
    $families[$row.family] = $true
    if ($row.match_kind -ne "regex") {
        throw "Unsupported match_kind '$($row.match_kind)' in '$($row.family)'"
    }
    if (-not ($row.symbol_pattern.StartsWith("^") -and $row.symbol_pattern.EndsWith('$'))) {
        throw "Contract family '$($row.family)' must use an anchored regular expression"
    }
}

function Get-PublicFunctionNames([string]$relativePath) {
    $path = Join-Path $RepositoryRoot $relativePath
    $text = Get-Content -Raw -LiteralPath $path
    $matches = [regex]::Matches(
        $text,
        '(?ms)^[A-Z_a-z][^;{}]*?\b(blst_[A-Za-z0-9_]+)\s*\([^;{}]*?\);'
    )
    return @($matches | ForEach-Object { $_.Groups[1].Value } | Sort-Object -Unique)
}

function Test-IsGroupDomainSymbol([string]$name) {
    return (
        $name -match '^blst_(avx_)?(p1|p2|fp12|gt_)' -or
        $name -match '^blst_(p1|p2|fp12)_from_avx$' -or
        $name -match '^blst_(final_exp|bendian_from_fp12|miller_loop|precompute_lines|pairing_|aggregate_in_g|aggregated_in_g|core_verify_pk|encode_to_g|hash_to_g|map_to_g|sk_to_pk|sign_pk|avx_enabled|avx_pairing)'
    )
}

$sources = @($rows | Select-Object -ExpandProperty source -Unique)
$coveredCount = 0
foreach ($source in $sources) {
    $sourceRows = @($rows | Where-Object { $_.source -eq $source })
    $symbols = @(Get-PublicFunctionNames $source | Where-Object { Test-IsGroupDomainSymbol $_ })
    if ($symbols.Count -eq 0) {
        throw "No group-domain public symbols found in '$source'"
    }

    foreach ($symbol in $symbols) {
        $matches = @($sourceRows | Where-Object { $symbol -match $_.symbol_pattern })
        if ($matches.Count -eq 0) {
            throw "Public group-domain symbol missing from contract inventory: ${source}::$symbol"
        }
        if ($matches.Count -ne 1) {
            $ids = ($matches | Select-Object -ExpandProperty family) -join ", "
            throw "Public group-domain symbol matches multiple contract families: ${source}::$symbol => $ids"
        }
        ++$coveredCount
    }

    foreach ($row in $sourceRows) {
        $rowMatches = @($symbols | Where-Object { $_ -match $row.symbol_pattern })
        if ($rowMatches.Count -eq 0) {
            throw "Contract family '$($row.family)' matches no public symbol in '$source'"
        }
    }
}

$forbiddenApiPatterns = @(
    'blst_avx_p[12]s?_unchecked_mult_pippenger',
    'blst_avx_fp12_mult_pippenger'
)
$sourceFiles = @(
    Get-ChildItem -LiteralPath (Join-Path $RepositoryRoot "bindings") -Recurse -File |
        Where-Object { $_.FullName -notmatch '[\\/]target[\\/]' -and $_.Extension -in '.h', '.cpp', '.rs' }
    Get-ChildItem -LiteralPath (Join-Path $RepositoryRoot "blst") -Recurse -File |
        Where-Object { $_.Extension -in '.h', '.cpp' }
    Get-ChildItem -LiteralPath (Join-Path $RepositoryRoot "test") -Recurse -File |
        Where-Object { $_.Extension -in '.h', '.cpp' }
)
foreach ($pattern in $forbiddenApiPatterns) {
    $matches = @($sourceFiles | Select-String -Pattern $pattern)
    if ($matches.Count -ne 0) {
        $locations = ($matches | ForEach-Object {
            "$($_.Path):$($_.LineNumber)"
        }) -join ", "
        throw "Removed group API spelling '$pattern' remains in source: $locations"
    }
}

# AUD-103: the selected scalar GT API promises secret-bit-safe conditional
# negation. Keep the implementation tied to the audited branchless selection
# primitive. The variable-time MSM path may retain its public-digit branch, but
# it must do so locally rather than weakening GT::cneg's contract.
$scalarGtText = Get-Content -Raw -LiteralPath (
    Join-Path $RepositoryRoot "blst/scalar/curve/GT.h")
$scalarCneg = [regex]::Match(
    $scalarGtText,
    '(?s)static\s+BLST_HD\s+void\s+cneg\s*\([^)]*\)\s*\{(?<body>.*?)\n\s*\}')
if (-not $scalarCneg.Success) {
    throw "Could not locate scalar GT::cneg for AUD-103 structural validation"
}
$scalarCnegBody = $scalarCneg.Groups['body'].Value
if ($scalarCnegBody -match '\bif\s*\(' -or
    $scalarCnegBody -notmatch '\bccopy\s*\(') {
    throw "Scalar GT::cneg must use branchless ccopy and contain no condition branch"
}

$scalarMsmText = Get-Content -Raw -LiteralPath (
    Join-Path $RepositoryRoot "blst/scalar/curve/multi_scalar.h")
$scalarGtMsmHelper = [regex]::Match(
    $scalarMsmText,
    '(?s)msm_add_or_sub\s*\(GT&\s+out,.*?\n\s*\}')
if (-not $scalarGtMsmHelper.Success -or
    $scalarGtMsmHelper.Value -notmatch '\bif\s*\(subtract\)' -or
    $scalarGtMsmHelper.Value -match 'GT::cneg') {
    throw "Scalar variable-time GT MSM must retain its local public-digit branch"
}

Write-Output "Group API contract validation passed: $coveredCount public C symbols covered by $($rows.Count) disjoint contract families."
