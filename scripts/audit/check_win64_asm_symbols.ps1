[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [string]$BuildDirectory,

    [Parameter(Mandatory = $true)]
    [ValidateSet("Portable", "NonPortable")]
    [string]$Mode,

    [string]$Dumpbin
)

$ErrorActionPreference = "Stop"

$resolvedBuildDirectory = (Resolve-Path -LiteralPath $BuildDirectory).Path
if (-not (Test-Path -LiteralPath $resolvedBuildDirectory -PathType Container)) {
    throw "Build directory is not a directory: $BuildDirectory"
}

if ([string]::IsNullOrWhiteSpace($Dumpbin)) {
    $dumpbinCommand = Get-Command dumpbin.exe -ErrorAction Stop
    $Dumpbin = $dumpbinCommand.Source
} else {
    $Dumpbin = (Resolve-Path -LiteralPath $Dumpbin).Path
}

$groups = @(
    @{
        Consumer = "ctq_inverse_mod_384-x86_64.asm"
        Provider = "ctx_inverse_mod_384-x86_64.asm"
        Symbols = @(
            'ct_inverse_mod_384$1'
        )
    },
    @{
        Consumer = "mulq_mont_256-x86_64.asm"
        Provider = "mulx_mont_256-x86_64.asm"
        Symbols = @(
            'mul_mont_sparse_256$1',
            'sqr_mont_sparse_256$1',
            'from_mont_256$1',
            'redc_mont_256$1'
        )
    },
    @{
        Consumer = "mulq_mont_384-x86_64.asm"
        Provider = "mulx_mont_384-x86_64.asm"
        Symbols = @(
            'mul_mont_384x$1',
            'sqr_mont_384x$1',
            'mul_382x$1',
            'sqr_382x$1',
            'mul_384$1',
            'sqr_384$1',
            'redc_mont_384$1',
            'from_mont_384$1',
            'sgn0_pty_mont_384$1',
            'sgn0_pty_mont_384x$1',
            'mul_mont_384$1',
            'sqr_mont_384$1',
            'sqr_n_mul_mont_384$1',
            'sqr_n_mul_mont_383$1',
            'sqr_mont_382x$1'
        )
    }
)

$objects = @(Get-ChildItem -LiteralPath $resolvedBuildDirectory -Recurse -Filter "*.asm.obj")
if ($objects.Count -eq 0) {
    throw "No MASM objects found below: $resolvedBuildDirectory"
}

$actual = [System.Collections.Generic.List[string]]::new()
foreach ($object in $objects) {
    $output = & $Dumpbin /symbols $object.FullName
    if ($LASTEXITCODE -ne 0) {
        throw "dumpbin failed for: $($object.FullName)"
    }

    foreach ($line in $output) {
        if ($line -match '^\s*[0-9A-F]+\s+[0-9A-F]+\s+(?<section>\S+).*\sExternal\s+\|\s+(?<symbol>\S+\$1)\s*$') {
            $definition = if ($Matches.section -eq "UNDEF") { "UNDEF" } else { "DEFINED" }
            $actual.Add("$($object.BaseName)|$definition|$($Matches.symbol)")
        }
    }
}

$expected = [System.Collections.Generic.List[string]]::new()
if ($Mode -eq "Portable") {
    foreach ($group in $groups) {
        foreach ($symbol in $group.Symbols) {
            $expected.Add("$($group.Consumer)|UNDEF|$symbol")
            $expected.Add("$($group.Provider)|DEFINED|$symbol")
        }
    }
}

$actualSorted = @($actual | Sort-Object)
$expectedSorted = @($expected | Sort-Object)
$difference = @(Compare-Object -ReferenceObject $expectedSorted -DifferenceObject $actualSorted)
if ($difference.Count -ne 0) {
    Write-Error "Unexpected Win64 fallback symbol surface for mode $Mode."
    $difference | Format-Table -AutoSize | Out-String | Write-Error
    exit 1
}

Write-Output (
    "Verified {0} Win64 fallback symbol entries for {1} mode across {2} MASM objects." -f
    $actualSorted.Count, $Mode, $objects.Count
)
