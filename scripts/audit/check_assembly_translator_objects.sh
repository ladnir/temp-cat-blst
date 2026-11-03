#!/bin/sh

set -eu

HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
ROOT=$(CDPATH= cd -- "$HERE/../.." && pwd)
CLANG=${CLANG:-clang}

if ! command -v "$CLANG" >/dev/null 2>&1; then
    echo "Clang executable not found: $CLANG" >&2
    exit 1
fi

if [ -n "${LLVM_READOBJ:-}" ]; then
    READOBJ=$LLVM_READOBJ
else
    READOBJ=
    for candidate in llvm-readobj llvm-readobj-20 llvm-readobj-19 \
                     llvm-readobj-18 llvm-readobj-17; do
        if command -v "$candidate" >/dev/null 2>&1; then
            READOBJ=$candidate
            break
        fi
    done
fi

if [ -z "$READOBJ" ] || ! command -v "$READOBJ" >/dev/null 2>&1; then
    echo "llvm-readobj executable not found" >&2
    exit 1
fi

if [ -n "${LLVM_NM:-}" ]; then
    NM=$LLVM_NM
else
    NM=
    for candidate in llvm-nm llvm-nm-20 llvm-nm-19 llvm-nm-18 llvm-nm-17; do
        if command -v "$candidate" >/dev/null 2>&1; then
            NM=$candidate
            break
        fi
    done
fi

if [ -z "$NM" ] || ! command -v "$NM" >/dev/null 2>&1; then
    echo "llvm-nm executable not found" >&2
    exit 1
fi

export LC_ALL=C

TEMP_ROOT=$(mktemp -d "${TMPDIR:-/tmp}/catblst-assembly-objects.XXXXXX")
cleanup()
{
    case "$TEMP_ROOT" in
        "${TMPDIR:-/tmp}"/catblst-assembly-objects.*)
            rm -rf -- "$TEMP_ROOT"
            ;;
        *)
            echo "refusing to remove unexpected temporary path: $TEMP_ROOT" >&2
            ;;
    esac
}
trap cleanup EXIT HUP INT TERM

printf '%s\n' '#define BLST_PORTABLE' '#define BLST_ASM' \
    '#define BLST_LIMB_SIZE 64' >"$TEMP_ROOT/config.h"

assemble()
{
    target=$1
    output=$2
    shift 2
    "$CLANG" --target="$target" -E -P "$@" \
        -I"$TEMP_ROOT" -I"$ROOT/build" \
        "$ROOT/blst/generated/assembly.S" >"$TEMP_ROOT/${output%.o}.s"
    "$CLANG" --target="$target" -Werror "$@" \
        -I"$TEMP_ROOT" -I"$ROOT/build" \
        -c "$ROOT/blst/generated/assembly.S" -o "$TEMP_ROOT/$output"
}

require()
{
    report=$1
    pattern=$2
    description=$3
    if ! grep -F "$pattern" "$report" >/dev/null; then
        echo "missing $description in $(basename -- "$report")" >&2
        exit 1
    fi
}

check_symbols()
{
    stem=$1
    sed -n -E \
        's/^\.(globl|global|comm)[[:space:]]+([^,[:space:]]+).*/\2/p' \
        "$TEMP_ROOT/$stem.s" | sort -u >"$TEMP_ROOT/$stem.expected"
    "$NM" --defined-only --extern-only --format=posix \
        "$TEMP_ROOT/$stem.o" | awk '{ print $1 }' | sort -u \
        >"$TEMP_ROOT/$stem.actual"
    if ! diff -u "$TEMP_ROOT/$stem.expected" "$TEMP_ROOT/$stem.actual"; then
        echo "object symbol surface differs from translated source: $stem" >&2
        exit 1
    fi
}

assemble x86_64-linux-gnu x86_64-elf.o
assemble x86_64-w64-windows-gnu x86_64-coff.o
assemble x86_64-apple-macos11 x86_64-macho.o
assemble arm64-apple-macos11 arm64-macho.o \
    -DBLST_APPLE_ARM64_SUPPORTED=1
assemble aarch64-linux-gnu arm64-elf.o \
    -DBLST_ALLOW_UNTESTED_ARM64=1
assemble aarch64-w64-windows-gnu arm64-coff.o \
    -DBLST_ALLOW_UNTESTED_ARM64=1

check_symbols x86_64-elf
check_symbols x86_64-coff
check_symbols x86_64-macho
check_symbols arm64-macho
check_symbols arm64-elf
check_symbols arm64-coff

"$READOBJ" --file-headers --sections --symbols --unwind \
    "$TEMP_ROOT/x86_64-elf.o" >"$TEMP_ROOT/x86_64-elf.txt"
require "$TEMP_ROOT/x86_64-elf.txt" "Format: elf64-x86-64" \
    "x86-64 ELF format"
require "$TEMP_ROOT/x86_64-elf.txt" "Name: .eh_frame" \
    "ELF unwind section"
require "$TEMP_ROOT/x86_64-elf.txt" "Name: add_mod_256" \
    "ELF arithmetic symbol"
require "$TEMP_ROOT/x86_64-elf.txt" \
    "Name: blst_sha256_block_data_order" "ELF SHA symbol"

"$READOBJ" --file-headers --sections --symbols --unwind \
    "$TEMP_ROOT/x86_64-coff.o" >"$TEMP_ROOT/x86_64-coff.txt"
require "$TEMP_ROOT/x86_64-coff.txt" "Format: COFF-x86-64" \
    "x86-64 COFF format"
require "$TEMP_ROOT/x86_64-coff.txt" "Name: .pdata" \
    "COFF procedure-data section"
require "$TEMP_ROOT/x86_64-coff.txt" "Name: .xdata" \
    "COFF unwind-data section"
require "$TEMP_ROOT/x86_64-coff.txt" "UnwindInformation [" \
    "decoded COFF unwind records"
require "$TEMP_ROOT/x86_64-coff.txt" "Name: add_mod_256" \
    "COFF arithmetic symbol"
require "$TEMP_ROOT/x86_64-coff.txt" \
    "Name: blst_sha256_block_data_order" "COFF SHA symbol"

"$READOBJ" --file-headers --sections --symbols \
    "$TEMP_ROOT/x86_64-macho.o" >"$TEMP_ROOT/x86_64-macho.txt"
require "$TEMP_ROOT/x86_64-macho.txt" "Format: Mach-O 64-bit x86-64" \
    "x86-64 Mach-O format"
require "$TEMP_ROOT/x86_64-macho.txt" "Name: __text" \
    "x86-64 Mach-O text section"
require "$TEMP_ROOT/x86_64-macho.txt" "Name: _add_mod_256" \
    "x86-64 Mach-O decorated arithmetic symbol"
require "$TEMP_ROOT/x86_64-macho.txt" \
    "Name: _blst_sha256_block_data_order" \
    "x86-64 Mach-O decorated SHA symbol"

"$READOBJ" --file-headers --sections --symbols \
    "$TEMP_ROOT/arm64-macho.o" >"$TEMP_ROOT/arm64-macho.txt"
require "$TEMP_ROOT/arm64-macho.txt" "Format: Mach-O arm64" \
    "ARM64 Mach-O format"
require "$TEMP_ROOT/arm64-macho.txt" "Name: __text" \
    "ARM64 Mach-O text section"
require "$TEMP_ROOT/arm64-macho.txt" "Name: _add_mod_256" \
    "ARM64 Mach-O decorated arithmetic symbol"
require "$TEMP_ROOT/arm64-macho.txt" \
    "Name: _blst_sha256_block_data_order" \
    "ARM64 Mach-O decorated SHA symbol"

"$READOBJ" --file-headers --sections --symbols \
    "$TEMP_ROOT/arm64-elf.o" >"$TEMP_ROOT/arm64-elf.txt"
require "$TEMP_ROOT/arm64-elf.txt" "Format: elf64-littleaarch64" \
    "ARM64 ELF format"
require "$TEMP_ROOT/arm64-elf.txt" "Name: .text" \
    "ARM64 ELF text section"
require "$TEMP_ROOT/arm64-elf.txt" "Name: add_mod_256" \
    "ARM64 ELF arithmetic symbol"
require "$TEMP_ROOT/arm64-elf.txt" \
    "Name: blst_sha256_block_data_order" "ARM64 ELF SHA symbol"

"$READOBJ" --file-headers --sections --symbols \
    "$TEMP_ROOT/arm64-coff.o" >"$TEMP_ROOT/arm64-coff.txt"
require "$TEMP_ROOT/arm64-coff.txt" "Format: COFF-ARM64" \
    "ARM64 COFF format"
require "$TEMP_ROOT/arm64-coff.txt" "Name: .text" \
    "ARM64 COFF text section"
require "$TEMP_ROOT/arm64-coff.txt" "Name: add_mod_256" \
    "ARM64 COFF arithmetic symbol"
require "$TEMP_ROOT/arm64-coff.txt" \
    "Name: blst_sha256_block_data_order" "ARM64 COFF SHA symbol"

echo "Assembly translator object checks passed for x86-64/ARM64 ELF, COFF, and Mach-O."
