#!/bin/sh

set -eu

if [ "$#" -ne 1 ]; then
    echo "usage: $0 OUTPUT_DIRECTORY" >&2
    exit 2
fi

HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
ASM_DIR="$HERE/../blst/asm"
MANIFEST="$HERE/assembly-output-manifest.txt"
OUTPUT=$1
PERL=${PERL:-perl}

export LC_ALL=C
export TZ=UTC
umask 022

if ! command -v "$PERL" >/dev/null 2>&1; then
    echo "Perl executable not found: $PERL" >&2
    exit 1
fi

if [ -e "$OUTPUT" ]; then
    if [ ! -d "$OUTPUT" ]; then
        echo "output exists and is not a directory: $OUTPUT" >&2
        exit 1
    fi
    if find "$OUTPUT" -mindepth 1 -print -quit | grep . >/dev/null 2>&1; then
        echo "output directory must be empty: $OUTPUT" >&2
        exit 1
    fi
else
    mkdir -p "$OUTPUT"
fi

for subdir in cheri coff elf mach-o win64; do
    mkdir "$OUTPUT/$subdir"
done

SOURCE_DIR="$OUTPUT/.generator-source"
mkdir "$SOURCE_DIR"
for source in "$ASM_DIR"/*.pl; do
    "$PERL" -0777 -pe 's/\r\n/\n/g' "$source" \
        >"$SOURCE_DIR/$(basename -- "$source")"
done

generate()
{
    generator=$1
    flavour=$2
    destination=$3
    raw="$OUTPUT/.raw-assembly-output"

    CC=false "$PERL" "$generator" "$flavour" >"$raw"
    if [ ! -s "$raw" ]; then
        echo "generator produced no output: $generator ($flavour)" >&2
        exit 1
    fi
    "$PERL" -0777 -pe 's/\r\n/\n/g' "$raw" >"$destination"
    rm -f "$raw"
}

for generator in "$SOURCE_DIR"/*-x86_64.pl; do
    name=$(basename -- "$generator" .pl)
    if [ "$name" = x86_64-xlate ]; then
        continue
    fi

    case "$name" in
        *portable*) ;;
        *) generate "$generator" masm "$OUTPUT/win64/$name.asm" ;;
    esac
    generate "$generator" elf "$OUTPUT/elf/$name.s"
    generate "$generator" mingw64 "$OUTPUT/coff/$name.s"
    generate "$generator" macosx "$OUTPUT/mach-o/$name.s"
done

for generator in "$SOURCE_DIR"/*-armv8.pl; do
    name=$(basename -- "$generator" .pl)
    generate "$generator" win64 "$OUTPUT/win64/$name.asm"
    generate "$generator" linux64 "$OUTPUT/elf/$name.S"
    generate "$generator" coff64 "$OUTPUT/coff/$name.S"
    generate "$generator" ios64 "$OUTPUT/mach-o/$name.S"
    generate "$generator" cheri64 "$OUTPUT/cheri/$name.S"
done

rm -rf -- "$SOURCE_DIR"

ACTUAL="$OUTPUT/.assembly-output-manifest.txt"
(
    cd "$OUTPUT"
    find cheri coff elf mach-o win64 -type f \
        \( -name '*.asm' -o -name '*.s' -o -name '*.S' \) -print |
        sort
) >"$ACTUAL"

if ! diff -u "$MANIFEST" "$ACTUAL"; then
    echo "generated assembly does not match the output allowlist" >&2
    exit 1
fi

rm -f "$ACTUAL"
