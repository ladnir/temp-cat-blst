#!/bin/sh

set -eu

if [ "$#" -ne 0 ]; then
    echo "usage: $0" >&2
    exit 2
fi

HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
MANIFEST="$HERE/assembly-output-manifest.txt"
PERL=${PERL:-perl}

export LC_ALL=C
export TZ=UTC

TEMP_ROOT=$(mktemp -d "${TMPDIR:-/tmp}/catblst-assembly.XXXXXX")
FIRST="$TEMP_ROOT/first"
SECOND="$TEMP_ROOT/second"

cleanup()
{
    case "$TEMP_ROOT" in
        "${TMPDIR:-/tmp}"/catblst-assembly.*) rm -rf -- "$TEMP_ROOT" ;;
        *) echo "refusing to remove unexpected temporary path: $TEMP_ROOT" >&2 ;;
    esac
}
trap cleanup EXIT HUP INT TERM

"$HERE/generate-assembly.sh" "$FIRST"
"$HERE/generate-assembly.sh" "$SECOND"

while IFS= read -r relative; do
    if ! cmp -s "$FIRST/$relative" "$SECOND/$relative"; then
        echo "nondeterministic generated output: $relative" >&2
        exit 1
    fi

    checked="$HERE/$relative"
    if [ ! -f "$checked" ]; then
        echo "missing checked-in generated output: build/$relative" >&2
        exit 1
    fi

    normalized="$TEMP_ROOT/checked/$relative"
    mkdir -p "$(dirname -- "$normalized")"
    "$PERL" -0777 -pe 's/\r\n/\n/g' "$checked" >"$normalized"
    if ! cmp -s "$FIRST/$relative" "$normalized"; then
        echo "generated output differs from build/$relative" >&2
        exit 1
    fi
done <"$MANIFEST"

CHECKED_LIST="$TEMP_ROOT/checked-manifest.txt"
(
    cd "$HERE"
    find cheri coff elf mach-o win64 -type f \
        \( -name '*.asm' -o -name '*.s' -o -name '*.S' \) -print |
        sort
) >"$CHECKED_LIST"

if ! diff -u "$MANIFEST" "$CHECKED_LIST"; then
    echo "checked-in generated assembly violates the output allowlist" >&2
    exit 1
fi

count=$(wc -l <"$MANIFEST")
echo "Verified $count deterministic assembly outputs against the checked-in tree."
