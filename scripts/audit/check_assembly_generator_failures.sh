#!/bin/sh

set -eu

HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
ROOT=$(CDPATH= cd -- "$HERE/../.." && pwd)
PERL=${PERL:-perl}
TEMP_ROOT=$(mktemp -d "${TMPDIR:-/tmp}/catblst-generator-failure.XXXXXX")

cleanup()
{
    case "$TEMP_ROOT" in
        "${TMPDIR:-/tmp}"/catblst-generator-failure.*) rm -rf -- "$TEMP_ROOT" ;;
        *) echo "refusing to remove unexpected temporary path: $TEMP_ROOT" >&2 ;;
    esac
}
trap cleanup EXIT HUP INT TERM

printf '%s\n' '#!/usr/bin/env perl' 'exit 7;' >"$TEMP_ROOT/x86_64-xlate.pl"
cp "$ROOT/blst/asm/add_mod_256-x86_64.pl" "$TEMP_ROOT/"
if "$PERL" "$TEMP_ROOT/add_mod_256-x86_64.pl" elf \
        "$TEMP_ROOT/x86-output.s" >/dev/null 2>&1; then
    echo "x86 generator ignored translator failure" >&2
    exit 1
fi

printf '%s\n' '#!/usr/bin/env perl' 'exit 7;' >"$TEMP_ROOT/arm-xlate.pl"
cp "$ROOT/blst/asm/add_mod_256-armv8.pl" "$TEMP_ROOT/"
if "$PERL" "$TEMP_ROOT/add_mod_256-armv8.pl" linux64 \
        "$TEMP_ROOT/arm-output.S" >/dev/null 2>&1; then
    echo "ARM generator ignored translator failure" >&2
    exit 1
fi

for generator in "$ROOT"/blst/asm/*-armv8.pl \
                 "$ROOT"/blst/asm/*-x86_64.pl; do
    case "$generator" in
        *-xlate.pl) continue ;;
    esac
    if ! grep -F 'close STDOUT or die "assembly output failed:' \
            "$generator" >/dev/null; then
        echo "generator does not check translator completion: $generator" >&2
        exit 1
    fi
done

echo "Assembly generator failure propagation checks passed."
