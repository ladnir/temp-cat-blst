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
    if ! grep -F 'open STDOUT, "|-", @xlate_args' \
            "$generator" >/dev/null; then
        echo "generator does not use an argument-list translator pipe: $generator" >&2
        exit 1
    fi
done

SPACED_ROOT="$TEMP_ROOT/path with spaces"
mkdir "$SPACED_ROOT"
cp "$ROOT/blst/asm/x86_64-xlate.pl" \
   "$ROOT/blst/asm/add_mod_256-x86_64.pl" "$SPACED_ROOT/"
"$PERL" "$SPACED_ROOT/add_mod_256-x86_64.pl" elf \
    "$SPACED_ROOT/x86 output.s"
test -s "$SPACED_ROOT/x86 output.s"

cp "$ROOT/blst/asm/arm-xlate.pl" \
   "$ROOT/blst/asm/add_mod_256-armv8.pl" "$SPACED_ROOT/"
"$PERL" "$SPACED_ROOT/add_mod_256-armv8.pl" linux64 \
    "$SPACED_ROOT/arm output.S"
test -s "$SPACED_ROOT/arm output.S"

if [ -c /dev/full ]; then
    if printf '.text\n' |
            "$PERL" "$ROOT/blst/asm/x86_64-xlate.pl" elf /dev/full \
                >/dev/null 2>&1; then
        echo "x86 translator ignored output failure" >&2
        exit 1
    fi

    if printf '.text\n' |
            "$PERL" "$ROOT/blst/asm/arm-xlate.pl" linux64 /dev/full \
                >/dev/null 2>&1; then
        echo "ARM translator ignored output failure" >&2
        exit 1
    fi
fi

for translator in "$ROOT/blst/asm/arm-xlate.pl" \
                  "$ROOT/blst/asm/x86_64-xlate.pl"; do
    if ! grep -F 'close STDOUT or die "assembly translation output failed:' \
            "$translator" >/dev/null; then
        echo "translator does not check output completion: $translator" >&2
        exit 1
    fi
done

echo "Assembly generator failure propagation checks passed."
