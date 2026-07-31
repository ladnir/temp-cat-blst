#!/bin/sh

set -eu

if [ "$#" -ne 1 ]; then
    echo "usage: $0 OUTPUT_FILE" >&2
    exit 2
fi

HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
BINDINGS="$HERE/../bindings"
OUTPUT=$1
CC=${CC:-cc}
PERL=${PERL:-perl}

case "$OUTPUT" in
    /*) ;;
    *) OUTPUT="$(pwd)/$OUTPUT" ;;
esac

export LC_ALL=C
export TZ=UTC
umask 022

if ! command -v "$CC" >/dev/null 2>&1; then
    echo "C preprocessor not found: $CC" >&2
    exit 1
fi
if ! command -v "$PERL" >/dev/null 2>&1; then
    echo "Perl executable not found: $PERL" >&2
    exit 1
fi
if [ ! -d "$(dirname -- "$OUTPUT")" ]; then
    echo "output parent directory does not exist: $(dirname -- "$OUTPUT")" >&2
    exit 1
fi

TEMP_ROOT=$(mktemp -d "${TMPDIR:-/tmp}/catblst-exports.XXXXXX")
PREPROCESSED="$TEMP_ROOT/blst.i"
GENERATED="$TEMP_ROOT/blst.def"

cleanup()
{
    case "$TEMP_ROOT" in
        "${TMPDIR:-/tmp}"/catblst-exports.*) rm -rf -- "$TEMP_ROOT" ;;
        *) echo "refusing to remove unexpected temporary path: $TEMP_ROOT" >&2 ;;
    esac
}
trap cleanup EXIT HUP INT TERM

(
    cd "$BINDINGS"
    "$CC" -E blst.h
) >"$PREPROCESSED"

(
    printf 'LIBRARY blst\n\nEXPORTS\n'
    "$PERL" -ne \
        '{ (/(blst_[\w]+)\s*\(/ || /(BLS12_[\w]+);/) && print "\t$1\n" }' \
        "$PREPROCESSED"
    printf '\n'
) >"$GENERATED"

mv -f -- "$GENERATED" "$OUTPUT"
