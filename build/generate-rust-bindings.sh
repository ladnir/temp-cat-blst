#!/bin/sh

set -eu

if [ "$#" -ne 1 ]; then
    echo "usage: $0 OUTPUT_FILE" >&2
    exit 2
fi

HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
BINDINGS="$HERE/../bindings"
OUTPUT=$1
BINDGEN=${BINDGEN:-bindgen}
PERL=${PERL:-perl}
BINDGEN_VERSION=0.65.1

case "$OUTPUT" in
    /*) ;;
    *) OUTPUT="$(pwd)/$OUTPUT" ;;
esac

export LC_ALL=C
export TZ=UTC
umask 022

if ! command -v "$BINDGEN" >/dev/null 2>&1; then
    echo "Rust bindgen executable not found: $BINDGEN" >&2
    exit 1
fi
ACTUAL_BINDGEN_VERSION=$($BINDGEN --version 2>/dev/null || true)
if [ "$ACTUAL_BINDGEN_VERSION" != "bindgen $BINDGEN_VERSION" ]; then
    echo "Rust bindgen $BINDGEN_VERSION is required; found: ${ACTUAL_BINDGEN_VERSION:-unknown}" >&2
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

OUTPUT_DIR=$(dirname -- "$OUTPUT")
TEMP_ROOT=$(mktemp -d "$OUTPUT_DIR/.catblst-bindings.XXXXXX")
RAW="$TEMP_ROOT/bindings.raw.rs"
GENERATED="$TEMP_ROOT/bindings.rs"

cleanup()
{
    case "$TEMP_ROOT" in
        "$OUTPUT_DIR"/.catblst-bindings.*) rm -rf -- "$TEMP_ROOT" ;;
        *) echo "refusing to remove unexpected temporary path: $TEMP_ROOT" >&2 ;;
    esac
}
trap cleanup EXIT HUP INT TERM

(
    cd "$BINDINGS"
    "$BINDGEN" \
        --opaque-type blst_pairing \
        --opaque-type blst_uniq \
        --with-derive-default \
        --with-derive-eq \
        --newtype-enum 'BLST.*' \
        blst_scalar.h -- -std=c99
) >"$RAW"

"$PERL" "$HERE/bindings_trim.pl" "$RAW" >"$GENERATED"
mv -f -- "$GENERATED" "$OUTPUT"
