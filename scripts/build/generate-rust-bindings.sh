#!/bin/sh

set -eu

if [ "$#" -ne 1 ]; then
    echo "usage: $0 OUTPUT_FILE" >&2
    exit 2
fi

HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
BINDINGS="$HERE/../../bindings"
OUTPUT=$1
BINDGEN=${BINDGEN:-bindgen}
RUSTFMT=${RUSTFMT:-rustfmt}
PERL=${PERL:-perl}
BINDGEN_VERSION=0.65.1
RUSTFMT_VERSION='rustfmt 1.5.1-stable (897e375 2022-11-02)'

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
if ! command -v "$RUSTFMT" >/dev/null 2>&1; then
    echo "Rust formatter executable not found: $RUSTFMT" >&2
    exit 1
fi
ACTUAL_RUSTFMT_VERSION=$($RUSTFMT --version 2>/dev/null || true)
if [ "$ACTUAL_RUSTFMT_VERSION" != "$RUSTFMT_VERSION" ]; then
    echo "Rust formatter $RUSTFMT_VERSION is required; found: ${ACTUAL_RUSTFMT_VERSION:-unknown}" >&2
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
UNFORMATTED="$TEMP_ROOT/bindings.unformatted.rs"
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
        --formatter none \
        --allowlist-function '^blst_.*' \
        --allowlist-type '^(blst_.*|BLST_.*)' \
        --blocklist-type '^blst_.*_must_be_4_bytes$' \
        --with-derive-default \
        --with-derive-eq \
        --newtype-enum 'BLST.*' \
        blst_scalar.h -- -std=c99 "-I$HERE/../.."
) >"$UNFORMATTED"

cp -- "$UNFORMATTED" "$RAW"
"$RUSTFMT" --edition 2018 --config-path "$BINDINGS/rust/rustfmt.toml" "$RAW"

"$PERL" "$HERE/bindings_trim.pl" "$RAW" >"$GENERATED"
mv -f -- "$GENERATED" "$OUTPUT"
