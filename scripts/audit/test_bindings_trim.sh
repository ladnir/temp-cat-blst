#!/bin/sh

set -eu

HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
TRIM="$HERE/../../build/bindings_trim.pl"
TEMP_ROOT=$(mktemp -d "${TMPDIR:-/tmp}/catblst-trim-test.XXXXXX")
INPUT="$TEMP_ROOT/input.rs"
OUTPUT="$TEMP_ROOT/output.rs"
BROKEN="$TEMP_ROOT/broken.rs"

cleanup()
{
    case "$TEMP_ROOT" in
        "${TMPDIR:-/tmp}"/catblst-trim-test.*) rm -rf -- "$TEMP_ROOT" ;;
        *) echo "refusing to remove unexpected temporary path: $TEMP_ROOT" >&2 ;;
    esac
}
trap cleanup EXIT HUP INT TERM

cat >"$INPUT" <<'EOF'
#[derive(Debug, Default, Copy, Clone, PartialEq, Eq)]
pub struct blst_scalar {}
#[derive(Debug, Default, Copy, Clone, PartialEq, Eq)]
pub struct blst_fp12 {}
#[derive(Debug, Default, Copy, Clone, PartialEq, Eq)]
pub struct blst_p1 {}
#[derive(Debug, Default, Copy, Clone, PartialEq, Eq)]
pub struct blst_p1_affine {}
#[derive(Debug, Default, Copy, Clone, PartialEq, Eq)]
pub struct blst_p2 {}
#[derive(Debug, Default, Copy, Clone, PartialEq, Eq)]
pub struct blst_p2_affine {}
#[derive(Debug, Default, Copy, Clone, PartialEq, Eq)]
pub struct blst_pairing {}
#[derive(Debug, Default, Copy, Clone, PartialEq, Eq)]
pub struct blst_uniq {}
pub struct BLST_ERROR(pub u32);
pub struct BLST_UNIQ_RESULT(pub u32);
pub struct blst_fr {}
pub struct blst_fp {}
pub struct blst_fp2 {}
pub struct blst_fp6 {}
pub fn marker() -> ::std::primitive::u8 { 0 }
EOF

perl "$TRIM" "$INPUT" >"$OUTPUT"
grep -q '#\[zeroize(drop)\]' "$OUTPUT"
grep -q 'Clone, PartialEq, Eq, Zeroize' "$OUTPUT"
grep -q 'pub fn marker() -> ::core::primitive::u8' "$OUTPUT"
if perl -0777 -ne 'exit 1 if /derive\([^\n]*PartialEq[^\n]*\)\]\npub struct blst_p[12](?:_affine)?\b/' "$OUTPUT"; then
    :
else
    echo "point PartialEq derive was not removed" >&2
    exit 1
fi

sed '0,/, Copy/s///' "$INPUT" >"$BROKEN"
if perl "$TRIM" "$BROKEN" >/dev/null 2>&1; then
    echo "binding trim accepted a missing required transformation" >&2
    exit 1
fi

echo "Rust binding postprocessor fixture passed."
