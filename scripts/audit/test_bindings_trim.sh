#!/bin/sh

set -eu

HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
TRIM="$HERE/../build/bindings_trim.pl"
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
#[doc = "Section contract.\n/\n/** Function contract."]
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
#[derive(Debug, Copy, Clone)]
pub struct blst_pairing_context {}
pub type blst_pairing = blst_pairing_context;
#[derive(Debug, Copy, Clone)]
pub struct blst_uniq_context {}
pub type blst_uniq = blst_uniq_context;
pub struct BLST_ERROR(pub u32);
pub struct BLST_UNIQ_RESULT(pub u32);
pub struct blst_fr {}
pub struct blst_fp {}
pub struct blst_fp2 {}
pub struct blst_fp6 {}
pub fn marker(p: *const ::std::os::raw::c_void) -> ::std::os::raw::c_uint { 0 }
EOF

perl "$TRIM" "$INPUT" >"$OUTPUT"
grep -q '#\[zeroize(drop)\]' "$OUTPUT"
grep -q 'Clone, PartialEq, Eq, Zeroize' "$OUTPUT"
grep -q 'pub fn marker(p: \*const ::core::ffi::c_void) -> u32' "$OUTPUT"
if grep -Fq '\n/\n/**' "$OUTPUT"; then
    echo "adjacent C documentation delimiters leaked into Rust docs" >&2
    exit 1
fi
for opaque in blst_pairing_context blst_uniq_context; do
    if grep -B1 "pub struct $opaque" "$OUTPUT" | grep -Eq 'Copy|Clone'; then
        echo "$opaque retained a duplicating derive" >&2
        exit 1
    fi
done
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
