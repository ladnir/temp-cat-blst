#!/bin/sh

HERE=`dirname "$0"`
cd "${HERE}"

CREATED_BLST_LINK=
PACKAGE_TARGET=

cleanup()
{
    if [ -n "${CREATED_BLST_LINK}" ] && [ -h blst ]; then
        rm -f blst
    fi
    if [ -n "${PACKAGE_TARGET}" ]; then
        rm -rf "${PACKAGE_TARGET}"
    fi
}

trap cleanup 0 2

if [ ! -d blst ]; then
    ln -s ../.. blst
    CREATED_BLST_LINK=1
fi

# Package verification relocates the native source tree. Reusing a CMake cache
# generated for the checkout source path makes Cargo verification fail, so use
# an isolated target unless the caller deliberately supplied one.
if [ -z "${CARGO_TARGET_DIR:-}" ]; then
    PACKAGE_TARGET=`mktemp -d "${TMPDIR:-/tmp}/blst-publish.XXXXXX"`
    CARGO_TARGET_DIR="${PACKAGE_TARGET}"
    export CARGO_TARGET_DIR
fi

# --allow-dirty because the temporary blst symbolic link is not committed
cargo +stable publish --allow-dirty "$@"
