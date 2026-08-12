#!/bin/sh

set -eu

CREATED_BLST_LINK=
PACKAGE_TARGET=
CARGO=${CARGO:-cargo}

cleanup()
{
    if [ -n "${CREATED_BLST_LINK}" ] && [ -h blst ]; then
        rm -f blst || :
    fi
    if [ -n "${PACKAGE_TARGET}" ]; then
        rm -rf "${PACKAGE_TARGET}" || :
    fi
}

trap cleanup 0
trap 'exit 129' 1
trap 'exit 130' 2
trap 'exit 143' 15

HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd -P)
cd -- "${HERE}"

REPOSITORY_ROOT=$(git rev-parse --show-toplevel)
if [ -n "$(git -C "${REPOSITORY_ROOT}" status --porcelain --untracked-files=no)" ]; then
    echo "refusing to publish with tracked or staged changes" >&2
    exit 1
fi

if [ ! -d blst ]; then
    ln -s ../.. blst
    CREATED_BLST_LINK=1
elif [ "$(CDPATH= cd -- blst && pwd -P)" != "${REPOSITORY_ROOT}" ]; then
    echo "refusing to publish with blst mapped outside the repository root" >&2
    exit 1
fi

# Package verification relocates the native source tree. Reusing a CMake cache
# generated for the checkout source path makes Cargo verification fail, so use
# an isolated target unless the caller deliberately supplied one.
if [ -z "${CARGO_TARGET_DIR:-}" ]; then
    PACKAGE_TARGET=$(mktemp -d "${TMPDIR:-/tmp}/blst-publish.XXXXXX")
    CARGO_TARGET_DIR="${PACKAGE_TARGET}"
    export CARGO_TARGET_DIR
fi

# --allow-dirty because the temporary blst symbolic link is not committed
"${CARGO}" +stable publish --allow-dirty "$@"
