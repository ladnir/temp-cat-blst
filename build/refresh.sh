#!/bin/sh

set -eu

HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)

usage()
{
    cat >&2 <<EOF
usage:
  $0 generate-assembly OUTPUT_DIRECTORY
  $0 verify-assembly
  $0 generate-win64-exports OUTPUT_FILE
  $0 generate-rust-bindings OUTPUT_FILE
EOF
    exit 2
}

if [ "$#" -lt 1 ]; then
    usage
fi

operation=$1
shift

case "$operation" in
    generate-assembly)
        [ "$#" -eq 1 ] || usage
        exec "$HERE/generate-assembly.sh" "$1"
        ;;
    verify-assembly)
        [ "$#" -eq 0 ] || usage
        exec "$HERE/verify-assembly.sh"
        ;;
    generate-win64-exports)
        [ "$#" -eq 1 ] || usage
        exec "$HERE/generate-win64-exports.sh" "$1"
        ;;
    generate-rust-bindings)
        [ "$#" -eq 1 ] || usage
        exec "$HERE/generate-rust-bindings.sh" "$1"
        ;;
    *)
        usage
        ;;
esac
