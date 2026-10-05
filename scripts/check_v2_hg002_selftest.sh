#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
RUNNER="$SCRIPT_DIR/check_v2_hg002.sh"

expect_error() {
    local expected=$1 output
    shift
    if output=$("$RUNNER" "$@" 2>&1); then
        printf 'expected failure containing %q\n' "$expected" >&2
        exit 1
    fi
    [[ $output == *"$expected"* ]] || {
        printf 'expected %q in error output, got:\n%s\n' "$expected" "$output" >&2
        exit 1
    }
}

"$RUNNER" --help >/dev/null
expect_error 'unknown option' --unexpected
expect_error 'positive integer' --threads 0
expect_error 'missing value' --v2-cli
printf 'check_v2_hg002 argument selftest passed\n'
