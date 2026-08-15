#!/bin/sh
set -eu

PROJECT_ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
CLI="$PROJECT_ROOT/graph-writing"

fail() {
    printf 'FAIL: %s\n' "$1" >&2
    exit 1
}

assert_contains() {
    output=$1
    expected=$2
    message=$3
    printf '%s\n' "$output" | grep -F "$expected" >/dev/null || fail "$message"
}

help_output=$("$CLI" --help)
assert_contains "$help_output" 'graph-writing <command>' 'global help does not show usage'

for command in barplot scatter bartonella; do
    assert_contains "$help_output" "$command" "global help does not list $command"
    command_help=$("$CLI" "$command" --help 2>&1)
    assert_contains "$command_help" 'Usage:' "$command help does not show usage"
done

set +e
unknown_output=$("$CLI" unknown 2>&1)
unknown_status=$?
set -e
[ "$unknown_status" -eq 2 ] || fail "unknown command exits with $unknown_status instead of 2"
assert_contains "$unknown_output" 'Unknown command: unknown' 'unknown command error is unclear'

outside_dir=$(mktemp -d)
trap 'rm -rf "$outside_dir"' EXIT
for command in barplot scatter bartonella; do
    set +e
    outside_output=$(cd "$outside_dir" && "$CLI" "$command" 2>&1)
    outside_status=$?
    set -e
    [ "$outside_status" -eq 1 ] || fail "$command argument error exits with $outside_status instead of 1"
    assert_contains "$outside_output" 'Usage: Rscript' "$command cannot load the project when called outside its directory"
done

printf 'PASS: graph-writing CLI\n'
