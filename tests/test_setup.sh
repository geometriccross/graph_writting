#!/bin/sh
set -eu

PROJECT_ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
SETUP="$PROJECT_ROOT/setup"
TEST_DIR=$(mktemp -d)
trap 'rm -rf "$TEST_DIR"' EXIT

fail() {
    printf 'FAIL: %s\n' "$1" >&2
    exit 1
}

assert_file_contains() {
    file=$1
    expected=$2
    message=$3
    grep -F -- "$expected" "$file" >/dev/null || fail "$message"
}

mkdir -p "$TEST_DIR/bin"
cat > "$TEST_DIR/bin/Rscript" <<'EOF'
#!/bin/sh
printf '%s\n' "$@" > "$SETUP_TEST_ARGS"
cat > "$SETUP_TEST_PROGRAM"
EOF
chmod +x "$TEST_DIR/bin/Rscript"

SETUP_TEST_ARGS="$TEST_DIR/args" \
SETUP_TEST_PROGRAM="$TEST_DIR/program.R" \
PATH="$TEST_DIR/bin:/bin:/usr/bin" \
    "$SETUP"

assert_file_contains "$TEST_DIR/args" '--vanilla' 'setup does not isolate bootstrap from user profiles'
assert_file_contains "$TEST_DIR/args" "$PROJECT_ROOT" 'setup does not pass its project root to R'
assert_file_contains "$TEST_DIR/program.R" 'renv::activate' 'setup does not activate renv'
assert_file_contains "$TEST_DIR/program.R" 'renv::restore' 'setup does not restore renv.lock'
assert_file_contains "$TEST_DIR/program.R" 'library = project_library' 'setup does not isolate restore to the project library'
assert_file_contains "$TEST_DIR/program.R" 'requireNamespace' 'setup does not verify runtime packages'

load_line=$(grep -n -m 1 'renv::load' "$TEST_DIR/program.R" | cut -d: -f1)
restore_line=$(grep -n -m 1 'renv::restore' "$TEST_DIR/program.R" | cut -d: -f1)
[ "$load_line" -lt "$restore_line" ] || fail 'setup restores before loading the project library'

printf 'PASS: setup\n'
