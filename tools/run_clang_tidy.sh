#!/bin/bash
#
# Run clang-tidy over the FANS sources. Arguments go to run-clang-tidy:
#
#   ./tools/run_clang_tidy.sh
#   ./tools/run_clang_tidy.sh -fix
#
set -euo pipefail

# cd rather than just resolve: clang-tidy looks for .clang-tidy relative to the
# working directory, and reports "No checks enabled" without it.
cd "$(dirname "$0")/.."
SRC=$PWD
BUILD="$SRC/build/clang-tidy"

# The test targets are fixtures, not code under review.
CC=clang CXX=clang++ cmake -S "$SRC" -B "$BUILD" \
    -D CMAKE_EXPORT_COMPILE_COMMANDS=ON \
    -D FANS_ENABLE_TESTING=OFF >/dev/null

# No -header-filter here: it would override HeaderFilterRegex in .clang-tidy.
run-clang-tidy -p "$BUILD" -quiet "$@"
