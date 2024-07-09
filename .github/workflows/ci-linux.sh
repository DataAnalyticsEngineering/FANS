#!/bin/bash

set -e

echo "Running CI linux"

if [ -z "$FANS_DIR" ]; then
    echo "FANS_DIR is not set"
    exit 1
fi

echo "---------------------------------------------------"
ARCH=$(uname -m)
NPROC=$(nproc)
echo "arch=$ARCH"
echo "Processors: ${NPROC}"
echo "---------------------------------------------------"

CMAKE_BUILD_PARALLEL_LEVEL=${NPROC}
CTEST_PARALLEL_LEVEL=${NPROC}

cd "$FANS_DIR" || exit

start=$(date +%s)

cmake --preset ci-linux
configure=$(date +%s)

# build and test in parallel, but retry single threaded in case of failure such
# that error messages are more readable.
cmake --build --preset ci-linux || cmake --build --preset ci-linux -j1
build=$(date +%s)

ctest --preset ci-linux || ctest --preset ci-linux -j1
test=$(date +%s)

cd build/release || exit
cpack -G "DEB"
packaging=$(date +%s)

echo "---------------------------------------------------"
echo "configure time: $((configure-start))s"
echo "build time:     $((build-configure))s"
echo "test time:      $((test-build))s"
echo "packaging time: $((packaging-test))s"
echo "total time:     $((test-start))s"
echo "---------------------------------------------------"

echo "CI linux finished"
