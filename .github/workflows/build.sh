#!/bin/bash

set -e

FANS_BUILD_DIR=${1:-build}

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

rm -rf ${FANS_BUILD_DIR}
mkdir ${FANS_BUILD_DIR}
cd ${FANS_BUILD_DIR}
cmake -DCMAKE_BUILD_TYPE=Release ..
configure=$(date +%s)

# run in parallel, but retry single threaded in case of failure such
# that error messages are more readable.
cmake --build . || cmake --build . -j1
build=$(date +%s)

echo "---------------------------------------------------"
echo "configure time: $((configure-start))s"
echo "build time:     $((build-configure))s"
echo "total time:     $((build-start))s"
echo "---------------------------------------------------"
