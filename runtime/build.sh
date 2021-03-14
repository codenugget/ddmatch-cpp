#!/bin/bash
set -x
pushd .
cd ..
if ! [ -d build ]; then
  mkdir -p build
  cd build
  cmake .. -DCMAKE_BUILD_TYPE=Debug
else
  cd build
fi
make
popd
