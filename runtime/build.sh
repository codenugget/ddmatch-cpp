#!/bin/bash
set -x
pushd .
cd ..
if ! [ -d build ]; then
  mkdir -p build
  cd build
  cmake ..
else
  cd build
fi
make
popd
