#!/usr/bin/env bash
if [ -e "build" ]; then
  rm -rf build
fi
mkdir build; cd build


CC=gcc CXX=g++ cmake -DGENERATOR_CXX_FLAGS="-O3;-march=native;-flto;-fno-strict-aliasing" ../
make -j16 VERBOSE=1