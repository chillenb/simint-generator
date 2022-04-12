#!/usr/bin/env bash
if [ -e "build" ]; then
  rm -r build
fi
mkdir build; cd build
_SF="-g -fopenmp -fopenmp-targets=nvptx64-nvidia-cuda -fopenmp-new-driver -Xopenmp-target -march=sm_86 --cuda-path=/opt/cuda/"
CC=clang CXX=clang++ cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=1 -DSIMINT_MAXAM=2 -DSIMINT_VECTOR=target -DSIMINT_OMP_TARGET_FLAGS="${_SF}" -DSIMINT_LIBRARY_TYPE="OBJECT" ../
make -j16 VERBOSE=1