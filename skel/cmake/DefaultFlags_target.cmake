# You may need to pass in extra flags when using OpenMP target offloading
# E.g. if you have clang and an NVIDIA Maxwell card, you might use the following flags:
# -DSIMINT_OMP_TARGET_FLAGS="-fopenmp -fopenmp-targets=nvptx64 -Xopenmp-target -march=sm_50 --cuda-path=/opt/cuda"
list(APPEND SIMINT_C_FLAGS "${SIMINT_OMP_TARGET_FLAGS}")
list(APPEND SIMINT_CXX_FLAGS "${SIMINT_OMP_TARGET_FLAGS}")
list(APPEND SIMINT_LINK_FLAGS "${SIMINT_OMP_TARGET_FLAGS}")

list(APPEND SIMINT_TESTS_CXX_FLAGS "${SIMINT_OMP_TARGET_FLAGS}")
list(APPEND SIMINT_TESTS_LINK_FLAGS "${SIMINT_OMP_TARGET_FLAGS}")
add_compile_definitions(SIMINT_TARGET)
