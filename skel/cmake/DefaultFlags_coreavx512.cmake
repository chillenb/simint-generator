if("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")

  list(APPEND SIMINT_C_FLAGS "-xCORE-AVX512;-fma")
  list(APPEND SIMINT_TESTS_CXX_FLAGS "-xCORE-AVX512;-fma")

elseif("${CMAKE_C_COMPILER_ID}" MATCHES "GNU" OR
       "${CMAKE_C_COMPILER_ID}" MATCHES "Clang" OR
       "${CMAKE_C_COMPILER_ID}" MATCHES "IntelLLVM")


  list(APPEND SIMINT_C_FLAGS "-mavx512f;-mavx512cd;-mavx512bw;-mavx512dq;-mavx512vl;-mfma")
  list(APPEND SIMINT_TESTS_CXX_FLAGS "-mavx512f;-mavx512cd;-mavx512bw;-mavx512dq;-mavx512vl;-mfma")

else()

  message(FATAL_ERROR "Unsupported compiler")

endif()
