#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "simint/simint.h"
#include "test/Common.hpp"

int main(int argc, char **argv) {
  int is_accelerator_available = 0;
  // Check to see if OpenMP target offloading works
#pragma omp target map(from : is_accelerator_available)
  {
    if (omp_is_initial_device() == 0)
      is_accelerator_available = 1;
  }
  if (is_accelerator_available)
    std::cout << "Accelerator is available!\n";
  else
    std::cout << "Accelerator is not available.\n";
  return 0;
}