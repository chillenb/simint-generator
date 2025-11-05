#include <cstdio>
#include <atomic>
#include <iostream>
#include <cmath>

#ifdef _OPENMP
  #include <omp.h>
#endif

#include "simint/simint.h"
#include "simint/boys/boys.h"
#include "test/Common.hpp"
#include "test/ValeevRef.hpp"


int main(int argc, char ** argv)
{
    // set up the function pointers
    simint_init();

    // parse command line
    if(argc != 2)
    {
        printf("Give me 1 argument! I got %d\n", argc-1);
        return 1;
    }

    // max n
    int nmax = atoi(argv[1]);
    if(nmax > 31)
    {
        printf("Max n too large! Got %d, max is 31\n", nmax);
        return 1;
    }

    const double step = 0.0001;
    const double end = BOYS_SHORTGRID_MAXX;
    double x = 0.0;

    long double *F_ref = new long double[nmax + 1];
    double *maxerr = new double[nmax + 1];
    double *maxrelerr = new double[nmax + 1];

    for(int n = 0; n <= nmax; n++)
    {
        maxerr[n] = 0.0;
        F_ref[n] = 0.0;
    }

    ValeevRef_Init();

    while(x < end)
    {
        Valeev_F(F_ref, nmax, (long double)x);

        for(int n = 0; n <= nmax - 1; n++)
        {
            double fn;
            if(x < BOYS_SHORTGRID_MAXX)
                fn = boys_F_taylor_single(x, n);
            else
                fn = boys_F_long_single(x, n);
            double err = (double) fabsl((long double) fn - F_ref[n]);
            double relerr = err / fabs((double)F_ref[n]);
            if(err > maxerr[n])
                maxerr[n] = err;
            if(relerr > maxrelerr[n])
                maxrelerr[n] = relerr;
        }
        x += step;
    }
    for(int n = 0; n < nmax; n++)
        printf("n=%2d max error = %.5e max rel error = %.5e\n", n, maxerr[n], maxrelerr[n]);

}