#include <cstdio>
#include <atomic>
#include <iostream>
#include <cmath>
#include <getopt.h>
#include <cstring>

#ifdef _OPENMP
  #include <omp.h>
#endif

#include "simint/simint.h"
#include "simint/boys/boys.h"
#include "test/Common.hpp"
#include "test/ValeevRef.hpp"

void usage()
{
    printf("Usage: test_boys [options]\n");
    printf("  -m [rational|taylor]    boys method (default: taylor)\n");
    printf("  -n nmax          max n (default: 20)\n");
    printf("  -h               display this help message\n");
}

int main(int argc, char ** argv)
{
    // set up the function pointers
    simint_init();

    int opt;
    int boys_method = 0;
    int nmax = 20;

    // parse command line
    while ((opt = getopt(argc, argv, "m:n:")) != -1)
    {
        switch (opt)
        {
        case 'm':
            if(strcmp(optarg, "rational") == 0)
                boys_method = 1;
            else if(strcmp(optarg, "taylor") == 0)
                boys_method = 0;
            else
            {
                printf("Option m must be either 'rational' or 'taylor'! Got '%s'\n", optarg);
                usage();
                return 1;
            }
            break;
        case 'n':
            nmax = atoi(optarg);
            break;
        case 'h':
            usage();
            return 0;
        case '?':
            printf("Unknown option: %c\n", optopt);
            usage();
            return 1;
        }
    }

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
            {
                if(boys_method == 1)
                    fn = boys_F_rational_single(x, n);
                else
                    fn = boys_F_taylor_single(x, n);
            }
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