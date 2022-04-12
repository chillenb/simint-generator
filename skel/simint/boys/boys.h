#pragma once

#include "simint/vectorization/vectorization.h"

#include "simint/boys/boys_taylor.h"
#include "simint/boys/boys_shortgrid.h"
#include "simint/boys/boys_long.h"

#ifdef __cplusplus
#include "simint/cpp_restrict.hpp"
extern "C" {
#endif

#ifdef SIMINT_TARGET
static inline
void boys_F_small_n_target(SIMINT_DBLTYPE * restrict F,
                          SIMINT_DBLTYPE x,
                          int n);
#pragma omp declare variant( boys_F_small_n_target ) match(device={kind(nohost)})
#endif // SIMINT_TARGET
static inline
void boys_F_split_small_n(SIMINT_DBLTYPE * restrict F,
                          SIMINT_DBLTYPE x,
                          int n)
{
    // n is small - just do it all of them via
    // lookup or longfac (no recursion)
    #ifndef SIMINT_BOYS_NOVECTOR
    if(vector_min(x) > BOYS_SHORTGRID_MAXX)
        boys_F_long_vec(F, x, n);
    else if(vector_max(x) < BOYS_SHORTGRID_MAXX)
        boys_F_taylor_vec(F, x, n);
    else
    #endif
    {
        double * restrict Fd = (double *)F;
        double const * restrict xd = (double *)(&x);

        for(int i = 0; i < SIMINT_SIMD_LEN; i++)
        {
            if(xd[i] < BOYS_SHORTGRID_MAXX)
                boys_F_taylor(Fd + i, xd[i], n);
            else
                boys_F_long(Fd + i, xd[i], n);
        }
    }
}


static inline
void boys_F_split_large_n(SIMINT_DBLTYPE * restrict F,
                          SIMINT_DBLTYPE x,
                          int n)
{
    // n is large - do only the highest, then recur down

    #ifndef SIMINT_BOYS_NOVECTOR
    if(vector_min(x) > BOYS_SHORTGRID_MAXX)
        F[n] = boys_F_long_single_vec(x, n);
    else if(vector_max(x) < BOYS_SHORTGRID_MAXX)
        F[n] = boys_F_taylor_single_vec(x, n);
    else
    #endif
    {
        double * restrict Fd = (double *)F;
        double const * restrict xd = (double *)(&x);
        for(int i = 0; i < SIMINT_SIMD_LEN; i++)
        {
            if(xd[i] < BOYS_SHORTGRID_MAXX)
                Fd[n*SIMINT_SIMD_LEN+i] = boys_F_taylor_single(xd[i], n);
            else
                Fd[n*SIMINT_SIMD_LEN+i] = boys_F_long_single(xd[i], n);
        }
    }

    // factors for the recursion
    const SIMINT_DBLTYPE x2 = SIMINT_MUL(SIMINT_DBLSET1(2.0), (x));
    const SIMINT_DBLTYPE ex = SIMINT_EXP(SIMINT_NEG(x));

    // now recur down
    for(int n2 = n-1; n2 >= 0; n2--)
    {
        const SIMINT_DBLTYPE den = SIMINT_DBLSET1(1.0 / (2.0 * n2 + 1));

        //F[n2] = den * (x2 * F[(n2+1)] + ex);
        F[n2] = SIMINT_MUL(den, ( SIMINT_FMADD(x2, F[(n2+1)], ex)));
    }
}

static inline
void boys_F_split(SIMINT_DBLTYPE * restrict F,
                  SIMINT_DBLTYPE x,
                  int n)
{
    if(n < 4)
        boys_F_split_small_n(F, x, n);
    else
        boys_F_split_large_n(F, x, n);
}

#ifdef SIMINT_TARGET
#pragma message "SIMINT_TARGET is defined"
// #pragma omp declare target to(boys_F_split,boys_F_split_small_n,boys_F_split_large_n)

static inline
void boys_F_small_n_target(double * restrict F,
                          double x,
                          int n)
{
    // n is small - just do it all of them via
    // lookup or longfac (no recursion)
    {
        if(x < BOYS_SHORTGRID_MAXX)
            boys_F_taylor(F, x, n);
        else
            boys_F_long(F, x, n);    
    }
}
#else
#pragma message "SIMINT_TARGET is not defined"

#endif // SIMINT_TARGET



#ifdef __cplusplus
}
#endif

