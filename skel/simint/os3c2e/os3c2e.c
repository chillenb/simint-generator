#include "simint/os3c2e/os3c2e.h"
#include "simint/os3c2e/os3c2e_config.h"

// This is the actual storage for this array
#define AMSIZE   SIMINT_OS3C2E_MAXAM+1
#define DERSIZE  SIMINT_OS3C2E_MAXDER+1
simint_os3c2efunc simint_os3c2efunc_array[DERSIZE][AMSIZE][AMSIZE][AMSIZE][AMSIZE];


int simint_compute_os3c2e(struct simint_multi_shellpair const * P,
                         struct simint_multi_shellpair const * Q,
                         double screen_tol,
                         double * restrict work,
                         double * restrict integrals)
{
    // don't forget that we don't include the square root in the screen values
    // stored in the shell pair
    double screen_tol2 = screen_tol * screen_tol;
    if(screen_tol > 0.0 && (P->screen_max * Q->screen_max) < screen_tol2 )
        return -1;

    return simint_os3c2efunc_array[0][P->am1][P->am2][Q->am1][Q->am2](*P, *Q,
                                                screen_tol2, work, integrals);
}


int simint_compute_os3c2e_deriv(int deriv,
                               struct simint_multi_shellpair const * P,
                               struct simint_multi_shellpair const * Q,
                               double screen_tol,
                               double * restrict work,
                               double * restrict integrals)
{
    // don't forget that we don't include the square root in the screen values
    // stored in the shell pair
    double screen_tol2 = screen_tol * screen_tol;
    if(screen_tol > 0.0 && (P->screen_max * Q->screen_max) < screen_tol2 )
        return -1;

    return simint_os3c2efunc_array[deriv][P->am1][P->am2][Q->am1][Q->am2](*P, *Q,
                                                  screen_tol2, work, integrals);
}

