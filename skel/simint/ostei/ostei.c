#include "simint/ostei/ostei.h"
#include "simint/ostei/ostei_config.h"
#include "simint/ostei/gen/ostei_generated.h"
#include <stdint.h>
#include <stdio.h>

// This is the actual storage for this array
#define AMSIZE   (SIMINT_OSTEI_MAXAM+1)
#define DERSIZE  (SIMINT_OSTEI_MAXDER+1)

#define SWITCHIDX(a, b, c, d)  ((d) + AMSIZE*((c) + AMSIZE*((b) + AMSIZE*(a))))

simint_osteifunc simint_osteifunc_array[DERSIZE][AMSIZE][AMSIZE][AMSIZE][AMSIZE];

#ifdef SIMINT_TARGET
simint_osteifunc simint_osteifunc_array_target[DERSIZE][AMSIZE][AMSIZE][AMSIZE][AMSIZE];
#pragma omp declare target to(simint_osteifunc_array_target) device_type(nohost)
#endif

#define uglymap map(to: screen_tol2)\
                map(from:retval)\
                map(alloc: work[:worksize])\
                map(tofrom:integrals[:nshell1234 * ncart1234])\
                map(to: PP,\
                PP.nshell12, \
                PP.nprim,\
                PP.am1,\
                PP.am2,\
                PP.nshell12_clip, \
                PP.nprim12[:PP.nshell12], \
                PP.AB_x[:PP.nshell12], \
                PP.AB_y[:PP.nshell12], \
                PP.AB_z[:PP.nshell12], \
                PP.x[:PP.nprim], \
                PP.y[:PP.nprim], \
                PP.z[:PP.nprim], \
                PP.PA_x[:PP.nprim], \
                PP.PA_y[:PP.nprim], \
                PP.PA_z[:PP.nprim], \
                PP.PB_x[:PP.nprim], \
                PP.PB_y[:PP.nprim], \
                PP.PB_z[:PP.nprim], \
                PP.alpha[:PP.nprim], \
                PP.prefac[:PP.nprim], \
                PP.screen_max,\
                QQ,\
                QQ.nshell12, \
                QQ.nprim,\
                QQ.am1,\
                QQ.am2,\
                QQ.nshell12_clip, \
                QQ.nprim12[:QQ.nshell12], \
                QQ.AB_x[:QQ.nshell12], \
                QQ.AB_y[:QQ.nshell12], \
                QQ.AB_z[:QQ.nshell12], \
                QQ.x[:QQ.nprim], \
                QQ.y[:QQ.nprim], \
                QQ.z[:QQ.nprim], \
                QQ.PA_x[:QQ.nprim], \
                QQ.PA_y[:QQ.nprim], \
                QQ.PA_z[:QQ.nprim], \
                QQ.PB_x[:QQ.nprim], \
                QQ.PB_y[:QQ.nprim], \
                QQ.PB_z[:QQ.nprim], \
                QQ.alpha[:QQ.nprim], \
                QQ.prefac[:QQ.nprim], \
                QQ.screen_max)

int simint_compute_ostei(struct simint_multi_shellpair const * P,
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

    #ifdef SIMINT_TARGET
    int retval = 0;

    struct simint_multi_shellpair const PP = *P;
    struct simint_multi_shellpair const QQ = *Q;
    const int nshell1234 = PP.nshell12 * QQ.nshell12;
    const int ncart1234 = NCART(PP.am1) * NCART(PP.am2) * NCART(QQ.am1) * NCART(QQ.am2);
    int maxam = MAX(PP.am1, PP.am2);
    maxam = MAX(maxam, QQ.am1);
    maxam = MAX(maxam, QQ.am2);
    const int worksize = simint_ostei_worksize(0, maxam);

    /*
    printf("P has %d primitive combinations in %d shell pair(s)\n", PP.nprim, PP.nshell12);
    printf("Q has %d primitive combinations in %d shell pair(s)\n", QQ.nprim, QQ.nshell12);
    printf("Work size: %d\n", worksize);
    printf("Integrals: %d\n", nshell1234*ncart1234);
    */

    int shellq_idx = SWITCHIDX(PP.am1, PP.am2, QQ.am1, QQ.am2);
    #include "simint/ostei/gen/ostei_switch.c"

    return retval;
    #else
    return simint_osteifunc_array[0][P->am1][P->am2][Q->am1][Q->am2](*P, *Q,
                                                screen_tol2, work, integrals);
    #endif
}


int simint_compute_ostei_deriv(int deriv,
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

    return simint_osteifunc_array[deriv][P->am1][P->am2][Q->am1][Q->am2](*P, *Q,
                                                  screen_tol2, work, integrals);
}

