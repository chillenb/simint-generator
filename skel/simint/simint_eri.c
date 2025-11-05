#include "simint/simint_eri.h"
#include "simint/ostei/ostei.h"

int simint_compute_eri(struct simint_multi_shellpair const * P,
                       struct simint_multi_shellpair const * Q,
                       double screen_tol,
                       double * restrict work,
                       double * restrict integrals)
{
    struct simint_eri_potential_data default_potential_data = { COULOMB_POTENTIAL, 0.0 };
    return simint_compute_ostei(P, Q, screen_tol, work, integrals, default_potential_data);
}

int simint_compute_eri_ex(struct simint_multi_shellpair const * P,
                           struct simint_multi_shellpair const * Q,
                           double screen_tol,
                           double * restrict work,
                           double * restrict integrals,
                           struct simint_eri_potential_data potential_data)
{
    return simint_compute_ostei(P, Q, screen_tol, work, integrals, potential_data);
}

int simint_compute_eri_deriv(int deriv,
                             struct simint_multi_shellpair const * P,
                             struct simint_multi_shellpair const * Q,
                             double screen_tol,
                             double * restrict work,
                             double * restrict integrals)
{
    struct simint_eri_potential_data default_potential_data = { COULOMB_POTENTIAL, 0.0 };
    return simint_compute_ostei_deriv(deriv, P, Q, screen_tol, work, integrals, default_potential_data);
}


size_t simint_eri_worksize(int derorder, int maxam)
{
    return simint_ostei_worksize(derorder, maxam);
}


size_t simint_eri_workmem(int derorder, int maxam)
{
    return simint_ostei_workmem(derorder, maxam);
}
