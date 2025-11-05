#pragma once

#include "simint/shell/shell.h"

#ifdef __cplusplus
#include "simint/cpp_restrict.hpp"
extern "C" {
#endif

enum simint_eri_potential_type {
    COULOMB_POTENTIAL = 0,
    ERF_COULOMB_POTENTIAL = 1,
    ERFC_COULOMB_POTENTIAL = 2
};

struct simint_eri_potential_data
{
    enum simint_eri_potential_type potential_type;
    double omega; // for erf and erfc potentials
};


//! A pointer to a function that calculates TEI utilizing a shared workspace
typedef int (*simint_osteifunc)(struct simint_multi_shellpair const,
                                struct simint_multi_shellpair const,
                                double,
                                double * restrict,
                                double * restrict,
                                struct simint_eri_potential_data const);


/*! \brief Compute an ostei given shell pair information
 *
 * \param [in] P The shell pairs for the bra side of the integral 
 * \param [in] Q The shell pairs for the ket side of the integral
 * \param [in] screen_tol Tolerance for screening (set to zero to disable)
 * \param [in] work Workspace to use in calculating the integrals
 * \param [inout] integrals Storage for the final integrals. Since size information
 *                          is not passed, you are expected to ensure that this buffer
 *                          is large enough
 */
int simint_compute_ostei(struct simint_multi_shellpair const * P,
                         struct simint_multi_shellpair const * Q,
                         double screen_tol,
                         double * restrict work,
                         double * restrict integrals,
                         struct simint_eri_potential_data const potential_data);

/*! \brief Compute an ostei derivative given shell pair information
 *
 * \param [in] deriv Order of the derivative to compute
 * \param [in] P The shell pairs for the bra side of the integral 
 * \param [in] Q The shell pairs for the ket side of the integral
 * \param [in] screen_tol Tolerance for screening (set to zero to disable)
 * \param [in] work Workspace to use in calculating the integrals
 * \param [inout] integrals Storage for the final integrals. Since size information
 *                          is not passed, you are expected to ensure that this buffer
 *                          is large enough
 */
int simint_compute_ostei_deriv(int deriv,
                               struct simint_multi_shellpair const * P,
                               struct simint_multi_shellpair const * Q,
                               double screen_tol,
                               double * restrict work,
                               double * restrict integrals,
                               struct simint_eri_potential_data const potential_data);



#ifdef __cplusplus
}
#endif

