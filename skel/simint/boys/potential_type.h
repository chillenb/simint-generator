#ifndef SIMINT_POTENTIAL_TYPE_H
#define SIMINT_POTENTIAL_TYPE_H

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

#endif // SIMINT_POTENTIAL_TYPE_H