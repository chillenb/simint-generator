#pragma once

#include "simint/boys/boys_ratinterp_data.h"
#include "simint/vectorization/vectorization.h"

#ifdef __cplusplus
#include "simint/cpp_restrict.hpp"
extern "C" {
#endif

static inline double boys_F_rational_single(double x, int L) {
    double num = 1.0, denom = 1.0;
    int nnum = boys_q2_numerator_nterms[L];
    int nden = boys_q2_denominator_nterms[L];
    const double *coefs = boys_q2_rationalfit_data[L];
    int n = boys_q2_N[L];
    int i = 1;
    int k = 0;
    denom *= coefs[0];
    for(k = 0; k < 2 * nden; k+=2) {
        double yi = coefs[i];
        double ai = coefs[i + 1];
        double term = (x + yi) * (x + yi) + ai;
        denom *= term;
        i += 2;
    }
    for(; k < n; k++) {
        denom *= (x + coefs[i]);
        i++;
    }
    for(k = 0; k < 2 * nnum; k+=2) {
        double zi = coefs[i];
        double bi = coefs[i + 1];
        double term = (x + zi) * (x + zi) + bi;
        num *= term;
        i += 2;
    }
    for(; k < n - 1; k++) {
        num *= (x + coefs[i]);
        i++;
    }
    double q2 = num / denom;
    double q1c = coefs[i];
    double q0c = coefs[i + 1];
    double q = (q2 * x + q1c) * x + q0c;
    double cn = coefs[i + 2];
    double qLhalf = pow(q, L + 0.5);
    double tmp = x + qLhalf * exp(-x);
    return pow(cn / tmp, L + 0.5);

}


static inline SIMINT_DBLTYPE pownhalf(SIMINT_DBLTYPE a, int n) {
    SIMINT_DBLTYPE result = SIMINT_SQRT(a);
    SIMINT_DBLTYPE z = a;
    if (n & 1) {
        result = SIMINT_MUL(result, z);
    }
    for (n >>= 1; n; n >>= 1) {
        z = SIMINT_MUL(z, z);
        if (n & 1) {
        result = SIMINT_MUL(result, z);
        }
    }
    return result;
}

static inline SIMINT_DBLTYPE recurdown(SIMINT_DBLTYPE x,
                                       SIMINT_DBLTYPE FLx,
                                       SIMINT_DBLTYPE expnx,
                                       int L) {
const double one_over_2Lplus1 = 1.0 / (2.0 * L - 1);
const SIMINT_DBLTYPE twox = SIMINT_MUL(SIMINT_DBLSET1(2.0), x);
const SIMINT_DBLTYPE numerator = SIMINT_FMADD(twox, FLx, expnx);
return SIMINT_MUL(SIMINT_DBLSET1(one_over_2Lplus1), numerator);
}

static inline int get_extrarecur(int L) {
    return 6;
}

static inline SIMINT_DBLTYPE boys_F_rational_single_vec(SIMINT_DBLTYPE x, int L) {
    SIMINT_DBLTYPE expnx = SIMINT_EXP(SIMINT_NEG(x));
    SIMINT_DBLTYPE num = SIMINT_DBLSET1(1.0);
    SIMINT_DBLTYPE denom = SIMINT_DBLSET1(1.0);

    int Lupper = L + get_extrarecur(L);
    int nnum = boys_q2_numerator_nterms[Lupper];
    int nden = boys_q2_denominator_nterms[Lupper];
    const double *coefs = boys_q2_rationalfit_data[Lupper];
    int n = boys_q2_N[Lupper];
    int i = 1;
    int k = 0;

    denom = SIMINT_MUL(SIMINT_DBLSET1(coefs[0]), denom);
    for(k = 0; k < 2 * nden; k+=2) {
        SIMINT_DBLTYPE yi = SIMINT_DBLSET1(coefs[i]);
        SIMINT_DBLTYPE ai = SIMINT_DBLSET1(coefs[i + 1]);
        SIMINT_DBLTYPE xpyi = SIMINT_ADD(x, yi);
        denom = SIMINT_MUL(denom, SIMINT_FMADD(xpyi, xpyi, ai));
        i += 2;
    }
    for(; k < n; k++) {
        denom = SIMINT_MUL(denom, SIMINT_ADD(x, SIMINT_DBLSET1(coefs[i])));
        i++;
    }
    for(k = 0; k < 2 * nnum; k+=2) {
        SIMINT_DBLTYPE zi = SIMINT_DBLSET1(coefs[i]);
        SIMINT_DBLTYPE bi = SIMINT_DBLSET1(coefs[i + 1]);
        SIMINT_DBLTYPE xpzi = SIMINT_ADD(x, zi);
        num = SIMINT_MUL(num, SIMINT_FMADD(xpzi, xpzi, bi));
        i += 2;
    }

    for(; k < n - 1; k++) {
        num = SIMINT_MUL(num, SIMINT_ADD(x, SIMINT_DBLSET1(coefs[i])));
        i++;
    }
    SIMINT_DBLTYPE q2 = SIMINT_DIV(num, denom);

    const double q1c = coefs[i];
    const double q0c = coefs[i + 1];
    SIMINT_DBLTYPE q = SIMINT_FMADD(q2, x, SIMINT_DBLSET1(q1c));
    q = SIMINT_FMADD(q, x, SIMINT_DBLSET1(q0c));
    const double cn = coefs[i + 2];
    SIMINT_DBLTYPE qLhalf = pownhalf(q, Lupper);
    SIMINT_DBLTYPE tmp = SIMINT_ADD(x, SIMINT_MUL(qLhalf, expnx));
    SIMINT_DBLTYPE FxL = pownhalf(SIMINT_DIV(SIMINT_DBLSET1(cn), tmp), Lupper);

    for(int ell = Lupper; ell > L; ell--) {
        FxL = recurdown(x, FxL, expnx, ell);
    }
    return FxL;

}


static inline void boys_F_rational_vec(SIMINT_DBLTYPE * restrict F, SIMINT_DBLTYPE x, int L) {
    SIMINT_DBLTYPE expnx = SIMINT_EXP(SIMINT_NEG(x));
    SIMINT_DBLTYPE num = SIMINT_DBLSET1(1.0);
    SIMINT_DBLTYPE denom = SIMINT_DBLSET1(1.0);

    int Lupper = L; // + get_extrarecur(L);
    int nnum = boys_q2_numerator_nterms[Lupper];
    int nden = boys_q2_denominator_nterms[Lupper];
    const double *coefs = boys_q2_rationalfit_data[Lupper];
    int n = boys_q2_N[Lupper];
    int i = 1;
    int k = 0;

    denom = SIMINT_MUL(SIMINT_DBLSET1(coefs[0]), denom);
    for(k = 0; k < 2 * nden; k+=2) {
        SIMINT_DBLTYPE yi = SIMINT_DBLSET1(coefs[i]);
        SIMINT_DBLTYPE ai = SIMINT_DBLSET1(coefs[i + 1]);
        SIMINT_DBLTYPE xpyi = SIMINT_ADD(x, yi);
        denom = SIMINT_MUL(denom, SIMINT_FMADD(xpyi, xpyi, ai));
        i += 2;
    }
    for(; k < n; k++) {
        denom = SIMINT_MUL(denom, SIMINT_ADD(x, SIMINT_DBLSET1(coefs[i])));
        i++;
    }
    for(k = 0; k < 2 * nnum; k+=2) {
        SIMINT_DBLTYPE zi = SIMINT_DBLSET1(coefs[i]);
        SIMINT_DBLTYPE bi = SIMINT_DBLSET1(coefs[i + 1]);
        SIMINT_DBLTYPE xpzi = SIMINT_ADD(x, zi);
        num = SIMINT_MUL(num, SIMINT_FMADD(xpzi, xpzi, bi));
        i += 2;
    }

    for(; k < n - 1; k++) {
        num = SIMINT_MUL(num, SIMINT_ADD(x, SIMINT_DBLSET1(coefs[i])));
        i++;
    }
    SIMINT_DBLTYPE q2 = SIMINT_DIV(num, denom);

    const double q1c = coefs[i];
    const double q0c = coefs[i + 1];
    SIMINT_DBLTYPE q = SIMINT_FMADD(q2, x, SIMINT_DBLSET1(q1c));
    q = SIMINT_FMADD(q, x, SIMINT_DBLSET1(q0c));
    const double cn = coefs[i + 2];
    SIMINT_DBLTYPE qLhalf = pownhalf(q, Lupper);
    SIMINT_DBLTYPE tmp = SIMINT_ADD(x, SIMINT_MUL(qLhalf, expnx));
    SIMINT_DBLTYPE FxL = pownhalf(SIMINT_DIV(SIMINT_DBLSET1(cn), tmp), Lupper);

    for(int ell = Lupper; ell > L; ell--) {
        FxL = recurdown(x, FxL, expnx, ell);
    }

    for(int ell = L; ell > 0; ell--) {
        F[ell] = FxL;
        FxL = recurdown(x, FxL, expnx, ell);
    }
    F[0] = FxL;
}

#ifdef __cplusplus
}
#endif
