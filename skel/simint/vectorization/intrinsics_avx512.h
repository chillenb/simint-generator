#pragma once

#include <immintrin.h>
#include <stdint.h>
#include <math.h>

#include "simint/vectorization/vector_config.h"
#include "simint/vectorization/intrinsics_avx.h"

#ifdef __cplusplus
#include "simint/cpp_restrict.hpp"
extern "C" {
#endif

union simint_double8
{
    __m512d v;
    double d[8];
};


#ifdef SIMINT_USE_SVML

// Use SVML if available


#if defined __GNUC__ || defined __INTEL_LLVM_COMPILER // GCC, Clang, Intel LLVM
#define SIMINT_ATTRIBUTE_CONST __attribute__((const))
#else
#define SIMINT_ATTRIBUTE_CONST
#endif

SIMINT_ATTRIBUTE_CONST __m512d __svml_exp8(__m512d x);
SIMINT_ATTRIBUTE_CONST __m512d __svml_erf8(__m512d x);
SIMINT_ATTRIBUTE_CONST __m512d __svml_pow8(__m512d a, __m512d p);

static inline __m512d simint_exp_vec8(__m512d x)
{
    return __svml_exp8(x);
}

static inline __m512d simint_erf_vec8(__m512d x)
{
    return __svml_erf8(x);
}

static inline __m512d simint_pow_vec8(__m512d a, __m512d p)
{
    return __svml_pow8(a, p);
}

#elif __GLIBC__ >= 2 && __GLIBC_MINOR__ >= 22

// Otherwise, use glibc vector math if available

    __m512d _ZGVeN8v_exp(__m512d x);
    static inline __m512d simint_exp_vec8(__m512d x) { return _ZGVeN8v_exp(x); }

    __m512d _ZGVeN8v_erf(__m512d x);
    static inline __m512d simint_erf_vec8(__m512d x) { return _ZGVeN8v_erf(x); }

    __m512d _ZGVeN8vv_pow(__m512d a, __m512d p);
    static inline __m512d simint_pow_vec8(__m512d a, __m512d p) { return _ZGVeN8vv_pow(a, p); }
#else

// Fallback to scalar implementations

    static inline __m512d simint_exp_vec8(__m512d x)
    {
        union simint_double8 u = { x };
        union simint_double8 res;
        for(int i = 0; i < 8; i++)
            res.d[i] = exp(u.d[i]);
        return res.v;
    }

    static inline __m512d simint_erf_vec8(__m512d x)
    {
        union simint_double8 u = { x };
        union simint_double8 res;
        for(int i = 0; i < 8; i++)
            res.d[i] = erf(u.d[i]);
        return res.v;
    }

    static inline __m512d simint_pow_vec8(__m512d a, __m512d p)
    {
        union simint_double8 ua = { a };
        union simint_double8 up = { p };
        union simint_double8 res;
        for(int i = 0; i < 8; i++)
            res.d[i] = pow(ua.d[i], up.d[i]);
        return res.v;
    }
#endif


#if defined SIMINT_COREAVX512 || defined SIMINT_MICAVX512

    #define SIMINT_SIMD_LEN 8

    #define SIMINT_DBLTYPE         __m512d
    #define SIMINT_I32VEC          __m256i
    #define SIMINT_DBLLOAD(p,i)    _mm512_load_pd((p) + (i))
    #define SIMINT_DBLSET1(a)      _mm512_set1_pd((a))
    #define SIMINT_NEG(a)          (SIMINT_MUL((a), (SIMINT_DBLSET1(-1.0)))) 
    #define SIMINT_ADD(a,b)        _mm512_add_pd((a), (b))
    #define SIMINT_SUB(a,b)        _mm512_sub_pd((a), (b))
    #define SIMINT_MUL(a,b)        _mm512_mul_pd((a), (b))
    #define SIMINT_DIV(a,b)        _mm512_div_pd((a), (b))
    #define SIMINT_SQRT(a)         _mm512_sqrt_pd((a))
    #define SIMINT_FMADD(a,b,c)    _mm512_fmadd_pd((a), (b), (c))
    #define SIMINT_FMSUB(a,b,c)    _mm512_fmsub_pd((a), (b), (c))

    #define SIMINT_ROUNDTO_I32(a)  _mm512_cvttpd_epi32((a))
    #define SIMINT_I32_TO_PD(a)    _mm512_cvtepi32_pd((a))
    #define SIMINT_MUL_I32(a,b)    _mm256_mullo_epi32((a), (b))
    #define SIMINT_I32SET1(a)      _mm256_set1_epi32((a))

    #define SIMINT_ALL_GREATER_THAN(v, t)  all_greater_than((v), (t))

    #define SIMINT_GATHER_DBL_BY_I32(vdx, base)  _mm512_i32gather_pd((vdx), (base), sizeof(double))

    #if defined __INTEL_COMPILER 
        #define SIMINT_EXP(a)       _mm512_exp_pd((a))
        #define SIMINT_ERF(a)       _mm512_erf_pd((a))
        #define SIMINT_POW(a,p)     _mm512_pow_pd((a), (p))
    #else
        #define SIMINT_EXP(a)       simint_exp_vec8((a))
        #define SIMINT_ERF(a)       simint_erf_vec8((a))
        #define SIMINT_POW(a,p)     simint_pow_vec8((a), (p))
    #endif



    ////////////////////////////////////////
    // Special functions
    ////////////////////////////////////////

    static inline
    void contract(int ncart,
                  int const * restrict offsets,
                  __m512d const * restrict src,
                  double * restrict dest)
    {
        for(int n = 0; n < SIMINT_SIMD_LEN; ++n)
        {
            double const * restrict src_tmp = (double *)src + n;
            double * restrict dest_tmp = dest + offsets[n]*ncart;

            for(int np = 0; np < ncart; ++np)
            {
                dest_tmp[np] += *src_tmp;
                src_tmp += SIMINT_SIMD_LEN;
            }
        }
    }


    static inline
    void contract_all(int ncart,
                      __m512d const * restrict src,
                      double * restrict dest)
    {
        #if defined __clang__ || defined __INTEL_COMPILER  || defined __GNUC__
        
        int ntrans   = ncart  / 8;
        int np_start = ntrans * 8;
        
        double tmp[64] __attribute__((aligned(64)));
        __m512d dst[8] __attribute__((aligned(64)));

        // Transpose-Add part
        double *src_ptr = (double*)src;
        for (int it = 0; it < ntrans; it++)
        {
            for (int i = 0; i < 8; i++)
            {
                for (int j = 0; j < 8; j++)
                    tmp[i * 8 + j] = src_ptr[j * 8 + i];
                dst[i] = _mm512_load_pd(tmp + i * 8);
            }
            src_ptr += 64;
            
            __m512d res = _mm512_loadu_pd(dest + it * 8);
            for (int i = 0; i < 8; i++)
                res = _mm512_add_pd(res, dst[i]);
            _mm512_storeu_pd(dest + it * 8, res);
        }
        
        // Remainder part
        for (int np = np_start; np < ncart; np++)
            dest[np] += _mm512_reduce_add_pd(src[np]);
        
        #else

        int offsets[8] = {0, 0, 0, 0, 0, 0, 0, 0};
        contract(ncart, offsets, src, dest);

        #endif
    }


    static inline
    void contract_fac(int ncart,
                      const __m512d factor,
                      int const * restrict offsets,
                      __m512d const * restrict src,
                      double * restrict dest)
    {
        for(int np = 0; np < ncart; ++np)
        {
            union simint_double8 vtmp = { SIMINT_MUL(src[np], factor) };

            for(int n = 0; n < SIMINT_SIMD_LEN; ++n)
                dest[offsets[n]*ncart+np] += vtmp.d[n]; 
        }
    }


    static inline
    void contract_all_fac(int ncart,
                          const __m512d factor,
                          __m512d const * restrict src,
                          double * restrict dest)
    {
        #if defined __clang__ || defined __INTEL_COMPILER || defined __GNUC__

        for(int np = 0; np < ncart; np++)
            dest[np] += _mm512_reduce_add_pd(_mm512_mul_pd(factor, src[np]));

        #else

        int offsets[8] = {0, 0, 0, 0, 0, 0, 0, 0};
        contract_fac(ncart, factor, offsets, src, dest);

        #endif
    }


    static inline
    double vector_min(__m512d v)
    {
        #if defined __clang__ || defined __INTEL_COMPILER  || defined __GNUC__
            return _mm512_reduce_min_pd(v);
        #else
            union simint_double8 u = { v };
            double min = u.d[0];
            for(int i = 1; i < 8; i++)
                min = (u.d[i] < min ? u.d[i] : min);
            return min;
        #endif
    }


    static inline
    double vector_max(__m512d v)
    {
        #if defined __clang__ || defined __INTEL_COMPILER  || defined __GNUC__
            return _mm512_reduce_max_pd(v);
        #else
            union simint_double8 u = { v };
            double max = u.d[0];
            for(int i = 1; i < 8; i++)
                max = (u.d[i] > max ? u.d[i] : max);
            return max;
        #endif
    }

    static inline
    unsigned char all_greater_than(__m512d v, __m512d threshold)
    {
        __mmask8 mask = _mm512_cmp_pd_mask(v, threshold, _CMP_GT_OQ);
        // we now have a mask of 1s where v > threshold
        unsigned char comp_result = _kortestc_mask8_u8(mask, mask);
        return comp_result;
    }


    static inline
    __m512d mask_load(int nlane, double * memaddr)
    {
        return _mm512_maskz_loadu_pd((1 << nlane) - 1, memaddr);
    }
    
    //#define SIMINT_PRIM_SCREEN_STAT
    static inline
    int count_prim_screen_survival(__m512d screen_val, const double screen_tol)
    {
        union simint_double8 u = {screen_val};
        int res = 0;
        for (int i = 0; i < SIMINT_SIMD_LEN; i++)
            if (u.d[i] >= screen_tol) res++;
        return res;
    }

#endif // defined SIMINT_COREAVX512 || defined SIMINT_MICAVX512

#ifdef __cplusplus
}
#endif

