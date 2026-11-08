// Copyright (C) 2011 CNRS - Ecole Polytechnique - INRIA.
//
// This file is part of TIFA.
//
// TIFA is free software; you can redistribute it and/or modify it under the 
// terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 2.1 of the License, or (at your option)
// any later version.
//
// TIFA is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
// more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation,
// Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.

/**
 * \file    tifa_factor.c
 * \author  Jerome Milan
 * \date    Tue Oct 11 2011
 * \version 2011-10-11
 */

#include <stdlib.h>
#include <stdbool.h>
#include <pthread.h>

#include "tifa_config.h"
#include "macros.h"
#include "funcs.h"
#include "factor_u64.h"
#include "tdiv.h"
#include "siqs.h"
#include "squfof.h"
#include "p_minus_1.h"
#include "p_plus_1.h"
#include "ecm.h"
#include "tifa_factor.h"
#include "ecm_schedule.h"
#include "stopping.h"

//-----------------------------------------------------------------------------
//
// Number of primes to use for trial division if the size of the composite to
// factor is below the threshold given by TFACTOR_TDIV_THRESH. The first 1028
// primes go up to 8191, the largest prime below sqrt(2^26), so trial division
// alone completely factors such numbers.
//
#define TFACTOR_NPRIMES_TDIV 1028
//
// Bitsize threshold for trial division.
//
#define TFACTOR_TDIV_THRESH 26
//
// Number of primes (all the primes below 1000) to use for trial division on
// larger numbers, before any other algorithm runs. SQUFOF, for one, fails on
// 2p and on factors shared with its multipliers.
//
#define TFACTOR_NPRIMES_PRE_TDIV 168
//
// Bitsize threshold for SQUFOF (only used without factor_u64), above which
// SIQS is faster. On 32 bit machines SQUFOF also starts to fail for 59-ish
// bit numbers.
//
#define TFACTOR_SQUFOF_THRESH 58
//
// Bitsize from which the p-1 and p+1 pre-stages are tried before SIQS, when
// the ECM pre-pass can not run (it finds their factors for less). Below it, a
// single run of both costs more than about 10% of the SIQS run.
//
#define TFACTOR_PRESTAGE_THRESH 90
//
// Stage-1 bound for the single Pollard p-1 run tried before SIQS on large
// numbers. It peels off a prime factor p whenever p-1 is smooth.
//
#define TFACTOR_PM1_B1 2000
//
// Stage-1 bound for the single Williams p+1 run tried right after the p-1
// one. It catches a prime p when p+1 (rather than p-1) is smooth.
//
#define TFACTOR_PP1_B1 2000
//
// Fraction of the expected SIQS time the ECM pre-pass may take. Balanced
// semiprimes pay it in vain, but random composites gain the most from 0.2,
// and numbers with a factor below a third of their size from 1.
//
#define TFACTOR_ECM_BUDGET 0.2
#define TFACTOR_ECM_BUDGET_UNBALANCED 1.0
//
// Size of the largest factors the ECM pre-pass looks for (bits): up to
// TFACTOR_WIDE_ECM_BITS bits, then above, where SIQS takes half a minute or
// more. The plans of the larger levels are built when such a number comes.
//
#define TFACTOR_ECM_FBITS 60
#define TFACTOR_WIDE_ECM_FBITS 85
#define TFACTOR_WIDE_ECM_BITS 256
//
// Above this size (that of mont_ecm), ecm() in GMP arithmetic searches for
// factors of up to 85 bits, without a budget: SIQS would take years.
//
#define TFACTOR_GMP_ECM_BITS 512
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//
// Adds the composite cofactor c left by trial division to the factors found.
//
static ecode_t add_cofactor(mpz_array_t* const factors,
                            uint32_array_t* const multis,
                            const mpz_t c,
                            const factoring_mode_t mode,
                            const tifa_factor_hint_t hint) {
    if (mode != FIND_COMPLETE_FACTORIZATION) {
        if (mode != FIND_SOME_PRIME_FACTORS) {
            append_mpz_to_array(factors, c);
        }
        return mode_to_outcome[mode];
    }
    mpz_array_t*    cfactors = alloc_mpz_array(ELONGATION);
    uint32_array_t* cmultis  = alloc_uint32_array(ELONGATION);

    ecode_t ecode = tifa_factor_with_hint(cfactors, cmultis, c, mode, hint);

    if (   (ecode == COMPLETE_FACTORIZATION_FOUND)
        || (ecode == PARTIAL_FACTORIZATION_FOUND)) {
        append_mpz_array(factors, cfactors);
        append_uint32_array(multis, cmultis);
    } else {
        append_mpz_to_array(factors, c);
        append_uint32_to_array(multis, 1);
        ecode = PARTIAL_FACTORIZATION_FOUND;
    }
    free_mpz_array(cfactors);
    free_uint32_array(cmultis);

    return ecode;
}
//-----------------------------------------------------------------------------
//
// Appends the factors of n built from those found (a proper factor and its
// cofactor, or the like), processed as the factoring mode requires.
//
static ecode_t add_found(mpz_array_t* const factors,
                         uint32_array_t* const multis,
                         mpz_array_t* const found, const mpz_t n,
                         const factoring_mode_t mode) {
    uint32_array_t* fmultis = alloc_uint32_array(4);

    ecode_t ecode = process_found_factors(found, fmultis, n, mode, tifa_factor);

    append_mpz_array(factors, found);
    if (mode == FIND_COMPLETE_FACTORIZATION) {
        append_uint32_array(multis, fmultis);
    }
    free_uint32_array(fmultis);

    return ecode;
}
//-----------------------------------------------------------------------------
#if TIFA_HAVE_MONT_ECM
static ecm_schedule_t* ecm_schedule           = NULL;
static pthread_once_t  ecm_schedule_once      = PTHREAD_ONCE_INIT;
static ecm_schedule_t* wide_ecm_schedule      = NULL;
static pthread_once_t  wide_ecm_schedule_once = PTHREAD_ONCE_INIT;

static void build_ecm_schedule(void) {
    ecm_schedule = alloc_ecm_schedule(TFACTOR_ECM_FBITS);
}

static void build_wide_ecm_schedule(void) {
    wide_ecm_schedule = alloc_ecm_schedule(TFACTOR_WIDE_ECM_FBITS);
}
//
// Searches n for a factor with ECM, for a fraction of the time SIQS would
// take: a factor found early leaves a remainder much faster to factor (most
// often below 64 bits) than n by SIQS.
//
static ecode_t run_ecm_prepass(mpz_array_t* const factors,
                               uint32_array_t* const multis,
                               const mpz_t n, const factoring_mode_t mode,
                               const double budget) {
    ecm_schedule_t* schedule;
    if (mpz_sizeinbase(n, 2) > TFACTOR_WIDE_ECM_BITS) {
        pthread_once(&wide_ecm_schedule_once, build_wide_ecm_schedule);
        schedule = wide_ecm_schedule;
    } else {
        pthread_once(&ecm_schedule_once, build_ecm_schedule);
        schedule = ecm_schedule;
    }
    ecode_t ecode = NO_FACTOR_FOUND;
    mpz_t g;
    mpz_init(g);
    if (ecm_schedule_prepass(g, n, schedule, budget)) {
        mpz_array_t* found = alloc_mpz_array(4);
        append_mpz_to_array(found, g);
        mpz_divexact(g, n, g);
        append_mpz_to_array(found, g);
        ecode = add_found(factors, multis, found, n, mode);
        free_mpz_array(found);
    }
    mpz_clear(g);

    return ecode;
}
#endif
//-----------------------------------------------------------------------------
//
// Single runs of p-1 and p+1. Without SINGLE_RUN, the factoring machine would
// retry each of them with growing bounds (and p+1 with several seeds), at
// many times the cost of the SIQS run they are meant to shortcut.
//
static ecode_t run_prestages(mpz_array_t* const factors,
                             uint32_array_t* const multis,
                             const mpz_t n,
                             const factoring_mode_t mode) {
    mpz_array_t* found = alloc_mpz_array(4);

    p_minus_1_params_t pm1_params;
    set_p_minus_1_params_to_default(n, &pm1_params);
    pm1_params.b1 = TFACTOR_PM1_B1;
    pm1_params.b2 = 0;
    ecode_t ecode = p_minus_1(found, NULL, n, &pm1_params, SINGLE_RUN);

    if (ecode != SOME_FACTORS_FOUND) {
        reset_mpz_array(found);

        p_plus_1_params_t pp1_params;
        set_p_plus_1_params_to_default(n, &pp1_params);
        pp1_params.b1 = TFACTOR_PP1_B1;
        pp1_params.b2 = 0;
        ecode = p_plus_1(found, NULL, n, &pp1_params, SINGLE_RUN);
    }
    if (ecode == SOME_FACTORS_FOUND) {
        ecode = add_found(factors, multis, found, n, mode);
    } else {
        ecode = NO_FACTOR_FOUND;
    }
    free_mpz_array(found);

    return ecode;
}
//-----------------------------------------------------------------------------
static ecode_t run_gmp_ecm(mpz_array_t* const factors,
                           uint32_array_t* const multis,
                           const mpz_t n,
                           const factoring_mode_t mode) {
    mpz_array_t* found = alloc_mpz_array(4);

    ecm_params_t params;
    set_ecm_params_to_default(n, &params);
    ecode_t ecode = ecm(found, NULL, n, &params, FIND_SOME_FACTORS);

    if (ecode == SOME_FACTORS_FOUND) {
        ecode = add_found(factors, multis, found, n, mode);
    } else {
        ecode = NO_FACTOR_FOUND;
    }
    free_mpz_array(found);

    return ecode;
}
//-----------------------------------------------------------------------------
//
// Factors an n > 2^TFACTOR_TDIV_THRESH with no prime factor below 1000.
//
static ecode_t factor_without_small_primes(mpz_array_t* const factors,
                                           uint32_array_t* const multis,
                                           const mpz_t n,
                                           const factoring_mode_t mode,
                                           const tifa_factor_hint_t hint) {
    //
    // Perfect powers can be problematic for some algorithms (i.e. SQUFOF), and
    // a congruence of squares never splits a prime power. Note that we can not
    // (conveniently) rely on the client code to perform such a check since
    // tifa_factor(...) is indirectly and recursively called with a
    // factoring_machine_t when using the FIND_*_FACTORIZATION modes.
    //
    if (mpz_perfect_power_p(n)) {
        mpz_t rootn;
        mpz_init(rootn);
        //
        // The largest k with n a k-th power, so that the root is no power.
        // With no prime factor below 1000 > 2^9, k is at most bits(n) / 9.
        //
        uint32_t k = (uint32_t) mpz_sizeinbase(n, 2) / 9 + 1;
        while (0 == mpz_root(rootn, n, k)) {
            k--;
        }

        if (MPZ_IS_PRIME(rootn)) {
            append_mpz_to_array(factors, rootn);
            if (mode == FIND_COMPLETE_FACTORIZATION) {
                append_uint32_to_array(multis, k);
            }
            mpz_clear(rootn);
            //
            // We have found the complete factorization but the exit code to
            // return actually depends on the factorization mode used.
            //
            return mode_to_outcome[mode];
        }
        const uint32_t old_flength = factors->length;
        const uint32_t old_length  = multis->length;

        ecode_t ecode = tifa_factor_with_hint(factors, multis, rootn, mode,
                                              hint);

        if (   (ecode != COMPLETE_FACTORIZATION_FOUND)
            && (ecode != PARTIAL_FACTORIZATION_FOUND)
            && (ecode != (ecode_t) mode_to_outcome[mode])) {
            //
            // Stopped: n = rootn^k is still a split, as add_cofactor gives.
            //
            factors->length = old_flength;
            multis->length  = old_length;
            if (mode == FIND_COMPLETE_FACTORIZATION) {
                append_mpz_to_array(factors, rootn);
                append_uint32_to_array(multis, 1);
                ecode = PARTIAL_FACTORIZATION_FOUND;
            } else if (mode != FIND_SOME_PRIME_FACTORS) {
                append_mpz_to_array(factors, rootn);
                ecode = (ecode_t) mode_to_outcome[mode];
            }
        }
        if (mode == FIND_COMPLETE_FACTORIZATION) {
            for (uint32_t i = old_length; i < multis->length; i++) {
                multis->data[i] *= k;
            }
        }
        mpz_clear(rootn);

        return ecode;
    }
    unsigned long int sizen = mpz_sizeinbase(n, 2);
    //
    // _NOTE_: The following choice of algorithm has been obtained by benching
    //         all TIFA's factoring algorithms with RSA moduli (i.e. numbers
    //         that factors in two primes of equal size). For more general
    //         numbers, things can be somewhat different, particularly for
    //         SQUFOF or ECM whose theoretical complexity directly depend on
    //         the number of prime factors or the size of the smaller one.
    //
    if (sizen <= TFACTOR_SQUFOF_THRESH) {
        squfof_params_t params;
        set_squfof_params_to_default(&params);
        ecode_t ecode = squfof(factors, multis, n, &params, mode);
        if (   (ecode == (ecode_t) mode_to_outcome[mode])
            || (ecode == PARTIAL_FACTORIZATION_FOUND)) {
            return ecode;
        }
        //
        // SQUFOF gave up (rarely): SIQS.
        //
    }
    ecode_t ecode = NO_FACTOR_FOUND;
    const bool prepass = (hint != TIFA_FACTOR_BALANCED);
#if TIFA_HAVE_MONT_ECM
    const bool ecm = mont_ecm_handles(n);
    if (ecm && prepass) {
        ecode = run_ecm_prepass(factors, multis, n, mode,
                                (hint == TIFA_FACTOR_UNBALANCED)
                                    ? TFACTOR_ECM_BUDGET_UNBALANCED
                                    : TFACTOR_ECM_BUDGET);
    }
#else
    const bool ecm = false;
#endif
    if (!ecm && prepass && (sizen >= TFACTOR_PRESTAGE_THRESH)) {
        ecode = run_prestages(factors, multis, n, mode);
    }
    if (   !ecm && prepass && (ecode == NO_FACTOR_FOUND)
        && (sizen > TFACTOR_GMP_ECM_BITS)) {
        ecode = run_gmp_ecm(factors, multis, n, mode);
    }
    if (ecode != NO_FACTOR_FOUND) {
        return ecode;
    }
    if (tifa_should_stop()) {
        return GIVING_UP;
    }
    siqs_params_t params;
    set_siqs_params_to_default(n, &params);
    return siqs(factors, multis, n, &params, mode);
}
//-----------------------------------------------------------------------------
#if TIFA_HAVE_FACTOR_U64
//
// Factors n < 2^64 entirely in machine words. Returns FAILURE if some
// composite factor could not be split, so that the caller falls back to the
// multi-precision algorithms.
//
static ecode_t factor_with_u64(mpz_array_t* const factors,
                               uint32_array_t* const multis,
                               const mpz_t n,
                               const factoring_mode_t mode) {
    uint64_t primes[FACTOR_U64_MAX_PRIMES];
    uint32_t pmultis[FACTOR_U64_MAX_PRIMES];

    const uint32_t nprimes = factor_u64(
                                 (uint64_t) mpz_get_ui(n), primes, pmultis
                             );
    if (nprimes == 0) {
        return FAILURE;
    }
    if ((nprimes == 1) && (pmultis[0] == 1)) {
        return NO_FACTOR_FOUND;
    }
    mpz_t p;
    mpz_init(p);
    for (uint32_t i = 0; i < nprimes; i++) {
        mpz_set_ui(p, (unsigned long int) primes[i]);
        append_mpz_to_array(factors, p);
        if (mode == FIND_COMPLETE_FACTORIZATION) {
            append_uint32_to_array(multis, pmultis[i]);
        }
    }
    mpz_clear(p);

    return mode_to_outcome[mode];
}
#endif
//-----------------------------------------------------------------------------
ecode_t tifa_factor(mpz_array_t* const factors, uint32_array_t* const multis,
                    const mpz_t n, const factoring_mode_t mode) {
    return tifa_factor_with_hint(factors, multis, n, mode, TIFA_FACTOR_RANDOM);
}
//-----------------------------------------------------------------------------
ecode_t tifa_factor_with_hint(mpz_array_t* const factors,
                              uint32_array_t* const multis,
                              const mpz_t n, const factoring_mode_t mode,
                              const tifa_factor_hint_t hint) {

    if (mode > FIND_COMPLETE_FACTORIZATION) {
        return UNKNOWN_FACTORING_MODE;
    }
    if (hint > TIFA_FACTOR_UNBALANCED) {
        return FAILURE;
    }
    if (mpz_cmp_ui(n, 1) <= 0) {
        return NO_FACTOR_FOUND;
    }
#if TIFA_HAVE_FACTOR_U64
    if (mpz_sizeinbase(n, 2) <= 64) {
        ecode_t ecode = factor_with_u64(factors, multis, n, mode);
        if (ecode != FAILURE) {
            return ecode;
        }
    }
#endif
    const bool tdiv_only = (mpz_sizeinbase(n, 2) <= TFACTOR_TDIV_THRESH);

    mpz_array_t*    tfactors = alloc_mpz_array(ELONGATION);
    uint32_array_t* tmultis  = alloc_uint32_array(ELONGATION);

    ecode_t ecode = tdiv(tfactors, tmultis, n,
                         tdiv_only ? TFACTOR_NPRIMES_TDIV
                                   : TFACTOR_NPRIMES_PRE_TDIV);
    //
    // A prime within the trial division bound "divides" itself.
    //
    if ((ecode == COMPLETE_FACTORIZATION_FOUND) && (tfactors->length == 1)
            && (tmultis->data[0] == 1)) {
        ecode = NO_FACTOR_FOUND;
    }
    if (ecode == NO_FACTOR_FOUND) {
        if (!tdiv_only && !MPZ_IS_PRIME(n)) {
            ecode = factor_without_small_primes(factors, multis, n, mode,
                                                hint);
        }
    } else {
        //
        // tdiv appends the cofactor last and reports SOME_FACTORS_FOUND only
        // when that cofactor is composite.
        //
        const bool composite_left = (ecode == SOME_FACTORS_FOUND);
        mpz_t cofactor;
        mpz_init(cofactor);

        if (composite_left) {
            tfactors->length--;
            tmultis->length--;
            mpz_set(cofactor, tfactors->data[tfactors->length]);
        }
        append_mpz_array(factors, tfactors);
        if (mode == FIND_COMPLETE_FACTORIZATION) {
            append_uint32_array(multis, tmultis);
        }
        if (composite_left) {
            ecode = add_cofactor(factors, multis, cofactor, mode, hint);
        } else {
            ecode = mode_to_outcome[mode];
        }
        mpz_clear(cofactor);
    }
    free_mpz_array(tfactors);
    free_uint32_array(tmultis);

    return ecode;
}
//-----------------------------------------------------------------------------
