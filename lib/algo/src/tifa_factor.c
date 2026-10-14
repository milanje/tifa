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

#include "tifa_config.h"
#include "macros.h"
#include "first_primes.h"
#include "siqs.h"
#include "squfof.h"
#include "p_minus_1.h"
#include "p_plus_1.h"
#include "tifa_factor.h"

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
// Number of small primes (all the primes below 1000) divided out of larger
// numbers before any other algorithm runs. SQUFOF, for one, fails on 2p and
// on factors shared with its multipliers.
//
#define TFACTOR_NPRIMES_STRIP 168
#if TIFA_WORDSIZE < 64
    //
    // Bitsize threshold for SQUFOF on 32 bit machines. On 32 bit machines
    // TIFA's single precision SQUFOF will start to fail for 59-ish bit
    // numbers. (Untested for less than 32 bit machines!)
    //
    #define TFACTOR_SQUFOF_THRESH 58
#else
    //
    // Bitsize threshold for SQUFOF on 64 bit machines.
    //
    #define TFACTOR_SQUFOF_THRESH 64
#endif
//
// Bitsize from which the p-1 and p+1 pre-stages are tried before SIQS. Below
// it, a single run of both costs more than about 10% of the SIQS run.
//
#define TFACTOR_PRESTAGE_THRESH 110
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
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//
// Divides the first nprimes primes out of n, records them with their
// multiplicities in primes and pmultis, and leaves what remains in cofactor.
// Returns true if cofactor is then known to be 1 or prime.
//
static bool strip_small_primes(mpz_array_t* const primes,
                               uint32_array_t* const pmultis,
                               mpz_t cofactor,
                               const mpz_t n,
                               const uint32_t nprimes) {
    mpz_t p;
    mpz_init(p);
    mpz_set(cofactor, n);

    unsigned long int ntwos = mpz_scan1(cofactor, 0);
    if (ntwos > 0) {
        mpz_tdiv_q_2exp(cofactor, cofactor, ntwos);
        mpz_set_ui(p, 2UL);
        append_mpz_to_array(primes, p);
        append_uint32_to_array(pmultis, (uint32_t) ntwos);
    }
    bool known_prime_or_one = false;

    for (uint32_t i = 1; i < nprimes; i++) {
        const unsigned long int pi = first_primes[i];

        if (mpz_cmp_ui(cofactor, pi * pi) < 0) {
            known_prime_or_one = true;
            break;
        }
        if (mpz_divisible_ui_p(cofactor, pi)) {
            uint32_t multi = 0;
            do {
                mpz_divexact_ui(cofactor, cofactor, pi);
                multi++;
            } while (mpz_divisible_ui_p(cofactor, pi));

            mpz_set_ui(p, pi);
            append_mpz_to_array(primes, p);
            append_uint32_to_array(pmultis, multi);
        }
    }
    if (!known_prime_or_one) {
        const unsigned long int pnext = first_primes[nprimes];
        known_prime_or_one = (mpz_cmp_ui(cofactor, pnext * pnext) < 0);
    }
    mpz_clear(p);

    return known_prime_or_one;
}
//-----------------------------------------------------------------------------
//
// Adds the cofactor c > 1 left by strip_small_primes to the factors found.
//
static ecode_t add_cofactor(mpz_array_t* const factors,
                            uint32_array_t* const multis,
                            const mpz_t c,
                            const bool c_is_prime,
                            const factoring_mode_t mode) {
    if (mode != FIND_COMPLETE_FACTORIZATION) {
        if (c_is_prime || (mode != FIND_SOME_PRIME_FACTORS)) {
            append_mpz_to_array(factors, c);
        }
        return mode_to_outcome[mode];
    }
    if (c_is_prime) {
        append_mpz_to_array(factors, c);
        append_uint32_to_array(multis, 1);
        return COMPLETE_FACTORIZATION_FOUND;
    }
    mpz_array_t*    cfactors = alloc_mpz_array(ELONGATION);
    uint32_array_t* cmultis  = alloc_uint32_array(ELONGATION);

    ecode_t ecode = tifa_factor(cfactors, cmultis, c, mode);

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
        uint32_array_t* fmultis = alloc_uint32_array(4);

        ecode = process_found_factors(found, fmultis, n, mode, tifa_factor);

        append_mpz_array(factors, found);
        if (mode == FIND_COMPLETE_FACTORIZATION) {
            append_uint32_array(multis, fmultis);
        }
        free_uint32_array(fmultis);
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
                                           const factoring_mode_t mode) {
    //
    // Perfect squares can be problematic for some algorithms (i.e. SQUFOF).
    // Note that we can not (conveniently) rely on the client code to perform
    // such a check since tifa_factor(...) is indirectly and recursively called
    // with a factoring_machine_t when using the FIND_*_FACTORIZATION modes.
    //
    if (MPZ_IS_SQUARE(n)) {
        mpz_t rootn;
        mpz_init(rootn);
        mpz_sqrt(rootn, n);

        if (MPZ_IS_PRIME(rootn)) {
            append_mpz_to_array(factors, rootn);
            if (mode == FIND_COMPLETE_FACTORIZATION) {
                append_uint32_to_array(multis, 2);
            }
            mpz_clear(rootn);
            //
            // We have found the complete factorization but the exit code to
            // return actually depends on the factorization mode used.
            //
            return mode_to_outcome[mode];
        }
        uint32_t old_length = 0;

        if (mode == FIND_COMPLETE_FACTORIZATION) {
            old_length = multis->length;
        }
        ecode_t ecode = tifa_factor(factors, multis, rootn, mode);

        if (mode == FIND_COMPLETE_FACTORIZATION) {
            for (uint32_t i = old_length; i < multis->length; i++) {
                multis->data[i] <<= 1;
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
        return squfof(factors, multis, n, &params, mode);
    }
    if (sizen >= TFACTOR_PRESTAGE_THRESH) {
        ecode_t ecode = run_prestages(factors, multis, n, mode);
        if (ecode != NO_FACTOR_FOUND) {
            return ecode;
        }
    }
    siqs_params_t params;
    set_siqs_params_to_default(n, &params);
    return siqs(factors, multis, n, &params, mode);
}
//-----------------------------------------------------------------------------
ecode_t tifa_factor(mpz_array_t* const factors, uint32_array_t* const multis,
                    const mpz_t n, const factoring_mode_t mode) {

    if (mode > FIND_COMPLETE_FACTORIZATION) {
        return UNKNOWN_FACTORING_MODE;
    }
    if (mpz_cmp_ui(n, 1) <= 0) {
        return NO_FACTOR_FOUND;
    }
    const bool tdiv_only = (mpz_sizeinbase(n, 2) <= TFACTOR_TDIV_THRESH);

    mpz_array_t*    primes  = alloc_mpz_array(ELONGATION);
    uint32_array_t* pmultis = alloc_uint32_array(ELONGATION);
    mpz_t cofactor;
    mpz_init(cofactor);

    bool cofactor_is_prime = strip_small_primes(
                                 primes, pmultis, cofactor, n,
                                 tdiv_only ? TFACTOR_NPRIMES_TDIV
                                           : TFACTOR_NPRIMES_STRIP
                             );
    ecode_t ecode;

    if (primes->length == 0) {
        if (cofactor_is_prime) {
            ecode = NO_FACTOR_FOUND;
        } else {
            ecode = factor_without_small_primes(factors, multis, n, mode);
        }
    } else {
        append_mpz_array(factors, primes);
        if (mode == FIND_COMPLETE_FACTORIZATION) {
            append_uint32_array(multis, pmultis);
        }
        if (mpz_cmp_ui(cofactor, 1) == 0) {
            ecode = mode_to_outcome[mode];
        } else {
            if (!cofactor_is_prime) {
                cofactor_is_prime = MPZ_IS_PRIME(cofactor);
            }
            ecode = add_cofactor(factors, multis, cofactor,
                                 cofactor_is_prime, mode);
        }
    }
    mpz_clear(cofactor);
    free_mpz_array(primes);
    free_uint32_array(pmultis);

    return ecode;
}
//-----------------------------------------------------------------------------
