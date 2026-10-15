// Copyright (C) 2026 Jerome Milan.
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
 * \file    p_minus_1.c
 * \author  Jerome Milan
 * \date    2026
 * \version 1.0
 */

#include <stdlib.h>
#include <inttypes.h>

#include "tifa_config.h"
#include "first_primes.h"
#include "macros.h"
#include "factoring_machine.h"
#include "p_minus_1.h"
#include "tifa_factor.h"
#include "tuning.h"
#include "mont_ecm.h"

#define __PREFIX__      "p-1: "
#define __VERBOSE__     TIFA_VERBOSE_PM1
#define __TIMING__      TIFA_TIMING_PM1

#include "messages.h"

//
// Maximum number of times the context is updated (i.e. B1 is grown and the
// stage is retried) before giving up.
//
#define MAX_NUPDATES     4
//
// Factor by which B1 grows on each update.
//
#define B1_GROWTH_FACTOR 2UL
//
// Default parameters.
//
#define DEFAULT_BASE     2
#define DEFAULT_B1       50000
//
// Default stage-2 bound, as a multiple of B1 (clamped to LARGEST_PRIME).
//
#define DEFAULT_B2_RATIO 20UL
//
// Half of the largest prime gap the stage-2 gap-power table covers. The
// maximal gap between consecutive primes below LARGEST_PRIME (821641) is 148,
// so gaps up to 200 are handled directly; any larger gap falls back to a
// modular exponentiation.
//
#define MAX_GAP_HALF     100
//
// Number of (a^q - 1) terms accumulated into a single product before one gcd
// with n in stage 2 (gcd dominates the per-prime cost, so batching it is a
// win). A batch whose gcd equals n is rescanned prime-by-prime.
//
#define STAGE2_BATCH     256

//
// An ad-hoc structure holding the data variables used by the Pollard's p-1
// implementation.
//
struct struct_p_minus_1_context_t {
    //
    // Number to factor.
    //
    mpz_t n;
    //
    // p-1 parameters.
    //
    p_minus_1_params_t* params;
    //
    // Number of times the context was updated.
    //
    int nupdates;
    //
    // Working variables.
    //
    mpz_t a;
    mpz_t g;
    mpz_t aq;
#if TIFA_HAVE_MONT_ECM
    //
    // Plan of the fixed-limb Montgomery implementation for the bounds
    // (mplan_b1, mplan_b2), used when n has at most MONT_ECM_MAX_LIMBS limbs.
    //
    mont_ecm_plan_t* mplan;
    uint32_t         mplan_b1;
    uint32_t         mplan_b2;
#endif
};
typedef struct struct_p_minus_1_context_t p_minus_1_context_t;

//------------------------------------------------------------------------------
//                 PROTOTYPES OF NON PUBLIC FUNCTION(S)
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
static ecode_t init_p_minus_1_context(factoring_machine_t* const machine);
//------------------------------------------------------------------------------
static ecode_t clear_p_minus_1_context(factoring_machine_t* const machine);
//------------------------------------------------------------------------------
static ecode_t update_p_minus_1_context(factoring_machine_t* const machine);
//------------------------------------------------------------------------------
static ecode_t perform_p_minus_1(factoring_machine_t* const machine);
static void backtrack_stage_1(mpz_t g, const mpz_t n,
                              const p_minus_1_params_t* const params);
//------------------------------------------------------------------------------
static ecode_t recurse(
    mpz_array_t* const,
    uint32_array_t* const,
    const mpz_t,
    factoring_mode_t
);
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//                 "PUBLIC" FUNCTION(S) IMPLEMENTATION
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
ecode_t p_minus_1(mpz_array_t* const factors, uint32_array_t* const multis,
                  const mpz_t n, const p_minus_1_params_t* const params,
                  const factoring_mode_t mode) {

    PRINT_INTRO_MSG("Pollard's p-1");
    INIT_TIMER;
    START_TIMER;

    ecode_t ecode = NO_FACTOR_FOUND;

    factoring_machine_t machine;

    mpz_init_set(machine.n, n);
    machine.mode                = mode;
    machine.params              = (void*) params;
    machine.init_context_func   = init_p_minus_1_context;
    machine.perform_algo_func   = perform_p_minus_1;
    machine.update_context_func = update_p_minus_1_context;
    machine.clear_context_func  = clear_p_minus_1_context;
    machine.recurse_func        = recurse;
    machine.factors             = factors;
    machine.multis              = multis;

    ecode = run_machine(&machine);

    mpz_clear(machine.n);

    STOP_TIMER;
    PRINT_STATUS(machine.success, ecode);

    return ecode;
}
//------------------------------------------------------------------------------
void set_p_minus_1_params_to_default(const mpz_t n,
                                     p_minus_1_params_t* const params) {
    params->base = DEFAULT_BASE;
    params->b1   = DEFAULT_B1;
    if (params->b1 > LARGEST_PRIME) {
        params->b1 = LARGEST_PRIME;
    }
    unsigned long b2 = DEFAULT_B2_RATIO * (unsigned long) params->b1;
    if (b2 > LARGEST_PRIME) {
        b2 = LARGEST_PRIME;
    }
    params->b2 = (uint32_t) b2;

    //
    // Override with machine-tuned values if a tuning file supplies a row for
    // this bit-size; otherwise the defaults above are kept as-is.
    //
    const tuning_row_t* row = tifa_tuning_row_for(
        "p_minus_1", TIFA_TUNING_NCOLS_P_MINUS_1,
        (uint32_t) mpz_sizeinbase(n, 2)
    );
    if (row != NULL) {
        params->b1   = row->cols[0];
        params->b2   = row->cols[1];
        params->base = row->cols[2];
        if (params->b1 > LARGEST_PRIME) {
            params->b1 = LARGEST_PRIME;
        }
        if (params->b2 > LARGEST_PRIME) {
            params->b2 = LARGEST_PRIME;
        }
    }
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//                  "PRIVATE" FUNCTION(S) IMPLEMENTATION
//------------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *                  Functions used by factoring machine
 *-----------------------------------------------------------------------------
 */

//------------------------------------------------------------------------------
static ecode_t init_p_minus_1_context(factoring_machine_t* const machine) {
    PRINT_INIT_MSG;
    INIT_TIMER;
    START_TIMER;

    machine->context = (void*) malloc(sizeof(p_minus_1_context_t));
    p_minus_1_context_t* context = (p_minus_1_context_t*) machine->context;

    mpz_init_set(context->n, machine->n);
    mpz_init_set_ui(context->a, 1);
    mpz_init_set_ui(context->g, 1);
    mpz_init_set_ui(context->aq, 1);

    context->params   = (p_minus_1_params_t*) machine->params;
    context->nupdates = 0;
#if TIFA_HAVE_MONT_ECM
    context->mplan    = NULL;
    context->mplan_b1 = 0;
    context->mplan_b2 = 0;
#endif

    STOP_TIMER;
    PRINT_TIMING;

    return SUCCESS;
}
//------------------------------------------------------------------------------
static ecode_t clear_p_minus_1_context(factoring_machine_t* const machine) {
    PRINT_CLEAN_MSG;
    INIT_TIMER;
    START_TIMER;

    p_minus_1_context_t* context = (p_minus_1_context_t*) machine->context;

    if (context != NULL) {
#if TIFA_HAVE_MONT_ECM
        free_mont_ecm_plan(context->mplan);
#endif
        mpz_clear(context->n);
        mpz_clear(context->a);
        mpz_clear(context->g);
        mpz_clear(context->aq);
        free(context);
    }

    STOP_TIMER;
    PRINT_TIMING;

    return SUCCESS;
}
//------------------------------------------------------------------------------
static ecode_t update_p_minus_1_context(factoring_machine_t* const machine) {
    INIT_TIMER;
    START_TIMER;

    p_minus_1_context_t* context = (p_minus_1_context_t*) machine->context;
    p_minus_1_params_t*  params  = context->params;
    ecode_t ecode = SUCCESS;

    context->nupdates++;

    if ((context->nupdates > MAX_NUPDATES) || (params->b1 >= LARGEST_PRIME)) {
        ecode = GIVING_UP;
        PRINT_UPDATE_GIVEUP_MSG;
    } else {
        //
        // Grow the smoothness bound and retry, clamping to the largest prime
        // available in the precomputed table.
        //
        unsigned long nb1 = (unsigned long) params->b1 * B1_GROWTH_FACTOR;
        if (nb1 > LARGEST_PRIME) {
            nb1 = LARGEST_PRIME;
        }
        params->b1 = (uint32_t) nb1;

        unsigned long nb2 = (unsigned long) params->b2 * B1_GROWTH_FACTOR;
        if (nb2 > LARGEST_PRIME) {
            nb2 = LARGEST_PRIME;
        }
        params->b2 = (uint32_t) nb2;
    }

    STOP_TIMER;
    PRINT_TIMING;

    return ecode;
}
//------------------------------------------------------------------------------
static void backtrack_stage_1(mpz_t g, const mpz_t n,
                              const p_minus_1_params_t* const params) {
    //
    // Stage 1 caught every prime factor of n at once (g = n). Redo it one
    // prime at a time, checking the gcd after each, and stop as soon as it
    // splits n. If a single prime takes g from 1 to n, g is left equal to n.
    //
    mpz_t b;
    mpz_init_set_ui(b, params->base);

    uint32_t ip = 0;
    while ((ip < NFIRST_PRIMES) && (first_primes[ip] <= params->b1)) {
        unsigned long pr = first_primes[ip];
        unsigned long pe = pr;
        while (pe <= params->b1 / pr) {
            pe *= pr;
        }
        for (unsigned long q = pr; ; q *= pr) {
            mpz_powm_ui(b, b, pr, n);
            mpz_sub_ui(g, b, 1);
            mpz_gcd(g, g, n);
            if (mpz_cmp_ui(g, 1) > 0) {
                mpz_clear(b);
                return;
            }
            if (q == pe) {
                break;
            }
        }
        ip++;
    }
    mpz_set(g, n);
    mpz_clear(b);
}
//------------------------------------------------------------------------------
static ecode_t perform_p_minus_1(factoring_machine_t* const machine) {
    p_minus_1_context_t* context = (p_minus_1_context_t*) machine->context;
    p_minus_1_params_t*  params  = (p_minus_1_params_t*)  machine->params;

    INIT_TIMER;
    START_TIMER;

    mpz_ptr n = context->n;
    mpz_ptr a = context->a;
    mpz_ptr g = context->g;

    ecode_t ecode = NO_FACTOR_FOUND;

    //
    // Cheap early check: a common factor between the base and n is itself a
    // factor of n.
    //
    mpz_set_ui(g, params->base);
    mpz_gcd(g, g, n);
    if ((mpz_cmp_ui(g, 1) > 0) && (mpz_cmp(g, n) < 0)) {
        append_mpz_to_array(machine->factors, g);
        mpz_divexact(g, n, g);
        append_mpz_to_array(machine->factors, g);

        STOP_TIMER;
        PRINT_TIMING;
        return SOME_FACTORS_FOUND;
    }

#if TIFA_HAVE_MONT_ECM
    //
    // Fixed-limb Montgomery arithmetic when n has at most 4 limbs. When a
    // gcd is n, the code below redoes the run and backtracks.
    //
    if (mont_ecm_handles(n) && (mpz_cmp_ui(n, params->base) > 0)) {
        if (   (context->mplan == NULL)
            || (context->mplan_b1 != params->b1)
            || (context->mplan_b2 != params->b2)) {
            free_mont_ecm_plan(context->mplan);
            context->mplan    = alloc_mont_ecm_plan(params->b1, params->b2);
            context->mplan_b1 = params->b1;
            context->mplan_b2 = params->b2;
        }
        const int found = mont_pm1(g, n, params->base, context->mplan);
        if (found != -1) {
            if (found == 1) {
                append_mpz_to_array(machine->factors, g);
                mpz_divexact(g, n, g);
                append_mpz_to_array(machine->factors, g);
                ecode = SOME_FACTORS_FOUND;
            }
            STOP_TIMER;
            PRINT_TIMING;
            return ecode;
        }
    }
#endif

    //
    // Stage 1: compute a <- base^E mod n, where
    //     E = prod_{p prime <= b1} p^floor(log_p(b1)),
    // accumulated one prime at a time (so no huge exponent is materialized).
    //
    mpz_set_ui(a, params->base);

    uint32_t ip = 0;
    while ((ip < NFIRST_PRIMES) && (first_primes[ip] <= params->b1)) {
        unsigned long pr = first_primes[ip];
        unsigned long pe = pr;
        //
        // pe <- largest power of pr that is <= b1 (no overflow: pe*pr <= b1).
        //
        while (pe <= params->b1 / pr) {
            pe *= pr;
        }
        mpz_powm_ui(a, a, pe, n);
        ip++;
    }

    //
    // Stage 1 result: g = gcd(base^E - 1, n). Compute it via g so that
    // a = base^E is preserved for stage 2. A nontrivial gcd is a factor.
    //
    mpz_sub_ui(g, a, 1);
    mpz_gcd(g, g, n);

    if (mpz_cmp(g, n) == 0) {
        backtrack_stage_1(g, n, params);
    }
    if ((mpz_cmp_ui(g, 1) > 0) && (mpz_cmp(g, n) < 0)) {
        append_mpz_to_array(machine->factors, g);
        mpz_divexact(g, n, g);
        append_mpz_to_array(machine->factors, g);

        STOP_TIMER;
        PRINT_TIMING;
        return SOME_FACTORS_FOUND;
    }

    //
    // Stage 2: if stage 1 did not overshoot (g != n) and a second bound is
    // set, look for a factor p such that p - 1 is a b1-smooth number times a
    // single prime q in (b1, b2]. Starting from a = base^E we have
    // ord_p(a) | q for such a p, hence gcd(a^q - 1, n) reveals p. We walk the
    // primes q in (b1, b2] updating a^q by a single multiplication per prime,
    // using a precomputed table of a raised to the (even) prime gaps.
    //
    // ip already indexes the first prime > b1 (left by the stage-1 loop).
    //
    if ((mpz_cmp(g, n) != 0) && (params->b2 > params->b1)
        && (ip < NFIRST_PRIMES) && (first_primes[ip] <= params->b2)) {

        mpz_ptr aq = context->aq;

        //
        // gap_pow[k] = a^(2k) mod n, for k = 1 .. MAX_GAP_HALF.
        //
        mpz_t gap_pow[MAX_GAP_HALF + 1];
        mpz_init(gap_pow[1]);
        mpz_mul(gap_pow[1], a, a);
        mpz_mod(gap_pow[1], gap_pow[1], n);
        for (int k = 2; k <= MAX_GAP_HALF; k++) {
            mpz_init(gap_pow[k]);
            mpz_mul(gap_pow[k], gap_pow[k - 1], gap_pow[1]);
            mpz_mod(gap_pow[k], gap_pow[k], n);
        }

        mpz_t t;
        mpz_t H;
        mpz_t aq0;
        mpz_init(t);
        mpz_init_set_ui(H, 1);
        mpz_init(aq0);

        //
        // Bootstrap: aq = a^{q0} mod n for the first prime q0 > b1.
        //
        unsigned long prev = first_primes[ip];
        mpz_powm_ui(aq, a, prev, n);

        //
        // Snapshot of the current batch's starting state, so a batch whose
        // product shares all of n can be rescanned prime-by-prime.
        //
        uint32_t      batch_ip   = ip;
        unsigned long batch_prev = prev;
        mpz_set(aq0, aq);
        uint32_t count = 0;

        while (true) {
            //
            // Accumulate (aq - 1) into the running product H (mod n); the gcd
            // is deferred to once per batch.
            //
            mpz_sub_ui(t, aq, 1);
            mpz_mul(H, H, t);
            mpz_mod(H, H, n);
            count++;

            //
            // Advance aq to the next prime <= b2 (aq <- aq * a^{gap} mod n),
            // if there is one.
            //
            bool more = false;
            if ((ip + 1 < NFIRST_PRIMES) && (first_primes[ip + 1] <= params->b2)) {
                unsigned long cur  = first_primes[ip + 1];
                unsigned long gap  = cur - prev;
                unsigned long half = gap >> 1;
                if (half <= (unsigned long) MAX_GAP_HALF) {
                    mpz_mul(aq, aq, gap_pow[half]);
                } else {
                    mpz_powm_ui(t, a, gap, n);
                    mpz_mul(aq, aq, t);
                }
                mpz_mod(aq, aq, n);
                prev = cur;
                ip++;
                more = true;
            }

            if ((count == STAGE2_BATCH) || !more) {
                mpz_gcd(g, H, n);
                if (mpz_cmp_ui(g, 1) > 0) {
                    if (mpz_cmp(g, n) < 0) {
                        append_mpz_to_array(machine->factors, g);
                        mpz_divexact(g, n, g);
                        append_mpz_to_array(machine->factors, g);
                        ecode = SOME_FACTORS_FOUND;
                    } else {
                        //
                        // The batch product is a multiple of n: rescan this
                        // batch prime-by-prime to isolate a proper factor (the
                        // rare case where >= 2 factors were caught in one batch).
                        //
                        mpz_set(aq, aq0);
                        unsigned long rprev = batch_prev;
                        uint32_t rip = batch_ip;
                        for (uint32_t j = 0; j < count; j++) {
                            mpz_sub_ui(t, aq, 1);
                            mpz_gcd(g, t, n);
                            if ((mpz_cmp_ui(g, 1) > 0) && (mpz_cmp(g, n) < 0)) {
                                append_mpz_to_array(machine->factors, g);
                                mpz_divexact(g, n, g);
                                append_mpz_to_array(machine->factors, g);
                                ecode = SOME_FACTORS_FOUND;
                                break;
                            }
                            if (j + 1 == count) {
                                break;
                            }
                            unsigned long cur  = first_primes[rip + 1];
                            unsigned long gap  = cur - rprev;
                            unsigned long half = gap >> 1;
                            if (half <= (unsigned long) MAX_GAP_HALF) {
                                mpz_mul(aq, aq, gap_pow[half]);
                            } else {
                                mpz_powm_ui(t, a, gap, n);
                                mpz_mul(aq, aq, t);
                            }
                            mpz_mod(aq, aq, n);
                            rprev = cur;
                            rip++;
                        }
                    }
                    break;
                }
                if (!more) {
                    break;
                }
                //
                // No factor in this batch: reset and snapshot the next one.
                //
                mpz_set_ui(H, 1);
                count      = 0;
                batch_ip   = ip;
                batch_prev = prev;
                mpz_set(aq0, aq);
            }
        }

        mpz_clear(t);
        mpz_clear(H);
        mpz_clear(aq0);
        for (int k = 1; k <= MAX_GAP_HALF; k++) {
            mpz_clear(gap_pow[k]);
        }
    }

    STOP_TIMER;
    PRINT_TIMING;

    return ecode;
}
//------------------------------------------------------------------------------
static ecode_t recurse(mpz_array_t* const factors, uint32_array_t* const multis,
                       const mpz_t n, factoring_mode_t mode) {

    return tifa_factor(factors, multis, n, mode);

    //
    // The following code should be used to factor a number using _only_
    // Pollard's p-1.
    //
    //p_minus_1_params_t params;
    //set_p_minus_1_params_to_default(n, &params);
    //return p_minus_1(factors, multis, n, &params, mode);
}
//------------------------------------------------------------------------------
