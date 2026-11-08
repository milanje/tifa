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
 * \file    squfof.c
 * \author  Jerome Milan
 * \date    Tue Oct 11 2011
 * \version 2011-10-11
 */

  //------------------------------------------------------------------------
  // References:
  //------------------------------------------------------------------------
  // - For a thorough explanation of the SQUFOF algorithm:
  //         "Square Form Factorization",
  //         Jason E. Gower, Samuel S. Wagstaff Jr.,
  //         Mathematics of Computation,
  //         S 0025-5718(07)02010-8,
  //         Article electronically published on May 14, 2007.
  //
  // - For an alternate explanation, including Shanks' very notes on SQUFOF,
  //   see David Joyner and Stephen McMath's work at:
  //       http://www.usna.edu/Users/math/wdj/mcmath/
  //   URL last retrieved on Tue 15 Mar 2011.
  //
  // - For a description of the "large step algorithm" used to quickly
  //   jump over several forms:
  //         "On the Parallel Generation of the Residues for the Continued
  //         Fraction Factoring Algorithm",
  //         Hugh C. Williams, Marvin C. Wunderlich,
  //         Mathematics of Computation, Volume 48, Number 177,
  //         January 1987, pages 405-423
  //------------------------------------------------------------------------

#include "tifa_config.h"

#include <stdlib.h>
#include <stdbool.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include <gmp.h>

#include "funcs.h"
#include "macros.h"
#include "factoring_machine.h"
#include "squfof.h"
#include "impl.h"
#include "tifa_factor.h"
#include "tifa_alloc.h"

#define SQUFOF_DEBUG 0
#if SQUFOF_DEBUG
  #include <stdio.h>
  #define SQUFOF_GMP_PRINT(...) gmp_printf(__VA_ARGS__)
  #define SQUFOF_PRINT(...)     printf(__VA_ARGS__)
#else
  #define SQUFOF_GMP_PRINT(...) /* intentionally left empty */
  #define SQUFOF_PRINT(...)     /* intentionally left empty */
#endif

#define __PREFIX__    "squfof: "
#define __VERBOSE__   TIFA_VERBOSE_SQUFOF
#define __TIMING__    TIFA_TIMING_SQUFOF
#define __EXTRA_MSG__ 1

#include "messages.h"

#include "stopwatch.h"

#define SUBTIMINGS 0

#if SUBTIMINGS
    static stopwatch_t timer_fwd;
    static int         timer_fwd_is_init = 0;
    static stopwatch_t timer_rev;
    static int         timer_rev_is_init = 0;
#endif

///init_stopwatch(&timer_std);

//-----------------------------------------------------------------------------
//                       GENERAL NON PUBLIC DEFINE(S)
//-----------------------------------------------------------------------------
//
// Maximum size of the queues.
//
#define MAX_QUEUE_ALLOC 64
//
// Number of multipliers raced together. Their forward cycles run in the same
// loop, so that the divisions and products of each overlap with the others'.
//
#define NRACERS 8
//
// Number of forward steps between two checks of the abort bounds (even).
//
#define NSTEPS_ROUND 32
//
// The loops over the racers must be unrolled for the racers' forms to stay
// in registers, and the rare queue and square checks kept out of them.
//
#define PRAGMA(X) _Pragma(#X)
#if defined(__GNUC__)
    #define UNROLLED(N) PRAGMA(GCC unroll N)
    #define COLD        __attribute__ ((noinline, cold))
#else
    #define UNROLLED(N) /* intentionally left empty */
    #define COLD        /* intentionally left empty */
#endif

//-----------------------------------------------------------------------------
//           SYMBOLS AND FUNCTIONS RELATED TO THE LARGE STEP ALGORITHM
//-----------------------------------------------------------------------------

//
// _NOTE_: We only implement here a simplified version of Williams &
//         Wunderlich's large step algorithm. This is used in step 4 of the
//         SQUFOF algorithm (using Gower & Wagstaff's description) to achieve
//         what Shanks calls the "fast return". Note that Gower & Wagstaff's
//         description does not include this "fast return".
//

//
// Should we perform the "fast return" variant?
//
// _NOTE_: The "fast return" variant is now enabled by default since the
//         mighty infinite loop bug has been solved (hopefully for good!).
//
#define PERFORM_FAST_RETURN 1

//
// Only use the large step algorithm if the number of forms to compute (from
// the inverse square root of a found square form) is larger than a certain
// threshold.
//
#if PERFORM_FAST_RETURN
    #define LARGE_STEP_THRESHOLD 4096
#else
    #define LARGE_STEP_THRESHOLD TIFA_LONG_MAX
#endif

//
// Mininum number of forms to step through sequentially before using the
// large step algorithm. If MIN_NSTEP_LS is too small, the large step's
// overhead will degrade performances. If set too high, we'll be less
// precise in the length's estimation of the large step.
//
#define MIN_NSTEP_LS 128
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//                           NON PUBLIC DATA
//-----------------------------------------------------------------------------

//
// Multipliers, in order of use: 1 and the squarefree divisors of
// 1155 = 3 * 5 * 7 * 11. The race starts with the first NRACERS, and the
// next one replaces each racer that stops.
//
static const unsigned int multipliers[] = {
       1,   3,   5,   7,  11,   15,
      21,  33,  35,  55,  77,  105,
     165, 231, 385, 1155
};
#define NMULTIPLIERS   (sizeof(multipliers) / sizeof(multipliers[0]))
#define MAX_MULTIPLIER 1155
//
// Squares modulo 64, 221 = 13 * 17 and 437 = 19 * 23, as bit sets. At the
// even steps, the Q are values of the principal form, hence squares modulo
// the primes dividing D: only moduli prime to the multipliers sift them.
// About 2.5% of them are squares modulo all three.
//
static const uint64_t squares_mod_64 = UINT64_C(0x0202021202030213);

static const uint64_t squares_mod_221[4] = {
    UINT64_C(0x00ba0c5846032213), UINT64_C(0x836001b040822035),
    UINT64_C(0x868c17402b011040), UINT64_C(0x0000000012113018)
};
static const uint64_t squares_mod_437[7] = {
    UINT64_C(0x44c2809807810253), UINT64_C(0x0ac80131b0a62201),
    UINT64_C(0x5802020a10415c28), UINT64_C(0x6100870211a68030),
    UINT64_C(0x0892605224208401), UINT64_C(0x32150242bc104518),
    UINT64_C(0x000008c58019a085)
};

//-----------------------------------------------------------------------------
//                   NON PUBLIC STRUCTURE(S) AND TYPEDEF(S)
//-----------------------------------------------------------------------------
//
// A SQUFOF racer: the cycle of the forms of discriminant D = k.n (or 2.k.n)
// for a multiplier k. Notations follow Gower and Wagstaff's paper.
//
struct struct_squfof_racer_t {
    unsigned long int multiplier;
    unsigned long int S;
    //
    // The first form is (1, P0, Q0), the current one (QQ, P, Q).
    //
    unsigned long int P0;
    unsigned long int Q0;
    unsigned long int P;
    unsigned long int Q;
    unsigned long int QQ;
    //
    // Pairs (x, y) only go to the queue for x <= L = floor(sqrt(2S + 1)):
    // since Q < 2 sqrt(D) < 2S + 2, larger x never match the root of a Q.
    // With S <= ULONG_MAX / 6, x and y < x fit in 32 bits.
    //
    unsigned long int L;
    //
    // Q / gcd(Q, 2k) <= L only if Q <= QL = 2k.L.
    //
    unsigned long int QL;
    //
    // Abort bound on the number of forward steps.
    //
    unsigned long int B;
    //
    // Number of forward steps of the race when this racer joined it.
    //
    unsigned long int start;

    uint32_t queue_x[MAX_QUEUE_ALLOC];
    uint32_t queue_y[MAX_QUEUE_ALLOC];
    uint32_t queue_head;
    uint32_t queue_length;

    bool running;
};
//-----------------------------------------------------------------------------
typedef struct struct_squfof_racer_t squfof_racer_t;
//-----------------------------------------------------------------------------
//
// An ad-hoc structure holding the data variables used by the SQUFOF
// implementation.
//
struct struct_squfof_context_t {
    squfof_params_t* params;
    mpz_ptr n;
    //
    // n, if 2.k.n fits in an unsigned long for all the multipliers k.
    //
    unsigned long int n_ui;

    squfof_racer_t racers[NRACERS];
    //
    // Number of racers set up so far, and index of the next multiplier.
    //
    uint32_t nracers;
    uint32_t next;
    //
    // Number of forward steps of the race.
    //
    unsigned long int nsteps;
    //
    // Temporaries: the discriminant of a racer, its square root and the
    // remainder, when 2.k.n does not fit in an unsigned long.
    //
    mpz_t D;
    mpz_t root;
    mpz_t rem;
};
//-----------------------------------------------------------------------------
typedef struct struct_squfof_context_t squfof_context_t;
//-----------------------------------------------------------------------------
//
// An ad-hoc structure holding temporary multiprecision integers used by
// the "fast return" variant.
//
struct struct_tmp_vars_t {
    mpz_t tmp_1;
    mpz_t tmp_2;
    mpz_t tmp_3;
    mpz_t tmp_4;
    mpz_t tmp_5;
    mpz_t tmp_6;
    mpz_t tmp_7;
    mpz_t tmp_8;

    mpz_ptr D;
    mpz_ptr N;

    unsigned long int S;
};
//-----------------------------------------------------------------------------
typedef struct struct_tmp_vars_t tmp_vars_t;
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//               PROTOTYPES OF SQUFOF'S NON PUBLIC FUNCTION(S)
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
static tmp_vars_t* alloc_init_tmp_vars(unsigned long int S, mpz_t D, mpz_t N);
//-----------------------------------------------------------------------------
static void free_tmp_vars(tmp_vars_t* vars);
//-----------------------------------------------------------------------------
static ecode_t init_squfof_context(factoring_machine_t* const machine);
//-----------------------------------------------------------------------------
static ecode_t clear_squfof_context(factoring_machine_t* const machine);
//-----------------------------------------------------------------------------
static ecode_t update_squfof_context(factoring_machine_t* const machine);
//-----------------------------------------------------------------------------
static ecode_t perform_squfof(factoring_machine_t* const machine);
//-----------------------------------------------------------------------------
static bool init_squfof_racer(squfof_racer_t* const racer,
                              squfof_context_t* const context,
                              unsigned long int multiplier);
//-----------------------------------------------------------------------------
static void next_racer(squfof_context_t* const context,
                       squfof_racer_t* const racer);
//-----------------------------------------------------------------------------
inline static unsigned long int gcd_2k(unsigned long int x,
                                       unsigned long int k);
//-----------------------------------------------------------------------------
static int race(squfof_context_t* const context);
//-----------------------------------------------------------------------------
static void enqueue(squfof_racer_t* const racer, unsigned long int P,
                    unsigned long int Q);
//-----------------------------------------------------------------------------
static bool proper_square(squfof_racer_t* const racer, unsigned long int P,
                          unsigned long int Q);
//-----------------------------------------------------------------------------
static unsigned long int reverse_cycle(squfof_racer_t* const racer,
                                       squfof_context_t* const context);
//-----------------------------------------------------------------------------
inline static void
single_step(unsigned long int* const P,  unsigned long int* const Q,
            unsigned long int* const QQ, const unsigned long int S
);
//-----------------------------------------------------------------------------
inline static bool
reverse_step(unsigned long int* const P,  unsigned long int* const Q,
             unsigned long int* const QQ, const unsigned long int S,
             unsigned long int* const factor
);
//-----------------------------------------------------------------------------
static unsigned long int
large_step(unsigned long int* const P,  unsigned long int* const Q,
           unsigned long int* const QQ, unsigned long int PA,
           unsigned long int QA,        unsigned long int QQA,
           unsigned long int PB,        unsigned long int QB,
           unsigned long int QQB,       tmp_vars_t* const vars
);
//-----------------------------------------------------------------------------
inline static unsigned long int
doubling_algo(unsigned long int* const P2,  unsigned long int* const Q2,
              unsigned long int* const QQ2, unsigned long int P,
              unsigned long int Q,          unsigned long int QQ,
              tmp_vars_t* const vars
);
//-----------------------------------------------------------------------------
inline static unsigned long int
jump_over_k_steps(unsigned long int* const P,  unsigned long int* const Q,
                  unsigned long int* const QQ, unsigned long int P0,
                  unsigned long int Q0,        unsigned long int QQ0,
                  unsigned long int K,         tmp_vars_t* const vars
);
//-----------------------------------------------------------------------------
static unsigned long int
jump_near_sym_point(unsigned long int* const P,  unsigned long int* const Q,
                    unsigned long int* const QQ, unsigned long int P0,
                    unsigned long int Q0,        unsigned long int QQ0,
                    unsigned long int S,         mpz_t D,
                    unsigned long int K,         mpz_t N
);
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//                 "PUBLIC" FUNCTION(S) IMPLEMENTATION
//-----------------------------------------------------------------------------

inline double squfof_get_elapsed_fwd() {
    return 0.0;
    //return get_stopwatch_elapsed(&timer_fwd);
}

inline double squfof_get_elapsed_rev() {
    return 0.0;
    //return get_stopwatch_elapsed(&timer_rev);
}


//-----------------------------------------------------------------------------
void set_squfof_params_to_default(squfof_params_t* const params
                                  MAYBE_UNUSED) {
    //
    // While there is no user parameters for this SQUFOF implementation, this
    // function is kept as a placeholder: it might be needed in future code
    // revisions.
    //
    return;
}
//-----------------------------------------------------------------------------
static ecode_t recurse(mpz_array_t* const factors, uint32_array_t* const multis,
                       const mpz_t n, factoring_mode_t mode) {

    return tifa_factor_rec(factors, multis, n, mode);

    //
    // The following code should be used to factor a number using _only_ SQUFOF.
    //
    //squfof_params_t params;
    //set_squfof_params_to_default(&params);
    //return squfof(factors, multis, n, &params, mode);
}
//-----------------------------------------------------------------------------
ecode_t squfof_impl(mpz_array_t* const factors, uint32_array_t* const multis,
                    const mpz_t n, const squfof_params_t* const params,
                    const factoring_mode_t mode) {

    PRINT_INTRO_MSG("square form");
    INIT_TIMER;
    START_TIMER;

    ecode_t ecode = NO_FACTOR_FOUND;

    factoring_machine_t machine;

    mpz_init_set(machine.n, n);

    machine.mode                = mode;
    machine.params              = (void*) params;
    machine.init_context_func   = init_squfof_context;
    machine.perform_algo_func   = perform_squfof;
    machine.update_context_func = update_squfof_context;
    machine.clear_context_func  = clear_squfof_context;
    machine.recurse_func        = recurse;
    machine.factors             = factors;
    machine.multis              = multis;

    ecode = run_machine(&machine);

    mpz_clear(machine.n);

    STOP_TIMER;
    PRINT_STATUS(machine.success, ecode);

    return ecode;
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//                  "PRIVATE" FUNCTION(S) IMPLEMENTATION
//-----------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *                  Functions used by factoring machine
 *-----------------------------------------------------------------------------
 */

//-----------------------------------------------------------------------------
static ecode_t init_squfof_context(factoring_machine_t* const machine) {
    //
    // Init the SQUFOF implementation specific variables
    //
    PRINT_INIT_MSG;
    INIT_TIMER;
    START_TIMER;

    machine->context          = (void*) tifa_malloc(sizeof(squfof_context_t));
    squfof_context_t* context = (squfof_context_t*) machine->context;

    context->params  = machine->params;
    context->n       = machine->n;
    context->nracers = 0;
    context->next    = 0;
    context->nsteps  = 0;
    context->n_ui    = 0;

    if (   mpz_fits_ulong_p(context->n)
        && (mpz_get_ui(context->n) <= ULONG_MAX / (2 * MAX_MULTIPLIER))) {
        context->n_ui = mpz_get_ui(context->n);
    }
    mpz_init(context->D);
    mpz_init(context->root);
    mpz_init(context->rem);

    for (uint32_t i = 0; i < NRACERS; i++) {
        next_racer(context, &context->racers[i]);
    }
    STOP_TIMER;
    PRINT_TIMING;

#if SUBTIMINGS
    if (!timer_fwd_is_init) {
        init_stopwatch(&timer_fwd);
        timer_fwd_is_init++;
    }
    if (!timer_rev_is_init) {
        init_stopwatch(&timer_rev);
        timer_rev_is_init++;
    }
#endif

    return SUCCESS;
}
//-----------------------------------------------------------------------------
static ecode_t clear_squfof_context(factoring_machine_t* const machine) {

    PRINT_CLEAN_MSG;
    INIT_TIMER;
    START_TIMER;

    squfof_context_t* context = (squfof_context_t*) machine->context;

    mpz_clear(context->D);
    mpz_clear(context->root);
    mpz_clear(context->rem);
    free(context);

    STOP_TIMER;
    PRINT_TIMING;

    return SUCCESS;
}
//-----------------------------------------------------------------------------
static ecode_t update_squfof_context(factoring_machine_t* const machine) {
    //
    // The race only ends once all the racers have stopped: there is nothing
    // left to try.
    //
    squfof_context_t* const context = (squfof_context_t*) machine->context;

    PRINT_UPDATE_GIVEUP_MSG;

    return (context->nracers == 0) ? INTEGER_TOO_LARGE : GIVING_UP;
}
//-----------------------------------------------------------------------------
static ecode_t perform_squfof(factoring_machine_t* const machine) {

    INIT_TIMER;

    squfof_context_t* const context = (squfof_context_t*) machine->context;

    if (context->nracers == 0) {
        return INTEGER_TOO_LARGE;
    }
    while (true) {
        PRINT_FWD_CYCL_MSG;
        START_TIMER;

        int winner = race(context);

        if (winner < 0) {
            STOP_TIMER;
            PRINT_ALL_RACERS_FAILED;
            PRINT_TIMING_ON_FAILURE;
            return GIVING_UP;
        }
        STOP_TIMER;
        PRINT_TIMING;
        PRINT_REV_CYCL_MSG;
        RESET_TIMER;
        START_TIMER;

        squfof_racer_t* const racer = &context->racers[winner];

        unsigned long int factor = reverse_cycle(racer, context);

        STOP_TIMER;
        PRINT_TIMING;
        //
        // Step 5: the factor found divides 2.k.n
        //
        factor /= gcd_2k(factor, racer->multiplier);

        if ((factor > 1) && (mpz_cmp_ui(context->n, factor) > 0)) {
            mpz_t factor_z;
            mpz_init_set_ui(factor_z, factor);
            append_mpz_to_array(machine->factors, factor_z);
            mpz_divexact(factor_z, context->n, factor_z);
            append_mpz_to_array(machine->factors, factor_z);
            mpz_clear(factor_z);

            return SOME_FACTORS_FOUND;
        }
        //
        // A trivial factor: let the other racers go on.
        //
        racer->running = false;
        RESET_TIMER;
    }
}
//-----------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *              Step 1) of the SQUFOF algorithm: the racers
 *-----------------------------------------------------------------------------
 */

//-----------------------------------------------------------------------------
static void set_discriminant(mpz_t D, const mpz_t n,
                             unsigned long int multiplier) {
    mpz_mul_ui(D, n, multiplier);
    if ((mpz_get_ui(D) & 3) == 1) {
        mpz_mul_2exp(D, D, 1);
    }
}
//-----------------------------------------------------------------------------
static unsigned long int floor_sqrt(unsigned long int x) {
    const unsigned long int rmax = ULONG_MAX >> (4 * sizeof(unsigned long int));

    unsigned long int r = (unsigned long int) sqrt((double) x);

    if (r > rmax) {
        r = rmax;
    }
    while (r * r > x) {
        r--;
    }
    while ((r < rmax) && ((r + 1) * (r + 1) <= x)) {
        r++;
    }
    return r;
}
//-----------------------------------------------------------------------------
inline static unsigned long int gcd_2k(unsigned long int x,
                                       unsigned long int k) {
    //
    // gcd(x, 2k) for a multiplier k: as 2k is a squarefree product of 2, 3,
    // 5, 7 and 11, the product of those of its primes that divide x.
    //
    unsigned long int d = IS_EVEN(x) ? 2 : 1;

    if ((k % 3 == 0) && (x % 3 == 0)) {
        d *= 3;
    }
    if ((k % 5 == 0) && (x % 5 == 0)) {
        d *= 5;
    }
    if ((k % 7 == 0) && (x % 7 == 0)) {
        d *= 7;
    }
    if ((k % 11 == 0) && (x % 11 == 0)) {
        d *= 11;
    }
    return d;
}
//-----------------------------------------------------------------------------
static bool init_squfof_racer(squfof_racer_t* const racer,
                              squfof_context_t* const context,
                              unsigned long int multiplier) {
    //
    // Sets up the racer of the given multiplier, unless k.n is a square or
    // too large for single precision.
    //
    unsigned long int S  = 0;
    unsigned long int Q0 = 0;

    if (context->n_ui != 0) {
        unsigned long int D = multiplier * context->n_ui;
        if ((D & 3) == 1) {
            D <<= 1;
        }
        S  = floor_sqrt(D);
        Q0 = D - S * S;
    } else {
        set_discriminant(context->D, context->n, multiplier);
        mpz_sqrtrem(context->root, context->rem, context->D);

        if (mpz_cmp_ui(context->root, TIFA_ULONG_MAX_DIVIDED_BY_6) > 0) {
            return false;
        }
        S  = mpz_get_ui(context->root);
        Q0 = mpz_get_ui(context->rem);
    }
    if (Q0 == 0) {
        return false;
    }
    racer->multiplier   = multiplier;
    racer->S            = S;
    racer->P0           = S;
    racer->Q0           = Q0;
    racer->P            = racer->P0;
    racer->Q            = racer->Q0;
    racer->QQ           = 1;
    racer->L            = floor_sqrt(2 * racer->S + 1);
    racer->QL           = 2 * multiplier * racer->L;
    //
    // Gower and Wagstaff's 2 floor(2 sqrt(2 sqrt(D))), give or take a few.
    //
    racer->B            = 4 * racer->L;
    racer->start        = context->nsteps;
    racer->queue_head   = 0;
    racer->queue_length = 0;
    racer->running      = true;

    return true;
}
//-----------------------------------------------------------------------------
static void next_racer(squfof_context_t* const context,
                       squfof_racer_t* const racer) {
    //
    // Sets up in the slot of a stopped racer the racer of the next usable
    // multiplier. With none left, the slot cycles through (1, 1, 1), a form
    // of discriminant 2, while the race goes on.
    //
    while (context->next < NMULTIPLIERS) {
        unsigned long int k = multipliers[context->next++];
        if (init_squfof_racer(racer, context, k)) {
            context->nracers++;
            return;
        }
    }
    racer->S       = 1;
    racer->P       = 1;
    racer->Q       = 1;
    racer->QQ      = 1;
    racer->QL      = 0;
    racer->running = false;
}
//-----------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *              Step 2) of the SQUFOF algorithm: forward cycle
 *-----------------------------------------------------------------------------
 */

//-----------------------------------------------------------------------------
inline static void forward_step(unsigned long int* const P,
                                unsigned long int* const Q,
                                unsigned long int* const QQ,
                                const unsigned long int S,
                                const unsigned long int QL,
                                squfof_racer_t* const racer) {
    unsigned long int q  = (S + *P) / *Q;
    unsigned long int PP = q * (*Q) - *P;

    if (IS_UNLIKELY(*Q <= QL)) {
        enqueue(racer, *P, *Q);
    }
    unsigned long int t = *QQ + q * (*P - PP);

    *QQ = *Q;
    *Q  = t;
    *P  = PP;
}
//-----------------------------------------------------------------------------
inline static bool maybe_square(const unsigned long int x) {
    unsigned long int r221 = x % 221;
    unsigned long int r437 = x % 437;

    return 0 != (  (squares_mod_64 >> (x & 63))
                 & (squares_mod_221[r221 >> 6] >> (r221 & 63))
                 & (squares_mod_437[r437 >> 6] >> (r437 & 63))
                 & 1);
}
//-----------------------------------------------------------------------------
static int race(squfof_context_t* const context) {
    //
    // Gower and Wagstaff's step 2 for all the racers at once, until one of
    // them finds a proper square form. Returns its index, or -1 once they
    // have all stopped.
    //
    squfof_racer_t* const racers = context->racers;

    unsigned long int S[NRACERS];
    unsigned long int P[NRACERS];
    unsigned long int Q[NRACERS];
    unsigned long int QQ[NRACERS];
    unsigned long int QL[NRACERS];

    for (uint32_t j = 0; j < NRACERS; j++) {
        S[j]  = racers[j].S;
        P[j]  = racers[j].P;
        Q[j]  = racers[j].Q;
        QQ[j] = racers[j].QQ;
        QL[j] = racers[j].QL;
    }
    unsigned long int nsteps = context->nsteps;
    int winner = -1;

    if (!IS_EVEN(nsteps)) {
        //
        // The last race stopped after the first step of a pair.
        //
        UNROLLED(NRACERS)
        for (uint32_t j = 0; j < NRACERS; j++) {
            forward_step(&P[j], &Q[j], &QQ[j], S[j], QL[j], &racers[j]);
        }
        nsteps++;
    }
    while (winner < 0) {
        bool running = false;
        for (uint32_t j = 0; j < NRACERS; j++) {
            squfof_racer_t* const racer = &racers[j];

            if (racer->running && (nsteps - racer->start > racer->B)) {
                racer->running = false;
                PRINT_UPPER_BOUND_REACHED_MULT(racer->multiplier);
            }
            if (!racer->running && (context->next < NMULTIPLIERS)) {
                context->nsteps = nsteps;
                next_racer(context, racer);
                S[j]  = racer->S;
                P[j]  = racer->P;
                Q[j]  = racer->Q;
                QQ[j] = racer->QQ;
                QL[j] = racer->QL;
            }
            running = running || racer->running;
        }
        if (!running) {
            break;
        }
        for (uint32_t s = 0; s < NSTEPS_ROUND; s += 2) {
            UNROLLED(NRACERS)
            for (uint32_t j = 0; j < NRACERS; j++) {
                forward_step(&P[j], &Q[j], &QQ[j], S[j], QL[j], &racers[j]);
            }
            nsteps++;
            //
            // Square forms only count after an odd number of steps.
            //
            UNROLLED(NRACERS)
            for (uint32_t j = 0; j < NRACERS; j++) {
                if (   IS_UNLIKELY(maybe_square(Q[j]))
                    && proper_square(&racers[j], P[j], Q[j])) {
                    winner = (int) j;
                    break;
                }
            }
            if (winner >= 0) {
                break;
            }
            UNROLLED(NRACERS)
            for (uint32_t j = 0; j < NRACERS; j++) {
                forward_step(&P[j], &Q[j], &QQ[j], S[j], QL[j], &racers[j]);
            }
            nsteps++;
        }
    }
    for (uint32_t j = 0; j < NRACERS; j++) {
        racers[j].P  = P[j];
        racers[j].Q  = Q[j];
        racers[j].QQ = QQ[j];
    }
    context->nsteps = nsteps;

    if (winner >= 0) {
        PRINT_SUCCESS_MULT(racers[winner].multiplier);
    }
    return winner;
}
//-----------------------------------------------------------------------------
COLD static void enqueue(squfof_racer_t* const racer, unsigned long int P,
                         unsigned long int Q) {
    //
    // Gower and Wagstaff's step 2b), for a Q <= QL: the pair
    // (g, P mod g) goes to the queue if g = Q / gcd(Q, 2k) <= L.
    //
    unsigned long int d = gcd_2k(Q, racer->multiplier);

    if ((Q > d * racer->L) || !racer->running) {
        return;
    }
    if (racer->queue_length == MAX_QUEUE_ALLOC) {
        racer->running = false;
        PRINT_QUEUE_OVERFLOW_MULT(racer->multiplier);
        return;
    }
    unsigned long int g    = Q / d;
    uint32_t          tail =   (racer->queue_head + racer->queue_length)
                             % MAX_QUEUE_ALLOC;

    racer->queue_x[tail] = (uint32_t) g;
    racer->queue_y[tail] = (uint32_t) (P % g);
    racer->queue_length++;
}
//-----------------------------------------------------------------------------
COLD static bool proper_square(squfof_racer_t* const racer,
                               unsigned long int P, unsigned long int Q) {
    //
    // Gower and Wagstaff's step 2d): whether Q = r^2 with no pair (r, y) in
    // the queue such that r divides P - y. Such a pair is removed, with all
    // the pairs before it.
    //
    if (!racer->running) {
        return false;
    }
    unsigned long int r = is_square(Q);

    if (r == 0) {
        return false;
    }
    for (uint32_t i = 0; i < racer->queue_length; i++) {
        uint32_t j = (racer->queue_head + i) % MAX_QUEUE_ALLOC;

        if (racer->queue_x[j] != r) {
            continue;
        }
        if (r == 1) {
            racer->running = false;
            PRINT_NO_PROPER_FORM_FOUND_MULT(racer->multiplier);
            return false;
        }
        if ((P % r) == racer->queue_y[j]) {
            racer->queue_head    = (j + 1) % MAX_QUEUE_ALLOC;
            racer->queue_length -= i + 1;
            return false;
        }
    }
    return true;
}
//-----------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *          Steps 3) and 4) of the SQUFOF algorithm: reverse cycle
 *-----------------------------------------------------------------------------
 */

//-----------------------------------------------------------------------------
static unsigned long int reverse_cycle(squfof_racer_t* const racer,
                                       squfof_context_t* const context) {
    //
    // Returns the factor of D given by the proper square form (QQ, P, Q)
    // that the racer just found.
    //
    unsigned long int S = racer->S;
    unsigned long int r = is_square(racer->Q);
    unsigned long int m = (S - racer->P) / r;
    //
    // Step 3: the inverse square root (r, P + m.r, Q) of the square form.
    // As D = P^2 + r^2.QQ, Q = (D - (P + m.r)^2) / r = QQ.r - m.(2P + m.r),
    // exact modulo the word size.
    //
    unsigned long int QQ = r;
    unsigned long int P  = racer->P + m * r;
    unsigned long int Q  = racer->QQ * r - m * (2 * racer->P + m * r);

    unsigned long int factor = 0;
    //
    // About as many forms separate the inverse square root from the point
    // of symmetry as the square form from the first form, half as many as
    // the steps taken.
    //
    unsigned long int half = (context->nsteps - racer->start - 1) / 2;

    if (half < LARGE_STEP_THRESHOLD) {
        //
        // Step 4: Cycle in the reverse direction to find a factor of N
        //
        while (!reverse_step(&P, &Q, &QQ, S, &factor)) { }
    } else {
        //
        // Use the large step algorithm to get closer to the point of symmetry.
        //
        set_discriminant(context->D, context->n, racer->multiplier);

        factor = jump_near_sym_point(&P, &Q, &QQ, racer->P0, racer->Q0, 1,
                                     S, context->D, half, context->n);
        if (factor == 0) {
            //
            // Step 4: Cycle in _both_ directions to find a factor of N
            //
            unsigned long int PR  = P;
            unsigned long int QR  = QQ;
            unsigned long int QQR = Q;

            while (true) {
                if (reverse_step(&P,  &Q,  &QQ,  S, &factor)) break;
                if (reverse_step(&PR, &QR, &QQR, S, &factor)) break;
            }
        }
    }
    return factor;
}
//-----------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *                 Functions used for the large step algorithm
 *-----------------------------------------------------------------------------
 */

//-----------------------------------------------------------------------------
static tmp_vars_t* alloc_init_tmp_vars(unsigned long int S, mpz_t D, mpz_t N) {
    tmp_vars_t* vars = tifa_malloc(sizeof(tmp_vars_t));

    mpz_init(vars->tmp_1);
    mpz_init(vars->tmp_2);
    mpz_init(vars->tmp_3);
    mpz_init(vars->tmp_4);
    mpz_init(vars->tmp_5);
    mpz_init(vars->tmp_6);
    mpz_init(vars->tmp_7);
    mpz_init(vars->tmp_8);

    vars->S = S;
    vars->D = D;
    vars->N = N;

    return vars;
}
//-----------------------------------------------------------------------------
static void free_tmp_vars(tmp_vars_t* vars) {
    mpz_clear(vars->tmp_1);
    mpz_clear(vars->tmp_2);
    mpz_clear(vars->tmp_3);
    mpz_clear(vars->tmp_4);
    mpz_clear(vars->tmp_5);
    mpz_clear(vars->tmp_6);
    mpz_clear(vars->tmp_7);
    mpz_clear(vars->tmp_8);
    free(vars);
}
//-----------------------------------------------------------------------------
inline static void single_step(unsigned long int* const P,
                               unsigned long int* const Q,
                               unsigned long int* const QQ,
                               const unsigned long int S) {
    //
    // Performs the single step algorithm. This is actually just a step in
    // a continued fraction expansion and is used to step from a binary
    // quadratic form to the next adjacent one (notations follow Gower's &
    // Wagstaff's).
    //
    unsigned long int  q  = 0;
    unsigned long int  PP = 0;
    unsigned long int  t  = 0;

    q = (S + *P) / *Q;
    PP   =  q;
    PP  *= *Q;
    PP  -= *P;
    t    = *P;
    t   -=  PP;
    t   *=  q;
    t   += *QQ;
    *QQ  = *Q;
    *Q   =  t;
    *P   =  PP;
}
//-----------------------------------------------------------------------------
inline static bool reverse_step(unsigned long int* const P,
                                unsigned long int* const Q,
                                unsigned long int* const QQ,
                                const unsigned long int S,
                                unsigned long int* const factor) {
    //
    // One step of the reverse cycle: this is single_step with a symmetry-point
    // check inserted before the state update. Returns true (and sets *factor)
    // when P == PP, i.e. the point of symmetry giving a factor of N.
    //
    unsigned long int q  = (S + *P) / *Q;
    unsigned long int PP = q;
    PP *= *Q;
    PP -= *P;
    if (IS_UNLIKELY(*P == PP)) {
        *factor = *Q;
        return true;
    }
    unsigned long int t = *P;
    t   -=  PP;
    t   *=   q;
    t   += *QQ;
    *QQ  = *Q;
    *Q   =   t;
    *P   =  PP;
    return false;
}
//-----------------------------------------------------------------------------
static unsigned long int large_step(unsigned long int* const P,
                                    unsigned long int* const Q,
                                    unsigned long int* const QQ,
                                    unsigned long int PA,
                                    unsigned long int QA,
                                    unsigned long int QQA,
                                    unsigned long int PB,
                                    unsigned long int QB,
                                    unsigned long int QQB,
                                    tmp_vars_t* const vars) {
    //
    // Performs the large step algorithm, "composing" the form given by
    // (-QQA, PA, QA) with the one given by (-QQB, PB, QB) (or a nearby form)
    // and stores the resulting form in (-QQ, P, Q).
    // Notations mostly follow Williams & Wunderlich's with S = floor(sqrt(D)).
    //
    unsigned long int gcd = 0;
    unsigned long int qk  = 0;
    unsigned long int PP  = 0;

    mpz_ptr P0   = vars->tmp_1;
    mpz_ptr Q0   = vars->tmp_2;
    mpz_ptr qa_z = vars->tmp_3;
    mpz_ptr qb_z = vars->tmp_4;
    mpz_ptr M    = vars->tmp_5;
    mpz_ptr u    = vars->tmp_6;
    mpz_ptr rem  = vars->tmp_7;
    mpz_ptr quot = vars->tmp_8;
    mpz_ptr D    = vars->D;
    mpz_ptr N    = vars->N;

    unsigned long int S      = vars->S;
    unsigned long int factor = 0;

    gcd = gcd_ulint(QA, QB);
    while (gcd != 1) {
        single_step(&PB, &QB, &QQB, S);
        gcd = gcd_ulint(QA, QB);
    }
    //
    // Step 1: Compute Q0 = QA * QB -- we won't keep its inverse mod D
    //
    mpz_set_ui(Q0, QA);
    mpz_mul_ui(Q0, Q0, QB);
    //
    // We require Q0 to be invertible in Z/DZ and gcd(QA, QB) = 1,
    // so we cycle though the nearby forms until we find a suitable
    // Q pairs...
    //
    while (mpz_invert(M, Q0, D) == 0) {
        mpz_gcd(M, Q0, N);
        //
        // _NOTE_: If the number to factor is not squarefree the current loop
        //         will (most of the time, but not always!) never terminates!
        //         (the gcd M will give a factor of n or k times that factor
        //         where k is the multiplier used: this has been determined
        //         empirically, theorists --or anyone using their brain, unlike
        //         yours truly--, please correct me!). The good news is that we
        //         can escape from that loop and extract that factor...
        //
        if (mpz_cmp_ui(M, 1) != 0) {
            factor = mpz_get_ui(M);
            goto clean_and_return;
        }
        single_step(&PA, &QA, &QQA, S);
        single_step(&PB, &QB, &QQB, S);
        gcd = gcd_ulint(QA, QB);
        while (gcd != 1) {
            single_step(&PB, &QB, &QQB, S);
            gcd = gcd_ulint(QA, QB);
        }
        mpz_set_ui(Q0, QA);
        mpz_mul_ui(Q0, Q0, QB);
    }
    //
    // Solve X*QA - Y*QB = PA - PB for X. Remember than gcd(QA, QB) = 1.
    //
    mpz_set_ui(qa_z, QA);
    mpz_set_ui(qb_z, QB);
    mpz_gcdext(M, u, NULL, qa_z, qb_z);
    mpz_set_ui(P0, PB);
    mpz_sub_ui(P0, P0, PA);
    mpz_mul(P0, u, P0);
    mpz_mul_ui(P0, P0, QA);
    mpz_add_ui(P0, P0, PA);
    mpz_mod(P0, P0, Q0);
    //
    // Step 2: Find Q and P' congruent to P mod Q
    //
    mpz_add_ui(quot, P0, S);
    mpz_fdiv_qr(quot, rem, quot, Q0);
    while (true) {
        mpz_ui_sub(P0, S, rem);
        mpz_set(M, D);
        mpz_submul(M, P0, P0);
        mpz_divexact(Q0, M, Q0);
        mpz_add_ui(quot, P0, S);
        mpz_fdiv_qr(quot, rem, quot, Q0);

        if ((mpz_sgn(Q0) > 0) && (mpz_cmp_ui(Q0, S) <= 0)) {
            break;
        }
    }
    //
    // Step 3: Keep Q and P mod Q
    //
    *Q = mpz_get_ui(Q0);
    //
    // Remark: P0 may be negative and/or not fit in a unsigned long
    // int, so the mpz_mod is necessary...
    //
    mpz_mod_ui(P0, P0, *Q);
    *P = mpz_get_ui(P0);
    //
    // While not stricly needed, we compute the actual value of P
    // (and not just P mod Q) together with QQ (which is the previous
    // or the next value of Q depending on whether or not we've gone
    // beyond the period of the form cycle)
    //
    qk = (S + *P) / *Q;
    PP = qk * (*Q) - *P;
    if (S >= *P) {
        qk += (S - *P) / *Q;
    } else {
        qk -= (*P - S + *Q - 1) / *Q;
    }
    *P = qk * (*Q) - PP;

    mpz_set_ui(M, *P);
    mpz_mul_ui(M, M, *P);
    mpz_sub(M, D, M);
    mpz_fdiv_q_ui(M, M, *Q);

    *QQ = mpz_get_ui(M);

  clean_and_return:

    return factor;
}
//-----------------------------------------------------------------------------
inline static unsigned long int doubling_algo(unsigned long int* const P2,
                                              unsigned long int* const Q2,
                                              unsigned long int* const QQ2,
                                              unsigned long int P,
                                              unsigned long int Q,
                                              unsigned long int QQ,
                                              tmp_vars_t* const vars) {
    //
    // Performs the doubling algorithm, composing the form given by (P, Q, QQ)
    // with the nearest suitable form and stores the result in (P2, Q2, QQ2).
    //
    unsigned long int PB  = P;
    unsigned long int QB  = Q;
    unsigned long int QQB = QQ;
    return large_step(P2, Q2, QQ2, P, Q, QQ, PB, QB, QQB, vars);
}
//-----------------------------------------------------------------------------
inline static unsigned long int jump_over_k_steps(unsigned long int* const P,
                                                  unsigned long int* const Q,
                                                  unsigned long int* const QQ,
                                                  unsigned long int P0,
                                                  unsigned long int Q0,
                                                  unsigned long int QQ0,
                                                  unsigned long int K,
                                                  tmp_vars_t* const vars) {
    //
    // Starting with the principal form (1, P0, Q0), jump over (approximately!)
    // K forms to get the form (QQ, P, Q) using Williams & Wunderlich's large
    // step algorithm.
    //
    unsigned long int S      = vars->S;
    unsigned long int factor = 0;
    if (K <= MIN_NSTEP_LS) {
        //
        // Just perform K single steps sequentially since we're under the
        // threshold for large step use.
        //
        *P  = P0;
        *Q  = Q0;
        *QQ = QQ0;
        for (uint32_t i = 0; i < K; i++) {
            single_step(P, Q, QQ, S);
        }
    } else {
        unsigned long int PMIN  = P0;
        unsigned long int QMIN  = Q0;
        unsigned long int QQMIN = QQ0;

        for (uint32_t i = 0; i < MIN_NSTEP_LS; i++) {
            single_step(&PMIN, &QMIN, &QQMIN, S);
        }
        uint32_t n = K / MIN_NSTEP_LS;

        *P  = PMIN;
        *Q  = QMIN;
        *QQ = QQMIN;

        uint32_t f = most_significant_bit(n);

        while (f != 0) {
            f--;
            factor = doubling_algo(P, Q, QQ, *P, *Q, *QQ, vars);
            if (factor != 0) {
                return factor;
            }
            if (BIT(n, f) == 1) {
                factor = large_step(P, Q, QQ, *P, *Q, *QQ,
                                    PMIN, QMIN, QQMIN, vars);
                if (factor != 0) {
                    return factor;
                }
            }
        }
    }
    return 0;
}
//-----------------------------------------------------------------------------
static unsigned long int jump_near_sym_point(unsigned long int* const P,
                                             unsigned long int* const Q,
                                             unsigned long int* const QQ,
                                             unsigned long int P0,
                                             unsigned long int Q0,
                                             unsigned long int QQ0,
                                             unsigned long int S,
                                             mpz_t D,
                                             unsigned long int K,
                                             mpz_t N) {
    //
    // Given the principal form (1, P0, Q0) and the inverse square root of a
    // found square form given by (QQ, P, Q), jump through K forms starting
    // from (QQ, P, Q) to get in the vicinity of a point of symmetry of the
    // form cycle of (QQ, P, Q). The resulting form is stored in (QQ, P, Q).
    //
    unsigned long int PH  = 0;
    unsigned long int QH  = 0;
    unsigned long int QQH = 0;

    unsigned long int factor = 0;

    tmp_vars_t* vars = alloc_init_tmp_vars(S, D, N);

    factor = jump_over_k_steps(&PH, &QH, &QQH, P0, Q0, QQ0, K, vars);
    if (factor != 0) {
        free_tmp_vars(vars);
        return factor;
    }
    factor = large_step(P, Q, QQ, *P, *Q, *QQ, PH, QH, QQH, vars);

    free_tmp_vars(vars);
    return factor;
}
//-----------------------------------------------------------------------------
