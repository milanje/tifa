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
 * \file    messages.h
 * \author  Jerome Milan
 * \date    Tue Oct 11 2011
 * \version 2011-10-11
 *
 * \brief Status / error messages output macros
 *
 * This file defines some macros used to output status or error messages if
 * some of the \c *_VERBOSE and/or \c *_TIMING symbols are set to non-zero.
 *
 * \warning The <tt>__VERBOSE__</tt>, <tt>__TIMING__</tt> and
 * <tt>__PREFIX__</tt> symbols should be defined in the file including this
 * header. It is under the responsability of the including file to check for
 * these symbol definitions or to define them, if needed.
 */

#if !defined(_TIFA_MESSAGES_H_)
   /**
    * \def _TIFA_MESSAGES_H_
    * Standard include guard.
    */
#define _TIFA_MESSAGES_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "timer.h"
#include "tifa_config.h"

//
// __VERBOSE__, __TIMING__ and __PREFIX__ should be defined in the file
// including this header
//
#if !defined(__VERBOSE__)
#error "__VERBOSE__ must be defined before including this header."
#endif
#if !defined(__TIMING__)
#error "__TIMING__ must be defined before including this header."
#endif
#if !defined(__PREFIX__)
#error "__PREFIX__ must be defined before including this header."
#endif

#if __VERBOSE__ || __TIMING__
    #include <stdio.h>
    #include "exit_codes.h"
#endif

#define __EXTRA_TIMING__ TIFA_EXTRA_TIMING
#if !defined(__EXTRA_MSG__)
	#define __EXTRA_MSG__ 0
#endif


//
// Decoration strings
//
#define TAB    "    "
#define SEPARATION_LINE "-----------------------------------" \
                        "------------------------------------"

//
// Strings pertaining to all factorization algorithms
//
#define INIT_STRING             "computing start up data..."
#define UPDATE_STRING           "updating context... "
#define UPDATE_GIVEUP_STRING    "updating context... give up"
#define CLEAN_STRING            "cleaning..."
#define FUNC_SUCCEEDED_STRING   "completed"
#define FUNC_FAILED_STRING      "failed"
//
// Strings pertaining to "congruences of squares" methods only
//
#define COLLECT_RELS_STRING      "collecting (X^2 = Y mod N) relations..."
#define COLLECT_RELS_DONE_STRING "collection of (X^2 = Y mod N) relations"
#define FACTOR_RES_STRING        "factoring residues on the factor base..."
#define LIN_ALG_STRING           "solving linear algebra system..."
#define DED_FACTORS_STRING       "deducing factors..."
#define UPDATE_MORE_RELS_STRING  "updating context to find more relations... "
#define UPDATE_NEW_MULT_STRING   "updating context to change multiplier... "
//
// Strings pertaining to Fermat's factorization only
//
#define FERMAT_FACT_STRING       "performing Fermat's factorization..."
#define FERMAT_FACT_DONE_STRING  "Fermat's factorization..."
#define UPDATE_MULTIPREC_STRING  "updating context to use multi-precision..."
//
// Strings pertaining to ECM factorization only
//
#define PERFORMING_ECM_STRING      "performing elliptic curve factorization..."
#define PERFORMING_P2_PRECOMP_STRING "performing phase 2 precomputations..."
#define CHOOSED_CURVES_IN_STRING     "choosing curves..."
#define PHASE_1_IN_STRING            "phase 1 (accumulated)..."
#define PHASE_2_IN_STRING            "phase 2 (accumulated)..."
#define UPDATE_MORE_CURVES_STRING    "updating context to try more curves... "
#define UPDATE_LARGER_BOUNDS_STRING  "updating context to use larger bounds... "
//
// Strings pertaining to SQUFOF only
//
#define FWD_CYCL_STRING    "forward cycling to find a proper form..."
#define INV_SQRT_STRING    "computing inverse square root of form..."
#define REV_CYCL_STRING    "reverse cycling to find a factor..."
#define FR_REV_CYCL_STRING "reverse cycling to find a factor (fast ret)..."
#define UPDATE_RACE_STRING "updating context to perform a race... "
#define ALL_RACERS_FAILED_STRING    "all racers failed"
#define QUEUE_OVERFLOW_STRING       "queue overflow"
#define NO_PROPER_FORM_FOUND_STRING "no proper form found"
#define UPPER_BOUND_REACHED_STRING  "reached upper bound on number of forms to try"
//
// Strings pertaining to SIQS only
//
#define NEXT_A_HACK_STRING "WARNING: running out of polynomials soon..."
//
// Strings pertaining to trial division only
//
#define DIVIDING_STRING_FORMAT "trial dividing by %lu primes..."
//
// Strings pertaining to Pollard's rho only
//
#define FLOYD_CYCLE_STRING "performing Floyd's cycle detection..."
#define BRENT_CYCLE_STRING "performing Brent's cycle detection..."
#define UPDATE_RHO_CONTEXT_STRING "updating iteration function..."

//
// No messages are printed if __VERBOSE__ and __TIMING__ are set to 0.
//
#if __VERBOSE__ || __TIMING__
    #define PRINTF(...)      printf(__VA_ARGS__); fflush(stdout)
    #define PRFX_PRINTF(...) printf(__PREFIX__ __VA_ARGS__); fflush(stdout)
#else
    #define PRINTF(...)          /* intentionally left empty */
    #define PRFX_PRINTF(...)     /* intentionally left empty */
#endif

//
// If needed, align the right margin by padding with white spaces
//
#define PRINT_LINE(STRING)    PRFX_PRINTF("%s\n",  STRING)
#define PADDED_PRINTF(STRING) PRFX_PRINTF("%-47s", STRING)

//
// No newline after message if __TIMING__ set to 1 to be able to write the
// timing results on the same line
//
#if !__TIMING__
    #define PRINT_MSG(STRING) PRINT_LINE(STRING)
#else
    #define PRINT_MSG(STRING) PADDED_PRINTF(STRING)
#endif

#if __VERBOSE__ || __TIMING__
    #define PRINT_SEPARATION_LINE PRINT_LINE(SEPARATION_LINE)
#else
    #define PRINT_SEPARATION_LINE /* intentionally left empty */
#endif

//
// Messages pertaining to all factorization algorithms
//
#define PRINT_FAILURE_NL           PRINTF("failure!\n")
#define PRINT_INIT_MSG             PRINT_MSG(INIT_STRING)
#define PRINT_UPDATE_MSG           PRINT_MSG(UPDATE_STRING)
#define PRINT_UPDATE_GIVEUP_MSG    PRINT_MSG(UPDATE_GIVEUP_STRING)
#define PRINT_CLEAN_MSG            PRINT_MSG(CLEAN_STRING)

#define PRINT_SUCCESS_MSG(ECODE) do {                         \
    char __tmp_str__[64];                                     \
    (void)sprintf(__tmp_str__, FUNC_SUCCEEDED_STRING " (%s)", \
                  ecode_to_str[ECODE]);                       \
    PRINT_MSG(__tmp_str__);                                   \
} while (0)

#define PRINT_FAILURE_MSG(ECODE) do {                         \
    char __tmp_str__[64];                                     \
    (void)sprintf(__tmp_str__, FUNC_FAILED_STRING " (%s)",    \
                  ecode_to_str[ECODE]);                       \
    PRINT_MSG(__tmp_str__);                                   \
} while (0)

#if __VERBOSE__ || __TIMING__
    //
    // Introductive messages for all algorithms
    //
    #define PRINT_INTRO_MSG(NAME)                                     \
        do {                                                          \
            PRINTF(__PREFIX__ "attempting %s factorization\n", NAME); \
            PRINT_SEPARATION_LINE;                                    \
        } while (0)
    //
    // Status message giving returned exit code
    //
    #define PRINT_STATUS(SUCCESS, ECODE)      \
        do {                                  \
            PRINT_SEPARATION_LINE;            \
            if (SUCCESS) {                    \
                PRINT_SUCCESS_MSG(ECODE);     \
            } else {                          \
                PRINT_FAILURE_MSG(ECODE);     \
            }                                 \
            PRINT_TOTAL;                      \
        } while (0)
#else
    #define PRINT_STATUS(SUCCESS, ECODE_STR) /* intentionally left empty */
    #define PRINT_INTRO_MSG(NAME)            /* intentionally left empty */
#endif
//
// Messages pertaining to "congruences of squares" methods only
//
#define PRINT_UPDATE_MORE_RELS_MSG PRINT_MSG(UPDATE_MORE_RELS_STRING)
#define PRINT_UPDATE_NEW_MULT_MSG  PRINT_MSG(UPDATE_NEW_MULT_STRING)

#define PRINT_FACTOR_RES_MSG     PRINT_MSG(FACTOR_RES_STRING)
#define PRINT_LIN_ALG_MSG        PRINT_MSG(LIN_ALG_STRING)
#define PRINT_DED_FACTORS_MSG    PRINT_MSG(DED_FACTORS_STRING)

#if __VERBOSE__
    #define PRINT_COLLECT_RELS_MSG PRINT_LINE(COLLECT_RELS_STRING)
#else
    #define PRINT_COLLECT_RELS_MSG PADDED_PRINTF(COLLECT_RELS_STRING)
#endif

#if (__VERBOSE__ && __TIMING__) || (__TIMING__ && __EXTRA_TIMING__)
    #define PRINT_COLLECT_RELS_DONE_MSG PRINT_MSG(COLLECT_RELS_DONE_STRING)
#else
    #define PRINT_COLLECT_RELS_DONE_MSG /* intentionally left empty */
#endif

#if __VERBOSE__
    #define PRINT_NRELS_FOUND(FOUND, TO_FIND)                          \
        PRFX_PRINTF(TAB TAB "(%"PRIu32"/%"PRIu32") relations found\n", \
                    FOUND, TO_FIND);
#else
    #define PRINT_NRELS_FOUND(FOUND, TO_FIND) /* intentionally left empty */
#endif

#if __VERBOSE__ && __TIMING__
    #define PRINT_RES_GENERATED_MSG(TIMING)                                   \
        PRFX_PRINTF(TAB "(relations generated in "TIMING_FORMAT" seconds)\n", \
                    TIMING);
    #define PRINT_RES_SELECTED_MSG(TIMING)                                    \
        PRFX_PRINTF(TAB "(relations selected in  "TIMING_FORMAT" seconds)\n", \
                    TIMING);
#else
    #define PRINT_RES_GENERATED_MSG(TIMING) /* intentionally left empty */
    #define PRINT_RES_SELECTED_MSG(TIMING)  /* intentionally left empty */
#endif

//
// Messages pertaining to Fermat's algorithm only
//
#define PRINT_UPDATE_MULTIPREC_MSG  PRINT_MSG(UPDATE_MULTIPREC_STRING)

#if __VERBOSE__
    #define PRINT_FERMAT_FACT_MSG PRINT_LINE(FERMAT_FACT_STRING)
#else
    #define PRINT_FERMAT_FACT_MSG PADDED_PRINTF(FERMAT_FACT_STRING)
#endif

#if __VERBOSE__ && __TIMING__
    #define PRINT_FERMAT_FACT_DONE_MSG PRINT_MSG(FERMAT_FACT_DONE_STRING)
    #define PRINT_NEXTPRIME_MSG(TIMING)                                     \
        PRFX_PRINTF(TAB "primes computed in   "TIMING_FORMAT" seconds\n",   \
                    TIMING);
    #define PRINT_SQRTM_MSG(TIMING)                                         \
        PRFX_PRINTF(TAB "sqrtm  computed in   "TIMING_FORMAT" seconds\n",   \
                    TIMING);
    #define PRINT_GREEDY_MSG(TIMING)                                        \
        PRFX_PRINTF(TAB "greedy phase done in "TIMING_FORMAT" seconds\n",   \
                    TIMING);
#else
    #define PRINT_NEXTPRIME_MSG(TIMING) /* intentionally left empty */
    #define PRINT_SQRTM_MSG(TIMING)     /* intentionally left empty */
    #define PRINT_GREEDY_MSG(TIMING)    /* intentionally left empty */
    #define PRINT_FERMAT_FACT_DONE_MSG  /* intentionally left empty */
#endif

//
// Messages pertaining to ECM only
//
#define PRINT_PERFORMING_ECM            PRINT_MSG(PERFORMING_ECM_STRING)
#define PRINT_PERFORMING_P2_PRECOMP     PRINT_MSG(PERFORMING_P2_PRECOMP_STRING)
#if __TIMING__
    #define PRINT_CHOOSED_CURVES_IN     PRINT_MSG(CHOOSED_CURVES_IN_STRING)
    #define PRINT_PHASE_1_IN            PRINT_MSG(PHASE_1_IN_STRING)
    #define PRINT_PHASE_2_IN            PRINT_MSG(PHASE_2_IN_STRING)
#else
    #define PRINT_CHOOSED_CURVES_IN     /* intentionally left empty */
    #define PRINT_PHASE_1_IN            /* intentionally left empty */
    #define PRINT_PHASE_2_IN            /* intentionally left empty */
#endif

#define PRINT_NCURVES_USED(NCURVES) \
            PRFX_PRINTF("used %"PRIu32" curve(s)\n", (uint32_t)NCURVES)

#define PRINT_UPDATE_MORE_CURVES_MSG   PRINT_MSG(UPDATE_MORE_CURVES_STRING)
#define PRINT_UPDATE_LARGER_BOUNDS_MSG PRINT_MSG(UPDATE_LARGER_BOUNDS_STRING)
//
// Messages pertaining to SQUFOF only
//
#define PRINT_UPDATE_RACE_MSG   PRINT_MSG(UPDATE_RACE_STRING)
#define PRINT_FWD_CYCL_MSG      PRINT_MSG(FWD_CYCL_STRING)
#define PRINT_INV_SQRT_MSG      PRINT_MSG(INV_SQRT_STRING)
#define PRINT_REV_CYCL_MSG      PRINT_MSG(REV_CYCL_STRING)
#define PRINT_FR_REV_CYCL_MSG   PRINT_MSG(FR_REV_CYCL_STRING)
#define PRINT_FR_REV_CYCL_MSG   PRINT_MSG(FR_REV_CYCL_STRING)
#if __EXTRA_MSG__
    #define PRINT_FWD_CYCL_EXTRA_MSG(MSG) PRINTF("\n");PRINT_LINE(TAB MSG);PRINT_FWD_CYCL_MSG
    #define PRINT_FWD_CYCL_EXTRA_MSG_MULT(MSG, MULT) PRINTF("\n");PRFX_PRINTF(TAB "mult. %lu: %s", MULT, MSG);//PRINT_FWD_CYCL_MSG
    #define PRINT_ALL_RACERS_FAILED       PRINT_FWD_CYCL_EXTRA_MSG(ALL_RACERS_FAILED_STRING)
    #define PRINT_QUEUE_OVERFLOW          PRINT_FWD_CYCL_EXTRA_MSG(QUEUE_OVERFLOW_STRING)
    #define PRINT_NO_PROPER_FORM_FOUND    PRINT_FWD_CYCL_EXTRA_MSG(NO_PROPER_FORM_FOUND_STRING)
    #define PRINT_UPPER_BOUND_REACHED     PRINT_FWD_CYCL_EXTRA_MSG(UPPER_BOUND_REACHED_STRING)
    #define PRINT_QUEUE_OVERFLOW_MULT(M)       PRINT_FWD_CYCL_EXTRA_MSG_MULT(QUEUE_OVERFLOW_STRING, M)
    #define PRINT_NO_PROPER_FORM_FOUND_MULT(M) PRINT_FWD_CYCL_EXTRA_MSG_MULT(NO_PROPER_FORM_FOUND_STRING, M)
    #define PRINT_UPPER_BOUND_REACHED_MULT(M)  PRINT_FWD_CYCL_EXTRA_MSG_MULT(UPPER_BOUND_REACHED_STRING, M)
    #define PRINT_SUCCESS_MULT(M)              PRINT_FWD_CYCL_EXTRA_MSG_MULT("success\n", M);PRINT_FWD_CYCL_MSG
#else
    #define PRINT_EXTRA_MSG(MSG)               /* intentionally left empty */
    #define PRINT_ALL_RACERS_FAILED            /* intentionally left empty */
    #define PRINT_QUEUE_OVERFLOW               /* intentionally left empty */
    #define PRINT_NO_PROPER_FORM_FOUND         /* intentionally left empty */
    #define PRINT_UPPER_BOUND_REACHED          /* intentionally left empty */
    #define PRINT_QUEUE_OVERFLOW_MULT(M)       /* intentionally left empty */
    #define PRINT_NO_PROPER_FORM_FOUND_MULT(M) /* intentionally left empty */
    #define PRINT_UPPER_BOUND_REACHED_MULT(M)  /* intentionally left empty */
    #define PRINT_SUCCESS_MULT(M)              /* intentionally left empty */
#endif
//
// Messages pertaining to SIQS only
//
#if __EXTRA_TIMING__
    #if __VERBOSE__
        #define PRINT_FILL_TIMING_MSG                                        \
            printf(__PREFIX__ "   Fill sieve:     done in %7.4f seconds\n",  \
                   GET_EXTRA_TIMING(fill));
    #else
        #define PRINT_FILL_TIMING_MSG                             \
            printf("\n" __PREFIX__                                \
                   "   Fill sieve:     done in %7.4f seconds\n",  \
                   GET_EXTRA_TIMING(fill));
    #endif

    #define PRINT_SCAN_TIMING_MSG                                        \
        printf(__PREFIX__ "   Scan sieve:     done in %7.4f seconds\n",  \
               GET_EXTRA_TIMING(scan));

    #define PRINT_INIT_POLY_TIMING_MSG                                   \
        printf(__PREFIX__ "    Init poly:     done in %7.4f seconds\n",  \
               GET_EXTRA_TIMING(init_poly));
#else
    #define PRINT_FILL_TIMING_MSG       /* intentionally left empty */
    #define PRINT_SCAN_TIMING_MSG       /* intentionally left empty */
    #define PRINT_INIT_POLY_TIMING_MSG  /* intentionally left empty */
    #define PRINT_BATCH_TIMING_MSG      /* intentionally left empty */
#endif

//
// Messages pertaining to trial division only
//
#if __TIMING__ || __VERBOSE__
    #define PRINT_DIVIDING_MSG(X)                              \
        do {                                                   \
            char __tmp_str__[64];                              \
            (void)sprintf(__tmp_str__, DIVIDING_STRING_FORMAT, \
                          (unsigned long int)X);               \
            PRINT_MSG(__tmp_str__);                            \
        } while (0)
#else
    #define PRINT_DIVIDING_MSG(X) /* intentionally left empty */
#endif

//
// Messages pertaining to Pollard's rho only
//
#if __TIMING__ || __VERBOSE__
    #define PRINT_PERFORMING_FLOYD_CYCLE PRINT_MSG(FLOYD_CYCLE_STRING)
    #define PRINT_PERFORMING_BRENT_CYCLE PRINT_MSG(BRENT_CYCLE_STRING)
    #define PRINT_UPDATING_RHO_CONTEXT   PRINT_MSG(UPDATE_RHO_CONTEXT_STRING)
#else
    #define PRINT_PERFORMING_FLOYD_CYCLE /* intentionally left empty */
    #define PRINT_PERFORMING_BRENT_CYCLE /* intentionally left empty */
    #define PRINT_UPDATING_RHO_CONTEXT   /* intentionally left empty */
#endif

//
// Messages pertaining to timings
//
#if __TIMING__
    #define PRINT_NAMED_TIMING(NAME)                      \
            PRINTF("done in " TIMING_FORMAT " seconds\n", \
                   GET_NAMED_TIMING(NAME));

    #define PRINT_NAMED_TIMING_ON_FAILURE(NAME)           \
            PRINTF("failed  " TIMING_FORMAT " seconds\n", \
                   GET_NAMED_TIMING(NAME));

    #define PRINT_NAMED_TOTAL(NAME)                       \
            PRINTF("total   " TIMING_FORMAT " seconds\n", \
                   GET_NAMED_TIMING(NAME));
#else
    #define PRINT_NAMED_TIMING(NAME)            /* intentionally left empty */
    #define PRINT_NAMED_TIMING_ON_FAILURE(NAME) /* intentionally left empty */
	#define PRINT_NAMED_TOTAL(NAME)             /* intentionally left empty */
#endif
#define PRINT_TIMING            PRINT_NAMED_TIMING()
#define PRINT_TIMING_ON_FAILURE PRINT_NAMED_TIMING_ON_FAILURE()
#define PRINT_TOTAL             PRINT_NAMED_TOTAL()

#ifdef __cplusplus
}
#endif

#endif
