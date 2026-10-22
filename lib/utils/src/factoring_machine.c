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
 * \file    factoring_machine.c
 * \author  Jerome Milan
 * \date    Tue Oct 11 2011
 * \version 2011-10-11
 */

#include <stdlib.h>
#include <inttypes.h>
#include <stdbool.h>
#include <stdio.h>

#include "macros.h"
#include "funcs.h"
#include "factoring_machine.h"

//------------------------------------------------------------------------------
//                 PROTOTYPES OF NON PUBLIC FUNCTION(S)
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
ecode_t single_run(factoring_machine_t*);
//------------------------------------------------------------------------------
ecode_t find_some_factors(factoring_machine_t*);
//------------------------------------------------------------------------------
ecode_t find_some_coprime_factors(factoring_machine_t*);
//------------------------------------------------------------------------------
ecode_t find_some_prime_factors(factoring_machine_t*);
//------------------------------------------------------------------------------
ecode_t find_complete_factorization(factoring_machine_t*);
//-----------------------------------------------------------------------------
bool are_all_prime(const mpz_array_t* const);
//------------------------------------------------------------------------------
ecode_t compute_multiplicities(
    uint32_array_t* const multis,
    const mpz_t n,
    const mpz_array_t* const factors
);
//-----------------------------------------------------------------------------

//------------------------------------------------------------------------------
//                 "PUBLIC" FUNCTION(S) IMPLEMENTATION
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
ecode_t run_machine(factoring_machine_t* machine) {
    if (mpz_cmp_ui(machine->n, 1) <= 0) {
        return NO_FACTOR_FOUND;
    }
    if (MPZ_IS_SQUARE(machine->n)) {
        return FAILURE;
    }
    //
    // The modes work on arrays of their own: the caller's may be full, or
    // already hold the factors of other numbers.
    //
    mpz_array_t*    const factors = machine->factors;
    uint32_array_t* const multis  = machine->multis;

    machine->factors = alloc_mpz_array(ELONGATION);
    machine->multis  = (multis != NULL) ? alloc_uint32_array(ELONGATION)
                                        : NULL;
    ecode_t ecode = UNKNOWN_FACTORING_MODE;

    switch (machine->mode) {
        case SINGLE_RUN:
            ecode = single_run(machine);
            break;
        case FIND_SOME_FACTORS:
            ecode = find_some_factors(machine);
            break;
        case FIND_SOME_PRIME_FACTORS:
            ecode = find_some_prime_factors(machine);
            break;
        case FIND_SOME_COPRIME_FACTORS:
            ecode = find_some_coprime_factors(machine);
            break;
        case FIND_COMPLETE_FACTORIZATION:
            ecode = find_complete_factorization(machine);
            break;
        default:
            break;
    }
    append_mpz_array(factors, machine->factors);
    free_mpz_array(machine->factors);
    machine->factors = factors;

    if (multis != NULL) {
        append_uint32_array(multis, machine->multis);
        free_uint32_array(machine->multis);
        machine->multis = multis;
    }
    return ecode;
}
//-----------------------------------------------------------------------------

//------------------------------------------------------------------------------
//                  "PRIVATE" FUNCTION(S) IMPLEMENTATION
//------------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *                  Functions directly used by run_machine
 *-----------------------------------------------------------------------------
 */

//-----------------------------------------------------------------------------
ecode_t single_run(factoring_machine_t* machine) {
    //
    // Perform a single run of the factoring algorithm.
    // Don't try to do anything "smartish" if the algorithm fails to find a
    // factor...
    //
    ecode_t exit_code;

    exit_code = machine->init_context_func(machine);
    if (exit_code != SUCCESS) {
        //
        // Context initialization failed: do not run the algorithm on a
        // half-built context. Clean up and report the failure.
        //
        machine->clear_context_func(machine);
        machine->success = false;
        return exit_code;
    }
    exit_code = machine->perform_algo_func(machine);
    
    switch (exit_code) {
    case SOME_FACTORS_FOUND:
    case COMPLETE_FACTORIZATION_FOUND:
        machine->success = true;
        break;
    default:
        machine->success = false;
    }
    (void) machine->clear_context_func(machine);

    return exit_code;
}
//-----------------------------------------------------------------------------
ecode_t find_some_factors(factoring_machine_t* machine) {
    //
    // Try to find a least "some factor(s)" even if they are not prime.
    //
    ecode_t exit_code = NO_FACTOR_FOUND;

    exit_code = machine->init_context_func(machine);
    if (exit_code != SUCCESS) {
        //
        // Context initialization failed: do not run the algorithm on a
        // half-built context. Clean up and report the failure.
        //
        machine->clear_context_func(machine);
        machine->success = false;
        return exit_code;
    }
    exit_code = machine->perform_algo_func(machine);

    while (true) {
        //
        // Try again and again until we either:
        //    - find some factor(s)
        //    - reach the abort limit
        //    - or end up with a really bad error...
        //
        if (exit_code == SOME_FACTORS_FOUND) {
            break;
        }
        if (exit_code == COMPLETE_FACTORIZATION_FOUND) {
            break;
        }
        if (exit_code == FATAL_INTERNAL_ERROR) {
            break;
        }
        exit_code = machine->update_context_func(machine);

        if (exit_code != SUCCESS) {
            break;
        }
        exit_code = machine->perform_algo_func(machine);
    }
    machine->clear_context_func(machine);

    switch (exit_code) {
    case SOME_FACTORS_FOUND:
    case COMPLETE_FACTORIZATION_FOUND:
        machine->success = true;
        break;
    
    default:
        machine->success = false;
    }
    return exit_code;
}
//-----------------------------------------------------------------------------
ecode_t find_some_coprime_factors(factoring_machine_t* machine) {
    //
    // Try to find "some coprime factors" but don't try to find _all_ prime
    // factors.
    //
    // Compared to find_some_factors(), this involves some processing on the
    // found factors to obtain a coprime base. Actually there is no guarantee
    // that the found factors are prime.
    //
    ecode_t exit_code = NO_FACTOR_FOUND;

    exit_code = machine->init_context_func(machine);
    if (exit_code != SUCCESS) {
        //
        // Context initialization failed: do not run the algorithm on a
        // half-built context. Clean up and report the failure.
        //
        machine->clear_context_func(machine);
        machine->success = false;
        return exit_code;
    }
    exit_code = machine->perform_algo_func(machine);

    while (true) {
        //
        // Try again and again until we either:
        //    - find some factor(s)
        //    - reach the abort limit
        //    - or end up with a really bad error...
        //
        if (exit_code == SOME_FACTORS_FOUND) {
            break;
        }
        if (exit_code == COMPLETE_FACTORIZATION_FOUND) {
            break;
        }
        if (exit_code == FATAL_INTERNAL_ERROR) {
            break;
        }
        exit_code = machine->update_context_func(machine);

        if (exit_code != SUCCESS) {
            break;
        }
        exit_code = machine->perform_algo_func(machine);
    }
    machine->clear_context_func(machine);

    switch (exit_code) {
    case SOME_FACTORS_FOUND:
        exit_code = process_found_factors(
                        machine->factors,
                        machine->multis,
                        machine->n,
                        FIND_SOME_COPRIME_FACTORS,
                        machine->recurse_func
                    );
        machine->success = true;
        break;

    case COMPLETE_FACTORIZATION_FOUND:
        machine->success = true;
        break;
    
    default:
        machine->success = false;
    }
    return exit_code;
}
//-----------------------------------------------------------------------------
ecode_t find_some_prime_factors(factoring_machine_t* machine) {
    //
    // Try to find "some prime factors" but don't try to find _all_ prime
    // factors.
    //
    // Compared to find_some_factors(), this involves some processing on the
    // found factors (to obtain a coprime base) and some compositeness checks.
    // Note that we return FAILURE as soon as one of the found factor is NOT
    // prime. SUCCESS is therefore restricted to the case where all found
    // factors (after processing) are prime. It is debatable whether or not
    // this is a suitable behaviour. Indeed we could return SUCCESS as soon as
    // _one_ of the found factors is prime, which would actually better match
    // the name of the function. But then, what would be the point of such a
    // factoring mode? Heck, I'm "not even sure" (to put it mildly) that it
    // is useful at all the way it is...
    //
    ecode_t exit_code = NO_FACTOR_FOUND;

    exit_code = machine->init_context_func(machine);
    if (exit_code != SUCCESS) {
        //
        // Context initialization failed: do not run the algorithm on a
        // half-built context. Clean up and report the failure.
        //
        machine->clear_context_func(machine);
        machine->success = false;
        return exit_code;
    }
    exit_code = machine->perform_algo_func(machine);

    while (true) {
        //
        // Try again and again until we either:
        //    - find some factor(s)
        //    - reach the abort limit
        //    - or end up with a really bad error...
        //
        if (exit_code == SOME_FACTORS_FOUND) {
            break;
        }
        if (exit_code == COMPLETE_FACTORIZATION_FOUND) {
            break;
        }
        if (exit_code == FATAL_INTERNAL_ERROR) {
            break;
        }
        exit_code = machine->update_context_func(machine);

        if (exit_code != SUCCESS) {
            break;
        }
        exit_code = machine->perform_algo_func(machine);
    }
    machine->clear_context_func(machine);


    switch (exit_code) {
    case SOME_FACTORS_FOUND:
        //
        // _NOTE_: We fail as soon as one element in the coprime base is
        //         not prime. Is that a sensible behavior?
        //
        exit_code = process_found_factors(
                        machine->factors,
                        machine->multis,
                        machine->n,
                        FIND_SOME_PRIME_FACTORS,
                        machine->recurse_func
                    );
        machine->success = (exit_code == SOME_PRIME_FACTORS_FOUND);
        break;

    case COMPLETE_FACTORIZATION_FOUND:
        machine->success = true;
        break;
    
    default:
        machine->success = false;
    }
    return exit_code;
}
//-----------------------------------------------------------------------------
ecode_t find_complete_factorization(factoring_machine_t* machine) {
    //
    // Try to find a "complete factorization", that is try to find _all_ prime
    // factors together with their multiplicities.
    //
    ecode_t exit_code = NO_FACTOR_FOUND;
    machine->success  = false;
    
    exit_code = machine->init_context_func(machine);
    if (exit_code != SUCCESS) {
        //
        // Context initialization failed: do not run the algorithm on a
        // half-built context. Clean up and report the failure.
        //
        machine->clear_context_func(machine);
        machine->success = false;
        return exit_code;
    }
    exit_code = machine->perform_algo_func(machine);
    
    while (true) {
        if (exit_code == SOME_FACTORS_FOUND) {
            break;
        }
        if (exit_code == COMPLETE_FACTORIZATION_FOUND) {
            break;
        }
        if (exit_code == FATAL_INTERNAL_ERROR) {
            break;
        }
        exit_code = machine->update_context_func(machine);

        if (exit_code != SUCCESS) {
            break;
        }
        exit_code = machine->perform_algo_func(machine);
    }
    machine->clear_context_func(machine);

    if (exit_code == COMPLETE_FACTORIZATION_FOUND) {
        //
        // The algorithm may report the complete factorization without the
        // multiplicities (e.g. the semiprime early exit of SIQS).
        //
        exit_code = compute_multiplicities(
                        machine->multis,
                        machine->n,
                        machine->factors
                    );
        if (exit_code != SUCCESS) {
            return exit_code;
        }
        machine->success = true;
        return COMPLETE_FACTORIZATION_FOUND;
    }

    if (exit_code != SOME_FACTORS_FOUND) {
        return exit_code;
    }
    exit_code = process_found_factors(
                    machine->factors,
                    machine->multis,
                    machine->n,
                    FIND_COMPLETE_FACTORIZATION,
                    machine->recurse_func
                );
    machine->success = (exit_code == COMPLETE_FACTORIZATION_FOUND);

    return exit_code;
}
//-----------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *                 Post-processing of the factors found by a run
 *-----------------------------------------------------------------------------
 */

//-----------------------------------------------------------------------------
static void replace_by_coprime_base(mpz_array_t* const factors,
                                    const mpz_t n) {
    uint32_t     flen = factors->length;
    mpz_array_t* base = alloc_mpz_array(flen * (flen + 1));

    find_coprime_base(base, n, factors);
    swap_mpz_array(factors, base);

    free_mpz_array(base);
}
//-----------------------------------------------------------------------------
static ecode_t complete_factorization(
    mpz_array_t* const factors,
    uint32_array_t* const multis,
    const mpz_t n,
    ecode_t (*recurse)(mpz_array_t* const, uint32_array_t* const,
                       const mpz_t, factoring_mode_t)) {

    //
    // Compute a coprime base for all the found factors, then compute the
    // multiplicity of each element in said base (even if they are not all
    // prime).
    //
    replace_by_coprime_base(factors, n);

    ecode_t exit_code = compute_multiplicities(multis, n, factors);

    if (exit_code != SUCCESS) {
        return exit_code;
    }
    
    //
    // We have a coprime base together with the multiplicity of each of its
    // element. Factor those which are not prime (each is tested once).
    //
    uint32_t alloced = 2 * factors->length;
    
    //
    // The all_* arrays will hold all the (hopefully prime) factors of
    // n and their multiplicities.
    //
    mpz_array_t*    all_factors = alloc_mpz_array(alloced);
    uint32_array_t* all_multis  = alloc_uint32_array(alloced);

    bool factorization_is_partial = false;

    for (uint32_t i = 0; i < factors->length; i++) {
        
        const mpz_ptr curfactor = factors->data[i];
        const uint32_t curmulti = multis->data[i];
        
        if (MPZ_IS_PRIME(curfactor)) {
            
            append_mpz_to_array(all_factors, curfactor);
            append_uint32_to_array(all_multis, curmulti);

        } else {
            //
            // Proceed to factorize this non prime factor...
            //
            mpz_array_t*   morefactors = alloc_mpz_array(ELONGATION);
            uint32_array_t* moremultis = alloc_uint32_array(ELONGATION);

            ecode_t retcode = recurse(
                                  morefactors,
                                  moremultis,
                                  curfactor,
                                  FIND_COMPLETE_FACTORIZATION
                              );

            switch (retcode) {

            case COMPLETE_FACTORIZATION_FOUND:
                for (uint32_t j = 0; j < moremultis->length; j++) {
                    moremultis->data[j] *= curmulti;
                }
                append_mpz_array(all_factors, morefactors);
                append_uint32_array(all_multis, moremultis);
                break;

            case PARTIAL_FACTORIZATION_FOUND:
                factorization_is_partial = true;
                
                for (uint32_t j = 0; j < moremultis->length; j++) {
                    moremultis->data[j] *= multis->data[i];
                }
                append_mpz_array(all_factors, morefactors);
                append_uint32_array(all_multis, moremultis);
                break;

            case FATAL_INTERNAL_ERROR:
            case NO_FACTOR_FOUND:
            case GIVING_UP:
            default:
                //
                // For the time being we just return with the found
                // factors and give up. Ideally, we could try different
                // factoring algorithms here before resigning.
                //
                factorization_is_partial = true;
                append_mpz_to_array(all_factors, curfactor);
                append_uint32_to_array(all_multis, curmulti);
                break;
            }
            free_mpz_array(morefactors);
            free_uint32_array(moremultis);
        }
    }
    swap_mpz_array(factors, all_factors);
    swap_uint32_array(multis, all_multis);

    free_mpz_array(all_factors);
    free_uint32_array(all_multis);

    if (factorization_is_partial) {
        return PARTIAL_FACTORIZATION_FOUND;
    } else {
        return COMPLETE_FACTORIZATION_FOUND;
    }
}
//-----------------------------------------------------------------------------
ecode_t process_found_factors(mpz_array_t* const factors,
                              uint32_array_t* const multis,
                              const mpz_t n,
                              const factoring_mode_t mode,
                              ecode_t (*recurse)(mpz_array_t* const,
                                                 uint32_array_t* const,
                                                 const mpz_t,
                                                 factoring_mode_t)) {
    switch (mode) {
    case SINGLE_RUN:
    case FIND_SOME_FACTORS:
        return SOME_FACTORS_FOUND;

    case FIND_SOME_COPRIME_FACTORS:
        replace_by_coprime_base(factors, n);
        return SOME_COPRIME_FACTORS_FOUND;

    case FIND_SOME_PRIME_FACTORS:
        replace_by_coprime_base(factors, n);
        if (are_all_prime(factors)) {
            return SOME_PRIME_FACTORS_FOUND;
        }
        return SOME_COPRIME_FACTORS_FOUND;

    case FIND_COMPLETE_FACTORIZATION:
        return complete_factorization(factors, multis, n, recurse);

    default:
        return UNKNOWN_FACTORING_MODE;
    }
}
//-----------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *                       Various utility functions
 *-----------------------------------------------------------------------------
 */

//-----------------------------------------------------------------------------
bool are_all_prime(const mpz_array_t* const array) {
    //
    // Returns true if all elements in array are (probably) prime.
    // Returns false otherwise.
    //
    bool all_prime = true;
    for (uint32_t i = 0; i < array->length; i++) {
        if (!MPZ_IS_PRIME(array->data[i])) {
            all_prime = false;
            break;
        }
    }
    return all_prime;
}
//-----------------------------------------------------------------------------
ecode_t compute_multiplicities(uint32_array_t* const multis,
                               const mpz_t n,
                               const mpz_array_t* const factors) {
    //
    // Computes the multiplicities of the factors of 'n' given by the coprime
    // base 'factors' and stores them in the 'multis' array.
    //
    mpz_t cofactor;
    mpz_init_set(cofactor, n);

    multis->length = 0;

    for (uint32_t i = 0; i < factors->length; i++) {
        //
        // _NOTE_: We already know that 'n' is divisible by each element in
        //         the array 'factors' _and_ we known that these factors are
        //         coprime with each other (i.e. 'factors' is a coprime base)
        //         so no divisibility test is required for the first division.
        //
        mpz_divexact(cofactor, cofactor, factors->data[i]);
        append_uint32_to_array(multis, 1);
    }
    if (mpz_cmp_ui(cofactor, 1) != 0) {
        for (uint32_t i = 0; i < factors->length; i++) {
            while (mpz_divisible_p(cofactor, factors->data[i])) {
                mpz_divexact(cofactor, cofactor, factors->data[i]);
                multis->data[i]++;
            }
            if (mpz_cmp_ui(cofactor, 1) == 0) {
                break;
            }
        }
    }
    if (mpz_cmp_ui(cofactor, 1) != 0) {
        //
        // This should not happen since 'n' is supposed to be
        // smooth on the coprime base given by 'factors'! We could return
        // the found factors and multiplities but that would mask a potential 
        // nasty bug in the code, so let's be convervative and return an error.
        //
        mpz_clear(cofactor);
        return FATAL_INTERNAL_ERROR;
    }
    mpz_clear(cofactor);
    return SUCCESS;
}
//-----------------------------------------------------------------------------
