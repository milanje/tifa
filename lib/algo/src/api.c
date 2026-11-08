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
 * \file    api.c
 * \brief   The public functions that take options: each runs its body
 *          under them.
 */

#include "tifa_config.h"
#include "tifa_factor.h"
#include "batch_gcd.h"
#include "stopping.h"
#include "impl.h"

//
// The body of a function returning an ecode_t, under the call's options.
//
#define UNDER_OPTIONS(OPTIONS, BODY) do {                                   \
    tifa_call_t call;                                                       \
    if (!tifa_call_begin(&call, (OPTIONS))) {                               \
        return FAILURE;                                                     \
    }                                                                       \
    const ecode_t ecode = (BODY);                                           \
    tifa_call_end(&call);                                                   \
    return ecode;                                                           \
} while (0)

//------------------------------------------------------------------------------
ecode_t tifa_factor(mpz_array_t* const factors, uint32_array_t* const multis,
                    const mpz_t n, const factoring_mode_t mode,
                    const tifa_options_t* const options) {
    const tifa_factor_hint_t hint = (options != NULL) ? options->hint
                                                      : TIFA_FACTOR_RANDOM;
    UNDER_OPTIONS(options, tifa_factor_impl(factors, multis, n, mode, hint));
}
//------------------------------------------------------------------------------
ecode_t tifa_factor_rec(mpz_array_t* const factors,
                        uint32_array_t* const multis,
                        const mpz_t n, factoring_mode_t mode) {
    return tifa_factor_impl(factors, multis, n, mode, TIFA_FACTOR_RANDOM);
}
//------------------------------------------------------------------------------
ecode_t cfrac(mpz_array_t* const factors, uint32_array_t* const multis,
              const mpz_t n, const cfrac_params_t* const params,
              const factoring_mode_t mode,
              const tifa_options_t* const options) {
    UNDER_OPTIONS(options, cfrac_impl(factors, multis, n, params, mode));
}
//------------------------------------------------------------------------------
ecode_t ecm(mpz_array_t* const factors, uint32_array_t* const multis,
            const mpz_t n, const ecm_params_t* const params,
            const factoring_mode_t mode, const tifa_options_t* const options) {
    UNDER_OPTIONS(options, ecm_impl(factors, multis, n, params, mode));
}
//------------------------------------------------------------------------------
ecode_t fermat(mpz_array_t* const factors, uint32_array_t* const multis,
               const mpz_t n, const fermat_params_t* const params,
               const factoring_mode_t mode,
               const tifa_options_t* const options) {
    UNDER_OPTIONS(options, fermat_impl(factors, multis, n, params, mode));
}
//------------------------------------------------------------------------------
ecode_t p_minus_1(mpz_array_t* const factors, uint32_array_t* const multis,
                  const mpz_t n, const p_minus_1_params_t* const params,
                  const factoring_mode_t mode,
                  const tifa_options_t* const options) {
    UNDER_OPTIONS(options, p_minus_1_impl(factors, multis, n, params, mode));
}
//------------------------------------------------------------------------------
ecode_t p_plus_1(mpz_array_t* const factors, uint32_array_t* const multis,
                 const mpz_t n, const p_plus_1_params_t* const params,
                 const factoring_mode_t mode,
                 const tifa_options_t* const options) {
    UNDER_OPTIONS(options, p_plus_1_impl(factors, multis, n, params, mode));
}
//------------------------------------------------------------------------------
ecode_t rho(mpz_array_t* const factors, uint32_array_t* const multis,
            const mpz_t n, const rho_params_t* const params,
            const factoring_mode_t mode, const tifa_options_t* const options) {
    UNDER_OPTIONS(options, rho_impl(factors, multis, n, params, mode));
}
//------------------------------------------------------------------------------
ecode_t siqs(mpz_array_t* const factors, uint32_array_t* const multis,
             const mpz_t n, const siqs_params_t* const params,
             const factoring_mode_t mode, const tifa_options_t* const options) {
    UNDER_OPTIONS(options, siqs_impl(factors, multis, n, params, mode));
}
//------------------------------------------------------------------------------
ecode_t squfof(mpz_array_t* const factors, uint32_array_t* const multis,
               const mpz_t n, const squfof_params_t* const params,
               const factoring_mode_t mode,
               const tifa_options_t* const options) {
    UNDER_OPTIONS(options, squfof_impl(factors, multis, n, params, mode));
}
//------------------------------------------------------------------------------
ecode_t tifa_batch_gcd(mpz_array_t* const gcds,
                       const mpz_array_t* const moduli,
                       const tifa_options_t* const options) {
    UNDER_OPTIONS(options, tifa_batch_gcd_impl(gcds, moduli));
}
//------------------------------------------------------------------------------
prime_cert_status_t tifa_prove_prime(prime_cert_t* cert, const mpz_t n,
                                     const tifa_options_t* const options) {
    tifa_call_t call;
    if (!tifa_call_begin(&call, options)) {
        cert->nsteps = 0;
        return PRIME_NOT_PROVEN;
    }
    const prime_cert_status_t status = tifa_prove_prime_impl(cert, n);
    tifa_call_end(&call);
    return status;
}
//------------------------------------------------------------------------------
cofactor_status_t tifa_cofactor(mpz_array_t* factors, uint32_array_t* multis,
                                const mpz_t n,
                                const cofactor_strategy_t* strategy,
                                const tifa_options_t* const options) {
    tifa_call_t call;
    if (!tifa_call_begin(&call, options)) {
        return COFACTOR_GAVE_UP;
    }
    const cofactor_status_t status = tifa_cofactor_impl(factors, multis, n,
                                                        strategy);
    tifa_call_end(&call);
    return status;
}
