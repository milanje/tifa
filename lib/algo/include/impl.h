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
 * \file    impl.h
 * \brief   The bodies of the public functions, without their options
 *          (internal): they run under those of the call they are made from.
 */

#if !defined(_TIFA_IMPL_H_)
#define _TIFA_IMPL_H_

#include <gmp.h>

#include "array.h"
#include "exit_codes.h"
#include "factoring_machine.h"
#include "tifa_options.h"
#include "cfrac.h"
#include "ecm.h"
#include "fermat.h"
#include "p_minus_1.h"
#include "p_plus_1.h"
#include "rho.h"
#include "siqs.h"
#include "squfof.h"
#include "cofactor.h"
#include "prime_cert.h"

ecode_t tifa_factor_impl(mpz_array_t* const factors,
                         uint32_array_t* const multis,
                         const mpz_t n, const factoring_mode_t mode,
                         const tifa_factor_hint_t hint);
//
// tifa_factor_impl without a hint, for the recursions of the factoring
// machine.
//
ecode_t tifa_factor_rec(mpz_array_t* const factors,
                        uint32_array_t* const multis,
                        const mpz_t n, factoring_mode_t mode);

ecode_t cfrac_impl(mpz_array_t* const factors, uint32_array_t* const multis,
                   const mpz_t n, const cfrac_params_t* const params,
                   const factoring_mode_t mode);
ecode_t ecm_impl(mpz_array_t* const factors, uint32_array_t* const multis,
                 const mpz_t n, const ecm_params_t* const params,
                 const factoring_mode_t mode);
ecode_t fermat_impl(mpz_array_t* const factors, uint32_array_t* const multis,
                    const mpz_t n, const fermat_params_t* const params,
                    const factoring_mode_t mode);
ecode_t p_minus_1_impl(mpz_array_t* const factors,
                       uint32_array_t* const multis, const mpz_t n,
                       const p_minus_1_params_t* const params,
                       const factoring_mode_t mode);
ecode_t p_plus_1_impl(mpz_array_t* const factors,
                      uint32_array_t* const multis, const mpz_t n,
                      const p_plus_1_params_t* const params,
                      const factoring_mode_t mode);
ecode_t rho_impl(mpz_array_t* const factors, uint32_array_t* const multis,
                 const mpz_t n, const rho_params_t* const params,
                 const factoring_mode_t mode);
ecode_t siqs_impl(mpz_array_t* const factors, uint32_array_t* const multis,
                  const mpz_t n, const siqs_params_t* const params,
                  const factoring_mode_t mode);
ecode_t squfof_impl(mpz_array_t* const factors, uint32_array_t* const multis,
                    const mpz_t n, const squfof_params_t* const params,
                    const factoring_mode_t mode);

prime_cert_status_t tifa_prove_prime_impl(prime_cert_t* cert, const mpz_t n);
//
// Drops the steps of a certificate.
//
void prime_cert_clear(prime_cert_t* const cert);
cofactor_status_t tifa_cofactor_impl(mpz_array_t* factors,
                                     uint32_array_t* multis, const mpz_t n,
                                     const cofactor_strategy_t* strategy);
ecode_t tifa_batch_gcd_impl(mpz_array_t* const gcds,
                            const mpz_array_t* const moduli);

#endif
