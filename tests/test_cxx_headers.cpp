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

// The public headers from C++: they declare C functions, and their includes
// (<gmp.h> among them) stay outside extern "C".

#include "tifa.h"
#include "bitstring_t.h"
#include "factor_u64.h"
#include "methods.h"
#include "tuning.h"

#include <cstdio>

int main() {
    mpz_t n;
    mpz_init_set_str(n, "1000000016000000063", 10);
    mpz_array_t* factors   = alloc_mpz_array(4);
    uint32_array_t* multis = alloc_uint32_array(4);
    tifa_options_t options = {};
    options.threads = 2;
    const ecode_t e = tifa_factor(factors, multis, n,
                                  FIND_COMPLETE_FACTORIZATION, &options);
    const bool ok = (e == COMPLETE_FACTORIZATION_FOUND)
                    && (factors->length == 2) && (first_primes[0] == 2);
    std::printf("== C++ headers ==\n  %s\n", ok ? "ok" : "FAIL");
    free_mpz_array(factors);
    free_uint32_array(multis);
    mpz_clear(n);
    return ok ? 0 : 1;
}
