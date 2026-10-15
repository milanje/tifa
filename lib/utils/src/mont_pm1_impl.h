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

//
// p-1 and p+1 kernels for moduli of exactly NL 64-bit limbs. Included once
// per limb count by mont_ecm.c after mont_ecm_impl.h (whose arithmetic they
// use), with NL and FN(name) defined.
//

//
// Bit i of the exponent E stored in plan->e (least significant limb first).
//
#define E_BIT(plan, i) (((plan)->e[(i) / 64] >> ((i) % 64)) & 1)

//
// x <- base^E. With base 2 the multiplications by the base are additions.
//
static void FN(pm1_stage1)(uint64_t* const x, const uint64_t* const base,
                           const bool base_is_two,
                           const mont_ecm_plan_t* const plan,
                           const mod_t* const m) {
    for (int j = 0; j < NL; j++) {
        x[j] = base[j];
    }
    for (long i = (long) plan->ebits - 2; i >= 0; i--) {
        FN(mul)(x, x, x, m);
        if (E_BIT(plan, i)) {
            if (base_is_two) {
                FN(add)(x, x, x, m);
            } else {
                FN(mul)(x, x, base, m);
            }
        }
    }
}

//
// (v0, v1) <- (V_k(p), V_{k+1}(p)) for k >= 1, where V is the Lucas
// sequence V_0 = 2, V_1 = p, V_{i+1} = p.V_i - V_{i-1}.
//
static void FN(lucas_ladder)(uint64_t* const v0, uint64_t* const v1,
                             const uint64_t* const p, const uint64_t k,
                             const mod_t* const m) {
    uint64_t two[NL];
    FN(add)(two, m->one, m->one, m);

    for (int j = 0; j < NL; j++) {
        v0[j] = p[j];
    }
    FN(mul)(v1, p, p, m);
    FN(sub)(v1, v1, two, m);

    for (int b = 62 - __builtin_clzll(k); b >= 0; b--) {
        uint64_t t[NL];
        FN(mul)(t, v0, v1, m);
        FN(sub)(t, t, p, m);
        if ((k >> b) & 1) {
            for (int j = 0; j < NL; j++) {
                v0[j] = t[j];
            }
            FN(mul)(v1, v1, v1, m);
            FN(sub)(v1, v1, two, m);
        } else {
            for (int j = 0; j < NL; j++) {
                v1[j] = t[j];
            }
            FN(mul)(v0, v0, v0, m);
            FN(sub)(v0, v0, two, m);
        }
    }
}

//
// v <- V_E(p): V_{2k} = V_k^2 - 2, V_{2k+1} = V_k.V_{k+1} - p.
//
static void FN(pp1_stage1)(uint64_t* const v, const uint64_t* const p,
                           const mont_ecm_plan_t* const plan,
                           const mod_t* const m) {
    uint64_t two[NL], v1[NL], t[NL];
    FN(add)(two, m->one, m->one, m);

    for (int j = 0; j < NL; j++) {
        v[j] = p[j];
    }
    FN(mul)(v1, p, p, m);
    FN(sub)(v1, v1, two, m);

    for (long i = (long) plan->ebits - 2; i >= 0; i--) {
        FN(mul)(t, v, v1, m);
        FN(sub)(t, t, p, m);
        if (E_BIT(plan, i)) {
            for (int j = 0; j < NL; j++) {
                v[j] = t[j];
            }
            FN(mul)(v1, v1, v1, m);
            FN(sub)(v1, v1, two, m);
        } else {
            for (int j = 0; j < NL; j++) {
                v1[j] = t[j];
            }
            FN(mul)(v, v, v, m);
            FN(sub)(v, v, two, m);
        }
    }
}

//
// Stage 2 shared by p-1 (q = a + 1/a, so that V_k(q) = a^k + a^-k) and p+1
// (q = V_E(P)): the product over the primes k.d +/- j of the plan of
// V_{k.d}(q) - V_j(q), which vanishes mod a prime factor p of n when
// a^(k.d - j) or a^(k.d + j) is 1 mod p. One product covers both.
//
static void FN(lucas_stage2)(uint64_t* const acc, const uint64_t* const q,
                             const mont_ecm_plan_t* const plan,
                             const mod_t* const m) {
    const uint32_t half = plan->d / 2;
    uint64_t (*const baby)[NL] = malloc(plan->nbaby * sizeof(*baby));

    uint64_t two[NL], v2[NL], prev[NL], cur[NL], next[NL];
    FN(add)(two, m->one, m->one, m);
    FN(mul)(v2, q, q, m);
    FN(sub)(v2, v2, two, m);
    //
    // Baby steps: V_{j+2} = V_j.V_2 - V_{j-2}, with V_{-1} = V_1.
    //
    if (plan->jslot[1] != UINT16_MAX) {
        for (int j = 0; j < NL; j++) {
            baby[plan->jslot[1]][j] = q[j];
        }
    }
    for (int j = 0; j < NL; j++) {
        prev[j] = q[j];
    }
    FN(mul)(cur, q, v2, m);
    FN(sub)(cur, cur, q, m);
    for (uint32_t i = 3; ; i += 2) {
        if (plan->jslot[i] != UINT16_MAX) {
            for (int j = 0; j < NL; j++) {
                baby[plan->jslot[i]][j] = cur[j];
            }
        }
        if (i == half) {
            break;
        }
        FN(mul)(next, cur, v2, m);
        FN(sub)(next, next, prev, m);
        for (int j = 0; j < NL; j++) {
            prev[j] = cur[j];
            cur[j]  = next[j];
        }
    }
    //
    // cur = V_{d/2}, so V_d = cur^2 - 2. Giant steps G_k = V_{k.d} from
    // k = plan->m_first on, by G_{k+1} = G_k.V_d - G_{k-1}.
    //
    uint64_t vd[NL], gprev[NL], gcur[NL], t[NL];
    FN(mul)(vd, cur, cur, m);
    FN(sub)(vd, vd, two, m);

    if (plan->m_first == 1) {
        for (int j = 0; j < NL; j++) {
            gprev[j] = two[j];
            gcur[j]  = vd[j];
        }
    } else {
        FN(lucas_ladder)(gprev, gcur, vd, plan->m_first - 1, m);
    }
    for (int j = 0; j < NL; j++) {
        acc[j] = m->one[j];
    }
    for (uint32_t k = 0; k < plan->nm; k++) {
        if (k > 0) {
            FN(mul)(next, gcur, vd, m);
            FN(sub)(next, next, gprev, m);
            for (int j = 0; j < NL; j++) {
                gprev[j] = gcur[j];
                gcur[j]  = next[j];
            }
        }
        for (uint32_t i = plan->offsets[k]; i < plan->offsets[k + 1]; i++) {
            FN(sub)(t, gcur, baby[plan->pairs[i]], m);
            FN(mul)(acc, acc, t, m);
        }
    }
    free(baby);
}

#undef E_BIT
