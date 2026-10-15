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
 * \file    factor_u64.c
 *
 * \brief Primality testing and factorization of integers below 2^64 with
 *        machine-word arithmetic.
 */

#include <math.h>
#include <pthread.h>

#include "factor_u64.h"

#if TIFA_HAVE_FACTOR_U64

#include "first_primes.h"

typedef unsigned __int128 u128;

//------------------------------------------------------------------------------
//                         Montgomery arithmetic mod n
//------------------------------------------------------------------------------

//
// Residues are kept in Montgomery form (a.2^64 mod n) and in [0, n). n odd.
//
typedef struct {
    uint64_t n;
    uint64_t ninv;   // n^-1 mod 2^64
    uint64_t one;    // 2^64 mod n
    uint64_t r2;     // 2^128 mod n
} mont_t;

static void mont_init(mont_t* const m, const uint64_t n) {
    uint64_t x = (3 * n) ^ 2;
    for (int i = 0; i < 4; i++) {
        x *= 2 - n * x;
    }
    m->n    = n;
    m->ninv = x;
    m->one  = (uint64_t) ((((u128) 1) << 64) % n);
    m->r2   = (uint64_t) (((u128) m->one * m->one) % n);
}

static inline uint64_t mont_mul(const uint64_t a, const uint64_t b,
                                const mont_t* const m) {
    const u128     t  = (u128) a * b;
    const uint64_t lo = (uint64_t) t;
    const uint64_t hi = (uint64_t) (t >> 64);
    const uint64_t q  = lo * m->ninv;
    const uint64_t qn = (uint64_t) (((u128) q * m->n) >> 64);
    return (hi < qn) ? (hi - qn + m->n) : (hi - qn);
}

static inline uint64_t mont_add(const uint64_t a, const uint64_t b,
                                const mont_t* const m) {
    const uint64_t s = a + b;
    return ((s < a) || (s >= m->n)) ? (s - m->n) : s;
}

static inline uint64_t mont_sub(const uint64_t a, const uint64_t b,
                                const mont_t* const m) {
    return (a < b) ? (a - b + m->n) : (a - b);
}

static inline uint64_t to_mont(const uint64_t a, const mont_t* const m) {
    return mont_mul(a % m->n, m->r2, m);
}

//------------------------------------------------------------------------------
//                              Integer helpers
//------------------------------------------------------------------------------

static uint64_t gcd_u64(uint64_t a, uint64_t b) {
    if (a == 0) {
        return b;
    }
    if (b == 0) {
        return a;
    }
    const int shift = __builtin_ctzll(a | b);
    a >>= __builtin_ctzll(a);
    do {
        b >>= __builtin_ctzll(b);
        if (a > b) {
            uint64_t t = a;
            a = b;
            b = t;
        }
        b -= a;
    } while (b != 0);
    return a << shift;
}

static uint64_t isqrt_u64(const uint64_t x) {
    uint64_t r = (uint64_t) sqrt((double) x);
    while ((u128) r * r > x) {
        r--;
    }
    while ((u128) (r + 1) * (r + 1) <= x) {
        r++;
    }
    return r;
}

static bool is_square_u64(const uint64_t x, uint64_t* const root) {
    //
    // Squares mod 64 and mod 63 first: rejects 95% of non-squares cheaply.
    //
    if (!((0x0202021202030213ULL >> (x & 63)) & 1)) {
        return false;
    }
    if (!((0x0402483012450293ULL >> (x % 63)) & 1)) {
        return false;
    }
    const uint64_t r = isqrt_u64(x);
    if (r * r != x) {
        return false;
    }
    *root = r;
    return true;
}

//
// a^-1 mod n, or 0 with *g = gcd(a, n) when a is not invertible.
//
static uint64_t inv_mod_u64(const uint64_t a, const uint64_t n,
                            uint64_t* const g) {
    __int128 t  = 0;
    __int128 nt = 1;
    uint64_t r  = n;
    uint64_t nr = a % n;
    while (nr != 0) {
        const uint64_t q   = r / nr;
        const uint64_t tmp = r - q * nr;
        r  = nr;
        nr = tmp;
        const __int128 tt = t - (__int128) q * nt;
        t  = nt;
        nt = tt;
    }
    *g = r;
    if (r != 1) {
        return 0;
    }
    return (uint64_t) ((t < 0) ? t + n : t);
}

//------------------------------------------------------------------------------
//                     Trial division by the primes below 1000
//------------------------------------------------------------------------------

//
// p divides n iff n.p^-1 mod 2^64 <= (2^64 - 1)/p (Granlund & Montgomery).
//
typedef struct {
    uint32_t p;
    uint64_t inv;
    uint64_t lim;
} tdiv_entry_t;

static const tdiv_entry_t tdiv_table[] = {
    {   3U, 0xaaaaaaaaaaaaaaabULL, 0x5555555555555555ULL },
    {   5U, 0xcccccccccccccccdULL, 0x3333333333333333ULL },
    {   7U, 0x6db6db6db6db6db7ULL, 0x2492492492492492ULL },
    {  11U, 0x2e8ba2e8ba2e8ba3ULL, 0x1745d1745d1745d1ULL },
    {  13U, 0x4ec4ec4ec4ec4ec5ULL, 0x13b13b13b13b13b1ULL },
    {  17U, 0xf0f0f0f0f0f0f0f1ULL, 0x0f0f0f0f0f0f0f0fULL },
    {  19U, 0x86bca1af286bca1bULL, 0x0d79435e50d79435ULL },
    {  23U, 0xd37a6f4de9bd37a7ULL, 0x0b21642c8590b216ULL },
    {  29U, 0x34f72c234f72c235ULL, 0x08d3dcb08d3dcb08ULL },
    {  31U, 0xef7bdef7bdef7bdfULL, 0x0842108421084210ULL },
    {  37U, 0x14c1bacf914c1badULL, 0x06eb3e45306eb3e4ULL },
    {  41U, 0x8f9c18f9c18f9c19ULL, 0x063e7063e7063e70ULL },
    {  43U, 0x82fa0be82fa0be83ULL, 0x05f417d05f417d05ULL },
    {  47U, 0x51b3bea3677d46cfULL, 0x0572620ae4c415c9ULL },
    {  53U, 0x21cfb2b78c13521dULL, 0x04d4873ecade304dULL },
    {  59U, 0xcbeea4e1a08ad8f3ULL, 0x0456c797dd49c341ULL },
    {  61U, 0x4fbcda3ac10c9715ULL, 0x04325c53ef368eb0ULL },
    {  67U, 0xf0b7672a07a44c6bULL, 0x03d226357e16ece5ULL },
    {  71U, 0x193d4bb7e327a977ULL, 0x039b0ad12073615aULL },
    {  73U, 0x7e3f1f8fc7e3f1f9ULL, 0x0381c0e070381c0eULL },
    {  79U, 0x9b8b577e613716afULL, 0x033d91d2a2067b23ULL },
    {  83U, 0xa3784a062b2e43dbULL, 0x03159721ed7e7534ULL },
    {  89U, 0xf47e8fd1fa3f47e9ULL, 0x02e05c0b81702e05ULL },
    {  97U, 0xa3a0fd5c5f02a3a1ULL, 0x02a3a0fd5c5f02a3ULL },
    { 101U, 0x3a4c0a237c32b16dULL, 0x0288df0cac5b3f5dULL },
    { 103U, 0xdab7ec1dd3431b57ULL, 0x027c45979c95204fULL },
    { 107U, 0x77a04c8f8d28ac43ULL, 0x02647c69456217ecULL },
    { 109U, 0xa6c0964fda6c0965ULL, 0x02593f69b02593f6ULL },
    { 113U, 0x90fdbc090fdbc091ULL, 0x0243f6f0243f6f02ULL },
    { 127U, 0x7efdfbf7efdfbf7fULL, 0x0204081020408102ULL },
    { 131U, 0x03e88cb3c9484e2bULL, 0x01f44659e4a42715ULL },
    { 137U, 0xe21a291c077975b9ULL, 0x01de5d6e3f8868a4ULL },
    { 139U, 0x3aef6ca970586723ULL, 0x01d77b654b82c339ULL },
    { 149U, 0xdf5b0f768ce2cabdULL, 0x01b7d6c3dda338b2ULL },
    { 151U, 0x6fe4dfc9bf937f27ULL, 0x01b2036406c80d90ULL },
    { 157U, 0x5b4fe5e92c0685b5ULL, 0x01a16d3f97a4b01aULL },
    { 163U, 0x1f693a1c451ab30bULL, 0x01920fb49d0e228dULL },
    { 167U, 0x8d07aa27db35a717ULL, 0x01886e5f0abb0499ULL },
    { 173U, 0x882383b30d516325ULL, 0x017ad2208e0ecc35ULL },
    { 179U, 0xed6866f8d962ae7bULL, 0x016e1f76b4337c6cULL },
    { 181U, 0x3454dca410f8ed9dULL, 0x016a13cd15372904ULL },
    { 191U, 0x1d7ca632ee936f3fULL, 0x01571ed3c506b39aULL },
    { 193U, 0x70bf015390948f41ULL, 0x015390948f40feacULL },
    { 197U, 0xc96bdb9d3d137e0dULL, 0x014cab88725af6e7ULL },
    { 199U, 0x2697cc8aef46c0f7ULL, 0x0149539e3b2d066eULL },
    { 211U, 0xc0e8f2a76e68575bULL, 0x013698df3de07479ULL },
    { 223U, 0x687763dfdb43bb1fULL, 0x0125e22708092f11ULL },
    { 227U, 0x1b10ea929ba144cbULL, 0x0120b470c67c0d88ULL },
    { 229U, 0x1d10c4c0478bbcedULL, 0x011e2ef3b3fb8744ULL },
    { 233U, 0x63fb9aeb1fdcd759ULL, 0x0119453808ca29c0ULL },
    { 239U, 0x64afaa4f437b2e0fULL, 0x0112358e75d30336ULL },
    { 241U, 0xf010fef010fef011ULL, 0x010fef010fef010fULL },
    { 251U, 0x28cbfbeb9a020a33ULL, 0x0105197f7d734041ULL },
    { 257U, 0xff00ff00ff00ff01ULL, 0x00ff00ff00ff00ffULL },
    { 263U, 0xd624fd1470e99cb7ULL, 0x00f92fb2211855a8ULL },
    { 269U, 0x8fb3ddbd6205b5c5ULL, 0x00f3a0d52cba8723ULL },
    { 271U, 0xd57da36ca27acdefULL, 0x00f1d48bcee0d399ULL },
    { 277U, 0xee70c03b25e4463dULL, 0x00ec979118f3fc4dULL },
    { 281U, 0xc5b1a6b80749cb29ULL, 0x00e939651fe2d8d3ULL },
    { 283U, 0x47768073c9b97113ULL, 0x00e79372e225fe30ULL },
    { 293U, 0x2591e94884ce32adULL, 0x00dfac1f74346c57ULL },
    { 307U, 0xf02806abc74be1fbULL, 0x00d578e97c3f5fe5ULL },
    { 311U, 0x7ec3e8f3a7198487ULL, 0x00d2ba083b445250ULL },
    { 313U, 0x58550f8a39409d09ULL, 0x00d161543e28e502ULL },
    { 317U, 0xec9e48ae6f71de15ULL, 0x00cebcf8bb5b4169ULL },
    { 331U, 0x2ff3a018bfce8063ULL, 0x00c5fe740317f9d0ULL },
    { 337U, 0x7f9ec3fcf61fe7b1ULL, 0x00c2780613c0309eULL },
    { 347U, 0x89f5abe570e046d3ULL, 0x00bcdd535db1cc5bULL },
    { 349U, 0xda971b23f1545af5ULL, 0x00bbc8408cd63069ULL },
    { 353U, 0x79d5f00b9a7862a1ULL, 0x00b9a7862a0ff465ULL },
    { 359U, 0x4dba1df32a128a57ULL, 0x00b68d31340e4307ULL },
    { 367U, 0x87530217b7747d8fULL, 0x00b2927c29da5519ULL },
    { 373U, 0x30baae53bb5e06ddULL, 0x00afb321a1496fdfULL },
    { 379U, 0xee70206c12e9b5b3ULL, 0x00aceb0f891e6551ULL },
    { 383U, 0xcdde9462ec9dbe7fULL, 0x00ab1cbdd3e2970fULL },
    { 389U, 0xafb64b05ec41cf4dULL, 0x00a87917088e262bULL },
    { 397U, 0x02944ff5aec02945ULL, 0x00a513fd6bb00a51ULL },
    { 401U, 0x2cb033128382df71ULL, 0x00a36e71a2cb0331ULL },
    { 409U, 0x1ccacc0c84b1c2a9ULL, 0x00a03c1688732b30ULL },
    { 419U, 0x19a93db575eb3a0bULL, 0x009c69169b30446dULL },
    { 421U, 0xcebeef94fa86fe2dULL, 0x009baade8e4a2f6eULL },
    { 431U, 0x6faa77fb3f8df54fULL, 0x00980e4156201301ULL },
    { 433U, 0x68a58af00975a751ULL, 0x00975a750ff68a58ULL },
    { 439U, 0xd56e36d0c3efac07ULL, 0x009548e4979e0829ULL },
    { 443U, 0xd8b44c47a8299b73ULL, 0x0093efd1c50e726bULL },
    { 449U, 0x02d9ccaf9ba70e41ULL, 0x0091f5bcb8bb02d9ULL },
    { 457U, 0x0985e1c023d9e879ULL, 0x008f67a1e3fdc261ULL },
    { 461U, 0x2a343316c494d305ULL, 0x008e2917e0e702c6ULL },
    { 463U, 0x70cb7916ab67652fULL, 0x008d8be33f95d715ULL },
    { 467U, 0xd398f132fb10fe5bULL, 0x008c55841c815ed5ULL },
    { 479U, 0x6f2a38a6bf54fa1fULL, 0x0088d180cd3a4133ULL },
    { 487U, 0x211df689b98f81d7ULL, 0x00869222b1acf1ceULL },
    { 491U, 0x0e994983e90f1ec3ULL, 0x0085797b917765abULL },
    { 499U, 0xad671e44bed87f3bULL, 0x008355ace3c897dbULL },
    { 503U, 0xf9623a0516e70fc7ULL, 0x00824a4e60b3262bULL },
    { 509U, 0x4b7129be9dece355ULL, 0x0080c121b28bd1baULL },
    { 521U, 0x190f3b7473f62c39ULL, 0x007dc9f3397d4c29ULL },
    { 523U, 0x63dacc9aad46f9a3ULL, 0x007d4ece8fe88139ULL },
    { 541U, 0xc1108fda24e8d035ULL, 0x0079237d65bcce50ULL },
    { 547U, 0xb77578472319bd8bULL, 0x0077cf53c5f7936cULL },
    { 557U, 0x473d20a1c7ed9da5ULL, 0x0075a8accfbdd11eULL },
    { 563U, 0xfbe85af0fea2c8fbULL, 0x007467ac557c228eULL },
    { 569U, 0x58a1f7e6ce0f4c09ULL, 0x00732d70ed8db8e9ULL },
    { 571U, 0x1a00e58c544986f3ULL, 0x0072c62a24c3797fULL },
    { 577U, 0x7194a17f55a10dc1ULL, 0x007194a17f55a10dULL },
    { 587U, 0x7084944785e33763ULL, 0x006fa549b41da7e7ULL },
    { 593U, 0xba10679bd84886b1ULL, 0x006e8419e6f61221ULL },
    { 599U, 0xebe9c6bb31260967ULL, 0x006d68b5356c207bULL },
    { 601U, 0x97a3fe4bd1ff25e9ULL, 0x006d0b803685c01bULL },
    { 607U, 0x6c6388395b84d99fULL, 0x006bf790a8b2d207ULL },
    { 613U, 0x8c51da6a1335df6dULL, 0x006ae907ef4b96c2ULL },
    { 617U, 0x46f3234475d5add9ULL, 0x006a37991a23aeadULL },
    { 619U, 0x905605ca3c619a43ULL, 0x0069dfbdd4295b66ULL },
    { 631U, 0xcee8dff304767747ULL, 0x0067dc4c45c8033eULL },
    { 641U, 0xff99c27f00663d81ULL, 0x00663d80ff99c27fULL },
    { 643U, 0xacca407f671ddc2bULL, 0x0065ec17e3559948ULL },
    { 647U, 0xe71298bac1e12337ULL, 0x00654ac835cfba5cULL },
    { 653U, 0xfa1e94309cd09045ULL, 0x00645c854ae10772ULL },
    { 659U, 0xbebccb8e91496b9bULL, 0x006372990e5f901fULL },
    { 661U, 0x312fa30cc7d7b8bdULL, 0x006325913c07beefULL },
    { 673U, 0x6160ff9e9f006161ULL, 0x006160ff9e9f0061ULL },
    { 677U, 0x6b03673b5e28152dULL, 0x0060cdb520e5e88eULL },
    { 683U, 0xfe802ffa00bfe803ULL, 0x005ff4017fd005ffULL },
    { 691U, 0xe66fe25c9e907c7bULL, 0x005ed79e31a4dccdULL },
    { 701U, 0x3f8b236c76528895ULL, 0x005d7d42d48ac5efULL },
    { 709U, 0xf6f923bf01ce2c0dULL, 0x005c6f35ccba5028ULL },
    { 719U, 0x6c3d3d98bed7c42fULL, 0x005b2618ec6ad0a5ULL },
    { 727U, 0x30981efcd4b010e7ULL, 0x005a2553748e42e7ULL },
    { 733U, 0x6f691fc81ebbe575ULL, 0x0059686cf744cd5bULL },
    { 739U, 0xb10480ddb47b52cbULL, 0x0058ae97bab79976ULL },
    { 743U, 0x74cd59ed64f3f0d7ULL, 0x0058345f1876865fULL },
    { 751U, 0x0105cb81316d6c0fULL, 0x005743d5bb24795aULL },
    { 757U, 0x9be64c6d91c1195dULL, 0x005692c4d1ab74abULL },
    { 761U, 0x71b3f945a27b1f49ULL, 0x00561e46a4d5f337ULL },
    { 769U, 0x77d80d50e508fd01ULL, 0x005538ed06533997ULL },
    { 773U, 0xa5eb778e133551cdULL, 0x0054c807f2c0bec2ULL },
    { 787U, 0x18657d3c2d8a3f1bULL, 0x005345efbc572d36ULL },
    { 797U, 0x2e40e220c34ad735ULL, 0x00523a758f941345ULL },
    { 809U, 0xa76593c70a714919ULL, 0x005102370f816c89ULL },
    { 811U, 0x1eef452124eea383ULL, 0x0050cf129fb94acfULL },
    { 821U, 0x38206dc242ba771dULL, 0x004fd31941cafdd1ULL },
    { 823U, 0x4cd4c35807772287ULL, 0x004fa1704aa75945ULL },
    { 827U, 0x83de917d5e69ddf3ULL, 0x004f3ed6d45a63adULL },
    { 829U, 0x882ef0403b4a6c15ULL, 0x004f0de57154ebedULL },
    { 839U, 0xf8fb6c51c606b677ULL, 0x004e1cae8815f811ULL },
    { 853U, 0xb4abaac446d3e1fdULL, 0x004cd47ba5f6ff19ULL },
    { 857U, 0xa9f83bbe484a14e9ULL, 0x004c78ae734df709ULL },
    { 859U, 0x0bebbc0d1ce874d3ULL, 0x004c4b19ed85cfb8ULL },
    { 863U, 0xbd418eaf0473189fULL, 0x004bf093221d1218ULL },
    { 877U, 0x44e3af6f372b7e65ULL, 0x004aba3c21dc633fULL },
    { 881U, 0xc87fdace4f9e5d91ULL, 0x004a6360c344de00ULL },
    { 883U, 0xec93479c446bd9bbULL, 0x004a383e9f74d68aULL },
    { 887U, 0xdac4d592e777c647ULL, 0x0049e28fbabb9940ULL },
    { 907U, 0xa63ea8c8f61f0c23ULL, 0x0048417b57c78cd7ULL },
    { 911U, 0xe476062ea5cbbb6fULL, 0x0047f043713f3a2bULL },
    { 919U, 0xdf68761c69daac27ULL, 0x00474ff2a10281cfULL },
    { 929U, 0xb813d737637aa061ULL, 0x00468b6f9a978f91ULL },
    { 937U, 0xa3a77aac1fb15099ULL, 0x0045f13f1caff2e2ULL },
    { 941U, 0x17f0c3e0712c5825ULL, 0x0045a5228cec23e9ULL },
    { 947U, 0xfd912a70ff30637bULL, 0x0045342c556c66b9ULL },
    { 953U, 0xfbb3b5dc01131289ULL, 0x0044c4a23feeced7ULL },
    { 967U, 0x856d560a0f5acdf7ULL, 0x0043c5c20d3c9fe6ULL },
    { 971U, 0x96472f314d3f89e3ULL, 0x00437e494b239798ULL },
    { 977U, 0xa76f5c7ed2253531ULL, 0x0043142d118e47cbULL },
    { 983U, 0x816eae7c7bf69fe7ULL, 0x0042ab5c73a13458ULL },
    { 991U, 0xb6a2bea4cfb1781fULL, 0x004221950db0f3dbULL },
    { 997U, 0xa3900c53318e81edULL, 0x0041bbb2f80a4553ULL },
};
#define NTDIV (sizeof(tdiv_table) / sizeof(tdiv_table[0]))
#define TDIV_BOUND 1000
//
// Primes tried by is_prime_u64 before its (conclusive) strong probable prime
// tests: more costs primes and the composites that reach them, as callers
// have usually trial-divided already.
//
#define IS_PRIME_NTDIV 12

//------------------------------------------------------------------------------
//                                 Primality
//------------------------------------------------------------------------------

static bool sprp_u64(const uint64_t a, const mont_t* const m) {
    const uint64_t n = m->n;
    uint64_t d = n - 1;
    const int s = __builtin_ctzll(d);
    d >>= s;

    const uint64_t am = to_mont(a, m);
    if (am == 0) {
        return true;
    }
    const uint64_t mone = n - m->one;
    uint64_t x = m->one;
    for (int b = 63 - __builtin_clzll(d); b >= 0; b--) {
        x = mont_mul(x, x, m);
        if ((d >> b) & 1) {
            x = mont_mul(x, am, m);
        }
    }
    if ((x == m->one) || (x == mone)) {
        return true;
    }
    for (int i = 1; i < s; i++) {
        x = mont_mul(x, x, m);
        if (x == mone) {
            return true;
        }
    }
    return false;
}

bool is_prime_u64(const uint64_t n) {
    if ((n & 1) == 0) {
        return n == 2;
    }
    if (n == 1) {
        return false;
    }
    for (uint32_t i = 0; i < IS_PRIME_NTDIV; i++) {
        const uint64_t p = tdiv_table[i].p;
        if (n == p) {
            return true;
        }
        if (n * tdiv_table[i].inv <= tdiv_table[i].lim) {
            return false;
        }
        if (p * p > n) {
            return true;
        }
    }
    //
    // Deterministic base sets: {2, 7, 61} below 2^32 (Jaeschke), and the
    // seven bases found by Jim Sinclair for 64-bit integers.
    //
    static const uint64_t bases32[] = { 2, 7, 61 };
    static const uint64_t bases64[] = {
        2, 325, 9375, 28178, 450775, 9780504, 1795265022
    };
    const bool      small  = (n < (1ULL << 32));
    const uint64_t* bases  = small ? bases32 : bases64;
    const uint32_t  nbases = small ? 3 : 7;

    mont_t m;
    mont_init(&m, n);
    for (uint32_t i = 0; i < nbases; i++) {
        if (!sprp_u64(bases[i], &m)) {
            return false;
        }
    }
    return true;
}

//------------------------------------------------------------------------------
//                          Hart's one line factoring
//------------------------------------------------------------------------------

#if !defined(HART_MAX_BITS)
#define HART_MAX_BITS   40
#endif
#if !defined(HART_MAX_ITERS)
#define HART_MAX_ITERS  4096
#endif
#define HART_MULTIPLIER 480

//
// W. B. Hart, "A one line factoring algorithm", J. Aust. Math. Soc. 92
// (2012). The multiplier 480 makes s^2 - k.n square much more often. Returns
// a proper factor of the odd composite n < 2^42, or 0.
//
static uint64_t hart_u64(const uint64_t n) {
    const uint64_t nm = n * HART_MULTIPLIER;
    //
    // Keep k.n below (2^32 - 2)^2 so that s^2 fits in 64 bits. Double
    // precision square roots are then within one of the exact ones, and
    // exact for the residues s^2 - k.n < 2^34 tested for squareness.
    //
    uint64_t imax = (0xFFFFFFFEULL * 0xFFFFFFFEULL) / nm;
    if (imax > HART_MAX_ITERS) {
        imax = HART_MAX_ITERS;
    }
    uint64_t ni = nm;
    for (uint64_t i = 1; i <= imax; i++, ni += nm) {
        const uint64_t s = (uint64_t) sqrt((double) ni) + 1;
        const uint64_t r = s * s - ni;

        if (!((0x0202021202030213ULL >> (r & 63)) & 1)) {
            continue;
        }
        const uint64_t t = (uint64_t) sqrt((double) r);
        if (t * t == r) {
            const uint64_t g = gcd_u64(n, s - t);
            if ((g > 1) && (g < n)) {
                return g;
            }
        }
    }
    return 0;
}

//------------------------------------------------------------------------------
//                      ECM with Montgomery curves (x, z)
//------------------------------------------------------------------------------

typedef struct {
    uint64_t x;
    uint64_t z;
} xz_t;

static inline xz_t xdbl(const xz_t p, const uint64_t a24,
                        const mont_t* const m) {
    const uint64_t s  = mont_add(p.x, p.z, m);
    const uint64_t d  = mont_sub(p.x, p.z, m);
    const uint64_t s2 = mont_mul(s, s, m);
    const uint64_t d2 = mont_mul(d, d, m);
    const uint64_t t  = mont_sub(s2, d2, m);
    xz_t r;
    r.x = mont_mul(s2, d2, m);
    r.z = mont_mul(t, mont_add(d2, mont_mul(a24, t, m), m), m);
    return r;
}

//
// p + q, given their difference diff = p - q.
//
static inline xz_t xadd(const xz_t p, const xz_t q, const xz_t diff,
                        const mont_t* const m) {
    const uint64_t u = mont_mul(mont_sub(p.x, p.z, m),
                                mont_add(q.x, q.z, m), m);
    const uint64_t v = mont_mul(mont_add(p.x, p.z, m),
                                mont_sub(q.x, q.z, m), m);
    const uint64_t s = mont_add(u, v, m);
    const uint64_t d = mont_sub(u, v, m);
    xz_t r;
    r.x = mont_mul(diff.z, mont_mul(s, s, m), m);
    r.z = mont_mul(diff.x, mont_mul(d, d, m), m);
    return r;
}

//
// xadd with a difference (xd : 1): one multiplication less.
//
static inline xz_t xadd_norm(const xz_t p, const xz_t q, const uint64_t xd,
                             const mont_t* const m) {
    const uint64_t u = mont_mul(mont_sub(p.x, p.z, m),
                                mont_add(q.x, q.z, m), m);
    const uint64_t v = mont_mul(mont_add(p.x, p.z, m),
                                mont_sub(q.x, q.z, m), m);
    const uint64_t s = mont_add(u, v, m);
    const uint64_t d = mont_sub(u, v, m);
    xz_t r;
    r.x = mont_mul(s, s, m);
    r.z = mont_mul(xd, mont_mul(d, d, m), m);
    return r;
}

//
// Montgomery form of the inverse of a (in Montgomery form), or 0 with
// *g = gcd(a, n) when a is not invertible. r3 = 2^192 mod n.
//
static uint64_t mont_inv(const uint64_t a, const uint64_t r3,
                         const mont_t* const m, uint64_t* const g) {
    const uint64_t x = inv_mod_u64(a, m->n, g);
    return (x == 0) ? 0 : mont_mul(x, r3, m);
}

//
// The work that depends only on the bounds: the stage 1 multiplier E (the
// prime powers up to b1, and the primes up to d/2, which stage 2 cannot
// reach), and the stage 2 schedule. Stage 2 covers the primes q in (b1, b2]
// as q = m.d +/- j, j odd and coprime to d, d = 2 mod 4: for each giant step
// m, the slots (among the j coprime to d) to multiply in, one per pair.
//
#define ECM_MAX_EWORDS 8
#define ECM_MAX_BABY   64
#define ECM_MAX_GIANTS 128
#define ECM_MAX_PAIRS  1024

typedef struct {
    uint32_t b1;
    uint32_t b2;
    uint32_t d;
    uint32_t ebits;
    uint64_t e[ECM_MAX_EWORDS];
    uint32_t nbaby;
    uint32_t nslots;
    uint8_t  slot_baby[ECM_MAX_BABY];
    uint32_t mfirst;
    uint32_t mlast;
    uint16_t offsets[ECM_MAX_GIANTS + 1];
    uint8_t  slots[ECM_MAX_PAIRS];
} ecm_plan_t;

//
// Bounds by size of n, for its smallest factor (at most half its bits),
// tuned on balanced semiprimes; d = 30 or 90 suits these b2. The ECM_MAX_*
// sizes above hold their plans.
//
static const struct {
    uint32_t bits;
    uint32_t b1;
    uint32_t b2;
    uint32_t d;
} ecm_bounds[] = {
    { 46,  40,  800, 30 },
    { 50,  70, 2100, 90 },
    { 54,  85, 2550, 90 },
    { 58, 100, 3000, 90 },
    { 62, 150, 3000, 90 },
    { 64, 200, 6000, 90 }
};
#define ECM_NPLANS (sizeof(ecm_bounds) / sizeof(ecm_bounds[0]))

static ecm_plan_t     ecm_plans[ECM_NPLANS];
static pthread_once_t ecm_plans_once = PTHREAD_ONCE_INIT;

static void build_ecm_plan(ecm_plan_t* const pl, const uint32_t b1,
                           const uint32_t b2, const uint32_t d) {
    pl->b1 = b1;
    pl->b2 = b2;
    pl->d  = d;
    const uint32_t s1 = (b1 > d / 2) ? b1 : d / 2;

    uint32_t nw = 1;
    pl->e[0] = 1;
    for (uint32_t i = 0; first_primes[i] <= s1; i++) {
        uint64_t q = first_primes[i];
        while ((first_primes[i] <= b1) && (q <= b1 / first_primes[i])) {
            q *= first_primes[i];
        }
        u128 carry = 0;
        for (uint32_t w = 0; w < nw; w++) {
            const u128 t = (u128) pl->e[w] * q + carry;
            pl->e[w] = (uint64_t) t;
            carry    = t >> 64;
        }
        if (carry != 0) {
            pl->e[nw++] = (uint64_t) carry;
        }
    }
    pl->ebits = 64 * nw - (uint32_t) __builtin_clzll(pl->e[nw - 1]);

    uint8_t slot_of[ECM_MAX_BABY];
    pl->nbaby  = (d / 2 + 1) / 2;
    pl->nslots = 0;
    for (uint32_t i = 0; i < pl->nbaby; i++) {
        if (gcd_u64(2 * i + 1, d) == 1) {
            slot_of[i] = (uint8_t) pl->nslots;
            pl->slot_baby[pl->nslots++] = (uint8_t) i;
        } else {
            slot_of[i] = UINT8_MAX;
        }
    }
    pl->mfirst = (s1 + 1 + d / 2) / d;
    pl->mlast  = (b2 + d / 2) / d;

    uint32_t ip = 0;
    uint32_t np = 0;
    while (first_primes[ip] <= s1) {
        ip++;
    }
    for (uint32_t mi = pl->mfirst; mi <= pl->mlast; mi++) {
        bool used[ECM_MAX_BABY] = { false };
        pl->offsets[mi - pl->mfirst] = (uint16_t) np;
        while (   (first_primes[ip] <= b2)
               && (first_primes[ip] <= mi * d + d / 2)) {
            const uint32_t q = first_primes[ip];
            const uint32_t j = (q > mi * d) ? (q - mi * d) : (mi * d - q);
            const uint8_t  s = slot_of[j / 2];
            if (!used[s]) {
                used[s] = true;
                pl->slots[np++] = s;
            }
            ip++;
        }
    }
    pl->offsets[pl->mlast - pl->mfirst + 1] = (uint16_t) np;
}

static void build_ecm_plans(void) {
    for (uint32_t i = 0; i < ECM_NPLANS; i++) {
        build_ecm_plan(&ecm_plans[i], ecm_bounds[i].b1, ecm_bounds[i].b2,
                       ecm_bounds[i].d);
    }
}

//
// Curves run ECM_LANES at a time: they share E and the stage 2 schedule, so
// they proceed in lockstep, and their independent multiplications overlap
// (1.5x to 2x faster than one at a time, even counting the curves run past
// the one that succeeds).
//
#define ECM_LANES 3

//
// The curves of Suyama's family with parameters sigma, sigma + 1...,
// normalized to a base point (x0 : 1) with a single inversion. Returns a
// proper factor of n, or 0.
//
static uint64_t ecm_curves(const mont_t* const m, const uint64_t r3,
                           const uint64_t sigma, const ecm_plan_t* const pl) {
    const uint64_t n = m->n;
    uint64_t g;
    uint64_t a24[ECM_LANES];
    uint64_t x0[ECM_LANES];
    uint64_t num[ECM_LANES];
    uint64_t den[ECM_LANES];
    uint64_t u3[ECM_LANES];
    uint64_t v3[ECM_LANES];
    uint64_t pre[ECM_LANES];

    const uint64_t five = to_mont(5, m);
    for (int c = 0; c < ECM_LANES; c++) {
        const uint64_t s   = to_mont(sigma + c, m);
        const uint64_t s2  = mont_add(s, s, m);
        const uint64_t u   = mont_sub(mont_mul(s, s, m), five, m);
        const uint64_t v   = mont_add(s2, s2, m);
        const uint64_t vmu = mont_sub(v, u, m);
        u3[c]  = mont_mul(mont_mul(u, u, m), u, m);
        v3[c]  = mont_mul(mont_mul(v, v, m), v, m);
        num[c] = mont_mul(mont_mul(mont_mul(vmu, vmu, m), vmu, m),
                          mont_add(mont_add(u, u, m), mont_add(u, v, m), m),
                          m);
        den[c] = mont_mul(u3[c], v, m);
        for (int i = 0; i < 4; i++) {
            den[c] = mont_add(den[c], den[c], m);
        }
        //
        // a24 = num / den and x0 = u^3 / v^3: both from 1 / (den.v^3).
        //
        const uint64_t dv = mont_mul(den[c], v3[c], m);
        pre[c] = (c == 0) ? dv : mont_mul(pre[c - 1], dv, m);
    }
    uint64_t inv = mont_inv(pre[ECM_LANES - 1], r3, m, &g);
    if (inv == 0) {
        return ((g > 1) && (g < n)) ? g : 0;
    }
    for (int c = ECM_LANES - 1; c >= 0; c--) {
        const uint64_t ic = (c == 0) ? inv : mont_mul(inv, pre[c - 1], m);
        if (c > 0) {
            inv = mont_mul(inv, mont_mul(den[c], v3[c], m), m);
        }
        a24[c] = mont_mul(mont_mul(num[c], v3[c], m), ic, m);
        x0[c]  = mont_mul(mont_mul(u3[c], den[c], m), ic, m);
    }
    //
    // Stage 1: Montgomery ladder over the bits of E.
    //
    xz_t r0[ECM_LANES];
    xz_t r1[ECM_LANES];
    for (int c = 0; c < ECM_LANES; c++) {
        r0[c].x = x0[c];
        r0[c].z = m->one;
        r1[c]   = xdbl(r0[c], a24[c], m);
    }
    for (int b = (int) pl->ebits - 2; b >= 0; b--) {
        if ((pl->e[b / 64] >> (b % 64)) & 1) {
            for (int c = 0; c < ECM_LANES; c++) {
                r0[c] = xadd_norm(r1[c], r0[c], x0[c], m);
                r1[c] = xdbl(r1[c], a24[c], m);
            }
        } else {
            for (int c = 0; c < ECM_LANES; c++) {
                r1[c] = xadd_norm(r1[c], r0[c], x0[c], m);
                r0[c] = xdbl(r0[c], a24[c], m);
            }
        }
    }
    for (int c = 0; c < ECM_LANES; c++) {
        g = gcd_u64(r0[c].z, n);
        if ((g > 1) && (g < n)) {
            return g;
        }
        if (g == n) {
            //
            // Both factors found at once: any point does for the other
            // lanes' sake (whatever gcd stage 2 finds still divides n).
            //
            r0[c].x = m->one;
            r0[c].z = m->one;
        }
    }
    //
    // Stage 2. Baby steps j.Q for odd j <= d/2; those coprime to d are
    // normalized (one inversion for all lanes), so that each pair costs
    // x(G) - x(jQ).z(G) and one multiplication into the accumulator.
    //
    xz_t     baby[ECM_LANES][ECM_MAX_BABY];
    uint64_t xs[ECM_LANES][ECM_MAX_BABY];
    uint64_t zpre[ECM_LANES * ECM_MAX_BABY];
    xz_t     dq[ECM_LANES];

    for (int c = 0; c < ECM_LANES; c++) {
        const xz_t q2 = xdbl(r0[c], a24[c], m);
        baby[c][0] = r0[c];
        baby[c][1] = xadd(q2, r0[c], r0[c], m);
        for (uint32_t i = 2; i < pl->nbaby; i++) {
            baby[c][i] = xadd(baby[c][i - 1], q2, baby[c][i - 2], m);
        }
        dq[c] = xdbl(baby[c][pl->nbaby - 1], a24[c], m);
    }
    uint32_t t = 0;
    for (int c = 0; c < ECM_LANES; c++) {
        for (uint32_t k = 0; k < pl->nslots; k++, t++) {
            const uint64_t z = baby[c][pl->slot_baby[k]].z;
            zpre[t] = (t == 0) ? z : mont_mul(zpre[t - 1], z, m);
        }
    }
    inv = mont_inv(zpre[t - 1], r3, m, &g);
    if (inv == 0) {
        return ((g > 1) && (g < n)) ? g : 0;
    }
    for (int c = ECM_LANES - 1; c >= 0; c--) {
        for (int k = (int) pl->nslots - 1; k >= 0; k--) {
            t--;
            const xz_t     b  = baby[c][pl->slot_baby[k]];
            const uint64_t iz = (t == 0) ? inv : mont_mul(inv, zpre[t - 1], m);
            if (t > 0) {
                inv = mont_mul(inv, b.z, m);
            }
            xs[c][k] = mont_mul(b.x, iz, m);
        }
    }
    xz_t     gprev[ECM_LANES];
    xz_t     gcur[ECM_LANES];
    uint64_t acc[ECM_LANES];
    for (int c = 0; c < ECM_LANES; c++) {
        gprev[c] = dq[c];
        gcur[c]  = dq[c];
        acc[c]   = m->one;
    }
    for (uint32_t mi = 1; mi <= pl->mlast; mi++) {
        if (mi == 2) {
            for (int c = 0; c < ECM_LANES; c++) {
                gcur[c] = xdbl(dq[c], a24[c], m);
            }
        } else if (mi > 2) {
            for (int c = 0; c < ECM_LANES; c++) {
                const xz_t next = xadd(gcur[c], dq[c], gprev[c], m);
                gprev[c] = gcur[c];
                gcur[c]  = next;
            }
        }
        if (mi < pl->mfirst) {
            continue;
        }
        const uint32_t k = mi - pl->mfirst;
        for (uint32_t i = pl->offsets[k]; i < pl->offsets[k + 1]; i++) {
            const uint8_t s = pl->slots[i];
            for (int c = 0; c < ECM_LANES; c++) {
                const uint64_t d = mont_sub(gcur[c].x,
                                            mont_mul(xs[c][s], gcur[c].z, m),
                                            m);
                acc[c] = mont_mul(acc[c], d, m);
            }
        }
    }
    for (int c = 0; c < ECM_LANES; c++) {
        g = gcd_u64(acc[c], n);
        if ((g > 1) && (g < n)) {
            return g;
        }
    }
    return 0;
}

static uint64_t ecm_u64(const uint64_t n, const uint32_t bits) {
    pthread_once(&ecm_plans_once, build_ecm_plans);

    uint32_t i = 0;
    while ((i < ECM_NPLANS - 1) && (bits > ecm_bounds[i].bits)) {
        i++;
    }
    mont_t m;
    mont_init(&m, n);
    const uint64_t r3 = mont_mul(m.r2, m.r2, &m);

    for (uint64_t sigma = 7; sigma < 7 + 300; sigma += ECM_LANES) {
        const uint64_t g = ecm_curves(&m, r3, sigma, &ecm_plans[i]);
        if (g != 0) {
            return g;
        }
    }
    return 0;
}

//------------------------------------------------------------------------------
//                                Brent's rho
//------------------------------------------------------------------------------

static uint64_t rho_u64(const uint64_t n, const uint64_t c) {
    mont_t m;
    mont_init(&m, n);

    const uint64_t cm    = to_mont(c, &m);
    const uint64_t batch = 128;

    uint64_t y  = m.one;
    uint64_t x  = y;
    uint64_t ys = y;
    uint64_t q  = m.one;
    uint64_t g  = 1;

    for (uint64_t r = 1; (g == 1) && (r <= (1ULL << 26)); r <<= 1) {
        x = y;
        for (uint64_t i = 0; i < r; i++) {
            y = mont_add(mont_mul(y, y, &m), cm, &m);
        }
        for (uint64_t k = 0; (k < r) && (g == 1); k += batch) {
            ys = y;
            const uint64_t lim = (batch < r - k) ? batch : (r - k);
            for (uint64_t i = 0; i < lim; i++) {
                y = mont_add(mont_mul(y, y, &m), cm, &m);
                q = mont_mul(q, (x > y) ? (x - y) : (y - x), &m);
            }
            g = gcd_u64(q, n);
        }
    }
    if (g == n) {
        //
        // The batch overshot: replay it one step at a time.
        //
        do {
            ys = mont_add(mont_mul(ys, ys, &m), cm, &m);
            g  = gcd_u64((x > ys) ? (x - ys) : (ys - x), n);
        } while (g == 1);
    }
    return ((g > 1) && (g < n)) ? g : 0;
}

//------------------------------------------------------------------------------
//                               Factorization
//------------------------------------------------------------------------------

//
// A proper factor of the odd composite n, not a perfect square, or 0.
//
static uint64_t find_factor_u64(const uint64_t n) {
    const uint32_t bits = 64 - __builtin_clzll(n);
    uint64_t g = 0;

    if (bits <= HART_MAX_BITS) {
        g = hart_u64(n);
    }
    if (g == 0) {
        g = ecm_u64(n, bits);
    }
    for (uint64_t c = 1; (g == 0) && (c < 32); c++) {
        g = rho_u64(n, c);
    }
    return g;
}

static void add_prime(const uint64_t p, const uint32_t multi,
                      uint64_t* const primes, uint32_t* const multis,
                      uint32_t* const nprimes) {
    for (uint32_t i = 0; i < *nprimes; i++) {
        if (primes[i] == p) {
            multis[i] += multi;
            return;
        }
    }
    uint32_t i = *nprimes;
    while ((i > 0) && (primes[i - 1] > p)) {
        primes[i] = primes[i - 1];
        multis[i] = multis[i - 1];
        i--;
    }
    primes[i] = p;
    multis[i] = multi;
    (*nprimes)++;
}

uint32_t factor_u64(uint64_t n, uint64_t* const primes,
                    uint32_t* const multis) {
    if (n < 2) {
        return 0;
    }
    uint32_t nprimes = 0;

    const int ntwos = __builtin_ctzll(n);
    if (ntwos > 0) {
        add_prime(2, (uint32_t) ntwos, primes, multis, &nprimes);
        n >>= ntwos;
    }
    for (uint32_t i = 0; (i < NTDIV) && (n > 1); i++) {
        const uint64_t p = tdiv_table[i].p;
        if (p * p > n) {
            break;
        }
        uint32_t multi = 0;
        while (n * tdiv_table[i].inv <= tdiv_table[i].lim) {
            n *= tdiv_table[i].inv;
            multi++;
        }
        if (multi > 0) {
            add_prime(p, multi, primes, multis, &nprimes);
        }
    }
    if (n == 1) {
        return nprimes;
    }
    //
    // Every prime factor left is above 1000 (or n is prime): split the
    // composites found until only primes remain. 40 entries hold all the
    // pending factors of a 64-bit integer.
    //
    uint64_t pending[64];
    uint32_t npending = 0;
    pending[npending++] = n;

    while (npending > 0) {
        const uint64_t c = pending[--npending];
        uint64_t r;

        if ((c < (uint64_t) TDIV_BOUND * TDIV_BOUND) || is_prime_u64(c)) {
            add_prime(c, 1, primes, multis, &nprimes);
        } else if (is_square_u64(c, &r)) {
            pending[npending++] = r;
            pending[npending++] = r;
        } else {
            const uint64_t g = find_factor_u64(c);
            if (g == 0) {
                return 0;
            }
            pending[npending++] = g;
            pending[npending++] = c / g;
        }
    }
    return nprimes;
}

#endif
