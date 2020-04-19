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
 * \file    funcs.h
 * \author  Jerome Milan
 * \date    Tue Oct 11 2011
 * \version 2011-10-11
 *
 * \brief Number theoretical, hash and comparison functions.
 *
 * Defines several number theoretical functions, hash functions and comparison
 * functions.
 */

#if !defined(_TIFA_FUNCS_H_)
   /**
    * \def _TIFA_FUNCS_H_
    * Standard include guard.
    */
#define _TIFA_FUNCS_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <limits.h>
#include <inttypes.h>
#include <math.h>
#include <gmp.h>

#include "macros.h"
#include "array.h"
#include "tifa_config.h"

   /**
    * \def LARGEST_MULTIPLIER
    * Largest multiplier allowed.
    */
#define LARGEST_MULTIPLIER 97

   /**
    * \def BITSIZE_LARGEST_MULTIPLIER
    * Size in bits of the largest multiplier allowed.
    */
#define BITSIZE_LARGEST_MULTIPLIER 7

   /**
    * \def MAX_IPRIME_IN_MULT_CALC
    * The <tt>MAX_IPRIME_IN_MULT_CALC</tt>-th smallest prime number is the
    * largest prime used in the determination of the best multiplier.
    */
#define MAX_IPRIME_IN_MULT_CALC 31

/*
 *-----------------------------------------------------------------------------
 *                      Structures and typedefs
 *-----------------------------------------------------------------------------
 */

   /**
    * \struct struct_mult_data_t funcs.h lib/utils/include/funcs.h
    * \brief  Ad hoc structure used in the computation of the multiplier
    * to use.
    *
    * This ad-hoc structure defines several variables needed in the
    * determination of the best multiplier, as described by M. A. Morrison
    * and J. Brillhart in the remark 5.3 of the paper "A Method of Factoring
    * and the Factorization of F_7" (Mathematics of Computation, vol 29,
    * #129, Jan 1975, pages 183-205).
    */
struct struct_mult_data_t {
       /**
        * The multiplier to use in the factoring algorithms.
        */
    uint32_t multiplier;
       /**
        * The number of primes p_i less than or equal to the
        * <tt>MAX_IPRIME_IN_MULT_CALC</tt>-th prime for which the legendre
        * symbol (<tt>k</tt>*<tt>N</tt>/p_i) is 0 or 1 \e and for which
        * either (<tt>k</tt>*<tt>N</tt>/3) or (<tt>k</tt>*<tt>N</tt>/5)
        * (or both) is 0 or 1.
        */
    uint32_t count;
       /**
        * The sum of 1/p_i where {p_i} is the set of primes previously
        * described.
        */
    double   sum_inv_pi;
};

   /**
    * \typedef mult_data_t
    * \brief Equivalent to <tt>struct struct_mult_data_t</tt>.
    */
typedef struct struct_mult_data_t mult_data_t;

/*
 *-----------------------------------------------------------------------------
 *                      Miscellaneous functions
 *-----------------------------------------------------------------------------
 */

   /**
    * \brief Most significant bit of a positive integer.
    *
    * Returns the value of the most significant bit of the integer \c n
    * in essentially constant time, or in other words, its logarithm in
    * base 2. The returned result is an integer from 0 (the least significant
    * bit) to 31 included (the most significant bit).
    *
    * \note This function is adapted from public domain code from the
    * Bit Twiddling Hacks web page:
    * http://graphics.stanford.edu/~seander/bithacks.html
    *
    * \param[in] n A positive integer.
    * \returns log(<tt>n</tt>) in base 2.
    */
uint32_t most_significant_bit(uint32_t n);

   /**
    * \brief Floor of logarithm in base 2 of a positive integer.
    *
    * Returns the value of the floor of the logarithm (in base 2) of a
    * positive integer in essentially constant time. In other words, it
    * returns the greatest natural \c i such that <tt>2^i <= n</tt>.
    *
    * \note This is actually just a call to <tt>most_significant_bit</tt>.
    *
    * \param[in] n A positive integer.
    * \returns Floor of log(<tt>n</tt>) in base 2.
    */
inline static uint32_t floor_log2(uint32_t n) {
    return most_significant_bit(n);
}

   /**
    * \brief Ceil of logarithm in base 2 of a positive integer.
    *
    * Returns the value of the ceil of the logarithm (in base 2) of a
    * positive integer in essentially constant time. In other words, it
    * returns the smallest natural \c i such that <tt>2^i >= n</tt>.
    *
    * \param[in] n A positive integer.
    * \returns Ceil of log(<tt>n</tt>) in base 2.
    */
inline static uint32_t ceil_log2(uint32_t n) {
    uint32_t e = most_significant_bit(n);
    if (! IS_POWER_OF_2_UI(n)) {
        e++;
    }
    return e;
}

   /**
    * \brief Ceil of logarithm in base 2 of a \c mp_limb_t.
    *
    * Returns the value of the ceil of the logarithm (in base 2) of a
    * \c mp_limb_t in essentially constant time. In other words, it
    * returns the smallest natural \c i such that <tt>2^i >= n</tt>.
    *
    * \param[in] n A positive integer as an \c mp_limb_t.
    * \returns Ceil of log(<tt>n</tt>) in base 2.
    */
inline static uint32_t ceil_log2_mp_limb(mp_limb_t limb) {
#if GMP_LIMB_BITS == 64
    uint64_t tmp = limb & GMP_NUMB_MASK;
    if (tmp >> 32) {
        return 32 + ceil_log2((uint32_t)(tmp >> 32));
    } else {
        return ceil_log2((uint32_t)tmp);
    }
#else
    uint32_t tmp = (uint32_t)limb & GMP_NUMB_MASK;
    return ceil_log2((uint32_t)tmp);
#endif
}

   /**
    * \brief Find a coprime base from a list of factors.
    *
    * Finds a coprime base for the list of factors of \c n given by the
    * array \c *factors and stores it in the allocated but \e uninitialized
    * array \c base. After invocation, we know that \c n is smooth on the
    * returned computed base and that all elements of the base are coprime
    * to each other.
    *
    * The resulting base is obtained:
    * <ol>
    *   <li>
    *     by completing the list of original factors with their cofactors,
    *   </li>
    *   <li>
    *     by keeping only factors (or non-trivial divisors of factors)
    *     coprime to all others.
    *   </li>
    * </ol>
    *
    * \warning There is absolutely no guarantee that the returned base elements
    * are prime.  If, by chance, the base only contains primes then it means
    * that we  have found the complete factorization of \c n (up to the prime
    * multiplicities).
    *
    * \note If the \c base array has not enough room to hold all the
    * coprimes found, it will be resized via a call to
    * <tt>resize_mpz_array</tt> with \c ELONGATION extra \c mpz_t slots to
    * avoid too frequent resizes.
    *
    * \note Consecutive invocation of this function with the same \c base and
    * \c n but for different \c factors arrays will build a coprime base for
    * all elements in all the aforementioned \c factors arrays.
    *
    * \param[in, out] base The found coprime base.
    * \param[in] n A positive integer.
    * \param[in] factors A pointer to an array holding some factors of \c n.
    * \param[in, out] A pointer to the \e unintialized \c mpz_array_t to hold
    *                 the coprime base.
    */
void find_coprime_base(mpz_array_t* const base, const mpz_t n,
                       const mpz_array_t* const factors);

/*
 *-----------------------------------------------------------------------------
 *                      Number theoretical functions
 *-----------------------------------------------------------------------------
 */

   /**
    * \def NO_SQRT_MOD_P
    * Value returned by the <tt>sqrtm(n, p)</tt> function if no modular square
    * root of \c n mod \c p exits.
    */
#define NO_SQRT_MOD_P (UINT32_MAX)

   /**
    * \def NO_SQRT_MOD_P2
    * Value returned by the <tt>sqrtm_p2(n, p)</tt> function if no modular
    * square  root of \c n mod <tt>p*p</tt> exits.
    */
#define NO_SQRT_MOD_P2 (ULONG_MAX)

   /**
    * \brief Kronecker symbol restricted to positive simple precision integers.
    *
    * Returns the value of the Kronecker symbol (<tt>a</tt>/<tt>b</tt>)
    * where \c a and \c b are positive integers.
    *
    * \param[in] a A positive integer.
    * \param[in] b A positive integer.
    * \returns The value of the kronecker symbol (<tt>a</tt>/<tt>b</tt>).
    */
int8_t kronecker_ui(uint32_t a, uint32_t b);

   /**
    * \brief Modular exponentiation restricted to positive simple precision
    * integers.
    *
    * Returns (<tt>base</tt>^<tt>power</tt>) mod <tt>modulus</tt> as an
    * unsigned integer.
    *
    * \param[in] base The base of the modular exponential.
    * \param[in] power The power of the modular exponential.
    * \param[in] modulus The modulus of the modular exponential.
    * \returns The modular exponential (<tt>base</tt>^<tt>power</tt>)
    *          mod <tt>modulus</tt>.
    */
uint32_t powm(uint32_t base, uint32_t power, uint32_t modulus);

   /**
    * \brief Shanks' algorithm for modular square roots computation.
    *
    * Returns the modular square root of <tt>a</tt> (mod <tt>p</tt>) (where
    * \c p is an \em odd \em prime) using Shanks' algorithm, that is, returns
    * the positive integer <tt>s</tt> such that <tt>s</tt>^2 = <tt>a</tt>
    * (mod <tt>p</tt>). If no such integer exists, returns
    * <tt>NO_SQRT_MOD_P</tt>.
    *
    * \warning The primality of \c p is not checked by <tt>sqrtm</tt>.
    * It is the responsability of the caller to check whether \c p is
    * indeed prime. Failure to assure such a precondition will lead to
    * an infinite loop.
    *
    * \param[in] a The modular square.
    * \param[in] p The modulus.
    * \returns The modular square root of <tt>a</tt> (mod <tt>p</tt>) if it
    *          exists. <tt>NO_SQRT_MOD_P</tt> otherwise.
    */
uint32_t sqrtm(uint32_t a, uint32_t p);

   /**
    * \brief Quadratic residues mod 221
    *
    * \c x is a square mod 221 if <tt>qres_mod_221[x % 221] == 1</tt>.
    */
extern const unsigned short qres_mod_221[221] MAYBE_UNUSED;

   /**
    * \brief Quadratic residues mod 256
    *
    * \c x is a square mod 256 if <tt>qres_mod_256[x % 256] == 1</tt>.
    */
extern const unsigned short qres_mod_256[256] MAYBE_UNUSED;

   /**
    * \brief Quadratic residues mod 315
    *
    * \c x is a square mod 315 if <tt>qres_mod_315[x % 315] == 1</tt>.
    */
extern const unsigned short qres_mod_315[315] MAYBE_UNUSED;

   /**
    * \brief Perfect square detection test
    *
    * Returns sqrt(\c x) if and only if \c x is a perfect square.
    * Returns 0 otherwise.
    *
    * \param[in] x The integer to test.
    * \returns sqrt(\c x) if \c x is a perfect square. 0 otherwise.
    */
inline static unsigned long int is_square(unsigned long int x) {
    //
    // Basic perfect square detection test.
    //
    // See for example algorithm 1.7.3 from the book "A Course in Computational
    // Algebraic Number Theory" by Henri Cohen, Springer-Verlag 1993.
    //
    // The description given in this book has been adapted to use larger
    // tables in exchange of a slight performance boost (about 30% on
    // Opteron 250).
    //
    if (qres_mod_256[x & 255] == 0) {
        //
        // Get rid of about 82.8% of non squares
        //
        return 0;
    }
    if (qres_mod_315[x % 315] == 0) {
        //
        // Get rid of about 84.8% of non squares
        //
        return 0;
    }
    if (qres_mod_221[x % 221] == 0) {
        //
        // Get rid of about 71.5% of non squares
        //
        return 0;
    }
    unsigned long int root = (unsigned long int)sqrt(x);
    if ((root * root) == x) {
        return root;
    }
    return 0;
}

   /**
    * \brief Composition test for \c uint32_t integers
    *
    * Returns \c false if \c n is definitely composite. Returns \c true if
    * \c n is \e probably prime.
    *
    * \note This is actually a basic Miller-Rabin composition test with
    * \c NMILLER_RABIN iterations preceded with some trial divisions if
    * \c n is sufficiently small.
    *
    * \param[in] n The \c uint32_t to be checked for composition.
    * \returns Returns \c false if \c n is found to be definitely composite.
    *          \c true otherwise.
    */
bool is_prime(uint32_t n);

   /**
    * \brief Greatest common divisor for unsigned long int
    *
    * Returns the greatest common divisor of \c a and \c b as an unsigned
    * long int.
    *
    * \param[in] a An unsigned long int.
    * \param[in] b An unsigned long int.
    * \returns The greatest common divisor of \c a and \c b.
    */
unsigned long int gcd_ulint(unsigned long int a, unsigned long int b);

   /**
    * \brief Modular inverse for unsigned long int
    *
    * Returns the modular inverse of \c n modulo the odd prime \c p as an
    * unsigned long int.
    *
    * \warning \c p must be a positive odd prime, strictly less than
    * <tt>LONG_MAX</tt> (yes, <tt>LONG_MAX</tt> and not <tt>ULONG_MAX</tt>!)
    * and, of course, <tt>n % p</tt> must be non-null.
    *
    * \param[in] n An unsigned long int.
    * \param[in] p An odd prime unsigned long int.
    * \returns The modular inverse of \c n mod <tt>p</tt>.
    */
unsigned long int modinv_ui(unsigned long int n, unsigned long int p);

   /**
    * \brief Modular square root modulo the square of a prime.
    *
    * Provided that \c a verifies 1 <= \c a < <tt>p*p</tt>, returns the
    * modular square root of <tt>a</tt> (mod <tt>p*p</tt>) (where
    * \c p is an \em odd \em prime) that is, returns a positive integer
    * <tt>s</tt> such that <tt>s</tt>^2 = <tt>a</tt> (mod <tt>p*p</tt>).
    * If no such integer exists, returns <tt>NO_SQRT_MOD_P2</tt>.
    *
    * \warning In order to use only single precision computation, the product
    * <tt>p*p</tt> should be strictly less than <tt>LONG_MAX</tt>.
    *
    * \warning The primality of \c p is not checked by <tt>sqrtm_p2</tt>.
    * It is the responsability of the caller to check whether \c p is
    * indeed prime. Failure to assure such a precondition will lead to
    * an infinite loop.
    *
    * \param[in] a The modular square.
    * \param[in] p The square root of the modulus.
    * \returns The modular square root of <tt>a</tt> (mod <tt>p*p</tt>) if it
    *          exists. <tt>NO_SQRT_MOD_P2</tt> otherwise.
    */
unsigned long int sqrtm_p2(uint32_t a, uint32_t p);

   /**
    * \brief Find best multiplier using the Knuth-Schroeppel function.
    *
    * Given the size of factor base \c size_base, returns the "best"
    * multiplier to factor \c n, using the modified version of the
    * Knuth-Schroeppel function described by Silverman in: "The Multiple
    * Quadratic Sieve".
    *
    * \note The greatest multiplier considered is given by
    * \c LARGEST_MULTIPLIER.
    *
    * \see "The Multiple Quadratic Sieve", Robert D. Silverman,
    * <i>Mathematics of Computation</i>, Volume 48, Number 177,
    * January 1987, pages 329-339.
    *
    * \param[in] n The number to factor
    * \param[in] size_base The desired size of the factor base
    * \returns The "best" multiplier to factor n.
    */
uint32_t ks_multiplier(const mpz_t n, const uint32_t size_base);

/*
 *-----------------------------------------------------------------------------
 *                          Hash functions
 *-----------------------------------------------------------------------------
 */

   /**
    * \brief Robert Jenkins' 32 bit mix function.
    *
    * Returns the hash of a uint32_t integer (passed as a pointer to \c void)
    * using Robert Jenkins' 32 bit mix function.
    *
    * \see \link http://www.concentric.net/~Ttwang/tech/inthash.htm \endlink
    *
    * \param[in] keyptr A pointer to the <tt>uint32_t</tt> to hash.
    * \returns The value of the hash function.
    */
uint32_t hash_rj_32(const void* const keyptr);

   /**
    * \brief An hash function for strings.
    *
    * Returns the hash of a C-style character string (passed as a pointer to
    * \c void) using an hash function attributed to P.J. Weinberger.
    *
    * \note This hash function and its implementation is extracted from
    * the famous Dragon book: "Compilers: Principles, Techniques and Tools",
    * Aho, Sethi, & Ullman.
    *
    * \param[in] keyptr A pointer to the character string to hash.
    * \returns The value of the hash function.
    */
uint32_t hash_pjw(const void* const keyptr);

   /**
    * \brief The "Super Fast Hash" function By Paul Hsieh.
    *
    * Returns the hash of a C-style character string (passed as a pointer to
    * \c void) using the so-called "SuperFastHash" function By Paul Hsieh.
    *
    * \see \link  http://www.azillionmonkeys.com/qed/hash.html \endlink
    *
    * \param[in] keyptr A pointer to the character string to hash.
    * \returns The value of the hash function.
    */
uint32_t hash_sfh_ph(const void* const keyptr);

/*
 *-----------------------------------------------------------------------------
 *                          Comparison functions
 *-----------------------------------------------------------------------------
 */

   /**
    * \brief Comparison function between two \c mpz_t.
    *
    * This is a natural order comparison function between two \c mpz_t
    * elements passed as pointers to \c void. It returns:
    * \li 1 if the first \c mpz_t is greater than the second one.
    * \li 0 if the first \c mpz_t is equal to the second one.
    * \li -1 if the first \c mpz_t is less than the second one.
    *
    * \note This function is actually nothing more than a wrapper for
    * \c mpz_cmp.
    *
    * \param[in] mpza A pointer to a \c mpz_t.
    * \param[in] mpzb A pointer to another \c mpz_t.
    * \returns The comparison between the two \c mpz_t.
    */
int mpz_cmp_func(const void* const mpza, const void* const mpzb);

   /**
    * \brief Comparison function between two \c uint32_t.
    *
    * This is a natural order comparison function between two \c uint32_t
    * elements passed as pointers to \c void. It returns:
    * \li 1 if the first \c uint32_t is greater than the second one.
    * \li 0 if the first \c uint32_t is equal to the second one.
    * \li -1 if the first \c uint32_t is less than the second one.
    *
    * \param[in] uinta A pointer to a \c uint32_t.
    * \param[in] uintb A pointer to another \c uint32_t.
    * \returns The comparison between the two \c uint32_t.
    */
int uint32_cmp_func(const void* const uinta, const void* const uintb);

   /**
    * \brief Comparison function between two strings.
    *
    * This is a lexicographical order comparison function between two
    * C-style character strings passed as pointers to \c void. It returns:
    * \li 1 if the first string is greater than the second one.
    * \li 0 if the first string is equal to the second one.
    * \li -1 if the first string is less than the second one.
    *
    * \note This function is actually nothing more than a wrapper for \c strcmp.
    *
    * \param[in] stra A pointer to a C-style character string.
    * \param[in] strb A pointer to another C-style character string.
    * \returns The lexicographical comparison between the two strings.
    */
int string_cmp_func(const void* const stra, const void* const strb);

   /**
    * \brief Comparison function between two <tt>mult_data_t</tt>.
    *
    * This is a comparison function between two \c mult_data_t structures
    * passed as pointers to <tt>void</tt>, according to the criteria
    * set forth in Morrison and Brillhart's paper "A Method of Factoring
    * and the Factorization of F_7" (Mathematics of Computation, vol 29,
    * #129, Jan 1975, pages 183-205).
    *
    * If \c a and \c b are the two underlying \c mult_data_t structures
    * to compare, it returns:
    *
    * <ul>
    *   <li>1 if \c a.count > <tt>b.count</tt></li>
    *   <li> -1 if \c a.count < <tt>b.count</tt></li>
    *   <li> If \c a.count == <tt>b.count</tt>, returns:</li>
    *   <ul>
    *     <li> 1 if \c a.sum_inv_pi > <tt>b.sum_inv_pi</tt></li>
    *     <li> -1 if \c a.sum_inv_pi > <tt>b.sum_inv_pi</tt></li>
    *     <li> If \c a.sum_inv_pi == <tt>b.sum_inv_pi</tt>, returns:<li>
    *     <ul>
    *       <li> 1 if \c a.multiplier < <tt>b.multiplier</tt>
    *         (Indeed, we prefer smaller multipliers)
    *       </li>
    *       <li> -1 if \c a.multiplier > <tt>b.multiplier</tt></li>
    *       <li> 0 if if \c a.multiplier > <tt>b.multiplier</tt></li>
    *     </ul>
    *   </ul>
    * </ul>
    *
    * \param[in] mda A pointer to the first <tt>mult_data_t</tt> to compare.
    * \param[in] mdb A pointer to the second <tt>mult_data_t</tt> to compare.
    * \returns The comparison between the two <tt>mult_data_t</tt>.
    */
int cmp_mult_data(const void* mda, const void* mdb);

/*
 *-----------------------------------------------------------------------------
 *                        Combinatorial functions
 *-----------------------------------------------------------------------------
 */

   /**
    * \brief Binomial coefficient C(n, k) (n choose k).
    * <p>
    * Returns the binomial coefficient C(<tt>n</tt>, <tt>k</tt>) (i.e
    * <tt>n</tt> choose <tt>k</tt>).
    * </p>
    * <p>
    * Note that this single precision function only returns correct results if
    * the actual value of the binomial coefficient fits in 32 bits.
    * </p>
    *
    * \param[in] n
    * \param[in] k
    *
    * \returns The binomial coefficient C(n, k).
    */
uint32_t n_choose_k(uint8_t n, uint8_t k);

   /**
    * \brief Generate the successor of a fixed cardinal subset from a
    * base set, in lexicographic order.
    * <p>
    * Starting with a subset of cardinal <tt>k</tt> of a base set of cardinal
    * <tt>n</tt>, generates the subset's successor in the lexicographic order.
    * The new subset is stored in <tt>subset</tt> and thus overrides the
    * previous one.
    * </p>
    * <p>
    * Subsets are decribed by an array of length <tt>k</tt> holding indexes
    * in the interval [1, <tt>n</tt>].
    * </p>
    * <p>
    * The first <tt>k</tt>-subset in the lexicographic order is given by
    * {1, 2, 3, ... <tt>k</tt>}. After a call to \c next_subset_lex
    * \c end is \c true if and only if the last <tt>k</tt>-subset
    * has been reached (i.e. the next one will be {1, 2, 3, ... <tt>k</tt>}).
    * </p>
    * <p>
    * This is actually algorithm 2.6 from the book <i>"Combinatorial
    * Algorithms - Generation, Enumeration, and Search"</i> by Donald L. Kreher
    * and Douglas Stinson.
    * </p>
    *
    * \param[in]     n      Cardinal of the base set.
    * \param[in]     k      Cardinal of the subset.
    * \param[in,out] subset Current subset to be replaced by its successor.
    * \param[in]     end    Have we reached the end of the cycle?
    */
void next_subset_lex(uint32_t n, uint32_t k, uint32_t* subset, bool* end);

   /**
    * \brief Generate a fixed cardinal subset from a base set, according
    * to a given rank.
    * <p>
    * Starting with a base set of cardinal <tt>n</tt>, constructs a subset
    * of cardinal <tt>k</tt> and rank <tt>r</tt> (in
    * [0, c(<tt>n</tt>, <tt>k</tt>)]) where the rank is given by the
    * lexicographic order. The constructed subset is stored in <tt>subset</tt>
    * and thus overrides the previous data.
    * </p>
    * <p>
    * Subsets are decribed by an array of length <tt>k</tt> holding indexes
    * in the interval [1, <tt>n</tt>].
    * </p>
    * <p>
    * This is actually algorithm 2.8 from the book <i>"Combinatorial
    * Algorithms - Generation, Enumeration, and Search"</i> by Donald L. Kreher
    * and Douglas Stinson.
    * </p>
    *
    * \param[in]  n      Cardinal of the base set.
    * \param[in]  k      Cardinal of the subset.
    * \param[in]  r      Rank of the subset (assuming lexicographic order).
    * \param[out] subset Subset to be returned.
    */
void unrank_subset_lex(uint32_t n, uint32_t k, uint32_t r, uint32_t* subset);

   /**
    * \brief Generate a fixed cardinal subset from a base set, according
    * to a given rank.
    * <p>
    * Starting with a base set of cardinal <tt>n</tt>, constructs a subset
    * of cardinal <tt>k</tt> and rank <tt>r</tt> (in
    * [0, c(<tt>n</tt>, <tt>k</tt>)]) where the rank is given by the
    * minimal change order. The constructed subset is stored in <tt>subset</tt>
    * and thus overrides the previous data.
    * </p>
    * <p>
    * Subsets are decribed by an array of length <tt>k</tt> holding indexes
    * in the interval [1, <tt>n</tt>].
    * </p>
    * <p>
    * This is actually algorithm 2.12 from the book <i>"Combinatorial
    * Algorithms - Generation, Enumeration, and Search"</i> by Donald L. Kreher
    * and Douglas Stinson.
    * </p>
    *
    * \param[in]  n      Cardinal of the base set.
    * \param[in]  k      Cardinal of the subset.
    * \param[in]  r      Rank of the subset (assuming minimal change order).
    * \param[out] subset Subset to be returned.
    */
void unrank_subset_revdoor(uint32_t n, uint32_t k, uint32_t r,
                           uint32_t* subset);

/*
 *-----------------------------------------------------------------------------
 *                      Random number generation functions
 *-----------------------------------------------------------------------------
 */

   /**
    * \brief Initializes TIFA's basic pseudo-random generator.
    * <p>
    *   Initializes TIFA's basic random number generator with a user
    *   defined seed.
    * </p>
    * \param[in] seed The seed as a \c uint32_t.
    */
void tifa_srand(uint32_t seed);

   /**
    * \brief Returns a pseudo-random integer.
    * <p>
    *   Returns a pseudo-random integer using TIFA's basic random number
    *   generator.
    * </p>
    * \param[in] seed The seed as a \c uint32_t.
    */
uint32_t tifa_rand();

#ifdef __cplusplus
}
#endif

#endif
