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
 * \file    macros.h
 * \author  Jerome Milan
 * \date    Tue Oct 11 2011
 * \version 2011-10-11
 *
 * \brief Various CPP macros.
 *
 * Defines some C preprocessor macros that should be kept internal to
 * the TIFA library to avoid poluting client code.
 */

#if !defined(_TIFA_MACROS_H_)
   /**
    * \def _TIFA_MACROS_H_
    * Standard include guard.
    */
#define _TIFA_MACROS_H_

#ifdef __cplusplus
extern "C" {
#endif

#if defined(__GMP_IMPL_H__)
    //
    // The following macros are defined in the gmp-impl.h header. Undefine
    // them to make sure we're always using TIFA's versions.
    //
    #undef MPN_NORMALIZE
    #undef SIZ
    #undef ABSIZ
    #undef PTR
    #undef ALLOC
    #undef MAX
    #undef MIN
    #undef ABS
#endif

    //
    // The following few macros are from the GMP library, in the gmp-impl.h
    // header. Since this header is not installed by "make install", chances
    // are that GMP users don't have access to it. Consequently, some macros
    // from gmp-impl.h are replicated here. These macros are only useful if
    // one wants to manipulate multi-precision integers with internal
    // functions from the mpn layer.
    //

   /**
    * \def MPN_NORMALIZE(dest, nlimbs)
    * Macro from the GMP library: Computes the effective size of an MPN number.
    *
    * Given \c dest, a pointer to an array of \c nlimbs \c mp_limbs_t
    * integers giving the representation of a multi-precision integer n,
    * computes the absolute value of the effective size of n, i.e the number
    * of significant \c mp_size_t integers needed to represent n and modifies 
    * the value of \c nlimbs accordingly.
    *
    * \note This macro is originally the MPN_NORMALIZE macro from the GMP
    * library. It has been slightly modified.
    */
#define MPN_NORMALIZE(dest, nlimbs)                          \
    do {                                                         \
        while (((nlimbs) > 0) && ((dest)[(nlimbs) - 1] == 0)) {  \
            (nlimbs)--;                                          \
        }                                                        \
    } while (0)

   /**
    * \def SIZ(x)
    * Macro from the GMP library: Returns the \c _mp_size field of an \c mpz_t
    * integer.
    *
    * Returns the \c _mp_size field of the variable \c x of type <tt>mpz_t</tt>,
    * that is to say the number of \c mp_limbs_t integers needed to represent
    * the value of <tt>x</tt>. The sign of the returned value is given by the
    * sign of <tt>x</tt>'s value.
    *
    * \note This macro is the SIZ macro from the GMP library. It is
    * redistributed under the GNU LGPL license.
    */
#if !defined(SIZ)
    #define SIZ(x) ((x)->_mp_size)
#endif

   /**
    * \def ABSIZ(x)
    * Macro from the GMP library: Returns the absolute value of
    * <tt>SIZ(x)</tt>.
    *
    * Returns the absolute value of <tt>SIZ(x)</tt> that is to say the number
    * of \c mp_limbs_t integers needed to represent the value of <tt>x</tt>.
    *
    * \note This macro is the ABSIZ macro from the GMP library. It is
    * redistributed under the GNU LGPL license.
    */
#if !defined(ABSIZ)
    #define ABSIZ(x) (ABS(SIZ(x)))
#endif

   /**
    * \def MPZ_TO_ABS(x)
    * Sets the \c mpz_t \c x to its absolute value.
    */
#if !defined(MPZ_TO_ABS)
    #define MPZ_TO_ABS(x) (SIZ(x) = ABSIZ(x))
#endif

   /**
    * \def PTR(x)
    * Macro from the GMP library: Returns the \c _mp_d field of an \c mpz_t
    * integer.
    *
    * Returns the \c _mp_d field of an \c mpz_t integer, that is to say a
    * pointer to an array of \c mp_limbs_t integers giving the representation
    * of the value of <tt>x</tt>.
    *
    * \note This macro is the PTR macro from the GMP library. It is
    * redistributed under the GNU LGPL license.
    */
#if !defined(PTR)
    #define PTR(x) ((x)->_mp_d)
#endif

   /**
    * \def ALLOC(x)
    * Macro from the GMP library: Returns the \c _mp_alloc field of an \c mpz_t
    * integer.
    *
    * Returns the \c _mp_alloc field of an \c mpz_t integer, that is to say
    * the size (in units of \c mp_limb_t) of the \c x->_mp_d array.
    *
    * \note This macro is the ALLOC macro from the GMP library. It is
    * redistributed under the GNU LGPL license.
    */
#if !defined(ALLOC)
    #define ALLOC(x) ((x)->_mp_alloc)
#endif

   /**
    * \def MPZ_LIMB_VALUE(x, i)
    * Returns the value of the <tt>i</tt>-th limb of the \c mpz_t integer
    * <tt>x</tt>. The least significant limb is given by <tt>i</tt> = 0.
    */
#define MPZ_LIMB_VALUE(x, i)   ( PTR(x)[(i)] & GMP_NUMB_MASK )

   /**
    * \def MPZ_LAST_LIMB_VALUE(x)
    * Returns the value of the most significant limb of the \c mpz_t integer
    * <tt>x</tt>. Is equivalent to <tt>MPZ_LIMB_VALUE(x, SIZ(x) - 1)</tt>.
    */
#define MPZ_LAST_LIMB_VALUE(x) ( PTR(x)[SIZ(x) - 1] & GMP_NUMB_MASK )

   /**
    * \def MAX(a, b)
    * Standard macro returning the maximum of a and b.
    *
    * \note As usual, be careful of possible side effects when using this kind
    * of macro. The standard disclaimers apply.
    *
    */
#if !defined(MAX)
    #define MAX(a,b) ( ((a) > (b)) ? (a) : (b) )
#endif

   /**
    * \def MIN(a, b)
    * Standard macro returning the minimum of a and b.
    *
    * \note As usual, be careful of possible side effects when using this kind
    * of macro. The standard disclaimers apply.
    *
    */
#if !defined(MIN)
  #define MIN(a,b) ( ((a) < (b)) ? (a) : (b) )
#endif

   /**
    * \def ABS(a)
    * Standard macro returning the absolute value of a.
    *
    * \note As usual, be careful of possible side effects when using this
    * kind of macro. The standard disclaimers apply.
    *
    */
#if !defined(ABS)
  #define ABS(a) ( ((a) < 0) ? (-(a)) : (a) )
#endif

   /**
    * \def IS_POWER_OF_2_UI(ui)
    * Macro returning a non-zero value if the unsigned integer \c ui is a
    * power of 2.
    *
    * \note As usual, be careful of possible side effects when using this
    * kind of macro. The standard disclaimers apply.
    *
    */
#define IS_POWER_OF_2_UI(ui) ( ((ui) & ((ui) - 1)) == 0 )

   /**
    * \def IS_EVEN(ui)
    * Macro returning True if the unsigned integer \c ui is even,
    * False otherwise.
    */
#define IS_EVEN(ui) (((ui) & 1) == 0)

   /**
    * \def IS_ODD(ui)
    * Macro returning True if the unsigned integer \c ui is odd,
    * False otherwise.
    */
#define IS_ODD(ui) (((ui) & 1) != 0)

   /**
    * \def ARE_EVEN(uia, uib)
    * Macro returning True if both of the unsigned integers \c uia and \c uib
    * are even, False otherwise.
    */
#define ARE_EVEN(uia, uib) ((((uia) | (uib)) & 1) == 0)

   /**
    * \def ARE_ODD(uia, uib)
    * Macro returning True if both of the unsigned integers \c uia and \c uib
    * are odd, False otherwise.
    */
#define ARE_ODD(uia, uib) ((((uia) | (uib)) & 1) != 0)

   /**
    * \def BIT(N, i)
    * Macro returning the value of the i-th least significant bit of N.
    * BIT(N, 0) returns the least significant bit of N.
    */
#define BIT(N, i) ( ((N) & (1<<(i))) ? 1 : 0 )

   /**
    * \def DUFF_DEVICE(COUNT, STATEMENTS)
    *
    * Implements the so-called "Duff device" which is a fairly well-known
    * (and so ugly looking!) loop unrolling technique. \c COUNT is the number
    * of times to perform the operations given by STATEMENTS.
    *
    * \warning \c COUNT should be strictly positive. Using \c COUNT equals to
    * zero will yield wrong results.
    *
    * \note This macro was actually inspired (borrowed? stolen?) from the
    * example given in the article "A Reusable Duff Device" written by Ralf 
    * Holly.
    *
    * \see Tom Duff's comments about this technique at
    * http://www.lysator.liu.se/c/duffs-device.html (URL last accessed on
    * Thu 17 Feb 2011)
    *
    * \see "A Reusable Duff Device", Ralf Holly, <i>Dr. Dobb's Journal</i>,
    * August 2005. Available online at:
    *    http://www.drdobbs.com/high-performance-computing/184406208
    * (URL last accessed on Thu 17 Feb 2011)
    *
    */
#define DUFF_DEVICE(COUNT, STATEMENT, ...)          \
    do {                                            \
        long int __count__ = (COUNT);               \
        long int __niter__ = (__count__ + 7) >> 3;  \
        switch (__count__ & 7) {                    \
            case 0: do { STATEMENT; __VA_ARGS__;    \
            case 7:      STATEMENT; __VA_ARGS__;    \
            case 6:      STATEMENT; __VA_ARGS__;    \
            case 5:      STATEMENT; __VA_ARGS__;    \
            case 4:      STATEMENT; __VA_ARGS__;    \
            case 3:      STATEMENT; __VA_ARGS__;    \
            case 2:      STATEMENT; __VA_ARGS__;    \
            case 1:      STATEMENT; __VA_ARGS__;    \
                         __niter__--;               \
                    } while (__niter__ > 0);        \
        }                                           \
    } while (0)

   /**
    * \def MPZ_IS_SQUARE(X)
    *
    * Syntaxic sugar macro wrapping a call to <tt>mpz_perfect_square_p </tt>.
    * "Returns" true if and only if the \c mpz_t X is a perfect square.
    *
    * Takes as parameter an <tt>mpz_t</tt>.
    */
#define MPZ_IS_SQUARE(X) (0 != mpz_perfect_square_p(X))

   /**
    * \def NMILLER_RABIN
    *
    * Number of Miller-Rabin iterations to perform for each compositeness test.
    */
#define NMILLER_RABIN    32

   /**
    * \def MPZ_IS_PRIME(X)
    *
    * Syntaxic sugar macro wrapping a call to <tt>mpz_probab_prime_p</tt>.
    * "Returns" true if and only if the \c mpz_t \c X is (probably) prime.
    *
    * Takes as parameter an <tt>mpz_t</tt>.
    */
#define MPZ_IS_PRIME(X) (0 != mpz_probab_prime_p((X), NMILLER_RABIN))

   /**
    * \def MPN_ADD(A, B, C)
    *
    * Syntaxic sugar macro wrapping a call to <tt>mpn_add</tt>. Performs
    * size normalization on the result and takes care of the possible carry
    * out. Does not perform any reallocation: the user should make sure
    * the result has enough space to accomodate the possible carry out.
    *
    * Takes as parameters three \c mpz_t (and not arrays of <tt>mp_limb_t</tt>).
    * \c SIZ(B) should be greater than or equal to <tt>SIZ(C)</tt>.
    *
    * \warning \c B and \c C should be both positive or the result will be
    * unpredictable.
    *
    * \see The GMP documentation for more information on the \c mpn_add
    * function.
    */
#define MPN_ADD(A, B, C)                                            \
    do {                                                            \
       if (mpn_add(PTR(A), PTR(B), SIZ(B), PTR(C), SIZ(C))) {       \
           SIZ(A) = SIZ(B);                                         \
           MPN_NORMALIZE(PTR(A), SIZ(A));                           \
           PTR(A)[SIZ(A)] = 1;                                      \
           SIZ(A)++;                                                \
       } else {                                                     \
           SIZ(A) = SIZ(B);                                         \
           MPN_NORMALIZE(PTR(A), SIZ(A));                           \
       }                                                            \
    } while (0)

   /**
    * \def MPN_ADD_CS(A, B, C)
    *
    * Syntaxic sugar macro wrapping a call to <tt>mpn_add</tt>. Performs
    * size normalization on the result, takes care of the possible carry
    * out, and Checks the Sizes of the operands to call \c mpn_add with the
    * proper parameters' order. However, does not perform any reallocation: the
    * user should make sure the result has enough space to accomodate the
    * possible carry out.
    *
    * Takes as parameters three \c mpz_t (and not arrays of <tt>mp_limb_t</tt>).
    * \c B and \c C can be used interchangeably.
    *
    * \warning \c B and \c C should be both positive or the result will be
    * unpredictable.
    *
    * \see The GMP documentation for more information on the \c mpn_add
    * function.
    */
#define MPN_ADD_CS(A, B, C)                                         \
    do {                                                            \
       if (SIZ(B) > SIZ(C)) {                                       \
           MPN_ADD(A, B, C);                                        \
       } else {                                                     \
           MPN_ADD(A, C, B);                                        \
       }                                                            \
    } while (0)

   /**
    * \def MPN_SUB(A, B, C)
    *
    * Syntaxic sugar macro wrapping a call to <tt>mpn_sub</tt>. Performs
    * size normalization on the result but does not take care of the possible
    * borrow out.
    *
    * Takes as parameters three \c mpz_t (and not arrays of <tt>mp_limb_t</tt>).
    * \c B should be greater than or equal to <tt>C</tt>.
    *
    * \warning \c B and \c C should be both positive or the result will be
    * unpredictable.
    *
    * \see The GMP documentation for more information on the \c mpn_sub
    * function.
    */
#define MPN_SUB(A, B, C)                                            \
    do {                                                            \
       mpn_sub(PTR(A), PTR(B), SIZ(B), PTR(C), SIZ(C));             \
       SIZ(A) = SIZ(B);                                             \
       MPN_NORMALIZE(PTR(A), SIZ(A));                               \
    } while (0)

   /**
    * \def MPN_SUB_N(A, B, C)
    *
    * Syntaxic sugar macro wrapping a call to <tt>mpn_sub_n</tt>. Performs
    * size normalization on the result but does not take care of the possible
    * borrow out.
    *
    * Takes as parameters three \c mpz_t (and not arrays of <tt>mp_limb_t</tt>).
    * \c B should be greater than or equal to <tt>C</tt>.
    *
    * \warning \c B and \c C should be both positive or the result will be
    * unpredictable.
    *
    * \see The GMP documentation for more information on the \c mpn_sub_n
    * function.
    */
#define MPN_SUB_N(A, B, C)                                          \
    do {                                                            \
       mpn_sub_n(PTR(A), PTR(B), PTR(C), SIZ(B));                   \
       SIZ(A) = SIZ(B);                                             \
       MPN_NORMALIZE(PTR(A), SIZ(A));                               \
    } while (0)

   /**
    * \def MPN_TDIV_QR(Q, R, N, D)
    *
    * Syntaxic sugar macro wrapping a call to <tt>mpn_tdiv_qr</tt>. Performs
    * size normalization on both the quotient and remainder.
    *
    * Takes as parameters four \c mpz_t (and not arrays of <tt>mp_limb_t</tt>).
    *
    * \warning \c N and \c D should be both positive or the result will be
    * unpredictable.
    *
    * \see The GMP documentation for more information on the \c mpn_tdiv_qr
    * function.
    */
#define MPN_TDIV_QR(Q, R, N, D)                                              \
    do {                                                                     \
        if (SIZ(N) >= SIZ(D)) {                                              \
            mpn_tdiv_qr(PTR(Q), PTR(R), 0,  PTR(N), SIZ(N), PTR(D), SIZ(D)); \
            SIZ(Q) = SIZ(N) - SIZ(D) + 1;                                    \
            MPN_NORMALIZE(PTR(Q), SIZ(Q));                                   \
            SIZ(R) = SIZ(D);                                                 \
            MPN_NORMALIZE(PTR(R), SIZ(R));                                   \
        }                                                                    \
    } while (0)

   /**
    * \def MPN_MUL(A, B, C)
    *
    * Syntaxic sugar macro wrapping a call to <tt>mpn_mul</tt>. Performs
    * size normalization on the result.
    *
    * Takes as parameters three \c mpz_t (and not arrays of <tt>mp_limb_t</tt>).
    * \c SIZ(B) should be greater than or equal to <tt>SIZ(C)</tt>.
    *
    * \warning \c B and \c C should be both positive or the result will be
    * unpredictable.
    *
    * \see The GMP documentation for more information on the \c mpn_mul
    * function.
    */
#define MPN_MUL(A, B, C)                                            \
    do {                                                            \
       mpn_mul(PTR(A), PTR(B), SIZ(B), PTR(C), SIZ(C));             \
       SIZ(A) = SIZ(B) + SIZ(C);                                    \
       MPN_NORMALIZE(PTR(A), SIZ(A));                               \
    } while (0)

   /**
    * \def MPN_MUL_N(A, B, C)
    *
    * Syntaxic sugar macro wrapping a call to <tt>mpn_mul_n</tt>. Performs
    * size normalization on the result.
    *
    * Takes as parameters three \c mpz_t (and not arrays of <tt>mp_limb_t</tt>).
    * \c SIZ(B) and \c SIZ(C) should be the same.
    *
    * \warning \c B and \c C should be both positive or the result will be
    * unpredictable.
    *
    * \see The GMP documentation for more information on the \c mpn_mul_n
    * function.
    */
#define MPN_MUL_N(A, B, C)                                          \
    do {                                                            \
       mpn_mul_n(PTR(A), PTR(B), PTR(C), SIZ(B));                   \
       SIZ(A) = SIZ(B) << 1;                                        \
       MPN_NORMALIZE(PTR(A), SIZ(A));                               \
    } while (0)

   /**
    * \def MPN_MUL_CS(A, B, C)
    *
    * Syntaxic sugar macro wrapping a call to <tt>mpn_mul</tt>. Performs
    * size normalization on the result, and Checks the Sizes of the operands to
    * call \c mpn_mul with the proper parameters' order.
    *
    * Takes as parameters three \c mpz_t (and not arrays of <tt>mp_limb_t</tt>).
    * \c B and \c C can be used interchangeably.
    *
    * \warning \c B and \c C should be both positive or the result will be
    * unpredictable.
    *
    * \see The GMP documentation for more information on the \c mpn_mul
    * function.
    */
#define MPN_MUL_CS(A, B, C)                                         \
    do {                                                            \
       if (SIZ(B) > SIZ(C)) {                                       \
           mpn_mul(PTR(A), PTR(B), SIZ(B), PTR(C), SIZ(C));         \
       } else {                                                     \
           mpn_mul(PTR(A), PTR(C), SIZ(C), PTR(B), SIZ(B));         \
       }                                                            \
       SIZ(A) = SIZ(B) + SIZ(C);                                    \
       MPN_NORMALIZE(PTR(A), SIZ(A));                               \
    } while (0)

   /**
    * \def MPN_MUL_CS_S(A, B, C)
    *
    * Syntaxic sugar macro wrapping a call to <tt>mpn_mul</tt>. Performs
    * size normalization on the result, and Checks the Sizes and the Signs
    * of the operands to call \c mpn_mul with the proper parameters' order.
    *
    * Takes as parameters three \c mpz_t (and not arrays of <tt>mp_limb_t</tt>).
    * \c B and \c C can be used interchangeably.
    *
    * \note \c B and \c C are allowed to be negative.
    *
    * \see The GMP documentation for more information on the \c mpn_mul
    * function.
    */
#define MPN_MUL_CS_S(A, B, C)                                       \
    do {                                                            \
       if (ABSIZ(B) > ABSIZ(C)) {                                   \
           mpn_mul(PTR(A), PTR(B), ABSIZ(B), PTR(C), ABSIZ(C));     \
       } else {                                                     \
           mpn_mul(PTR(A), PTR(C), ABSIZ(C), PTR(B), ABSIZ(B));     \
       }                                                            \
       SIZ(A) = ABSIZ(B) + ABSIZ(C);                                \
       MPN_NORMALIZE(PTR(A), SIZ(A));                               \
       if ((SIZ(B) ^ SIZ(C)) < 0) {                                 \
            SIZ(A) = -SIZ(A);                                       \
       }                                                            \
    } while (0)

   /**
    * \def DECLARE_MPZ_SWAP_VARS
    *
    * Macro declaring local variables needed by the \c MPZ_SWAP macro.
    * Should be called \e once prior to any use of the \c MPZ_SWAP macro.
    *
    * \warning
    * Declares the variables \c __TMPPTR__MACROS_H__a9b3c01__ and
    * \c __TMPSIZ__MACROS_H__a9b3c01__ . Hopefully their names are fancy
    * enough to avoid any local conflict.
    */
#define DECLARE_MPZ_SWAP_VARS                                       \
    mp_ptr __TMPPTR__MACROS_H__a9b3c01__;                           \
    mp_size_t __TMPSIZ__MACROS_H__a9b3c01__;

   /**
    * \def MPZ_SWAP(A, B)
    *
    * Macro swapping the values of the two \c mpz_t \c A and <tt>B</tt>.
    *
    * \warning
    * The macro \c DECLARE_MPZ_SWAP_VARS should be called \e once before
    * using <tt>MPZ_SWAP</tt>.
    */
#define MPZ_SWAP(A, B)                                              \
    do {                                                            \
        __TMPPTR__MACROS_H__a9b3c01__ = PTR(A);                     \
        __TMPSIZ__MACROS_H__a9b3c01__ = SIZ(A);                     \
        PTR(A) = PTR(B);                                            \
        SIZ(A) = SIZ(B);                                            \
        PTR(B) = __TMPPTR__MACROS_H__a9b3c01__;                     \
        SIZ(B) = __TMPSIZ__MACROS_H__a9b3c01__;                     \
    } while (0)

   /**
    * \def TIFA_DEBUG_MSG(...)
    *
    * Macro printing debug message with filename and line number.
    *
    * \warning
    * The symbol \c ENABLE_TIFA_DEBUG_MSG should be defined to non zero
    * \e before including this file.
    */
#if (defined(ENABLE_TIFA_DEBUG_MSG) && (ENABLE_TIFA_DEBUG_MSG != 0))
    #define TIFA_DEBUG_MSG(...)                                     \
        do {                                                        \
            printf("\nDEBUG (%s:%i) ", __FILE__, __LINE__);         \
            printf(__VA_ARGS__);fflush(stdout);                     \
          } while (0)
#else
    #define TIFA_DEBUG_MSG(...) /* Intentionally left empty */
#endif

#ifdef __cplusplus
}
#endif

#endif
