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
 * \file    methods.h
 * \brief   The methods the parameters of the algorithms choose from: linear
 *          algebra and smoothness tests.
 */

#if !defined(_TIFA_METHODS_H_)
   /**
    * \def _TIFA_METHODS_H_
    * Standard include guard.
    */
#define _TIFA_METHODS_H_

#ifdef __cplusplus
extern "C" {
#endif

   /**
    * \enum linalg_method_enum
    *
    * Enumeration listing the different linear system resolution method
    * implemented.
    *
    * For the time being, only one method is available.
    */
enum linalg_method_enum {
       /**
        * "Smart" gaussian elimination described in:
        * "A compact algorithm for Gaussian elimination over GF(2) implemented
        * on highly parallel computers", by D. Parkinson and M. Wunderlich
        * (Parallel Computing 1 (1984) 65-73).
        */
    SMART_GAUSS_ELIM = 0
};

   /**
    * \typedef linalg_method_t
    * \brief Equivalent to <tt>enum linalg_method_enum</tt>.
    */
typedef enum linalg_method_enum linalg_method_t;

   /**
    * Global constant array mapping linalg methods to their string 
    * representations.
    */
static const char* const linalg_method_to_str[1] = {
    "smart gaussian elimination"
};

   /**
    * \enum smooth_filter_method_enum
    *
    * An enumeration of the possible methods used to test residue smoothness.
    */
enum smooth_filter_method_enum {
        /**
         * Simple trial division.
         */
    TDIV = 0,
        /**
         * Simple trial division with (multi-step) early abort.
         */
    TDIV_EARLY_ABORT,
        /**
         * D. Bernstein's batch method described in "How to find smooth parts
         * of integers", http://cr.yp.to/factorization/smoothparts-20040510.pdf.
         */
    DJB_BATCH,
};

   /**
    * \typedef smooth_filter_method_t
    * \brief Equivalent to <tt>enum smooth_filter_method_enum</tt>.
    */
typedef enum smooth_filter_method_enum smooth_filter_method_t;

   /**
    * Global constant array mapping filter methods to their string 
    * representations.
    */
static const char* const filter_method_to_str[3] = {
    "trial division",
    "trial division + early abort",
    "batch",
};

#ifdef __cplusplus
}
#endif

#endif
