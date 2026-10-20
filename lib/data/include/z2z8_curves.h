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
 * \file    z2z8_curves.h
 *
 * \brief Elliptic curves with torsion Z/2 x Z/8 over Q, for ECM.
 *
 * Their order modulo every prime is divisible by 16 (12 for Suyama's
 * curves), hence 1-6% fewer curves to find a factor. They come from Edwards
 * curves x^2 + y^2 = 1 + d.x^2.y^2 with a point (x8, x8) of order 8,
 * d = (2.x8^2 - 1)/x8^4, x8 = (k^2 - 2k + 2)/(k^2 - 2) for rational k (d
 * square), on which a search found a point (x0, y0) of infinite order with
 * small coordinates (tools/z2z8 regenerates and checks this table). In
 * Montgomery form: a24 = 1/(1 - d) and starting point u0 = (1 + y0)/(1 - y0).
 */

#if !defined(_TIFA_Z2Z8_CURVES_H_)
   /**
    * \def _TIFA_Z2Z8_CURVES_H_
    * Standard include guard.
    */
#define _TIFA_Z2Z8_CURVES_H_

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

   /**
    * \def NZ2Z8_CURVES
    * Number of curves in \c z2z8_curves.
    */
#define NZ2Z8_CURVES 27

   /**
    * \var z2z8_curves
    * One curve per row: d (numerator, denominator), y0 (idem).
    */
extern const uint64_t z2z8_curves[NZ2Z8_CURVES][4];

#ifdef __cplusplus
}
#endif

#endif
