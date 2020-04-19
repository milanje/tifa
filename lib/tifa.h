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
 * \file    tifa.h
 * \author  Jerome Milan
 * \date    Tue Oct 11 2011
 * \version 2011-10-11
 *
 * \brief Library wide public include file.
 *
 * Includes only TIFA's structures and functions needed from client code
 * perspective.
 */

#if !defined(_TIFA_TIFA_H_)
   /**
    * \def _TIFA_TIFA_H_
    * Standard include guard.
    */
#define _TIFA_TIFA_H_

//
// The configuration file is not strictly needed but nice to have...
//
#include "tifa_config.h"
//
// Includes from lib/algo
//
#include "cfrac.h"
#include "ecm.h"
#include "fermat.h"
#include "qs.h"
#include "siqs.h"
#include "squfof.h"
#include "tdiv.h"
#include "tifa_factor.h"
//
// Includes from lib/data
//
#include "first_primes.h"
//
// Includes from lib/utils
//
#include "array.h"
#include "exit_codes.h"
#include "factoring_machine.h"

#endif
