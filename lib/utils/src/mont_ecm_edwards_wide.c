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
 * \file    mont_ecm_edwards_wide.c
 *
 * \brief The Edwards stage 1 of mont_ecm_edwards.c on 5 to 8 limbs, checking
 *        the time limit and the interrupt now and then (run by
 *        mont_ecm_wide.c; the long curves on fewer limbs keep the ladder).
 */

#include "mont_ecm.h"

#if TIFA_HAVE_MONT_ECM

#define MONT_MOD_LIMBS MONT_ECM_MAX_LIMBS
#include "mont_ecm_internal.h"
#include "stopping.h"

#define MONT_ECM_STOP_CHECKS
#define MONT_ECM_EDWARDS
#define EDW_NAME stage1_edwards_long

#define NL 5
#include "mont_ecm_impl.h"
#undef NL
#define NL 6
#include "mont_ecm_impl.h"
#undef NL
#define NL 7
#include "mont_ecm_impl.h"
#undef NL
#define NL 8
#include "mont_ecm_impl.h"
#undef NL

#endif
