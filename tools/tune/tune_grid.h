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
 * \file    tune_grid.h
 * \brief   A small cartesian-product grid over integer parameter axes.
 *
 * A \c grid_t holds one \c grid_axis_t per tunable parameter column, each axis
 * being an explicit list of candidate values. The grid enumerates the cartesian
 * product of the axes by linear index, so callers can iterate every parameter
 * combination without maintaining an odometer.
 */

#if !defined(_TIFA_TUNE_GRID_H_)
#define _TIFA_TUNE_GRID_H_

#include <inttypes.h>

#include "tuning.h"   // TIFA_TUNING_MAX_COLS

   /**
    * \def TUNE_MAX_AXIS_VALUES
    * Maximum number of candidate values on a single axis.
    */
#define TUNE_MAX_AXIS_VALUES 16

typedef struct {
    uint32_t nvalues;
    uint32_t values[TUNE_MAX_AXIS_VALUES];
} grid_axis_t;

typedef struct {
    uint32_t    naxes;
    grid_axis_t axes[TIFA_TUNING_MAX_COLS];
} grid_t;

   /**
    * \brief Resets a grid to zero axes.
    */
void grid_clear(grid_t* g);

   /**
    * \brief Appends an axis built from a center value and multiplicative ratios.
    *
    * Each candidate value is <tt>round(center * ratio)</tt>, clamped to the
    * inclusive range <tt>[minval, maxval]</tt>; adjacent duplicates are
    * dropped so repeated clamped values do not inflate the grid. Adds nothing (and does not fail)
    * if the grid already has \c TIFA_TUNING_MAX_COLS axes.
    *
    * \param[in,out] g       The grid to extend.
    * \param[in]     center  The center value.
    * \param[in]     ratios  Array of multiplicative ratios.
    * \param[in]     nratios Number of ratios.
    * \param[in]     minval  Lower clamp (inclusive).
    * \param[in]     maxval  Upper clamp (inclusive).
    */
void grid_add_axis_mul(grid_t* g, double center,
                       const double* ratios, uint32_t nratios,
                       uint32_t minval, uint32_t maxval);

   /**
    * \brief Appends an axis holding a single fixed value.
    */
void grid_add_axis_fixed(grid_t* g, uint32_t value);

   /**
    * \brief Appends an axis holding an explicit list of values.
    *
    * At most \c TUNE_MAX_AXIS_VALUES values are stored. Adds nothing if the
    * grid is already full or the list is empty.
    */
void grid_add_axis_list(grid_t* g, const uint32_t* values, uint32_t nvalues);

   /**
    * \brief Returns the number of points in the grid (product of axis sizes).
    */
uint64_t grid_size(const grid_t* g);

   /**
    * \brief Writes the grid point at linear index \c idx into \c out_cols.
    *
    * \param[in]  g        The grid.
    * \param[in]  idx      A point index in <tt>[0, grid_size(g))</tt>.
    * \param[out] out_cols Receives one value per axis (length \c g->naxes).
    */
void grid_point(const grid_t* g, uint64_t idx, uint32_t out_cols[]);

#endif
