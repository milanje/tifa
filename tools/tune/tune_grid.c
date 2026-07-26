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
 * \file    tune_grid.c
 */

#include <math.h>

#include "tune_grid.h"

//------------------------------------------------------------------------------
void grid_clear(grid_t* g) {
    g->naxes = 0;
}

//
// Appends v to the axis if it is within [0, TUNE_MAX_AXIS_VALUES) and not
// already present (values are added in nondecreasing order so a duplicate can
// only be the previous entry).
//
static void axis_push_unique(grid_axis_t* axis, uint32_t v) {
    if (axis->nvalues == TUNE_MAX_AXIS_VALUES) {
        return;
    }
    if ((axis->nvalues > 0) && (axis->values[axis->nvalues - 1] == v)) {
        return;
    }
    axis->values[axis->nvalues++] = v;
}

//------------------------------------------------------------------------------
void grid_add_axis_mul(grid_t* g, double center,
                       const double* ratios, uint32_t nratios,
                       uint32_t minval, uint32_t maxval) {
    if (g->naxes == TIFA_TUNING_MAX_COLS) {
        return;
    }
    grid_axis_t* axis = &g->axes[g->naxes];
    axis->nvalues = 0;

    for (uint32_t i = 0; i < nratios; i++) {
        double d = center * ratios[i];
        if (d < 0.0) {
            d = 0.0;
        }
        double r = floor(d + 0.5);
        uint32_t v;
        if (r >= (double) maxval) {
            v = maxval;
        } else {
            v = (uint32_t) r;
        }
        if (v < minval) {
            v = minval;
        }
        axis_push_unique(axis, v);
    }
    if (axis->nvalues > 0) {
        g->naxes++;
    }
}

//------------------------------------------------------------------------------
void grid_add_axis_fixed(grid_t* g, uint32_t value) {
    if (g->naxes == TIFA_TUNING_MAX_COLS) {
        return;
    }
    grid_axis_t* axis = &g->axes[g->naxes];
    axis->nvalues  = 1;
    axis->values[0] = value;
    g->naxes++;
}

//------------------------------------------------------------------------------
void grid_add_axis_list(grid_t* g, const uint32_t* values, uint32_t nvalues) {
    if ((g->naxes == TIFA_TUNING_MAX_COLS) || (nvalues == 0)) {
        return;
    }
    grid_axis_t* axis = &g->axes[g->naxes];
    axis->nvalues = 0;
    for (uint32_t i = 0; (i < nvalues) && (i < TUNE_MAX_AXIS_VALUES); i++) {
        axis->values[axis->nvalues++] = values[i];
    }
    g->naxes++;
}

//------------------------------------------------------------------------------
uint64_t grid_size(const grid_t* g) {
    uint64_t n = 1;
    for (uint32_t i = 0; i < g->naxes; i++) {
        n *= (g->axes[i].nvalues == 0) ? 1 : g->axes[i].nvalues;
    }
    return n;
}

//------------------------------------------------------------------------------
void grid_point(const grid_t* g, uint64_t idx, uint32_t out_cols[]) {
    // Mixed-radix decomposition: axis 0 is the fastest-varying digit.
    for (uint32_t i = 0; i < g->naxes; i++) {
        uint32_t base = (g->axes[i].nvalues == 0) ? 1 : g->axes[i].nvalues;
        uint32_t digit = (uint32_t) (idx % base);
        out_cols[i] = g->axes[i].values[digit];
        idx /= base;
    }
}
