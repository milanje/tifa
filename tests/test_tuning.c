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

// Unit tests for the runtime tuning-file loader (lib/algo/include/tuning.h,
// lib/algo/src/tuning.c).

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tifa_config.h"
#include "tifa_test.h"
#include "tuning.h"
#include "factoring_machine.h"
#include "p_minus_1.h"

#if !defined(TIFA_TEST_FIXTURE_DIR)
    #error "TIFA_TEST_FIXTURE_DIR must be defined by the build"
#endif

// When the baked tables are compiled in, "no runtime file" no longer means
// "hardcoded defaults": the baked middle layer supplies values. The assertions
// that check the hardcoded default therefore only apply when baking is off.
#if defined(TIFA_USE_BAKED_TUNING) && TIFA_USE_BAKED_TUNING
    #define TUNING_BAKED_IN 1
#else
    #define TUNING_BAKED_IN 0
#endif

static void make_path(char* buf, size_t buflen, const char* name) {
    snprintf(buf, buflen, "%s/%s", TIFA_TEST_FIXTURE_DIR, name);
}

int main(void) {
    TIFA_TEST_BEGIN("tuning loader");

    // The resolver falls back to it once the override is cleared.
    unsetenv("TIFA_TUNING_DIR");

    char path[1024];

    // --- Load a well-formed file. ---------------------------------------
    make_path(path, sizeof(path), "p_minus_1.tuning");
    tuning_table_t* table = NULL;
    int rc = tifa_tuning_load_file(path, TIFA_TUNING_NCOLS_P_MINUS_1, &table);
    CHECK(rc == 0, "well-formed file should load");
    CHECK(table != NULL, "loaded table should be non-NULL");

    if (table != NULL) {
        CHECK_EQ_UI(table->ncols, TIFA_TUNING_NCOLS_P_MINUS_1, "ncols");
        CHECK_EQ_UI(table->nrows, 3, "nrows");

        // --- Lookup semantics. ------------------------------------------
        const tuning_row_t* r;

        r = tifa_tuning_lookup(table, 30);   // below the first row -> first
        CHECK(r != NULL, "lookup(30) non-NULL");
        if (r != NULL) {
            CHECK_EQ_UI(r->min_bitsize, 40, "lookup(30) clamps to first row");
            CHECK_EQ_UI(r->cols[0], 2000, "lookup(30) b1");
        }

        r = tifa_tuning_lookup(table, 50);   // between rows -> lower row
        CHECK(r != NULL, "lookup(50) non-NULL");
        if (r != NULL) {
            CHECK_EQ_UI(r->min_bitsize, 40, "lookup(50) picks lower row");
        }

        r = tifa_tuning_lookup(table, 56);   // exact boundary -> that row
        CHECK(r != NULL, "lookup(56) non-NULL");
        if (r != NULL) {
            CHECK_EQ_UI(r->min_bitsize, 56, "lookup(56) exact");
            CHECK_EQ_UI(r->cols[1], 240000, "lookup(56) b2");
        }

        r = tifa_tuning_lookup(table, 100);  // above the last row -> last
        CHECK(r != NULL, "lookup(100) non-NULL");
        if (r != NULL) {
            CHECK_EQ_UI(r->min_bitsize, 72, "lookup(100) clamps to last row");
            CHECK_EQ_UI(r->cols[2], 2, "lookup(100) base");
        }

        // --- Interpolation between rows. --------------------------------
        tuning_row_t row;
        CHECK(tifa_tuning_interpolate(table, 60, TIFA_TUNING_INTERP_P_MINUS_1,
                                      &row), "interpolate(60) found");
        CHECK_EQ_UI(row.cols[0], 24000, "interpolate(60) b1, a quarter on");
        CHECK_EQ_UI(row.cols[1], 480000, "interpolate(60) b2");
        CHECK_EQ_UI(row.cols[2], 2, "interpolate(60) base");

        tifa_tuning_interpolate(table, 60, 0, &row);
        CHECK_EQ_UI(row.cols[0], 12000, "interpolate(60) without a mask");

        tifa_tuning_interpolate(table, 30, TIFA_TUNING_INTERP_P_MINUS_1, &row);
        CHECK_EQ_UI(row.cols[0], 2000, "interpolate(30) clamps to first row");
        tifa_tuning_interpolate(table, 56, TIFA_TUNING_INTERP_P_MINUS_1, &row);
        CHECK_EQ_UI(row.cols[0], 12000, "interpolate(56) exact");
        tifa_tuning_interpolate(table, 100, TIFA_TUNING_INTERP_P_MINUS_1, &row);
        CHECK_EQ_UI(row.cols[0], 60000, "interpolate(100) clamps to last row");
    }
    tifa_tuning_free(table);

    // --- A malformed file must be rejected without crashing. ------------
    make_path(path, sizeof(path), "bad.tuning");
    tuning_table_t* bad = NULL;
    rc = tifa_tuning_load_file(path, TIFA_TUNING_NCOLS_P_MINUS_1, &bad);
    CHECK(rc != 0, "malformed file should fail to load");
    CHECK(bad == NULL, "malformed file yields NULL table");

    // --- A missing file must fail cleanly. ------------------------------
    make_path(path, sizeof(path), "does_not_exist.tuning");
    tuning_table_t* missing = NULL;
    rc = tifa_tuning_load_file(path, TIFA_TUNING_NCOLS_P_MINUS_1, &missing);
    CHECK(rc != 0, "missing file should fail to load");
    CHECK(missing == NULL, "missing file yields NULL table");

    // --- The directory-based resolver. ----------------------------------
    tifa_tuning_set_dir(TIFA_TEST_FIXTURE_DIR);

    const tuning_row_t* row =
        tifa_tuning_row_for("p_minus_1", TIFA_TUNING_NCOLS_P_MINUS_1, 60);
    CHECK(row != NULL, "resolver finds p_minus_1 at bitsize 60");
    if (row != NULL) {
        CHECK_EQ_UI(row->min_bitsize, 56, "resolver picks row 56 for bitsize 60");
    }

    // An algorithm with no fixture file resolves to NULL (defaults kept).
    const tuning_row_t* none =
        tifa_tuning_row_for("nonexistent_algo", TIFA_TUNING_NCOLS_P_MINUS_1, 60);
    CHECK(none == NULL, "resolver returns NULL for a missing algorithm file");

    // Clearing the override must make the resolver miss again (unless baked
    // tables are compiled in, in which case the baked layer answers, or the
    // compiled-in tuning directory has the file).
    tifa_tuning_set_dir(NULL);
#if !TUNING_BAKED_IN
    const tuning_row_t* cleared =
        tifa_tuning_row_for("p_minus_1", TIFA_TUNING_NCOLS_P_MINUS_1, 60);
    CHECK(cleared == NULL, "resolver misses after the directory is cleared");
#endif

    // --- Fallback proof via set_p_minus_1_params_to_default. ------------
    // A 60-bit number; only its bit-size matters here.
    mpz_t n;
    mpz_init(n);
    mpz_ui_pow_ui(n, 2, 59);       // bit-size 60
    CHECK_EQ_UI(mpz_sizeinbase(n, 2), 60, "n is 60-bit");

    p_minus_1_params_t params;

    // With the override cleared (and no p_minus_1.tuning in the compiled-in
    // tuning directory), the compiled-in defaults must be used.
    // (Only meaningful when the baked tables are not compiled in.)
#if !TUNING_BAKED_IN
    tifa_tuning_set_dir(NULL);
    set_p_minus_1_params_to_default(n, &params);
    CHECK_EQ_UI(params.b1, 50000, "default b1 when no tuning file");
    CHECK_EQ_UI(params.base, 2, "default base when no tuning file");
#endif

    // With the fixture directory, rows 56 and 72 interpolated at bit-size 60
    // must win. A runtime tuning file takes precedence in every build mode.
    tifa_tuning_set_dir(TIFA_TEST_FIXTURE_DIR);
    set_p_minus_1_params_to_default(n, &params);
    CHECK_EQ_UI(params.b1, 24000, "tuned b1 overrides default");
    CHECK_EQ_UI(params.b2, 480000, "tuned b2 overrides default");
    CHECK_EQ_UI(params.base, 2, "tuned base");

    // Clearing the override restores the default behavior (same assumption).
#if !TUNING_BAKED_IN
    tifa_tuning_set_dir(NULL);
    set_p_minus_1_params_to_default(n, &params);
    CHECK_EQ_UI(params.b1, 50000, "default b1 restored after clearing dir");
#endif

    mpz_clear(n);

    TIFA_TEST_END();
}
