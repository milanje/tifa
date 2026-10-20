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
 * \file    tune.c
 * \brief   tifa_tune: measure and record good algorithm parameters per bit-size.
 *
 * For each selected algorithm and bit-size, tifa_tune generates a shared set of
 * test composites, sweeps the algorithm's parameter grid (evaluating every
 * combination on the same composites for a fair comparison), times each run,
 * and selects the parameter combination minimizing the expected time to find a
 * factor. The winning combination for each bit-size is written as one row of a
 * per-algorithm <tt>&lt;algo&gt;.tuning</tt> file that the library can load at
 * runtime (see lib/algo/include/tuning.h).
 */

// This tool relies on POSIX facilities (fork/pipe/select/waitpid/kill for the
// per-run watchdog, isatty/fileno for the progress bar, gettimeofday for wall
// time, strtok_r for argument parsing). Request them explicitly before any
// include, since the library is built in strict C99 mode where a conforming
// libc otherwise hides these declarations. _GNU_SOURCE is a no-op where they
// are already visible (e.g. macOS).
#if !defined(_GNU_SOURCE)
#define _GNU_SOURCE
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>
#include <stdbool.h>
#include <errno.h>
#include <signal.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/select.h>
#include <sys/wait.h>
#include <unistd.h>

#include <gmp.h>

#include "tifa_config.h"
#include "tuning.h"

#include "tune_algos.h"
#include "tune_grid.h"
#include "tune_numgen.h"

#define TUNE_EPS         1e-9
#define TUNE_MAX_ROWS    256
#define TUNE_MAX_SELECTED 32   // cap on algorithms named in a --algo list
// Consecutive sizes a "sticky-on" column must win before it is locked on, so a
// single noisy win cannot prematurely stop the off case from being tested.
#define TUNE_STICKY_HYSTERESIS 2

typedef struct {
    const char*   algo;         // "all" or a specific algorithm name
    unsigned int  min_bits;
    unsigned int  max_bits;
    unsigned int  step;
    unsigned int  trials;       // test composites per (bit-size, factor-count)
    double        threshold;    // required success rate for deterministic algos
    double        min_success;  // required success rate for probabilistic algos
    double        timeout_s;    // per-run wall-clock budget (0 disables watchdog)
    unsigned long seed;
    const char*   outdir;
    const char*   emit_c;       // if set, code-generate this .c from tuning files
    const char*   emit_h;       // matching .h path for --emit-c
    const char*   indir;        // input dir of tuning files for --emit-c
} options_t;

//
// Standard TIFA license header, written verbatim to generated source files.
//
static void write_license(FILE* f) {
    fprintf(f,
        "// Copyright (C) 2026 Jerome Milan.\n"
        "//\n"
        "// This file is part of TIFA.\n"
        "//\n"
        "// TIFA is free software; you can redistribute it and/or modify it under the\n"
        "// terms of the GNU Lesser General Public License as published by the Free\n"
        "// Software Foundation; either version 2.1 of the License, or (at your option)\n"
        "// any later version.\n"
        "//\n"
        "// TIFA is distributed in the hope that it will be useful, but WITHOUT ANY\n"
        "// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS\n"
        "// FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for\n"
        "// more details.\n"
        "//\n"
        "// You should have received a copy of the GNU Lesser General Public License\n"
        "// along with this library; if not, write to the Free Software Foundation,\n"
        "// Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.\n");
}

//
// Code-generates the baked tuning tables: reads <indir>/<algo>.tuning for every
// known algorithm and emits a committed C source (cpath) + header (hpath) that
// expose tifa_tuning_baked_row_for(). Mirrors the offline genprimes.pl model.
//
static int emit_baked(const char* cpath, const char* hpath, const char* indir) {
    FILE* h = fopen(hpath, "w");
    if (h == NULL) {
        fprintf(stderr, "tifa_tune: cannot open '%s' for writing\n", hpath);
        return 1;
    }
    write_license(h);
    fprintf(h,
        "\n// Automatically generated by tifa_tune --emit-c. Do not edit by hand.\n\n"
        "#if !defined(_TIFA_TUNING_BAKED_H_)\n"
        "#define _TIFA_TUNING_BAKED_H_\n\n"
        "#ifdef __cplusplus\nextern \"C\" {\n#endif\n\n"
        "#include <inttypes.h>\n\n"
        "#include \"tuning.h\"\n\n"
        "// Returns the compiled-in tuning row for an algorithm and bit-size, or\n"
        "// NULL if none applies. Used as a fallback by tifa_tuning_row_for().\n"
        "const tuning_row_t* tifa_tuning_baked_row_for(const char* algo,\n"
        "                                              uint32_t expect_ncols,\n"
        "                                              uint32_t bitsize);\n\n"
        "#ifdef __cplusplus\n}\n#endif\n\n#endif\n");
    fclose(h);

    FILE* c = fopen(cpath, "w");
    if (c == NULL) {
        fprintf(stderr, "tifa_tune: cannot open '%s' for writing\n", cpath);
        return 1;
    }
    write_license(c);
    fprintf(c,
        "\n// Automatically generated by tifa_tune --emit-c. Do not edit by hand.\n\n"
        "#include <string.h>\n#include <stdio.h>\n\n"
        "#include \"tuning_baked.h\"\n\n");

    const char* loaded_names[TUNE_MAX_ROWS];
    uint32_t    loaded_nrows[TUNE_MAX_ROWS];
    uint32_t    loaded_ncols[TUNE_MAX_ROWS];
    unsigned    nloaded = 0;

    for (unsigned int i = 0; i < TUNE_NALGOS; i++) {
        const tune_algo_desc_t* desc = &TUNE_ALGOS[i];
        char path[1024];
        snprintf(path, sizeof(path), "%s/%s.tuning", indir, desc->name);

        tuning_table_t* t = NULL;
        if (tifa_tuning_load_file(path, desc->ncols, &t) != 0) {
            continue;   // no file for this algorithm; skip
        }
        fprintf(c, "static const tuning_row_t %s_baked_rows[] = {\n", desc->name);
        for (uint32_t r = 0; r < t->nrows; r++) {
            fprintf(c, "    { %" PRIu32 ", {", t->rows[r].min_bitsize);
            for (uint32_t col = 0; col < desc->ncols; col++) {
                fprintf(c, " %" PRIu32 "%s", t->rows[r].cols[col],
                        (col + 1 < desc->ncols) ? "," : "");
            }
            fprintf(c, " } },\n");
        }
        fprintf(c, "};\n\n");

        loaded_names[nloaded] = desc->name;
        loaded_nrows[nloaded] = t->nrows;
        loaded_ncols[nloaded] = desc->ncols;
        nloaded++;
        tifa_tuning_free(t);
    }

    if (nloaded > 0) {
        fprintf(c,
            "typedef struct {\n"
            "    const char*         name;\n"
            "    const tuning_row_t* rows;\n"
            "    uint32_t            nrows;\n"
            "    uint32_t            ncols;\n"
            "} baked_entry_t;\n\n"
            "static const baked_entry_t baked_tables[] = {\n");
        for (unsigned int i = 0; i < nloaded; i++) {
            fprintf(c, "    { \"%s\", %s_baked_rows, %" PRIu32 ", %" PRIu32 " },\n",
                    loaded_names[i], loaded_names[i],
                    loaded_nrows[i], loaded_ncols[i]);
        }
        fprintf(c,
            "};\n\n"
            "static const unsigned int nbaked_tables =\n"
            "    sizeof(baked_tables) / sizeof(baked_tables[0]);\n\n");
    }

    fprintf(c,
        "const tuning_row_t* tifa_tuning_baked_row_for(const char* algo,\n"
        "                                              uint32_t expect_ncols,\n"
        "                                              uint32_t bitsize) {\n");
    if (nloaded == 0) {
        fprintf(c,
            "    (void) algo; (void) expect_ncols; (void) bitsize;\n"
            "    return NULL;\n");
    } else {
        fprintf(c,
            "    for (unsigned int i = 0; i < nbaked_tables; i++) {\n"
            "        if ((strcmp(baked_tables[i].name, algo) == 0)\n"
            "            && (baked_tables[i].ncols == expect_ncols)) {\n"
            "            tuning_table_t t;\n"
            "            memset(&t, 0, sizeof(t));\n"
            "            snprintf(t.algo, sizeof(t.algo), \"%%s\", baked_tables[i].name);\n"
            "            t.ncols = baked_tables[i].ncols;\n"
            "            t.nrows = baked_tables[i].nrows;\n"
            "            t.rows  = (tuning_row_t*) baked_tables[i].rows;\n"
            "            return tifa_tuning_lookup(&t, bitsize);\n"
            "        }\n"
            "    }\n"
            "    return NULL;\n");
    }
    fprintf(c, "}\n");
    fclose(c);

    printf("wrote %s and %s (%u baked table(s))\n", cpath, hpath, nloaded);
    return 0;
}

//
// Creates `path` and any missing parent directories (like "mkdir -p"). Returns
// 0 on success (including when the directory already exists).
//
static int make_dirs(const char* path) {
    char tmp[1024];
    snprintf(tmp, sizeof(tmp), "%s", path);

    size_t len = strlen(tmp);
    if (len == 0) {
        return 0;
    }
    if (tmp[len - 1] == '/') {
        tmp[len - 1] = '\0';
    }

    // Create each intermediate component, then the leaf.
    for (char* p = tmp + 1; *p != '\0'; p++) {
        if (*p == '/') {
            *p = '\0';
            if ((mkdir(tmp, 0755) != 0) && (errno != EEXIST)) {
                return -1;
            }
            *p = '/';
        }
    }
    if ((mkdir(tmp, 0755) != 0) && (errno != EEXIST)) {
        return -1;
    }
    return 0;
}

//------------------------------------------------------------------------------
static void usage(FILE* out, const char* prog) {
    fprintf(out,
        "Usage: %s [options]\n"
        "\n"
        "  --algo <names|all>  Comma-separated algorithm names, or 'all'\n"
        "                      (default: all)\n"
        "  --min <bits>        Smallest composite bit-size (default: 40)\n"
        "  --max <bits>        Largest composite bit-size (default: 200)\n"
        "  --step <bits>       Bit-size increment (default: 16)\n"
        "  --trials <n>        Test composites per bit-size & factor-count (25)\n"
        "  --threshold <p>     Min success rate for sieve algos (default: 0.90)\n"
        "  --min-success <p>   Min success rate for probabilistic algos (0.05)\n"
        "  --seed <u>          RNG seed for reproducible composites (default: 1)\n"
        "  --timeout <sec>     Per-run wall-clock budget; a run that exceeds it\n"
        "                      is killed and scored as a failure (0 disables;\n"
        "                      default: 30)\n"
        "  --out <dir>         Output directory (default: .)\n"
        "  --quick             Fast preset (min=48 max=64 step=16 trials=3)\n"
        "  --emit-c <file>     Code-generate baked C tables (no tuning is run)\n"
        "  --emit-h <file>     Header path to write alongside --emit-c\n"
        "  --in <dir>          Directory of tuning files to bake (default: --out)\n"
        "  --help              Print this message\n"
        "\n"
        "The requested [min,max] range is intersected with each algorithm's own\n"
        "supported range. Composites are products of similar-size primes; the\n"
        "factor-finding methods (p-1, p+1, ecm) also test 3-prime products.\n",
        prog);

    fprintf(out, "\nValid algorithm names for --algo:\n  all");
    for (unsigned int i = 0; i < TUNE_NALGOS; i++) {
        fprintf(out, " %s", TUNE_ALGOS[i].name);
    }
    fprintf(out, "\n");
}

//------------------------------------------------------------------------------
static int write_tuning_file(const options_t* opt, const tune_algo_desc_t* desc,
                             const tuning_row_t* rows, uint32_t nrows) {
    char path[1024];
    snprintf(path, sizeof(path), "%s/%s.tuning", opt->outdir, desc->name);

    FILE* f = fopen(path, "w");
    if (f == NULL) {
        fprintf(stderr, "tifa_tune: cannot open '%s' for writing\n", path);
        return 1;
    }

    double thr = desc->probabilistic ? opt->min_success : opt->threshold;

    fprintf(f, "# TIFA tuning file (schema v1)\n");
    fprintf(f, "# algo:      %s\n", desc->name);
    fprintf(f, "# generated: by tifa_tune (%s)\n", TIFA_STRING);
    fprintf(f, "# host:      %s\n", TIFA_COMPILED_ON);
    fprintf(f, "# objective: min mean_cpu / max(p_success, eps)"
               " s.t. p_success >= %.2f\n", thr);
    fprintf(f, "# trials:    %u base (x2 below 90 bits, x4 below 64 bits)\n",
            opt->trials);
    if (desc->max_factors > 2) {
        fprintf(f, "# composites: 2 balanced primes below %u bits,"
                   " %u similar-size primes at/above\n",
                desc->multi_factor_min_bits, desc->max_factors);
    } else {
        fprintf(f, "# composites: 2 balanced primes\n");
    }
    fprintf(f, "# key:       bitsize(n)\n");
    fprintf(f, "columns: min_bitsize %s\n", desc->colnames);

    for (uint32_t i = 0; i < nrows; i++) {
        fprintf(f, "%" PRIu32, rows[i].min_bitsize);
        for (uint32_t c = 0; c < desc->ncols; c++) {
            fprintf(f, " %" PRIu32, rows[i].cols[c]);
        }
        fprintf(f, "\n");
    }

    fclose(f);
    printf("  wrote %s (%" PRIu32 " rows)\n", path, nrows);
    return 0;
}

//------------------------------------------------------------------------------
//                          Progress reporting
//------------------------------------------------------------------------------
//
// A dependency-free, self-erasing progress bar for the tuning of one bit-size.
// It appears only when stdout is a terminal AND the work for that size is slow
// (measured or projected to exceed PROGRESS_THRESHOLD_S), so fast sizes and
// redirected output stay clean.
//
#define PROGRESS_THRESHOLD_S 5.0
#define PROGRESS_BAR_WIDTH   24
#define PROGRESS_LINE_WIDTH  72   // fixed field width so redraws fully overwrite

typedef struct {
    unsigned int bits;
    uint64_t     total;      // total factorization runs for this bit-size
    uint64_t     done;
    double       start;      // wall-clock start
    double       last_draw;  // wall-clock time of the last redraw
    bool         enabled;    // stdout is a terminal
    bool         active;     // the bar is currently being shown
} progress_t;

static double wall_now(void) {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double) tv.tv_sec + (double) tv.tv_usec / 1.0e6;
}

static void progress_init(progress_t* p, unsigned int bits, uint64_t total) {
    p->bits      = bits;
    p->total     = total;
    p->done      = 0;
    p->start     = wall_now();
    p->last_draw = 0.0;
    p->enabled   = (isatty(fileno(stdout)) != 0);
    p->active    = false;
}

static void progress_draw(const progress_t* p) {
    double elapsed = wall_now() - p->start;
    double frac    = (p->total > 0) ? (double) p->done / (double) p->total : 1.0;
    if (frac > 1.0) {
        frac = 1.0;
    }
    double eta = (p->done > 0)
               ? elapsed * (double) (p->total - p->done) / (double) p->done
               : 0.0;

    char bar[PROGRESS_BAR_WIDTH + 1];
    int filled = (int) (frac * PROGRESS_BAR_WIDTH);
    for (int i = 0; i < PROGRESS_BAR_WIDTH; i++) {
        bar[i] = (i < filled) ? '#' : '-';
    }
    bar[PROGRESS_BAR_WIDTH] = '\0';

    char eta_str[16];
    if (eta < 60.0) {
        snprintf(eta_str, sizeof(eta_str), "%.0fs", eta);
    } else {
        snprintf(eta_str, sizeof(eta_str), "%dm%02ds",
                 (int) (eta / 60.0), (int) eta % 60);
    }

    char line[128];
    int n = snprintf(line, sizeof(line),
        "  %u bits: [%s] %3.0f%%  %llu/%llu runs, ETA %s",
        p->bits, bar, frac * 100.0,
        (unsigned long long) p->done, (unsigned long long) p->total, eta_str);
    while ((n < PROGRESS_LINE_WIDTH) && (n < (int) sizeof(line) - 1)) {
        line[n++] = ' ';
    }
    line[n] = '\0';

    // Carriage return + no newline redraws in place; stdout is line-buffered,
    // so an explicit flush is required.
    printf("\r%s", line);
    fflush(stdout);
}

static void progress_tick(progress_t* p) {
    p->done++;
    if (!p->enabled) {
        return;
    }
    double now     = wall_now();
    double elapsed = now - p->start;

    if (!p->active) {
        bool slow = (elapsed > PROGRESS_THRESHOLD_S);
        if (!slow && (p->done >= 4)) {
            // Project total time from the work done so far.
            double projected = (elapsed / (double) p->done) * (double) p->total;
            slow = (projected > PROGRESS_THRESHOLD_S);
        }
        if (!slow) {
            return;
        }
        p->active = true;
    }
    if (now - p->last_draw >= 0.25) {
        progress_draw(p);
        p->last_draw = now;
    }
}

static void progress_finish(const progress_t* p) {
    if (p->enabled && p->active) {
        // Erase the bar so the caller's result line starts on a clean line.
        printf("\r%*s\r", PROGRESS_LINE_WIDTH, "");
        fflush(stdout);
    }
}

//
// Number of test composites to use per factor-count at a given bit-size.
// Small numbers run in microseconds, so we use proportionally more of them to
// beat down timing noise and sharpen the success-rate estimate; larger numbers
// are expensive, so we stay at the requested base count.
//
static unsigned int trials_for_bits(unsigned int base, unsigned int bits) {
    if (bits < 64) {
        return base * 4;
    }
    if (bits < 90) {
        return base * 2;
    }
    return base;
}

//------------------------------------------------------------------------------
//                          Per-run watchdog
//------------------------------------------------------------------------------
//
// Result a child process reports back to the parent.
//
typedef struct {
    int    success;
    double cpu;
} run_result_t;

// How a watched run ended (for reporting; both non-OK outcomes are scored as
// failures with the full budget as their cost, so the tuner rejects them).
typedef enum { RUN_OK = 0, RUN_TIMEOUT, RUN_CRASH } run_outcome_t;

//
// Runs one factorization under a wall-clock budget. The factorization executes
// in a forked child; if it exceeds budget_s the parent kills it and scores the
// run as a failure (so pathological parameter/input combinations can never
// stall the whole tune). A killed child's memory is discarded by the OS, so no
// GMP allocation leaks back into the tuner -- the reason fork is used instead
// of interrupting the call in-process.
//
// When budget_s <= 0 the watchdog is disabled and the run executes in-process.
// On any fork/pipe failure it also falls back to running in-process.
//
static bool run_watched(const tune_algo_desc_t* desc,
                        const mpz_t n, const mpz_t p, const uint32_t cols[],
                        double budget_s, double* secs_out,
                        run_outcome_t* outcome) {
    if (outcome != NULL) {
        *outcome = RUN_OK;
    }

    // Watchdog disabled or unavailable: run in-process (a crash here would take
    // the tuner down with it -- the trade-off of turning the watchdog off).
    if (budget_s <= 0.0) {
        return desc->run_once(n, p, cols, secs_out);
    }

    int fds[2];
    if (pipe(fds) != 0) {
        return desc->run_once(n, p, cols, secs_out);
    }

    pid_t pid = fork();
    if (pid < 0) {
        close(fds[0]);
        close(fds[1]);
        return desc->run_once(n, p, cols, secs_out);
    }

    if (pid == 0) {
        // Child: compute, report, and _exit WITHOUT flushing the parent's
        // stdio buffers (so a half-written progress line is not duplicated).
        close(fds[0]);
        double secs = 0.0;
        bool ok = desc->run_once(n, p, cols, &secs);
        run_result_t r = { ok ? 1 : 0, secs };
        ssize_t w = write(fds[1], &r, sizeof(r));
        (void) w;
        close(fds[1]);
        _exit(0);
    }

    // Parent: wait for a result on the pipe, up to the budget.
    close(fds[1]);
    double start = wall_now();
    run_result_t  r  = { 0, budget_s };
    bool          got = false;
    run_outcome_t oc  = RUN_CRASH;   // default if the child never reports

    for (;;) {
        double remaining = budget_s - (wall_now() - start);
        if (remaining <= 0.0) {
            oc = RUN_TIMEOUT;
            break;
        }
        fd_set rfds;
        FD_ZERO(&rfds);
        FD_SET(fds[0], &rfds);
        struct timeval tv;
        tv.tv_sec  = (time_t) remaining;
        tv.tv_usec = (suseconds_t) ((remaining - (double) tv.tv_sec) * 1.0e6);

        int sel = select(fds[0] + 1, &rfds, NULL, NULL, &tv);
        if (sel > 0) {
            // Readable: a full result means success; EOF/short read means the
            // child died without reporting (e.g. it crashed).
            got = (read(fds[0], &r, sizeof(r)) == (ssize_t) sizeof(r));
            oc  = got ? RUN_OK : RUN_CRASH;
            break;
        }
        if (sel == 0) {
            oc = RUN_TIMEOUT;            // budget exhausted
            break;
        }
        if (errno == EINTR) {
            continue;
        }
        oc = RUN_CRASH;                  // select error: give up on this run
        break;
    }
    close(fds[0]);

    if (!got) {
        kill(pid, SIGKILL);             // harmless no-op if already dead
        r.success = 0;
        r.cpu     = budget_s;           // penalize equally, so it is never chosen
    }
    int status;
    while ((waitpid(pid, &status, 0) < 0) && (errno == EINTR)) {
        /* retry */
    }

    if (outcome != NULL) {
        *outcome = oc;
    }
    *secs_out = r.cpu;
    return r.success != 0;
}

//
// Tunes a single algorithm across the configured bit-size range and writes its
// tuning file. Returns 0 on success.
//
static int tune_one_algo(const options_t* opt, const tune_algo_desc_t* desc,
                         gmp_randstate_t rng) {
    // Intersect the requested range with the algorithm's supported range.
    unsigned int lo = (opt->min_bits > desc->min_bits) ? opt->min_bits
                                                        : desc->min_bits;
    unsigned int hi = (opt->max_bits < desc->max_bits) ? opt->max_bits
                                                        : desc->max_bits;

    printf("Tuning %s (%s, %u..%u bits, ",
           desc->name, desc->probabilistic ? "probabilistic" : "deterministic",
           lo, hi);
    if (desc->max_factors > 2) {
        printf("2-prime <%u bits, %u-prime >=%u bits)\n",
               desc->multi_factor_min_bits, desc->max_factors,
               desc->multi_factor_min_bits);
    } else {
        printf("2 balanced primes)\n");
    }

    if (lo > hi) {
        printf("  (no bit-sizes in the requested range; skipping)\n");
        return 0;
    }

    double thr = desc->probabilistic ? opt->min_success : opt->threshold;

    tuning_row_t rows[TUNE_MAX_ROWS];
    uint32_t     nrows = 0;
    // Once a "sticky-on" column (the large-prime variation) has won at
    // TUNE_STICKY_HYSTERESIS consecutive sizes, stop re-testing its "off" value
    // at larger bit-sizes -- the crossover is one-way. Hysteresis guards against
    // a single noisy win locking off out prematurely; an "off" win resets it.
    unsigned sticky_wins   = 0;
    bool     sticky_locked = false;

    for (unsigned int bits = lo; bits <= hi && nrows < TUNE_MAX_ROWS;
         bits += opt->step) {

        // Use more test composites for small (cheap) numbers to beat down
        // timing noise.
        unsigned int trials = trials_for_bits(opt->trials, bits);
        unsigned int setlen = trials;

        // Choose the factor count for this size: the requested max_factors at
        // or above the crossover, else 2 balanced primes.
        unsigned int nfactors =
            (bits >= desc->multi_factor_min_bits) ? desc->max_factors : 2;

        // Build the shared test set once, so every grid point is compared on
        // identical composites.
        mpz_t* nset = malloc(setlen * sizeof(mpz_t));
        mpz_t* pset = malloc(setlen * sizeof(mpz_t));
        for (unsigned int t = 0; t < setlen; t++) {
            mpz_init(nset[t]);
            mpz_init(pset[t]);
            tune_gen_composite(nset[t], pset[t], nfactors, bits, rng);
        }

        grid_t grid;
        desc->build_grid(bits, &grid);
        // If the large-prime variation already won at a smaller size, collapse
        // its axis to "on" so we don't re-benchmark the losing "off" case.
        if (sticky_locked && (desc->sticky_col >= 0)
            && ((uint32_t) desc->sticky_col < grid.naxes)) {
            grid.axes[desc->sticky_col].nvalues  = 1;
            grid.axes[desc->sticky_col].values[0] = 1;
        }
        uint64_t npoints = grid_size(&grid);

        progress_t pr;
        progress_init(&pr, bits, npoints * (uint64_t) setlen);

        bool     have_best = false;
        bool     best_qualified = false;
        double   best_score = 0.0;
        uint32_t best_cols[TIFA_TUNING_MAX_COLS];
        unsigned timeouts = 0;
        unsigned crashes  = 0;

        for (uint64_t g = 0; g < npoints; g++) {
            uint32_t cols[TIFA_TUNING_MAX_COLS];
            grid_point(&grid, g, cols);

            // Once a point has spent more than the qualified best's total, its
            // score cannot beat it: stop timing it, and cut a run off before
            // it alone would. The slack keeps fork and wall-clock jitter from
            // cutting short runs.
            const double bound = best_qualified ? best_score * (double) setlen
                                                : INFINITY;
            double   sum_cpu = 0.0;
            unsigned succ    = 0;
            bool     pruned  = false;
            for (unsigned int t = 0; t < setlen; t++) {
                double budget = opt->timeout_s;
                bool   cut    = false;
                if ((budget > 0.0) && isfinite(bound)) {
                    const double left = 2.0 * (bound - sum_cpu) + 0.1;
                    if (left < budget) {
                        budget = left;
                        cut    = true;
                    }
                }
                double secs = 0.0;
                run_outcome_t oc = RUN_OK;
                bool ok = run_watched(desc, nset[t], pset[t], cols,
                                      budget, &secs, &oc);
                sum_cpu += secs;
                if (ok) {
                    succ++;
                }
                if ((oc == RUN_TIMEOUT) && !cut) {
                    timeouts++;
                } else if (oc == RUN_CRASH) {
                    crashes++;
                }
                progress_tick(&pr);
                if (sum_cpu > bound) {
                    pr.done += setlen - t - 1;
                    pruned   = true;
                    break;
                }
            }
            if (pruned) {
                continue;
            }
            double p_success = (double) succ / (double) setlen;
            double mean_cpu  = sum_cpu / (double) setlen;
            double score     = mean_cpu / fmax(p_success, TUNE_EPS);
            bool   qualified = (p_success >= thr);

            // Prefer a qualifying point; among equally (un)qualified points,
            // prefer the lower score.
            bool better;
            if (!have_best) {
                better = true;
            } else if (qualified != best_qualified) {
                better = qualified;                 // qualifying always wins
            } else {
                better = (score < best_score);
            }
            if (better) {
                have_best      = true;
                best_qualified = qualified;
                best_score     = score;
                memcpy(best_cols, cols, sizeof(best_cols));
            }
        }

        progress_finish(&pr);

        if (have_best) {
            rows[nrows].min_bitsize = bits;
            memcpy(rows[nrows].cols, best_cols, sizeof(best_cols));
            printf("  %3u bits [%u nums, %u-prime]:", bits, setlen, nfactors);
            for (uint32_t c = 0; c < desc->ncols; c++) {
                printf(" %" PRIu32, rows[nrows].cols[c]);
            }
            printf("   (score %.6g%s", best_score,
                   best_qualified ? "" : ", below threshold");
            if (timeouts > 0) {
                printf(", %u run%s hit the %gs limit", timeouts,
                       (timeouts == 1) ? "" : "s", opt->timeout_s);
            }
            if (crashes > 0) {
                printf(", %u run%s crashed", crashes,
                       (crashes == 1) ? "" : "s");
            }
            printf(")\n");
            nrows++;

            // Track consecutive "on" wins; lock only after enough of them.
            if ((desc->sticky_col >= 0) && !sticky_locked) {
                if (best_cols[desc->sticky_col] == 1) {
                    sticky_wins++;
                    if (sticky_wins >= TUNE_STICKY_HYSTERESIS) {
                        sticky_locked = true;
                    }
                } else {
                    sticky_wins = 0;   // an "off" win breaks the streak
                }
            }
        }

        for (unsigned int t = 0; t < setlen; t++) {
            mpz_clear(nset[t]);
            mpz_clear(pset[t]);
        }
        free(nset);
        free(pset);
    }

    return write_tuning_file(opt, desc, rows, nrows);
}

//------------------------------------------------------------------------------
int main(int argc, char** argv) {
    // Line-buffer stdout so per-bit-size progress is visible even when the
    // output is piped or redirected to a file.
    setvbuf(stdout, NULL, _IOLBF, 0);

    options_t opt = {
        .algo        = "all",
        .min_bits    = 40,
        .max_bits    = 200,
        .step        = 16,
        .trials      = 25,
        .threshold   = 0.90,
        .min_success = 0.05,
        .timeout_s   = 30.0,
        .seed        = 1,
        .outdir      = "."
    };

    // Manual argument parsing keeps the tool portable and dependency-free.
    for (int i = 1; i < argc; i++) {
        const char* a = argv[i];

        #define NEXT_VAL(dstptr)                                              \
            do {                                                              \
                if (i + 1 >= argc) {                                          \
                    fprintf(stderr, "tifa_tune: %s needs a value\n", a);      \
                    return 2;                                                 \
                }                                                             \
                *(dstptr) = argv[++i];                                        \
            } while (0)

        if (strcmp(a, "--help") == 0) {
            usage(stdout, argv[0]);
            return 0;
        } else if (strcmp(a, "--quick") == 0) {
            // A deliberately coarse, fast preset for a first pass. Widen the
            // range/trials (or drop --quick) for a more thorough search.
            opt.min_bits = 48;
            opt.max_bits = 64;
            opt.step     = 16;
            opt.trials   = 3;
        } else if (strcmp(a, "--algo") == 0) {
            NEXT_VAL(&opt.algo);
        } else if (strcmp(a, "--out") == 0) {
            NEXT_VAL(&opt.outdir);
        } else if (strcmp(a, "--min") == 0) {
            const char* v; NEXT_VAL(&v); opt.min_bits = (unsigned) strtoul(v, NULL, 10);
        } else if (strcmp(a, "--max") == 0) {
            const char* v; NEXT_VAL(&v); opt.max_bits = (unsigned) strtoul(v, NULL, 10);
        } else if (strcmp(a, "--step") == 0) {
            const char* v; NEXT_VAL(&v); opt.step = (unsigned) strtoul(v, NULL, 10);
        } else if (strcmp(a, "--trials") == 0) {
            const char* v; NEXT_VAL(&v); opt.trials = (unsigned) strtoul(v, NULL, 10);
        } else if (strcmp(a, "--threshold") == 0) {
            const char* v; NEXT_VAL(&v); opt.threshold = strtod(v, NULL);
        } else if (strcmp(a, "--min-success") == 0) {
            const char* v; NEXT_VAL(&v); opt.min_success = strtod(v, NULL);
        } else if (strcmp(a, "--timeout") == 0) {
            const char* v; NEXT_VAL(&v); opt.timeout_s = strtod(v, NULL);
        } else if (strcmp(a, "--seed") == 0) {
            const char* v; NEXT_VAL(&v); opt.seed = strtoul(v, NULL, 10);
        } else if (strcmp(a, "--emit-c") == 0) {
            NEXT_VAL(&opt.emit_c);
        } else if (strcmp(a, "--emit-h") == 0) {
            NEXT_VAL(&opt.emit_h);
        } else if (strcmp(a, "--in") == 0) {
            NEXT_VAL(&opt.indir);
        } else {
            fprintf(stderr, "tifa_tune: unknown option '%s'\n", a);
            usage(stderr, argv[0]);
            return 2;
        }
        #undef NEXT_VAL
    }

    // Code-generation mode: read existing tuning files and emit baked tables.
    // No benchmarking is performed.
    if (opt.emit_c != NULL) {
        if (opt.emit_h == NULL) {
            fprintf(stderr, "tifa_tune: --emit-c requires --emit-h\n");
            return 2;
        }
        const char* indir = (opt.indir != NULL) ? opt.indir : opt.outdir;
        return emit_baked(opt.emit_c, opt.emit_h, indir);
    }

    if ((opt.step == 0) || (opt.min_bits < 4) || (opt.max_bits < opt.min_bits)
        || (opt.trials == 0)) {
        fprintf(stderr, "tifa_tune: invalid bit-size range or trial count\n");
        return 2;
    }

    if (make_dirs(opt.outdir) != 0) {
        fprintf(stderr, "tifa_tune: cannot create output directory '%s': %s\n",
                opt.outdir, strerror(errno));
        return 1;
    }

    gmp_randstate_t rng;
    gmp_randinit_default(rng);
    gmp_randseed_ui(rng, opt.seed);

    int rc = 0;
    if (strcmp(opt.algo, "all") == 0) {
        for (unsigned int i = 0; i < TUNE_NALGOS; i++) {
            rc |= tune_one_algo(&opt, &TUNE_ALGOS[i], rng);
        }
    } else {
        // Parse a comma-separated list of algorithm names. Validate every name
        // first so a typo fails immediately, before any (possibly long) run.
        char buf[256];
        snprintf(buf, sizeof(buf), "%s", opt.algo);

        const tune_algo_desc_t* selected[TUNE_MAX_SELECTED];
        unsigned int nsel = 0;

        char* save = NULL;
        for (char* tok = strtok_r(buf, ",", &save);
             tok != NULL;
             tok = strtok_r(NULL, ",", &save)) {
            // Trim surrounding whitespace so "a, b" works.
            while ((*tok == ' ') || (*tok == '\t')) {
                tok++;
            }
            char* end = tok + strlen(tok);
            while ((end > tok)
                   && ((end[-1] == ' ') || (end[-1] == '\t') || (end[-1] == '\n'))) {
                *--end = '\0';
            }
            if (*tok == '\0') {
                continue;   // empty field (e.g. trailing comma)
            }
            const tune_algo_desc_t* desc = tune_find_algo(tok);
            if (desc == NULL) {
                fprintf(stderr, "tifa_tune: unknown algorithm '%s'\n", tok);
                gmp_randclear(rng);
                return 2;
            }
            if (nsel < TUNE_MAX_SELECTED) {
                selected[nsel++] = desc;
            }
        }
        if (nsel == 0) {
            fprintf(stderr, "tifa_tune: no algorithm named in --algo\n");
            gmp_randclear(rng);
            return 2;
        }
        for (unsigned int i = 0; i < nsel; i++) {
            rc |= tune_one_algo(&opt, selected[i], rng);
        }
    }

    gmp_randclear(rng);
    return rc;
}
