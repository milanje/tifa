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
 * \brief   tifa-tune: measure and record good algorithm parameters per bit-size.
 *
 * For each selected algorithm and bit-size, tifa-tune generates a shared set of
 * test composites, searches the algorithm's parameter grid around its defaults
 * for the size (every combination, or one parameter at a time with --search
 * descent, all on the same composites for a fair comparison, and dropping one
 * as soon as it has taken longer than the best), times each run, and selects
 * the parameter combination minimizing the expected time to find a factor. The
 * winning combination for each bit-size is written as one row of a
 * per-algorithm <tt>&lt;algo&gt;.tuning</tt> file that the library can load at
 * runtime (see lib/algo/include/tuning.h), or smoothed across the sizes with
 * --fit.
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
#include <limits.h>
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
#include "tifa_alloc.h"

#define TUNE_EPS         1e-9
#define TUNE_MAX_ROWS    256
#define TUNE_MAX_SELECTED 32   // cap on algorithms named in a --algo list
// Consecutive sizes a "sticky-on" column must win before it is locked on, so a
// single noisy win cannot prematurely stop the off case from being tested.
#define TUNE_STICKY_HYSTERESIS 2
// --fit smooths each row over this many tuned sizes at most, the nearest.
#define TUNE_FIT_WINDOW 5
// The automatic watchdog: this many times the defaults' longest run, at least
// TUNE_MIN_TIMEOUT seconds.
#define TUNE_TIMEOUT_FACTOR 10.0
#define TUNE_MIN_TIMEOUT    30.0

typedef struct {
    const char*   algo;         // "all" or a specific algorithm name
    unsigned int  min_bits;
    unsigned int  max_bits;
    unsigned int  step;
    unsigned int  sizes[TUNE_MAX_ROWS];   // --sizes, sorted, else none
    unsigned int  nsizes;
    bool          descent;      // coordinate descent rather than the whole grid
    unsigned int  fit_step;     // --fit: rows every fit_step bits (0: no fit)
    unsigned int  trials;       // test composites per (bit-size, factor-count)
    double        threshold;    // required success rate for deterministic algos
    double        min_success;  // required success rate for probabilistic algos
    double        timeout_s;    // per-run wall-clock budget (0 disables
                                // watchdog, < 0 automatic)
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
// expose tifa_tuning_baked_table(). Mirrors the offline genprimes.py model.
//
static int emit_baked(const char* cpath, const char* hpath, const char* indir) {
    // The output may overwrite committed tables: refuse a mistyped --in.
    struct stat st;
    if ((stat(indir, &st) != 0) || !S_ISDIR(st.st_mode)) {
        fprintf(stderr, "tifa-tune: '%s' is not a directory\n", indir);
        return 1;
    }
    FILE* h = fopen(hpath, "w");
    if (h == NULL) {
        fprintf(stderr, "tifa-tune: cannot open '%s' for writing\n", hpath);
        return 1;
    }
    write_license(h);
    fprintf(h,
        "\n// Automatically generated by tifa-tune --emit-c. Do not edit by hand.\n\n"
        "#if !defined(_TIFA_TUNING_BAKED_H_)\n"
        "#define _TIFA_TUNING_BAKED_H_\n\n"
        "#ifdef __cplusplus\nextern \"C\" {\n#endif\n\n"
        "#include <inttypes.h>\n\n"
        "#include \"tuning.h\"\n\n"
        "// Returns the compiled-in tuning table for an algorithm, or NULL if there\n"
        "// is none. Used as a fallback when no runtime tuning file is found.\n"
        "const tuning_table_t* tifa_tuning_baked_table(const char* algo,\n"
        "                                              uint32_t expect_ncols);\n\n"
        "#ifdef __cplusplus\n}\n#endif\n\n#endif\n");
    fclose(h);

    FILE* c = fopen(cpath, "w");
    if (c == NULL) {
        fprintf(stderr, "tifa-tune: cannot open '%s' for writing\n", cpath);
        return 1;
    }
    write_license(c);
    fprintf(c,
        "\n// Automatically generated by tifa-tune --emit-c. Do not edit by hand.\n\n"
        "#include <string.h>\n\n"
        "#include \"tuning_baked.h\"\n\n");

    const char* loaded_names[TUNE_MAX_ROWS];
    uint32_t    loaded_nrows[TUNE_MAX_ROWS];
    uint32_t    loaded_ncols[TUNE_MAX_ROWS];
    unsigned    nloaded = 0;

    for (unsigned int i = 0; i < TUNE_NALGOS; i++) {
        const tune_algo_desc_t* desc = &TUNE_ALGOS[i];
        char path[1024];
        const int length = snprintf(path, sizeof(path), "%s/%s.tuning", indir,
                                    desc->name);
        if ((length < 0) || ((size_t) length >= sizeof(path))) {
            fprintf(stderr, "tifa-tune: path too long in '%s'\n", indir);
            fclose(c);
            return 1;
        }

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
        fprintf(c, "static const tuning_table_t baked_tables[] = {\n");
        for (unsigned int i = 0; i < nloaded; i++) {
            fprintf(c, "    { \"%s\", %" PRIu32 ", %" PRIu32 ", "
                       "(tuning_row_t*) %s_baked_rows },\n",
                    loaded_names[i], loaded_ncols[i], loaded_nrows[i],
                    loaded_names[i]);
        }
        fprintf(c,
            "};\n\n"
            "static const unsigned int nbaked_tables =\n"
            "    sizeof(baked_tables) / sizeof(baked_tables[0]);\n\n");
    }

    fprintf(c,
        "const tuning_table_t* tifa_tuning_baked_table(const char* algo,\n"
        "                                              uint32_t expect_ncols) {\n");
    if (nloaded == 0) {
        fprintf(c,
            "    (void) algo; (void) expect_ncols;\n"
            "    return NULL;\n");
    } else {
        fprintf(c,
            "    for (unsigned int i = 0; i < nbaked_tables; i++) {\n"
            "        if ((strcmp(baked_tables[i].algo, algo) == 0)\n"
            "            && (baked_tables[i].ncols == expect_ncols)) {\n"
            "            return &baked_tables[i];\n"
            "        }\n"
            "    }\n"
            "    return NULL;\n");
    }
    fprintf(c, "}\n");
    fclose(c);

    if (nloaded == 0) {
        fprintf(stderr, "tifa-tune: no tuning file in '%s': no baked table\n",
                indir);
    }
    printf("wrote %s and %s (%u baked table(s))\n", cpath, hpath, nloaded);
    return 0;
}

//
// Creates `path` and any missing parent directories (like "mkdir -p"). Returns
// 0 on success (including when the directory already exists).
//
static int make_dirs(const char* path) {
    char tmp[1024];
    const int length = snprintf(tmp, sizeof(tmp), "%s", path);
    if ((length < 0) || ((size_t) length >= sizeof(tmp))) {
        fprintf(stderr, "tifa-tune: path too long: '%s'\n", path);
        return 1;
    }

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
//
// Parses an unsigned decimal option value: strtoul alone would take a sign
// (wrapping "-1" to UINT_MAX) and ignore trailing junk.
//
static bool parse_unsigned(const char* v, unsigned int* out) {
    if ((v[0] < '0') || (v[0] > '9')) {
        return false;
    }
    char* end = NULL;
    errno = 0;
    unsigned long x = strtoul(v, &end, 10);
    if ((*end != '\0') || (errno == ERANGE) || (x > UINT_MAX)) {
        return false;
    }
    *out = (unsigned int) x;
    return true;
}

//------------------------------------------------------------------------------
static bool parse_double_in(const char* v, double lo, double hi, double* out) {
    char* end = NULL;
    double x = strtod(v, &end);
    if ((end == v) || (*end != '\0') || !isfinite(x) || (x < lo) || (x > hi)) {
        return false;
    }
    *out = x;
    return true;
}

//------------------------------------------------------------------------------
//
// Parses a comma-separated list of bit-sizes into opt->sizes, sorted, without
// duplicates.
//
static bool parse_sizes(const char* v, options_t* opt) {
    opt->nsizes = 0;
    const char* p = v;
    for (;;) {
        char field[16];
        const size_t length = strcspn(p, ",");
        if ((length == 0) || (length >= sizeof(field))
            || (opt->nsizes == TUNE_MAX_ROWS)) {
            return false;
        }
        memcpy(field, p, length);
        field[length] = '\0';
        unsigned int bits;
        if (!parse_unsigned(field, &bits) || (bits < 4)) {
            return false;
        }
        uint32_t i = 0;
        while ((i < opt->nsizes) && (opt->sizes[i] < bits)) {
            i++;
        }
        if ((i == opt->nsizes) || (opt->sizes[i] != bits)) {
            memmove(&opt->sizes[i + 1], &opt->sizes[i],
                    (opt->nsizes - i) * sizeof(opt->sizes[0]));
            opt->sizes[i] = bits;
            opt->nsizes++;
        }
        if (p[length] == '\0') {
            return true;
        }
        p += length + 1;
    }
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
        "  --step <bits>       Bit-size increment; --max is always tuned\n"
        "                      (default: 16)\n"
        "  --sizes <list>      Comma-separated bit-sizes to tune, instead of\n"
        "                      --min, --max and --step\n"
        "  --search <how>      'grid' tries every combination of parameters\n"
        "                      (default); 'descent' varies one at a time from\n"
        "                      the defaults, much faster for large sizes, and\n"
        "                      varies the SIQS sieve threshold too\n"
        "  --fit <bits>        Write rows every <bits> bits from the smallest\n"
        "                      tuned size to the largest, each fitted on the\n"
        "                      5 nearest tuned sizes (a least-squares quadratic)\n"
        "  --trials <n>        Test composites per bit-size, x2 below 90 bits\n"
        "                      and x4 below 64 (default: 25)\n"
        "  --threshold <p>     Min success rate for sieve algos (default: 0.90)\n"
        "  --min-success <p>   Min success rate for probabilistic algos (0.05)\n"
        "  --seed <u>          RNG seed for reproducible composites (default: 1)\n"
        "  --timeout <sec>     Per-run wall-clock budget; a run that exceeds it\n"
        "                      is killed and scored as a failure (0 disables;\n"
        "                      default: none for the defaults, then 10 times\n"
        "                      their longest run, at least 30)\n"
        "  --out <dir>         Output directory (default: .)\n"
        "  --quick             Fast preset (min=48 max=64 step=16 trials=3)\n"
        "  --emit-c <file>     Code-generate baked C tables (no tuning is run)\n"
        "  --emit-h <file>     Header path to write alongside --emit-c\n"
        "  --in <dir>          Directory of tuning files to bake (default: --out)\n"
        "  --help              Print this message\n"
        "\n"
        "The requested [min,max] range is intersected with each algorithm's own\n"
        "supported range. Composites are products of similar-size primes; the\n"
        "factor-finding methods (p-1, p+1, rho, ecm) use 3-prime products from\n"
        "64 bits.\n",
        prog);

    fprintf(out, "\nValid algorithm names for --algo:\n  all");
    for (unsigned int i = 0; i < TUNE_NALGOS; i++) {
        fprintf(out, " %s", TUNE_ALGOS[i].name);
    }
    fprintf(out, "\n");
}

//------------------------------------------------------------------------------
//
// Writes rows, fitted on the tuned rows if those are given (else NULL).
//
static int write_tuning_file(const options_t* opt, const tune_algo_desc_t* desc,
                             const tuning_row_t* rows, uint32_t nrows,
                             const tuning_row_t* tuned, uint32_t ntuned) {
    char path[1024];
    const int length = snprintf(path, sizeof(path), "%s/%s.tuning",
                                opt->outdir, desc->name);
    if ((length < 0) || ((size_t) length >= sizeof(path))) {
        fprintf(stderr, "tifa-tune: path too long in '%s'\n", opt->outdir);
        return 1;
    }

    FILE* f = fopen(path, "w");
    if (f == NULL) {
        fprintf(stderr, "tifa-tune: cannot open '%s' for writing\n", path);
        return 1;
    }

    double thr = desc->probabilistic ? opt->min_success : opt->threshold;

    fprintf(f, "# TIFA tuning file (schema v1)\n");
    fprintf(f, "# algo:      %s\n", desc->name);
    fprintf(f, "# generated: by tifa-tune (%s)\n", TIFA_STRING);
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
    fprintf(f, "# search:    %s\n", opt->descent ? "descent" : "grid");
    fprintf(f, "# key:       bitsize(n)\n");
    if (tuned != NULL) {
        fprintf(f, "# fit:       rows every %u bits, each a least-squares"
                   " quadratic on the %u nearest tuned rows:\n",
                opt->fit_step, TUNE_FIT_WINDOW);
        for (uint32_t i = 0; i < ntuned; i++) {
            fprintf(f, "#   %" PRIu32, tuned[i].min_bitsize);
            for (uint32_t c = 0; c < desc->ncols; c++) {
                fprintf(f, " %" PRIu32, tuned[i].cols[c]);
            }
            fprintf(f, "\n");
        }
    }
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
// When budget_s <= 0 the watchdog is disabled and the run executes in-process;
// an infinite budget_s only guards against crashes. On any fork/pipe failure it
// also falls back to running in-process.
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
        struct timeval* timeout = NULL;
        if (isfinite(remaining)) {
            tv.tv_sec  = (time_t) remaining;
            tv.tv_usec = (suseconds_t) ((remaining - (double) tv.tv_sec) * 1.0e6);
            timeout    = &tv;
        }

        int sel = select(fds[0] + 1, &rfds, NULL, NULL, timeout);
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
// The search for the best parameters at one bit-size.
//
typedef struct {
    const options_t*        opt;
    const tune_algo_desc_t* desc;
    const grid_t*           grid;
    mpz_t*       nset;
    mpz_t*       pset;
    unsigned int setlen;
    double       thr;
    double       budget;     // per-run watchdog budget (see run_watched)
    progress_t   pr;
    bool         have_best;
    bool         best_qualified;
    double       best_score;
    uint32_t     best_cols[TIFA_TUNING_MAX_COLS];
    unsigned     timeouts;
    unsigned     crashes;
    uint64_t     visited[TUNE_MAX_ROWS];
    unsigned     nvisited;
} search_t;

//
// Times the point with the value of index idx[i] on each axis i, keeps it if it
// is the best yet, and returns whether it is. Under the automatic watchdog, the
// first point timed has no budget and sets the others'.
//
static bool evaluate(search_t* s, const uint32_t idx[]) {
    uint64_t linear = 0;
    for (uint32_t i = s->grid->naxes; i-- > 0;) {
        linear = linear * s->grid->axes[i].nvalues + idx[i];
    }
    for (unsigned i = 0; i < s->nvisited; i++) {
        if (s->visited[i] == linear) {
            return false;
        }
    }
    if (s->nvisited < TUNE_MAX_ROWS) {
        s->visited[s->nvisited++] = linear;
    }

    uint32_t cols[TIFA_TUNING_MAX_COLS] = { 0 };
    grid_point_at(s->grid, idx, cols);
    if (s->desc->tie != NULL) {
        s->desc->tie(cols);
    }
    // Once a point has spent more than the qualified best's total, its score
    // cannot beat it: stop timing it, and cut a run off before it alone
    // would. The slack keeps fork and wall-clock jitter from cutting short
    // runs.
    const double bound = s->best_qualified ? s->best_score * (double) s->setlen
                                           : INFINITY;
    double   sum_cpu = 0.0;
    double   longest = 0.0;
    unsigned succ    = 0;
    for (unsigned int t = 0; t < s->setlen; t++) {
        double budget = s->budget;
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
        bool ok = run_watched(s->desc, s->nset[t], s->pset[t], cols,
                              budget, &secs, &oc);
        sum_cpu += secs;
        if (ok) {
            succ++;
        }
        if (oc == RUN_OK) {
            longest = fmax(longest, secs);
        } else if ((oc == RUN_TIMEOUT) && !cut) {
            s->timeouts++;
        } else if (oc == RUN_CRASH) {
            s->crashes++;
        }
        progress_tick(&s->pr);
        if (sum_cpu > bound) {
            s->pr.done += s->setlen - t - 1;
            return false;
        }
    }
    if (!isfinite(s->budget)) {
        s->budget = fmax(TUNE_MIN_TIMEOUT, TUNE_TIMEOUT_FACTOR * longest);
    }
    if (succ == 0) {
        return false;
    }
    double p_success = (double) succ / (double) s->setlen;
    double mean_cpu  = sum_cpu / (double) s->setlen;
    double score     = mean_cpu / fmax(p_success, TUNE_EPS);
    bool   qualified = (p_success >= s->thr);

    // Prefer a qualifying point; among equally (un)qualified points, prefer
    // the lower score.
    bool better;
    if (!s->have_best) {
        better = true;
    } else if (qualified != s->best_qualified) {
        better = qualified;
    } else {
        better = (score < s->best_score);
    }
    if (better) {
        s->have_best      = true;
        s->best_qualified = qualified;
        s->best_score     = score;
        memcpy(s->best_cols, cols, sizeof(s->best_cols));
    }
    return better;
}

//
// Times every point of the grid, the defaults first.
//
static void search_grid(search_t* s) {
    const grid_t* g = s->grid;
    uint32_t idx[TIFA_TUNING_MAX_COLS];
    for (uint32_t i = 0; i < g->naxes; i++) {
        idx[i] = g->axes[i].start;
    }
    evaluate(s, idx);
    memset(idx, 0, sizeof(idx));
    for (;;) {
        evaluate(s, idx);
        uint32_t i = 0;
        while ((i < g->naxes) && (++idx[i] == g->axes[i].nvalues)) {
            idx[i++] = 0;
        }
        if (i == g->naxes) {
            return;
        }
    }
}

//
// From the defaults, tries each value of one axis at a time, the others at their
// best yet, until a pass over the axes improves nothing (two passes at most).
//
static void search_descent(search_t* s) {
    const grid_t* g = s->grid;
    uint32_t idx[TIFA_TUNING_MAX_COLS];
    for (uint32_t i = 0; i < g->naxes; i++) {
        idx[i] = g->axes[i].start;
    }
    evaluate(s, idx);
    for (unsigned pass = 0; pass < 2; pass++) {
        bool improved = false;
        for (uint32_t a = 0; a < g->naxes; a++) {
            const uint32_t from = idx[a];
            uint32_t best = from;
            for (uint32_t v = 0; v < g->axes[a].nvalues; v++) {
                if (v == from) {
                    continue;
                }
                idx[a] = v;
                if (evaluate(s, idx)) {
                    best     = v;
                    improved = true;
                }
            }
            idx[a] = best;
        }
        if (!improved) {
            return;
        }
        uint64_t more = 0;
        for (uint32_t a = 0; a < g->naxes; a++) {
            more += g->axes[a].nvalues - 1;
        }
        s->pr.total += more * s->setlen;
    }
}

//
// Fits the column values of the tuned rows at the sizes (x) of the rows, with a
// least-squares polynomial of degree at most 2 on x - at, and returns its value
// at x = at: the fit's constant coefficient.
//
static double fit_at(const double* xs, const double* ys, uint32_t n, double at) {
    double sx[5] = { 0 };
    double sy[3] = { 0 };
    for (uint32_t i = 0; i < n; i++) {
        const double x = xs[i] - at;
        double p = 1.0;
        for (int k = 0; k < 5; k++) {
            if (k < 3) {
                sy[k] += p * ys[i];
            }
            sx[k] += p;
            p *= x;
        }
    }
    for (uint32_t degree = (n > 2) ? 2 : n - 1; degree > 0; degree--) {
        // Gaussian elimination on the normal equations.
        const uint32_t m = degree + 1;
        double a[3][4];
        for (uint32_t r = 0; r < m; r++) {
            for (uint32_t c = 0; c < m; c++) {
                a[r][c] = sx[r + c];
            }
            a[r][m] = sy[r];
        }
        bool singular = false;
        for (uint32_t c = 0; (c < m) && !singular; c++) {
            uint32_t pivot = c;
            for (uint32_t r = c + 1; r < m; r++) {
                if (fabs(a[r][c]) > fabs(a[pivot][c])) {
                    pivot = r;
                }
            }
            if (fabs(a[pivot][c]) < 1e-12 * fabs(sx[0])) {
                singular = true;
                break;
            }
            for (uint32_t k = 0; k <= m; k++) {
                const double tmp = a[c][k];
                a[c][k] = a[pivot][k];
                a[pivot][k] = tmp;
            }
            for (uint32_t r = 0; r < m; r++) {
                if (r != c) {
                    const double f = a[r][c] / a[c][c];
                    for (uint32_t k = c; k <= m; k++) {
                        a[r][k] -= f * a[c][k];
                    }
                }
            }
        }
        if (!singular) {
            return a[0][m] / a[0][0];
        }
    }
    return sy[0] / sx[0];
}

//
// Writes into out a row every opt->fit_step bits from the first tuned size to
// the last, each fitted on the TUNE_FIT_WINDOW nearest tuned rows and kept
// within their values. Returns the number of rows, 0 if there are too many.
//
static uint32_t fit_rows(const options_t* opt, const tune_algo_desc_t* desc,
                         const tuning_row_t* tuned, uint32_t ntuned,
                         tuning_row_t* out) {
    const uint32_t first = tuned[0].min_bitsize;
    const uint32_t last  = tuned[ntuned - 1].min_bitsize;
    const uint32_t k     = (ntuned < TUNE_FIT_WINDOW) ? ntuned : TUNE_FIT_WINDOW;
    uint32_t nout = 0;
    for (uint64_t at = first; at <= last;
         at = (at == last) ? at + 1
                           : ((last - at < opt->fit_step) ? last
                                                          : at + opt->fit_step)) {
        if (nout == TUNE_MAX_ROWS) {
            return 0;
        }
        uint32_t lo = 0;
        while ((lo + k < ntuned)
               && ((int64_t) at - tuned[lo].min_bitsize
                   > (int64_t) tuned[lo + k].min_bitsize - (int64_t) at)) {
            lo++;
        }
        uint32_t below = 0;
        while ((below + 1 < ntuned) && (tuned[below + 1].min_bitsize <= at)) {
            below++;
        }
        tuning_row_t* row = &out[nout++];
        memset(row, 0, sizeof(*row));
        row->min_bitsize = (uint32_t) at;
        for (uint32_t c = 0; c < desc->ncols; c++) {
            const char kind = desc->fit[c];
            if (kind == 'b') {
                row->cols[c] = tuned[below].cols[c];
                continue;
            }
            double xs[TUNE_FIT_WINDOW];
            double ys[TUNE_FIT_WINDOW];
            double lowest  = INFINITY;
            double highest = 0.0;
            for (uint32_t i = 0; i < k; i++) {
                const double v = (double) tuned[lo + i].cols[c];
                xs[i]   = (double) tuned[lo + i].min_bitsize;
                ys[i]   = (kind == 'g') ? log2(fmax(v, 1.0)) : v;
                lowest  = fmin(lowest, v);
                highest = fmax(highest, v);
            }
            double v = lowest;
            if (lowest < highest) {
                v = fit_at(xs, ys, k, (double) at);
                if (kind == 'g') {
                    v = exp2(v);
                }
                v = fmin(fmax(v, lowest), highest);
            }
            row->cols[c] = (uint32_t) floor(v + 0.5);
        }
        if (desc->tie != NULL) {
            desc->tie(row->cols);
        }
    }
    return nout;
}

//
// Tunes a single algorithm across the configured bit-sizes and writes its
// tuning file. Returns 0 on success.
//
static int tune_one_algo(const options_t* opt, const tune_algo_desc_t* desc,
                         gmp_randstate_t rng) {
    // The requested sizes within the algorithm's supported range; the last
    // step of a range stops at its end, so the table covers the whole range.
    unsigned int sizes[TUNE_MAX_ROWS];
    unsigned int nsizes = 0;
    if (opt->nsizes > 0) {
        for (unsigned int i = 0; i < opt->nsizes; i++) {
            if ((opt->sizes[i] >= desc->min_bits)
                && (opt->sizes[i] <= desc->max_bits)) {
                sizes[nsizes++] = opt->sizes[i];
            } else {
                printf("%s: skipping %u bits, outside %u..%u\n", desc->name,
                       opt->sizes[i], desc->min_bits, desc->max_bits);
            }
        }
    } else {
        unsigned int lo = (opt->min_bits > desc->min_bits) ? opt->min_bits
                                                            : desc->min_bits;
        unsigned int hi = (opt->max_bits < desc->max_bits) ? opt->max_bits
                                                            : desc->max_bits;
        for (unsigned int bits = lo; bits <= hi && nsizes < TUNE_MAX_ROWS;
             bits = (bits == hi) ? hi + 1
                                 : ((hi - bits < opt->step) ? hi
                                                            : bits + opt->step)) {
            sizes[nsizes++] = bits;
        }
    }

    if (nsizes == 0) {
        printf("Tuning %s: no bit-sizes in the requested range; skipping\n",
               desc->name);
        return 0;
    }
    printf("Tuning %s (%s, %u..%u bits, ",
           desc->name, desc->probabilistic ? "probabilistic" : "deterministic",
           sizes[0], sizes[nsizes - 1]);
    if (desc->max_factors > 2) {
        printf("2-prime <%u bits, %u-prime >=%u bits)\n",
               desc->multi_factor_min_bits, desc->max_factors,
               desc->multi_factor_min_bits);
    } else {
        printf("2 balanced primes)\n");
    }

    tuning_row_t rows[TUNE_MAX_ROWS];
    uint32_t     nrows = 0;
    // Once a "sticky-on" column (the large-prime variation) has won at
    // TUNE_STICKY_HYSTERESIS consecutive sizes, stop re-testing its "off" value
    // at larger bit-sizes -- the crossover is one-way. Hysteresis guards against
    // a single noisy win locking off out prematurely; an "off" win resets it.
    unsigned sticky_wins   = 0;
    bool     sticky_locked = false;

    for (unsigned int is = 0; is < nsizes; is++) {
        const unsigned int bits = sizes[is];

        // Use more test composites for small (cheap) numbers to beat down
        // timing noise.
        unsigned int setlen = trials_for_bits(opt->trials, bits);

        // Choose the factor count for this size: the requested max_factors at
        // or above the crossover, else 2 balanced primes.
        unsigned int nfactors =
            (bits >= desc->multi_factor_min_bits) ? desc->max_factors : 2;

        // Build the shared test set once, so every grid point is compared on
        // identical composites.
        mpz_t* nset = tifa_malloc(setlen * sizeof(mpz_t));
        mpz_t* pset = tifa_malloc(setlen * sizeof(mpz_t));
        for (unsigned int t = 0; t < setlen; t++) {
            mpz_init(nset[t]);
            mpz_init(pset[t]);
            tune_gen_composite(nset[t], pset[t], nfactors, bits, rng);
        }

        grid_t grid;
        desc->build_grid(bits, opt->descent, &grid);
        // If the large-prime variation won at TUNE_STICKY_HYSTERESIS consecutive
        // smaller sizes, collapse its axis to "on" so we don't re-benchmark the
        // losing "off" case.
        if (sticky_locked && (desc->sticky_col >= 0)
            && ((uint32_t) desc->sticky_col < grid.naxes)) {
            grid.axes[desc->sticky_col].nvalues   = 1;
            grid.axes[desc->sticky_col].values[0] = 1;
            grid.axes[desc->sticky_col].start     = 0;
        }

        search_t* s = tifa_malloc(sizeof(search_t));
        memset(s, 0, sizeof(*s));
        s->opt    = opt;
        s->desc   = desc;
        s->grid   = &grid;
        s->nset   = nset;
        s->pset   = pset;
        s->setlen = setlen;
        s->thr    = desc->probabilistic ? opt->min_success : opt->threshold;
        s->budget = (opt->timeout_s < 0.0) ? INFINITY : opt->timeout_s;

        uint64_t npoints = grid_size(&grid);
        if (opt->descent) {
            npoints = 1;
            for (uint32_t a = 0; a < grid.naxes; a++) {
                npoints += grid.axes[a].nvalues - 1;
            }
        }
        progress_init(&s->pr, bits, npoints * (uint64_t) setlen);
        if (opt->descent) {
            search_descent(s);
        } else {
            search_grid(s);
        }
        progress_finish(&s->pr);

        if (s->have_best) {
            rows[nrows].min_bitsize = bits;
            memcpy(rows[nrows].cols, s->best_cols, sizeof(s->best_cols));
            printf("  %3u bits [%u nums, %u-prime]:", bits, setlen, nfactors);
            for (uint32_t c = 0; c < desc->ncols; c++) {
                printf(" %" PRIu32, rows[nrows].cols[c]);
            }
            printf("   (score %.6g%s", s->best_score,
                   s->best_qualified ? "" : ", below threshold");
            if (s->timeouts > 0) {
                printf(", %u run%s hit the %gs limit", s->timeouts,
                       (s->timeouts == 1) ? "" : "s", s->budget);
            }
            if (s->crashes > 0) {
                printf(", %u run%s crashed", s->crashes,
                       (s->crashes == 1) ? "" : "s");
            }
            printf(")\n");
            nrows++;

            // Track consecutive "on" wins; lock only after enough of them.
            if ((desc->sticky_col >= 0) && !sticky_locked) {
                if (s->best_cols[desc->sticky_col] == 1) {
                    sticky_wins++;
                    if (sticky_wins >= TUNE_STICKY_HYSTERESIS) {
                        sticky_locked = true;
                    }
                } else {
                    sticky_wins = 0;   // an "off" win breaks the streak
                }
            }
        } else {
            printf("  %3u bits: no point factored any number, no row\n", bits);
        }

        free(s);
        for (unsigned int t = 0; t < setlen; t++) {
            mpz_clear(nset[t]);
            mpz_clear(pset[t]);
        }
        free(nset);
        free(pset);
    }

    if ((opt->fit_step == 0) || (nrows < 2)) {
        return write_tuning_file(opt, desc, rows, nrows, NULL, 0);
    }
    tuning_row_t* fitted = tifa_malloc(TUNE_MAX_ROWS * sizeof(tuning_row_t));
    const uint32_t nfitted = fit_rows(opt, desc, rows, nrows, fitted);
    int rc;
    if (nfitted == 0) {
        fprintf(stderr, "tifa-tune: --fit %u makes more than %d rows\n",
                opt->fit_step, TUNE_MAX_ROWS);
        rc = 1;
    } else {
        rc = write_tuning_file(opt, desc, fitted, nfitted, rows, nrows);
    }
    free(fitted);
    return rc;
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
        .timeout_s   = -1.0,
        .seed        = 1,
        .outdir      = "."
    };

    // Manual argument parsing keeps the tool portable and dependency-free.
    for (int i = 1; i < argc; i++) {
        const char* a = argv[i];

        #define NEXT_VAL(dstptr)                                              \
            do {                                                              \
                if (i + 1 >= argc) {                                          \
                    fprintf(stderr, "tifa-tune: %s needs a value\n", a);      \
                    return 2;                                                 \
                }                                                             \
                *(dstptr) = argv[++i];                                        \
            } while (0)
        #define BAD_VAL() do {                                                \
                fprintf(stderr, "tifa-tune: bad value for %s\n", a);         \
                return 2;                                                     \
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
            const char* v; NEXT_VAL(&v);
            if (!parse_unsigned(v, &opt.min_bits)) { BAD_VAL(); }
        } else if (strcmp(a, "--max") == 0) {
            const char* v; NEXT_VAL(&v);
            if (!parse_unsigned(v, &opt.max_bits)) { BAD_VAL(); }
        } else if (strcmp(a, "--step") == 0) {
            const char* v; NEXT_VAL(&v);
            if (!parse_unsigned(v, &opt.step)) { BAD_VAL(); }
        } else if (strcmp(a, "--sizes") == 0) {
            const char* v; NEXT_VAL(&v);
            if (!parse_sizes(v, &opt)) { BAD_VAL(); }
        } else if (strcmp(a, "--search") == 0) {
            const char* v; NEXT_VAL(&v);
            if (strcmp(v, "grid") == 0) {
                opt.descent = false;
            } else if (strcmp(v, "descent") == 0) {
                opt.descent = true;
            } else {
                BAD_VAL();
            }
        } else if (strcmp(a, "--fit") == 0) {
            const char* v; NEXT_VAL(&v);
            if (!parse_unsigned(v, &opt.fit_step) || (opt.fit_step == 0)) {
                BAD_VAL();
            }
        } else if (strcmp(a, "--trials") == 0) {
            const char* v; NEXT_VAL(&v);
            if (!parse_unsigned(v, &opt.trials)) { BAD_VAL(); }
        } else if (strcmp(a, "--threshold") == 0) {
            const char* v; NEXT_VAL(&v);
            if (!parse_double_in(v, 0.0, 1.0, &opt.threshold)) { BAD_VAL(); }
        } else if (strcmp(a, "--min-success") == 0) {
            const char* v; NEXT_VAL(&v);
            if (!parse_double_in(v, 0.0, 1.0, &opt.min_success)) { BAD_VAL(); }
        } else if (strcmp(a, "--timeout") == 0) {
            // select() refuses timeouts above about 1e8 s.
            const char* v; NEXT_VAL(&v);
            if (!parse_double_in(v, 0.0, 1.0e6, &opt.timeout_s)) { BAD_VAL(); }
        } else if (strcmp(a, "--seed") == 0) {
            const char* v; NEXT_VAL(&v);
            unsigned int seed;
            if (!parse_unsigned(v, &seed)) { BAD_VAL(); }
            opt.seed = seed;
        } else if (strcmp(a, "--emit-c") == 0) {
            NEXT_VAL(&opt.emit_c);
        } else if (strcmp(a, "--emit-h") == 0) {
            NEXT_VAL(&opt.emit_h);
        } else if (strcmp(a, "--in") == 0) {
            NEXT_VAL(&opt.indir);
        } else {
            fprintf(stderr, "tifa-tune: unknown option '%s'\n", a);
            usage(stderr, argv[0]);
            return 2;
        }
        #undef NEXT_VAL
        #undef BAD_VAL
    }

    // Code-generation mode: read existing tuning files and emit baked tables.
    // No benchmarking is performed.
    if (opt.emit_c != NULL) {
        if (opt.emit_h == NULL) {
            fprintf(stderr, "tifa-tune: --emit-c requires --emit-h\n");
            return 2;
        }
        const char* indir = (opt.indir != NULL) ? opt.indir : opt.outdir;
        return emit_baked(opt.emit_c, opt.emit_h, indir);
    }

    if ((opt.step == 0) || (opt.min_bits < 4) || (opt.max_bits < opt.min_bits)
        || (opt.trials == 0)) {
        fprintf(stderr, "tifa-tune: invalid bit-size range or trial count\n");
        return 2;
    }

    if (make_dirs(opt.outdir) != 0) {
        fprintf(stderr, "tifa-tune: cannot create output directory '%s': %s\n",
                opt.outdir, strerror(errno));
        return 1;
    }
    struct stat st;
    if ((stat(opt.outdir, &st) != 0) || !S_ISDIR(st.st_mode)
            || (access(opt.outdir, W_OK) != 0)) {
        fprintf(stderr, "tifa-tune: '%s' is not a writable directory\n",
                opt.outdir);
        return 1;
    }
    // Center the grids on the built-in defaults, not on tuning files or baked
    // tables from an earlier run.
    tifa_tuning_set_dir("");

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
                fprintf(stderr, "tifa-tune: unknown algorithm '%s'\n", tok);
                gmp_randclear(rng);
                return 2;
            }
            if (nsel < TUNE_MAX_SELECTED) {
                selected[nsel++] = desc;
            }
        }
        if (nsel == 0) {
            fprintf(stderr, "tifa-tune: no algorithm named in --algo\n");
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
