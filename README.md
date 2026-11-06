# TIFA &mdash; Tools for Integer FActorization

<img src="readme/tifacube.png" align="left" width="160" hspace="12" alt="TIFA logo">

`IDDN.FR.001.220019.000.S.A.2011.000.31235`

`TIFA` (_Tools for Integer FActorization_) is a small C library that implements
utilities and algorithms to perform integer factorizations, with an emphasis on
small to medium-size composites, say from 40 to 280 bits. It can be used on its
own, or inside a larger program such as a number field sieve, to factor its
cofactors one per thread.

This library was written around 2010 when I was working at
[LIX](https://www.lix.polytechnique.fr), the Computer Science Laboratory
of [Ecole Polytechnique](https://www.polytechnique.edu/en). The source code was
uploaded "as is" to [GitHub](https://github.com) in April 2020. It has been
completely revamped in the fall of 2026, for correctness, thread safety and
above all speed: its build system moved from [SCons](https://scons.org) to
[CMake](https://cmake.org), and on the same balanced semiprimes, its SIQS is
3 to 5 times as fast as the 2020 code on an Apple M5 Pro and 2.3 to 3 times
on x86-64 (an AMD EPYC 7763) from 100 to 220 bits, `tifa_factor()` 2.2 to
3.3 times from 100 to 180 bits. On x86-64, its SIQS takes 0.6 to 0.75 of the
time of [yafu](https://github.com/bbuhrow/yafu) 3.1.9's AVX2 build at 200 and
220 bits, and about as long as its AVX-512 build.

<br clear="left">

## What's new

Version 0.4.7 (November 2026) makes SIQS 1.1 to 4.1 times as fast from 232
to 280 bits on an Apple M5 Pro (a 280-bit semiprime in 98 s rather than
407): its factor base, capped at 16384 primes, goes to 65536, the exponents
of its relations are kept sparse beyond 16384, and its parameters from 232
bits are new. `tifa_prove_prime()` proves that an integer is prime, with a
certificate anyone can check; `tifa_set_time_limit()` limits the time of a
thread's factorizations; `tifa-tune` tunes large sizes; and the Python
package proves primes, takes a `timeout` and reads tuning files. A review
for untrusted input fixed the certificate checker, the tuning files, the
factoring core and the Python wheels (see
[readme/news.md](./readme/news.md#4-november-2026--security-fixes-047)).
The ABI only grows (`libtifa.so.0.4`).

Version 0.4.6 (November 2026) brings a Python package: `pip install
tifa-factor`, then `tifa.factor(n)` returns the prime factorization of `n`
as `{prime: exponent}`, in 0.24 to 0.51 of the time of PARI/GP's
`factorint` on balanced semiprimes from 100 to 220 bits (see
[Python](#python)). Where the processor has AVX-512, the bucket fill of
SIQS takes sixteen primes at a time: `siqs()` takes 0.84 to 0.87 of the
time of 0.4.5 at 220 bits (EPYC 9V45 and 9V74, Xeon 6975P-C). The ABI is
unchanged (`libtifa.so.0.4`).

Version 0.4.5 (November 2026) makes SIQS 0.94 to 0.98 of the time of 0.4.4
from 140 to 220 bits on an Apple M5 Pro, and 0.96 to 0.99 on x86-64 (EPYC
7763 and 9V74, Xeon 8370C): the candidate pre-check finds the small primes
of each candidate with doubles and a mask, rather than 128-bit reductions,
random branches and GMP calls, and a new linear-algebra method,
`FILTERED_GAUSS_ELIM`, filters the matrix before the elimination (0.29 to
0.46 of the time of
`SMART_GAUSS_ELIM` from 180 to 220 bits). The table of primes is generated
by a Python 3 script: TIFA no longer has any Perl. The ABI is unchanged
(`libtifa.so.0.4`).

Version 0.4.4 (November 2026) makes SIQS 1.25 to 1.35 times faster at 200
and 220 bits, on an Apple M5 Pro as on x86-64 (EPYC 7763 and Xeon 8370C):
the candidates with no chance of a relation are dropped before the batch
smoothness test, 40 small primes are left out of the sieve, the polynomial
switch takes no divisions, and on x86-64 the scan uses SSE2 and the bucket
fill no longer the branches GCC made. `tifa_factor()` is 2% to 18% faster on
balanced semiprimes from 80 to 180 bits. The ABI is unchanged
(`libtifa.so.0.4`).

Version 0.4.3 (November 2026) makes SIQS 1.5 to 1.6 times faster at 200
and 220 bits on x86-64 (EPYC 7763 and 9V74) and up to 1.3 times on an Apple
M5 Pro: GCC no longer turns the sieve's branch-free code into branches, the
primes from 8192 are sieved eight at a time with SSE2 on x86-64, Q2
polynomials are used when kn = 1 mod 8, and the parameters from 160 to 220
bits are retuned. `tifa_factor()` is 5% to 30% faster on balanced semiprimes
from 80 to 180 bits. The ABI is unchanged (`libtifa.so.0.4`).

Version 0.4.2 (November 2026) makes `tifa_factor()` search 74- to 83-bit
numbers with ECM before SIQS again: those with a factor of 20 to 26 bits take
0.50 to 0.86 of 0.4.1's time on an Apple M5 Pro and on x86-64 (an EPYC 9V74,
Zen 4), while balanced semiprimes of these sizes take 1.17 to 1.28. `ecm()`
keeps the plans of its levels from one call to the next: 2% to 6% faster from
70 to 120 bits on the M5, 4% to 10% on x86-64. The ABI is unchanged
(`libtifa.so.0.4`).

Version 0.4.1 (November 2026) makes SQUFOF 4 to 8 times faster on an Apple
M5 Pro (2.5 to 4 times on x86-64, an EPYC 7763, Zen 3) and p+1 12% to 32%
faster; `tifa_factor()` is 5% to 19% faster on random composites from 100 to
160 bits, and 3 to 10 times faster on numbers with a factor of 26 to 40 bits
from 120 to 160 bits. The new `tifa_factor_with_hint()` takes what is known
of the number to factor (two factors of about the same size, or a small
one). The ABI only grows: the shared library is still `libtifa.so.0.4`.

Version 0.4.0 (October 2026) makes ECM 19% to 33% faster on an Apple M5 Pro
and 9% to 19% on x86-64 (an EPYC 7763, Zen 3), and SIQS 10% to 14% faster
from 110 to 130 bits on the M5 (7% to 11% on x86-64) and 3% to 9% above;
`tifa_factor()` is 7% to 11% faster from 128 to 160 bits on the M5 (3% to 8%
on x86-64). The ABI changes (`libtifa.so.0.4`): the shared library exports
the public API alone, and only its headers are installed.

Version 0.3.4 (October 2026) makes SIQS 15% to 20% faster from 135 to 190
bits on x86-64 (an EPYC 9V74, Zen 4), whose 32 KiB sieve chunks the default
half widths no longer overflow, and 3% to 6% faster on x86-64 from 100 to 130
bits; on an Apple M5 Pro, 4% to 8% faster from 150 to 200 bits. `tifa_factor()`
is 17% to 18% faster at 140 and 160 bits on x86-64. The ABI is unchanged
(`libtifa.so.0.3`).

Version 0.3.3 (October 2026) makes SIQS 14% to 18% faster from 66 to 90 bits
and 22% to 30% from 95 to 200, against 0.3.2 on one thread of an Apple M5 Pro;
`tifa_factor()` is 14% to 15% faster at 66 and 80 bits, 26% to 28% from 100 to
160. SIQS tests its candidates for smoothness eight at a time in machine
integers and sieves without unpredictable branches; the benchmark programs are
in `bench/`. The ABI is unchanged (`libtifa.so.0.3`).

Version 0.3.2 (October 2026) makes SIQS 32% faster at 66 and 70 bits, 21% to
25% from 75 to 85, 11% to 18% from 90 to 110 and 1% to 10% above, against
0.3.1 on one thread of an Apple M5 Pro; `tifa_factor()` is 12% to 26% faster
from 66 to 100 bits. SIQS keeps its largest arrays from one run to the next
in each thread, uses fewer GMP calls outside its sieve, and sieves faster on
small factor bases. The ABI is unchanged (`libtifa.so.0.3`).

Version 0.3.1 (October 2026) makes SIQS 31% faster from 66 to 75 bits, 18% to
28% from 80 to 100 and 3% to 13% above, against 0.3.0 on one thread of an
Apple M5 Pro; `tifa_factor()` is 14% to 28% faster from 66 to 100 bits. SIQS
spends less time outside its sieve (Gaussian elimination, trial division, the
setup of each new `a`) and collects 8 extra relations instead of 24, more
only when none splits the number. The programs take `-h` and `--help`. The ABI
is unchanged (`libtifa.so.0.3`).

Version 0.3.0 (October 2026) makes SIQS faster and fixes latent defects found
in a review of the whole library. Against 0.2.0, on one thread of an Apple M5
Pro:

- SIQS is 7% to 20% faster from 66 to 110 bits (a faster setup, retuned
  defaults and number of primes in `a`), 1% to 10% above; `tifa_factor()`
  is 10% to 15% faster from 66 to 100 bits.
- `tifa_factor()` takes apart any perfect power (p^3 gave up), `cfrac()` and
  `squfof()` no longer hang on some inputs, and Fermat's square roots modulo
  p^2 are right for p above 65535.
- ECM, p+1 and p-1 split small composites they gave up on (up to 25% of them
  near 10^5) when a stage caught every prime factor at once.
- The programs are named `tifa-<name>` (`tifa-factor`, `tifa-ecm`, ...) and
  exit with 0, 1, 2 or 3 for a complete, partial, failed factorization or
  invalid input. The ABI changes: `libtifa.so.0.3`.

Version 0.2.0 (October 2026) reworked most of the library for speed,
correctness and thread safety: `tifa_factor()` became 6x to 100x faster from
32 to 110 bits, `tifa_cofactor()` appeared, and the factorization functions
can be called from several threads at once.

[readme/news.md](./readme/news.md) has all the measurements and the list of
changes.

## What's in there

- `tifa_factor()` (`tifa_factor.h`) factors an integer completely, choosing the
  methods by size: machine-word arithmetic below 2^64; above, trial division by
  the primes below 1000, a short ECM search for factors of up to 60 bits (85
  above 256 bits; without 128-bit integers, single p-1 and p+1 runs up to 512
  bits), then SIQS. Above 512 bits, where SIQS would take years, single p-1
  and p+1 runs and ECM in GMP arithmetic look for factors of up to 85 bits,
  for as long as that takes (about 12 minutes at 520 bits on an Apple M5 Pro).
  `tifa_factor_with_hint()` does the same, told what is known of the number:
  two prime factors of about the same size (SIQS runs right away), or one
  below a third of its size (ECM searches for longer).
- `tifa_set_time_limit()` (`time_limit.h`) limits the time of the calling
  thread's factorizations: SIQS, its linear algebra and ECM give up soon
  after it, or when the function given to `tifa_set_interrupt_check()` says
  so (Ctrl-C in Python).
- `tifa_cofactor()` (`cofactor.h`) decides whether an integer splits into primes
  below a bound, as the cofactorization step of a sieve asks (see
  [Cofactorization](#cofactorization) below).
- `tifa_prove_prime()` (`prime_cert.h`) proves that an integer is prime, and
  returns a certificate that `tifa_check_prime_cert()` checks, as anyone can:
  steps from primes of p - 1 and p + 1 (Pocklington, Morrison, and
  Brillhart, Lehmer and Selfridge), which `tifa_prime_cert_to_str()` writes
  as JSON and `tifa_prime_cert_from_str()` reads back. `tifa_factor()` tells
  primes by a Baillie-PSW test, which proves nothing above 2^64;
  `tifa-factor --prove` proves the factors it prints.
- `tifa_batch_gcd()` (`batch_gcd.h`) finds, among many integers, those that
  share a prime with another: the gcd of each with the product of the
  others, from product and remainder trees, as for the RSA keys of the
  Internet (Heninger et al., "Mining Your Ps and Qs", 2012). 100,000 moduli
  of 2048 bits take 30 s on an Apple M5 Pro. `tifa-batch-gcd` does it on
  files of integers.
- The algorithms themselves, each with its own header and parameters: SIQS
  (`siqs.h`), CFRAC (`cfrac.h`), ECM (`ecm.h`), Pollard's p-1 (`p_minus_1.h`),
  Williams' p+1 (`p_plus_1.h`), Pollard's rho (`rho.h`), SQUFOF (`squfof.h`),
  McKee's variant of Fermat's method (`fermat.h`) and trial division
  (`tdiv.h`).
- Below 2^64, deterministic primality testing and factorization in machine
  words (`factor_u64.h`); ECM, p-1 and p+1 in Montgomery arithmetic on one to
  eight words, up to 512 bits (`mont_ecm.h`).

- A Python package (`python/`): `tifa.factor(n)` is `tifa_factor()`, with
  proofs of the primes on request; `tifa.batch_gcd()` is `tifa_batch_gcd()`,
  and `tifa.Cofactorizer` gives `tifa_cofactor()` for one integer or a whole
  list at a time (see [Python](#python) below).

All of them can run in several threads at once. SIQS keeps its largest arrays
from one run to the next in each thread (170 KiB at most, freed when the thread
exits).

The [readme.txt](./readme/readme.txt) file in the `readme` directory has more
detail on the algorithms and the source layout; its build and installation
steps use CMake as well, matching the instructions below.

## Building

Requirements:

- a C99 compiler,
- [GMP](https://gmplib.org) (>= 6.0), the sole external dependency,
- POSIX threads,
- [CMake](https://cmake.org) >= 3.20.

Optional: [Doxygen](https://www.doxygen.nl) (+ LaTeX for the PDF) for the API
documentation.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build           # the test suite
cmake --install build            # honors -DCMAKE_INSTALL_PREFIX=...
```

The machine-word code (below 2^64), the Montgomery ECM, p-1 and p+1, and the
values, smoothness test and square roots of SIQS's candidates need a compiler
with 128-bit integers (GCC and Clang have them on 64-bit platforms); without
them, the library falls back on its slower GMP code.

If GMP lives in a non-standard location, pass `-DGMP_INCDIR=<dir>` and
`-DGMP_LIBDIR=<dir>`. Build options are set with `-D<OPTION>=<value>` (run
`cmake -LH build` to list them); for example `-DTIFA_BITSTRING_T=...`,
`-DTIFA_ALLOW_VERBOSE=ON`, `-DTIFA_ALLOW_TIMING=ON`.

The build produces both a static (`libtifa.a`) and a shared (`libtifa.so` /
`libtifa.dylib`) library, plus the command-line factorization programs and the
`tifa-tune` parameter-tuning tool (see [Parameter tuning](#parameter-tuning)
below). `tifa-factor` factors a number with the methods the `tifa_factor()`
function chooses, given a hint (`random`, the default, `balanced` or
`unbalanced`) if there is one:

```sh
./build/tools/tifa-factor 19409887069306060631113
./build/tools/tifa-factor balanced 206010655745875572166197184221848916071
./build/tools/tifa-factor --prove 340282366920938463463374607431768211457
```

With `--prove`, `tifa-factor` also proves the prime factors and prints
their certificates in JSON, as `tifa.check_certificate` takes them (see
[Python](#python)); a factor it cannot prove prime makes the factorization
partial.

The other programs run one algorithm each (`tifa-cfrac`, `tifa-ecm`,
`tifa-fermat`, `tifa-siqs`, `tifa-squfof`, `tifa-rho`, `tifa-p-minus-1`,
`tifa-p-plus-1`, `tifa-tdiv`); given only the number, they
use the algorithm's default parameters (`-h` or `--help` lists the others,
which come before the number). All the factorization programs exit with
status 0 for a complete factorization (or a prime, or 1), 1 for a partial one,
2 if no factor was found, and 3 for invalid input.

`tifa-batch-gcd` reads positive integers (RSA moduli, say) from files, or
from its standard input, and prints the position (from 1) and the gcd of
each one that shares a prime with another, from `tifa_batch_gcd()`:

```sh
./build/tools/tifa-batch-gcd moduli.txt          # one decimal integer per line
./build/tools/tifa-batch-gcd --hex moduli.txt    # hexadecimal, as are the gcds
./build/tools/tifa-batch-gcd --all moduli.txt    # every gcd, 1 included
```

A gcd that is neither 1 nor the integer is a factor of it; the integer
itself means that all its primes are shared, as with a duplicate. It exits
with status 0 if some integers share primes, 2 if none do, and 3 for invalid
input. 20,000 integers of 2048 bits take 3.9 s on an Apple M5 Pro.

Installation also
provides a `tifa.pc` pkg-config file and a CMake package, so downstream projects
can consume TIFA with either:

```cmake
find_package(Tifa REQUIRED)
target_link_libraries(your_app PRIVATE Tifa::tifa)          # or Tifa::tifa_static
```

```sh
cc app.c $(pkg-config --cflags --libs tifa) -o app
```

Other targets: `cmake --build build --target doc` (Doxygen documentation, if
available) and `cpack --config build/CPackSourceConfig.cmake` (source tarball).

The [2011 user's guide](https://www.lix.polytechnique.fr/Labo/Jerome.Milan/tifa/doc/tifa_user.pdf)
(PDF, on the
[project's LIX webpage](https://www.lix.polytechnique.fr/Labo/Jerome.Milan/tifa/tifa.xhtml))
describes the library before the 2026 rework; the Doxygen documentation
(`cmake --build build --target doc`) describes the current one.

## Python

The `tifa` package gives Python 3.10 and later `tifa_factor()`:

```python
>>> import tifa
>>> tifa.factor(360)
{2: 3, 3: 2, 5: 1}
>>> tifa.factor(2**64 + 1, hint="unbalanced")
{274177: 1, 67280421310721: 1}
```

`factor(n)` returns the prime factorization of a positive integer as
`{prime: exponent}`, and releases the GIL, so that threads factor in parallel.
On balanced semiprimes it takes 0.24 to 0.51 of the time of PARI/GP's
`factorint` from 100 to 220 bits. `factor(n, proof=True)` also proves the
primes, and `prime_certificate(p)` returns a proof that
`check_certificate(p, cert)` checks in plain Python; `batch_gcd(moduli)` is
`tifa_batch_gcd()`, and `Cofactorizer(lpb, ...)` gives `tifa_cofactor()`,
for one integer or a whole list at a time. `factor(n, timeout=60)`
gives up after about a minute (`tifa.FactorizationTimeout`, a
`TimeoutError`), and Ctrl-C interrupts it. The package installs
`tifa-tune` as well, and `set_tuning_dir(path)` takes the parameters from
the files it writes (see [Parameter tuning](#parameter-tuning)).

From version 0.4.6 on, `pip install tifa-factor` installs it from PyPI, with
wheels for Linux and macOS. From this tree, at its top:

```sh
pip install ".[test]"    # builds the package (CMake >= 3.26 and GMP, as above)
python -m pytest         # its tests
```

`GMP_ROOT=<dir> pip install .` takes GMP from `<dir>/include` and `<dir>/lib`.
[python/README.md](./python/README.md), the package's page on PyPI, has the
rest.

### Wheels for PyPI

`.github/workflows/wheels.yml` builds the wheels with
[cibuildwheel](https://cibuildwheel.pypa.io) (its settings are in
`pyproject.toml`) on each `v*` tag and on demand: Linux x86-64 and arm64
(manylinux and musllinux) and macOS arm64 and x86-64, one abi3 wheel each,
tested with Python 3.10 and 3.15, and the sdist. GMP 6.3.0 is built for them
by `python/build_gmp.sh` (its checksum checked; on x86-64, with the code of
each processor chosen at run time) and linked statically. The package's
version is the library's, from `CMakeLists.txt`.

The workflow publishes only on demand, by PyPI's trusted publishing (no
token). Once, before the first upload:

1. on [PyPI](https://pypi.org/manage/account/publishing/) and on
   [TestPyPI](https://test.pypi.org/manage/account/publishing/), add a pending
   trusted publisher for the project `tifa-factor`: owner `milanje`,
   repository `tifa`, workflow `wheels.yml`, environment `pypi` (on PyPI) or
   `testpypi` (on TestPyPI);
2. in the repository's settings (Environments), create the environments
   `pypi` and `testpypi`, with required reviewers if wanted.

Then, for each release, once its tag is pushed: Actions, Wheels, Run
workflow, choose the tag and `testpypi`; check
`pip install -i https://test.pypi.org/simple/ tifa-factor`; then the same with
`pypi`. The workflow publishes from a tag only, and stops when the sdist's
version is not the tag's.

## Cofactorization

For the cofactors left by the sieve of the number field sieve,
`tifa_cofactor()` (`cofactor.h`) answers whether every prime factor is below
2^lpb, and stops as soon as the answer is known instead of factoring
completely. The bounds go into a strategy built once, which threads may share:

```c
cofactor_params_t params = {
    .lpb              = 31,         // accept prime factors below 2^31
    .fb_bound         = 1UL << 22,  // the cofactors have no smaller factor
    .max_large_primes = 3,
    .max_siqs_bits    = 62,         // above 2 lpb bits, ECM only: give up
                                    // rather than run SIQS (0: always exact)
};
cofactor_strategy_t* strategy = alloc_cofactor_strategy(&params);

if (tifa_cofactor(factors, multis, n, strategy) == COFACTOR_SMOOTH) {
    // factors and multis hold the prime factorization of n
}
free_cofactor_strategy(strategy);
```

In Python, `tifa.Cofactorizer(31, fb_bound=2**22, max_large_primes=3,
max_siqs_bits=62)` holds the strategy, and its `many(numbers)` answers for a
list at 3% more than the time of the C function on 64-bit integers.

Above `max_siqs_bits`, a composite that ECM does not split within the time SIQS
would take is given up (`COFACTOR_GAVE_UP`). With three large primes and
`max_siqs_bits = 2 lpb`, that takes 60% to 65% of the time of exact answers and
misses 0.25% to 0.30% of the smooth cofactors (see
[readme/news.md](./readme/news.md#cofactorization)).

## Parameter tuning

Each algorithm depends on numeric parameters (smoothness bounds `b1`/`b2`, curve
counts, sieve widths, factor-base sizes, ...) whose best values depend on the
bit-size of the number being factored. `set_<algo>_params_to_default(n, params)`
fills these; by default it uses built-in, size-indexed values. The `tifa-tune`
tool lets you measure better values on your own machine and feed them back to
the library.

There are three ways parameters reach the library, tried in order, with the
built-in defaults always as the final fallback:

1. **A runtime tuning file** — `<dir>/<algo>.tuning`, where `<dir>` is the one
   set by `tifa_tuning_set_dir()` (`tuning.h`, for embedders), else the
   `TIFA_TUNING_DIR` environment variable (ignored by setuid and setgid
   programs), else the directory compiled in with `-DTIFA_TUNING_DIR=<dir>`
   (default `<prefix>/share/tifa/tuning`; none in the Python package). You
   can swap files without recompiling.
2. **Compiled-in ("baked") tables** — enabled with `-DTIFA_USE_BAKED_TUNING=ON`;
   consulted when no runtime file is found (see below).
3. **The built-in defaults** — used whenever neither of the above applies, so an
   untouched build behaves exactly as before this feature existed.
   `tifa_tuning_set_dir("")` selects them alone, skipping 1 and 2.

A table applies only to sizes up to 8 bits beyond its first and last rows:
outside, the built-in defaults do, rather than the parameters of numbers much
smaller or larger.

A tuning file lists rows of parameters for increasing bit-sizes. Between two
rows, most quantities (bounds, sizes, counts) are interpolated; the flags,
seeds, choices, the SIQS sieve width (whose cost jumps at each multiple of the
sieve chunk) and its extra primes in `a` are taken from the row below.

### Generating tuning files

```sh
cmake --build build --target tifa_tune     # the CMake target of tifa-tune
./build/tools/tune/tifa-tune --algo all --out ~/.tifa/tuning
```

`tifa-tune` builds test composites (products of similar-size primes), sweeps a
parameter grid around the built-in defaults for each bit-size, times each
run, and writes the combination minimizing the expected time to find a factor.
A combination is dropped as soon as its runs have taken longer in total than
the best one's. Useful options (`--help` lists
all): `--algo <names|all>` (a comma-separated list of algorithm names, or `all`;
`--help` lists the valid names), `--min`/`--max`/`--step` (bit-size range; `--max` is always tuned),
`--sizes <list>` (comma-separated bit-sizes instead of a range), `--trials`,
`--search descent` and `--fit <bits>` (see [Large SIQS sizes](#large-siqs-sizes)),
`--seed` (reproducible composites), `--timeout <sec>` (per-run wall-clock budget;
a run that exceeds it is killed in its child process and scored as a failure, so
a pathological parameter/input combination can never stall the tune — by
default none for the defaults, then 10 times their longest run and at least
30 s; `0` disables), and `--quick` (a fast, coarse preset). More
test composites are used for small numbers (twice as many below 90 bits, four
times below 64), where runs are cheap and timing noise is otherwise larger. When
a given size takes more than a few seconds, a progress bar with an ETA is shown —
only on an interactive terminal, so redirected output stays plain. A
convenience target runs the quick preset into `<build>/tuning`:

```sh
cmake --build build --target tune-quick
```

Point the library at the generated files at runtime:

```sh
TIFA_TUNING_DIR=~/.tifa/tuning ./build/tools/tifa-siqs <number_to_factor>
```

#### Large SIQS sizes

From about 240 bits, a SIQS run takes from seconds to minutes, and the best
parameters differ between processors: the built-in sieve half widths, chosen
on an Apple M5 Pro and recent x86-64 cores, took up to 9% longer than twice as
wide ones on an AMD EPYC 7763 (Zen 3) from 264 to 280 bits. To tune a few large
sizes:

```sh
./build/tools/tune/tifa-tune --algo siqs --sizes 232,248,264,280 \
    --search descent --trials 2 --fit 8 --out ~/.tifa/tuning
```

`--search descent` starts from the built-in defaults and tries each value of
one parameter at a time — the factor base, the sieve half width, the sieve
threshold, the large-prime variation and the extra primes in `a` — keeping the
others at their best yet, for at most two passes: 14 combinations per pass
rather than the grid's 100. `--fit 8` then writes a row every 8 bits from
the smallest tuned size to the largest, each a least-squares quadratic on the 5
nearest tuned sizes (on the logarithm of the factor base and of the half
width), kept within their values: this smooths the noise of a few numbers per
size and fills in the sizes between. The tuned rows are kept as comments at the
top of the file.

### Baking tuned values into the library

To compile the tuned tables in (so no runtime file or environment variable is
needed), regenerate the committed baked tables and enable them:

```sh
# One-off: run a full tuning pass and regenerate lib/algo/{src,include}/tuning_baked.*
cmake -S . -B build -DTIFA_REGEN_TUNING=ON
cmake --build build --target regen-tuning

# Then build with the baked tables compiled in
cmake -S . -B build -DTIFA_USE_BAKED_TUNING=ON
cmake --build build
```

The committed `tuning_baked.c` ships empty (the accessor returns nothing, so the
built-in defaults are used); `regen-tuning` overwrites it with your machine's
measured values, which you may then commit. Tuning data is machine-specific, so
the shipped baseline is intentionally empty.

## Measuring performance

The timings of [news.md](./readme/news.md) come from the programs of `bench/`,
built with `-DTIFA_BUILD_BENCH=ON` (they are never installed):

- `siqs-bench <bits> [<count> [<seconds>]]` times `siqs()` on `<count>`
  balanced semiprimes of `<bits>` bits (`<lo>-<hi>` for random sizes), factored
  again and again for at least `<seconds>` of CPU time. Its options scale or
  shift the default parameters: `--fb-mul`, `--hw-mul`, `--threshold-delta`,
  `--large-primes`, `--nrelations`.
- `factor-bench <algo> <bits> <pbits> [<count> [<seconds>]]` times
  `tifa_factor()` (`factor`, or `factor-balanced` and `factor-unbalanced`
  with a hint), ECM (`ecm`), p-1 (`pm1`), p+1 (`pp1`) or SQUFOF (`squfof`) on
  numbers n = p.q, p of `<pbits>` bits (random composites if `<pbits>` is 0),
  and counts the runs that found no factor.
- `siqs-check <lo> <hi> <count>` factors semiprimes of random sizes with SIQS
  and counts the incomplete and wrong results; the test suite runs it from 50
  to 120 bits.
- `bench/ab.sh <base> <new> <args...>` compares two builds, or two sets of
  options: it runs both alternately over four rounds, each order twice, and
  prints the ratio of their times.

```sh
cmake -S . -B new -DCMAKE_BUILD_TYPE=Release -DTIFA_BUILD_BENCH=ON
cmake --build new
bench/ab.sh old/bench/siqs-bench new/bench/siqs-bench 70 100 3
```

where `old` is a build of the reference version (from a `git worktree`, say).
The older tables of news.md use 100 numbers per size up to 120 bits (3
seconds), 30 up to 150 (6 seconds) and 8 above (12 seconds), on one thread. A
build compared with itself stays within 1% of 1 on an idle machine, a few
percent under load: smaller differences are noise.

Since 0.4.3, the `siqs()` tables time each number with both builds instead,
one after the other (`siqs-bench --seed <s> <bits> 1 <seconds>`, the order
alternating from one number to the next), and give the geometric mean of the
ratios over 8 to 16 numbers per size: less noisy from 160 bits, where a
number takes long enough to time alone. Those of x86-64 come from one-off
workflows on GitHub's runners, whose processors vary from one job to the
next.

## Licence

The TIFA library is Copyright (C) 2011 CNRS, Ecole Polytechnique and INRIA.
The files added in 2026 are Copyright (C) 2026 Jerome Milan.

The TIFA library is free software; you can redistribute it and/or modify it
under the terms of the
[GNU Lesser General Public License](https://www.gnu.org/licenses/lgpl.html)
as published by the Free Software Foundation; either
[version 2.1](https://www.gnu.org/licenses/old-licenses/lgpl-2.1.html)
of the License, or (at your option) any later version.

The TIFA library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

The texts of the licenses are in the `readme` directory:
[LGPL 2.1](./readme/license_v2.1.txt), [LGPL 3](./readme/license_v3.txt) and
the [GPL 3](./readme/gpl_v3.txt) it refers to.

The TIFA "puzzle cube" logo was created by Jerome Milan and is
Copyright (C) 2011 [CNRS](https://www.cnrs.fr),
[INRIA](https://www.inria.fr/en) and
[Ecole Polytechnique](https://www.polytechnique.edu/en). All rights reserved.

Finally, note that TIFA has been registered in June 2011 at the French
agency for software protection ([APP](https://www.app.asso.fr/en)) with the
Inter Deposit Digital Number `IDDN.FR.001.220019.000.S.A.2011.000.31235`.
