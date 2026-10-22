# TIFA news

## 22 October 2026 &mdash; command-line tools and documentation

A `tifa_factor` program factors a number with the methods the `tifa_factor()`
function chooses, and the `fermat` and `tdiv` programs now take the number
alone, as the others do. The Perl scripts of `tools/scripts` (a wrapper for the
programs, benchmark and plotting scripts, and their modules) are removed: they
had not worked as installed for a long time, and `tifa_tune` measures the
algorithms' parameters instead.

The documentation (the readme files and the headers) was brought in line with
the code. The unused `SIQS_DFLT_*` and `CFRAC_DFLT_*` macros for the sizes of
the factor base, of the trial division and of the sieve are removed: these
defaults depend on the size of the number.

## October 2026 &mdash; faster, correct, thread-safe

A series of 56 commits (`aac2494..6e6deaf`) reworked most of the library for
speed, correctness and use inside a sieve (one factorization per thread, as in
the cofactorization step of the number field sieve). The tables below compare
the code before (`aac2494`) and after (`6e6deaf`) on the same inputs.

All times are means per call on one thread of an Apple M5 Pro (Apple clang 21,
GMP 6.3, release builds, `-O3`). "Failed" counts the outputs that are not a
complete factorization into primes with multiplicities.

### `tifa_factor`

Balanced semiprimes:

| Bits | Before | After | Speed-up | Failed before |
|-----:|-------:|------:|---------:|--------------:|
|   32 |   10.7 µs |    1.80 µs |   5.9x | 0 / 300 |
|   40 |   34.6 µs |    3.84 µs |   9.0x | 0 / 300 |
|   48 |   75.9 µs |    5.56 µs |  13.7x | 0 / 300 |
|   56 |  249 µs   |    9.22 µs |  27x   | 0 / 300 |
|   64 |  856 µs   |   20.9 µs  |  41x   | 0 / 300 |
|   72 |  17.1 ms  |  166 µs    | 103x   | 133 / 300 |
|   80 |  23.0 ms  |  289 µs    |  80x   | 196 / 300 |
|   90 |  24.8 ms  |  511 µs    |  48x   | 30 / 40 |
|  100 |  40.4 ms  |  1.05 ms   |  38x   | 38 / 40 |
|  110 |  37.8 ms  |  1.94 ms   |  19x   | 36 / 40 |
|  128 |  44.7 ms  |  5.53 ms   |  8.1x  | 12 / 12 |
|  140 |  68.5 ms  |  12.1 ms   |  5.7x  | 12 / 12 |
|  160 |   139 ms  |  62.2 ms   |  2.2x  | 4 / 4 |
|  180 |   313 ms  |   173 ms   |  1.8x  | 2 / 2 |
|  200 |  1.39 s   |   963 ms   |  1.4x  | 2 / 2 |

Other inputs (the "after" code had no failure on any of them):

| Input | Bits | Before | After | Speed-up | Failed before |
|-------|-----:|-------:|------:|---------:|--------------:|
| Composites without factors below 1000 |  64 |   847 µs |  9.95 µs |  85x   | 0 / 200 |
|                                        |  80 |  6.06 ms |  99.3 µs |  61x   | 30 / 200 |
|                                        | 100 |  12.2 ms |   339 µs |  36x   | 10 / 40 |
|                                        | 128 |  11.9 ms |  1.27 ms |  9.4x  | 3 / 12 |
|                                        | 160 |  45.3 ms |  11.9 ms |  3.8x  | 1 / 4 |
| Random composites                      |  64 |  17.2 ms |   3.3 µs | 5221x  | 11 / 40 |
|                                        | 100 |  9.41 ms |   174 µs |  54x   | 9 / 40 |
|                                        | 128 |  17.7 ms |  1.01 ms |  17x   | 16 / 40 |
| Primes                                 |  64 |   383 ms |   1.7 µs | about 230000x | 0 / 3 |
|                                        | 100 |  40.3 ms |   8.0 µs | 5036x  | 0 / 3 |
|                                        | 128 |  46.0 ms |  11.3 µs | 4056x  | 0 / 3 |
|                                        | 160 |   128 ms |  17.7 µs | 7260x  | 0 / 3 |

### The algorithms

| Algorithm | Input | Before | After | Speed-up |
|-----------|-------|-------:|------:|---------:|
| SIQS | 80-bit semiprimes  |   492 µs |   262 µs | 1.9x |
|      | 100-bit semiprimes |  1.39 ms |   902 µs | 1.5x |
|      | 128-bit semiprimes |  7.83 ms |  5.27 ms | 1.5x |
|      | 160-bit semiprimes |  83.6 ms |  57.9 ms | 1.4x |
|      | 200-bit semiprimes |  1.32 s  |   950 ms | 1.4x |
| ECM, default bounds | 25-bit factor, 90-bit n  |   875 µs |  54.4 µs |  16x |
|                     | 30-bit factor, 100-bit n |  2.37 ms |   168 µs |  14x |
|                     | 35-bit factor, 120-bit n |  9.96 ms |   499 µs |  20x |
|                     | 40-bit factor, 128-bit n |  21.6 ms |  1.47 ms |  15x |
|                     | 40-bit factor, 190-bit n |  3.37 s  |  1.82 ms | 1851x |
| p-1, B1 = 10000 | 128 bits | 2.45 ms |   417 µs | 5.9x |
|                 | 200 bits | 4.92 ms |   794 µs | 6.2x |
| p+1, B1 = 10000 | 128 bits | 2.24 ms |   464 µs | 4.8x |
|                 | 200 bits | 3.80 ms |  1.03 ms | 3.7x |

The SIQS "before" outputs all failed the check: `siqs()` did not fill in the
multiplicities.

### Cofactorization

Sets of integers as a sieve leaves them: windows of consecutive integers
without prime factors below the factor base bound, of lpb + 1 bits up to the
cofactor bound, each with its true answer. A cofactor is smooth when its prime
factors are all below 2^lpb, with at most the given number of large primes.
Before, there was no cofactorization function: the test is a size and
primality filter, then `tifa_factor`. After, `tifa_cofactor` answers exactly
(`max_siqs_bits = 0`), or gives up on the composites above 2 lpb bits that ECM
does not split within the time SIQS would take (`max_siqs_bits = 2 lpb`).

| Set | Before | After, same filter | `tifa_cofactor`, exact | `tifa_cofactor`, giving up |
|-----|-------:|-------------------:|-----------------------:|---------------------------:|
| lpb 32, 2 large primes, base bound 2^23, 33-64 bits (3000, 176 smooth)  |  46.5 µs | 2.80 µs | 1.25 µs (37x) | 1.33 µs (35x) |
| lpb 31, 3 large primes, base bound 2^22, 32-92 bits (1500, 72 smooth)   | 4.02 ms  | 63.8 µs | 56.1 µs (72x) | 33.4 µs (120x) |
| lpb 36, 3 large primes, base bound 2^25, 38-108 bits (800, 27 smooth)   | 6.06 ms  |  157 µs |  133 µs (46x) | 85.3 µs (71x) |

No answer was wrong and no smooth cofactor was missed on these sets, in either
mode. On sets of smooth cofactors only, the mode giving up misses 0.25% to 0.30%
of them.

### What changed

Correctness:

- `tifa_factor` did no trial division above 26 bits, ran its p-1 and p+1
  pre-stages in the caller's mode, and left the multiplicities unset when an
  algorithm reported a complete factorization.
- SIQS seeded its polynomials from the clock, could use a prime dividing kN in
  'a', gave up on duplicate relations, kept or rejected the last positions of
  each negative-side block on the wrong sieve values, could run forever on
  unlucky parameters (it now detects a stalled collection and retries wider),
  read out of bounds with small factor bases, and wrapped around when
  interpolating decreasing table columns.
- The batch smoothness filter (CFRAC) dropped or retested candidates when its
  output filled up, and cleared integers still in use; ECM could step past the
  prime table; Fermat took its
  split modulo kN instead of n; p-1 failed when stage 1 caught every factor at
  once; `modinv_ui` shifted signed values.

Thread safety: no global state, a locked tuning cache, the caller's parameters
left untouched, per-call random curves, and a test running factorizations
concurrently (`test_threads`). Throughput scales to 14-15x on 18 threads.

New code:

- `factor_u64`: primality (deterministic, the Miller-Rabin bases run in
  lockstep) and factorization (Hart's method, then ECM three curves at a time)
  in machine words below 2^64.
- Montgomery arithmetic on 1 to 4 limbs for ECM, p-1 and p+1, up to 256 bits;
  ECM curves with torsion Z/2 x Z/8, a normalized stage 2, and bounds chosen by
  the size of the factor sought.
- `tifa_cofactor` (`cofactor.h`): the cofactorization of sieves, with a bounded
  ECM search and a strategy shared by the threads.
- Primality: `is_prime_u64` below 2^64, Baillie-PSW alone above.

SIQS: the standard Knuth-Schroeppel multiplier, polynomial switches without
divisions, sieve chunks sized to the L1 data cache, one sieve kernel for all
primes chosen by hits per chunk, a faster scan, native trial division of the
relations, a direct batch smoothness test for small factor bases, a cheaper
setup, and the defaults retuned (the large prime variation from 103 bits).

Parameter tuning (`tifa_tune`): searches both ways around the defaults with
pruning, tuning files interpolated between rows, `--max` always tuned, SIQS
tuned from 65 bits.

## 10 June 2011

TIFA has been registered in June 2011 at the French agency for software
protection ([APP](https://www.app.asso.fr/en)) with the Inter Deposit Digital
Number `IDDN.FR.001.220019.000.S.A.2011.000.31235`.
