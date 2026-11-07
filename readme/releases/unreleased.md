# TIFA, unreleased

The changes on master since the last release, newest first. At the next
release, this file becomes that version's notes. [All releases](../news.md).

## 7 November 2026 &mdash; SIQS in several threads

`tifa_set_threads(n)` (`time_limit.h`) lets the calling thread's SIQS runs
collect their relations in up to `n` threads, the calling one included, at
most one per 100 primes of the factor base: none below about 100 bits (200
primes), where starting threads takes longer than the run. Each thread
sieves its own polynomials, from its own place among the subsets of primes
of 'a', with its own sieve and candidates; their candidates go through the
batch smoothness test and the large prime table of the run, under a lock.
The time limit and the interrupt check stay the calling thread's: it checks
them and stops the others. `tifa_factor()` and `tifa_prove_prime()` follow
it, `tifa-factor` and `python -m tifa` take `--threads <n>`, and Python's
`factor()` and `prime_certificate()` take `threads=`.

On an Apple M5 Pro (six performance cores), the wall-clock time of SIQS on
balanced semiprimes, against one thread (at 100 bits, 200 primes, two
threads at most):

| bits | 100 | 120 | 140 | 160 | 180 | 200 | 220 | 240 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 2 threads | 0.88 | 0.68 | 0.61 | 0.59 | 0.57 | 0.55 | 0.53 | 0.51 |
| 4 threads | 0.88 | 0.57 | 0.43 | 0.38 | 0.35 | 0.33 | 0.31 | 0.28 |
| 6 threads | 0.92 | 0.57 | 0.42 | 0.33 | 0.29 | 0.26 | 0.24 | 0.23 |

The threads take about 8% more CPU time in all. With one thread, the
default, SIQS runs as before. Against 0.4.9: the same instructions within
0.05% from 100 to 180 bits, and for `tifa_factor()` on random composites
within 0.02% (2000 per size from 100 to 180 bits, three seeds); 0.99 to
1.01 of the time for SIQS from 70 to 220 bits, `tifa_factor()` from 70 to
100 bits, `ecm()` and p-1 at 160 bits, in both orders (`bench/ab.sh`).
On x86-64 (GitHub's AMD EPYC 9V45 and Intel Xeon Platinum 8573C runners,
GCC 15), the same runs took 0.99 to 1.01 of the time as well.

Six SIQS runs at once, one per thread, each took as long as one alone (0.42
to 0.45 s at 200 bits, 1.58 to 1.70 s at 220): the shared L2 cache does
not slow them. To factor many numbers, as a sieve's cofactorization does,
one factorization per thread stays the better use of the cores.
