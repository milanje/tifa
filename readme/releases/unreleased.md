# TIFA, unreleased

The changes on master since the last release, newest first. At the next
release, this file becomes that version's notes. [All releases](../news.md).

## 8 November 2026 &mdash; Fixes from a bug hunt over the whole code

Five reviews, one per part of the code, and randomized checks against
known factorizations under ASan, UBSan and TSan, with and without 128-bit
integers, found the following.

- Block Lanczos could not drop empty columns (relations that are squares):
  with more than 100 of them, it looped forever, past any time limit. CFRAC
  met them on some numbers. After a breakdown near the end, on small
  matrices, it now keeps the dependencies its check accepts rather than
  losing the attempt.
- CFRAC put the primes the batch test found in a second run's relations
  into the first run's rows, a FATAL_INTERNAL_ERROR whenever `nprimes_tdiv`
  was below the base; with `nprimes_tdiv` 0, the batch test never found 2.
  It now checks the time limit while collecting relations too.
- With threads, the ECM search's other threads ran their curves to the end
  once the calling thread had stopped, and `tifa_factor()` then started
  SIQS anyway: 0.2 to 0.6 s past a time limit at 500 bits, 0.01 s now.
- Stopped on the composite root of a perfect power, `tifa_factor()` gave
  up with nothing; it now returns the power, as for other cofactors.
- With GMP's limbs of 32 bits, SIQS dropped every candidate from 140 bits:
  its pre-check reads 64-bit limbs, and is now off there.
- `tifa_prove_prime()` ran SIQS on composites that ECM had not searched
  (above 512 bits, or without 128-bit integers); its documentation now
  says what a proof can take above 256 bits.
- A tuning name of 32 characters or more, looked up 16 times, turned the
  tuning files off. CFRAC's trial division could clear a relation twice.
- Python: a call made from a signal handler during another call removed
  the outer call's time limit, threads and Ctrl-C. Integers of more than
  4300 digits turned the errors meant into a `ValueError`, and stopped
  `python -m tifa`.
- The algorithms' programs said "not complete" (exit status 1) for
  complete factorizations, and lost the exponent of a power; numbers of
  4095 digits are read from the standard input too. `tifa-tune` reported
  success over an empty file when a write failed; it now writes through a
  temporary file. `tifa-batch-gcd` dropped what followed a NUL byte.

In one thread, against the polish before it: the same instructions within
0.03% for SIQS, `tifa_factor()` and `tifa_cofactor()`, and 0.996 to 1.006
of the time in both orders.

## 7 November 2026 &mdash; Fixes to the threads, and polish

- In several threads, SIQS built the other threads' sieves on its half
  width rounded up to whole chunks, from which their 'a' could take
  another number of primes than the calling thread's: the primes of 'a'
  were then read and written out of bounds (by default, on a few percent
  of the numbers of some sizes from 230 bits). They now take the calling
  thread's sieve width.
- In several threads, a `nrelations` of some hundreds could fill the set of
  the relations seen, whose search then never ended. The set now holds as
  many as the rounds of more relations can ask for.
- The time limit's functions read whether their thread key exists with an
  acquire, as it is written, rather than plainly.
- `siqs-bench` and `factor-bench` take `--threads <n>`, timed by the wall
  clock; the documentation names the ECM search wherever threads are.

## 7 November 2026 &mdash; The ECM search in several threads

With `tifa_set_threads(n)`, the ECM search that `tifa_factor()` and
`tifa_prove_prime()` run before SIQS shares its curves among up to `n`
threads as well, at most one per 0.1 ms of curves (from about 120 bits
without a hint, 100 with the unbalanced one, 135 when proving a prime):
each thread takes the next curve, the first factor found stops them all,
and the calling thread alone checks the time limit and the interrupt. Run
in one thread, the search had come to take as long as SIQS in six: 4.6 s,
then 6.0 s of SIQS, at 264 bits.

On an Apple M5 Pro, the wall-clock time of `tifa_factor()` on balanced
semiprimes, without a hint, in 6 threads against one:

| bits | 100 | 120 | 140 | 160 | 180 | 200 | 220 | 240 | 264 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| SIQS in threads | 0.87 | 0.62 | 0.46 | 0.43 | 0.44 | 0.39 | 0.33 | 0.26 | 0.32 |
| and the ECM search | 0.88 | 0.61 | 0.41 | 0.32 | 0.27 | 0.24 | 0.23 | 0.22 | 0.22 |

The threads take more CPU time in all: 24% more at 140 bits, 12% at 180, 9%
at 200 and 7% at 220.

In one thread, the same instructions within 0.04% for SIQS, `tifa_factor()`
and `tifa_cofactor()` (SIQS at 70 bits aside, whose count varies by 0.2%
from run to run), and 0.99 to 1.01 of the time in both orders.

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
