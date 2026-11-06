# TIFA, unreleased

The changes on master since the last release, newest first. At the next
release, this file becomes that version's notes. [All releases](../news.md).

## 6 November 2026 &mdash; ECM above 512 bits

Above 512 bits, where the Montgomery ECM stops (eight limbs), `tifa_factor()`
ran single p-1 and p+1 runs, then SIQS, which takes years there, whatever the
size of the factors. After those runs, it now searches with `ecm()` in GMP
arithmetic, up its levels to factors of 85 bits, without a budget, then runs
SIQS as before. On an Apple M5 Pro, with n = p.q, q prime:

| `tifa_factor()`, bits of n, of p | 520, 40 | 600, 50 | 1024, 40 |
|---|---:|---:|---:|
| before | endless | endless | endless |
| after | 31 ms | 0.57 s | 92 ms |

A curve takes up to 2.1 s at 520 bits and 5.2 s at 1024 (B1 = 576000), and
when nothing is found the search takes about 12 minutes at 520 bits and 30 at
1024, from the time of its curves. The stage 2 table of the GMP curves came
from a primality test of each number up to B2, 0.8 s at the largest B2; it
now comes from a sieve, in milliseconds.

The time limit and the interrupt check (Ctrl-C in Python) stopped a
factorization only between curves, which take up to 0.7 s on eight limbs
since ECM went up to 512 bits. The curves in GMP arithmetic and on five to
eight limbs now check them every 4096 bits of the stage 1 product, and
every 64 (GMP) or 256 (Montgomery) giant steps of stage 2: `ecm()` on 600
bits stops 6 ms after its limit, inside a curve of 2 s. The kernels on one
to four limbs are compiled without these checks.

The kernels' loops are now aligned on 64 bytes with Clang on arm64: without
that, the 412 bytes of code added before them took p+1 at 160 bits to 1.04
of its time on the M5, with the same instructions. Against the code before,
on the M5: SIQS 0.99 to 1.00 from 70 to 140 bits, `tifa_factor()` 0.99 to
1.00 from 70 to 300 bits, `ecm()` 1.00 from 100 to 400 bits, p-1 and p+1 at
160 bits 0.99 to 1.03 (noise: the same instructions, and cycles 0.98 of
those before on average).

## 6 November 2026 &mdash; tifa-batch-gcd

The `tifa-batch-gcd` program runs `tifa_batch_gcd()` on files of integers
(RSA moduli, say), or on its standard input: decimal integers, or
hexadecimal ones with `--hex` (with or without `0x`), separated by spaces or
newlines, `#` starting a comment. It prints the position and the gcd of each
integer that shares a prime with another (of every one with `--all`), and
exits with status 0 if some do, 2 if none do, 3 for invalid input. 20,000
integers of 2048 bits take 3.9 s and 137 MB on an Apple M5 Pro, the reading
of the integers a small part of it.

## 6 November 2026 &mdash; Cofactorization in Python

`tifa.Cofactorizer` gives the Python package `tifa_cofactor()`, the
cofactorization of sieves: built once from the bounds (`lpb`, `fb_bound`,
`max_large_primes`, `max_siqs_bits`), it answers `{prime: exponent}` for an
integer whose prime factors are all below 2^lpb, `None` otherwise, for one
integer or, with `many`, a list at once, in C with the GIL released. On an
Apple M5 Pro, on 3000 integers of 33 to 64 bits without prime factors below
2^23 (lpb 32, two large primes), the C function takes 1.48 us per integer,
`many` 1.53 us and calls one by one 1.92 us. Integers that fit in a long now
cross between Python and C without a hexadecimal string.

## 6 November 2026 &mdash; Small sizes as fast as in 0.4.6

0.4.7 took 1% to 2% more time than 0.4.6 in `tifa_factor()` at 80 bits and
about 1% more in SIQS at 90, which looked like a matter of code layout. The
instructions retired (`/usr/bin/time -l` on fixed runs) said otherwise: 1.3%
more at 80 bits, 0.6% at 90. The sparse decompositions of the relations
(for factor bases beyond 16384 primes) tested which kind they were for each
exponent added, even with the dense ones of small sizes, and the coprime
base of the factors found allocated its stack before its early exits, at
each of the four calls of a two-factor result. The trial division of the
relations is now compiled once per kind, and the early exits come first:
0.16% more instructions than 0.4.6 at 70 bits, 0.05% at 90, and against
0.4.7, `tifa_factor()` takes 0.98 to 0.99 of its time at 70 and 80 bits on an
Apple M5 Pro, the rest unchanged.

## 6 November 2026 &mdash; Batch GCD

`tifa_batch_gcd()` (`batch_gcd.h`) gives the gcd of each integer of a set
with the product of the others, and the Python package's
`tifa.batch_gcd()` the same: the RSA moduli that share a prime with
another, as Heninger, Durumeric, Wustrow and Halderman found among the keys
of the Internet ("Mining Your Ps and Qs", USENIX Security 2012). A product
tree gives the product P of the moduli; then, down the tree, each node N
gives (P / N) mod N to its children, until gcd(n, (P / n) mod n) at each
modulus n. Heninger et al.'s P mod N^2, computed the same way, took 1.1 times
as long and 1.5 times the memory. On an Apple M5 Pro, one thread:

| moduli | 10,000 of 1024 bits | 10,000 of 2048 | 100,000 of 1024 | 100,000 of 2048 |
|---|---:|---:|---:|---:|
| time | 0.71 s | 1.6 s | 12.5 s | 30 s |
| memory | 30 MB | 73 MB | 0.39 GB | 0.77 GB |

With P mod N^2, 87% of the time went to the divisions of the remainder
tree; Bernstein's scaled remainder tree, which multiplies instead, took 0.77
to 0.91 of its time with GMP's full products (it wants middle products,
which GMP does not export), no better than (P / N) mod N. The time limit and the interrupt check of the
calling thread stop it (Ctrl-C in Python). The code is last in the
library, and the benchmark programs are the same, byte for byte.

## 6 November 2026 &mdash; ECM up to 512 bits

ECM, p-1 and p+1 in Montgomery arithmetic (`mont_ecm.h`) now take moduli of
five to eight 64-bit limbs, up to 512 bits: the kernels of one to four
limbs, compiled for each size, in a file of their own at the end of the
library. Above 256 bits, `tifa_factor()` ran single p-1 and p+1 runs, then
SIQS, whatever the size of the factors: 31 s at 264 bits, 98 s at 280,
minutes to hours above, without end from about 330 bits. It now searches with ECM
first, as below 256 bits, for factors of up to 85 bits rather than 60, for
at most 0.2 of the time SIQS would take (version 0.4.7's from 264 bits). On
an Apple M5 Pro, with n = p.q, q prime:

| `tifa_factor()`, bits of n, of p | 264, 40 | 280, 55 | 300, 65 | 400, 50 | 500, 60 |
|---|---:|---:|---:|---:|---:|
| before | 31 s | 98 s | minutes | no end | no end |
| now | 11 ms | 0.14 s | 0.44 s | 50 ms | 0.69 s |

The new levels double B1 from one to the next (51200 to 819200), their
curve counts measured on 20 to 80 factors each; their plans, 7 to 70 ms
each to build (about one curve), are built when a search first reaches
them. `ecm()`, `p_minus_1()`, `p_plus_1()` and the proofs of
`tifa_prove_prime()` use the new arithmetic too: `ecm()` finds 40-bit
factors of 300 to 500 bits 3 to 8 times as fast, and the proofs of five
288-bit primes took 4.2 s on average rather than 24.5 s.

Below 257 bits nothing changes but the layout of the code. Time of the new
version over that of the previous one, `bench/ab.sh`: SIQS 0.995 to 1.001
from 70 to 220 bits, `tifa_factor()` 0.994 to 1.003 from 80 to 140 bits
(0.98 to 1.005 with a factor of 40 to 50 bits, at 160 and 200), `ecm()`
0.990 to 1.006 from 100 to 220 bits, p-1 and p+1 0.975 and 0.995 at 160.

## 5 November 2026 &mdash; Python: Ctrl-C

In the main thread, Ctrl-C interrupts `tifa.factor()` and
`tifa.prime_certificate()`: they raised `KeyboardInterrupt` only once done,
a minute later and more from 272 bits. Now, where the factorizations check their
time limit (every 16 chunks of SIQS's sieve, after each curve of ECM, every
256 columns of the elimination), they also call a function that
`tifa_set_interrupt_check()` (`time_limit.h`) sets, for the calling thread;
once it returns true, they give up as at the time limit, and
`tifa_interrupted()` tells so. The Python package's function takes the GIL
at most every 50 ms to run the signal handlers (`PyErr_CheckSignals()`), which
Python runs in the main thread only: in other threads, a `timeout` bounds a
factorization. On an Apple M5 Pro, `KeyboardInterrupt` (or the exception
of another handler) came 7 to 67 ms after the signal, from 240 to 320 bits,
in SIQS as in ECM. The signals are first polled 50 ms into a call, so that
shorter ones never take the GIL: a call takes as long as before, but for
70 ns at 64 bits (0.3%), and the C library 0.99 to 1.00 of its time. The ABI
only grows.
