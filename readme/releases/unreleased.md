# TIFA, unreleased

The changes on master since the last release, newest first. At the next
release, this file becomes that version's notes. [All releases](../news.md).

## 9 November 2026 &mdash; ECM: stage 1 on Edwards curves

The curves of `z2z8_curves.h`, with which every ECM search starts, are
Edwards curves x^2 + y^2 = 1 + d.x^2.y^2; the table now holds their x0 too,
so that their stage 1 runs on that form rather than on the Montgomery
ladder: a doubling costs 7 products against 10 for a step of the ladder,
and a window over E (signed, recoded once per plan, below 2^64; sliding,
of 2 to 8 bits by the size of E, from 65 to 256 bits) adds a multiple of
the starting point every few doublings (Hisil, Wong, Carter and Dawson's
formulas). The point found then passes to the Montgomery form of the same
group, u = (1 + y)/(1 - y), for stage 2, which does not change. Suyama's
curves, the long curves of the levels from 90 bits, and moduli of 5 to 8
limbs keep the ladder.

Each curve reaches the same point as before, and finds the same factors,
but for an exceptional point of the addition law (d is a square): meeting
one modulo a prime factor p derails the computation modulo p alone, and in
every case seen it revealed p: of 270000 curves on numbers with a prime of
18 to 40 bits, 6 more found it, none fewer. In one thread on an Apple M5 Pro, a curve takes
0.83 of the time on two limbs and 0.91 on three at B1 = 300, 0.75 on two at
B1 = 2000, 0.81 on four at B1 = 6400; below 2^64, `split_u64` on numbers
of two primes of 25 to 32 bits 0.94. On the cofactors of CADO-NFS's
siever, `tifa_cofactor()` takes 0.94 to 0.96 of the time at c100 to c160
(two large primes, below 2^64), and 0.97 at c180 (three large primes,
two limbs, mostly SIQS), with the same relations.

Against 0.5.1, in one thread (`bench/ab.sh`, both orders): `tifa_factor()`
0.97 of the time at 70 bits, 0.99 at 100 to 180 bits and on unbalanced
numbers at 160 bits (1.9% fewer instructions there); the cofactor bench of
two large primes 0.87 (6% fewer instructions); SIQS 0.984 to 1.001, the
same instructions within 0.3%. The code of `mont_ecm.c` does not change
but for the choice of the stage 1: the Edwards one is in
`mont_ecm_edwards.c`, last in `TIFA_SOURCES`.

## 9 November 2026 &mdash; Cofactorization on a sieve's own cofactors

`tifa_cofactor()` had only been measured on synthetic cofactors. On the
pairs that CADO-NFS's siever hands to its own cofactorization (CADO-NFS's
polynomials and parameters, its `las` patched locally to write each pair out
with CADO's verdict), in one thread on an Apple M5 Pro, it finds the same
relations as CADO, in 0.50 to 0.69 of the time:

| Set | Pairs | Relations | TIFA, per pair | CADO-NFS, per pair |
|---|---:|---:|---:|---:|
| c100 (lpb 25 and 26) | 16910 | 12892 | 3.77 µs | 7.55 µs |
| c120 (27 and 28) | 4460 | 3221 | 5.91 µs | 10.23 µs |
| c140 (30 and 30) | 15946 | 11417 | 8.90 µs | 14.30 µs |
| c160 (31 and 32) | 5249 | 3958 | 10.97 µs | 17.92 µs |
| c180 (31 and 32, three large primes on side 1) | 117979 | 1717 | 195 µs | 283 µs |

TIFA answers exactly, side 1 first and side 0 only when side 1 is smooth;
CADO's time is its siever's own timer of its cofactorization. At c180 CADO
leaves 26196 pairs undecided; TIFA decides them, and none is a relation.

With `max_siqs_bits = 2 lpb`, as the documentation suggested for three large
primes, 145 of the 1717 relations at c180 were lost (8.4%), where synthetic
sets had lost 0.25% to 0.30% of their smooth numbers. On a sieve's cofactors
the large primes all lie near the large-prime bound, and the first ECM level
within reach, for factors of 30 bits, took longer than SIQS would on an 82-
to 88-bit cofactor: no curve ran before giving up. Without SIQS, ECM now
spends all of that time, the last level cut to fit: 6 relations lost
(0.35%), in 0.62 of the time of exact answers (121 µs a pair, 0.43 of
CADO's). On synthetic cofactors (sets B and C of 0.2.0, 20000 and 12000
numbers), giving up now takes 0.67 and 0.68 of the time of exact answers,
rather than 0.48 and 0.39, and misses 0 of 785 and 3 of 575 smooth numbers,
rather than 3 and 17: close to what 0.2.0 measured on smaller sets (0.60
and 0.64, none missed). Exact answers do not change. Against 0.5.1, in one
thread (`bench/ab.sh`, both orders): SIQS 0.986 to 1.002 of the time,
`tifa_factor()` 0.985 to 1.008, the same instructions within 0.25% (1.4%
fewer for `tifa_factor()` at 70 bits); the cofactor bench of two large
primes, which this change does not reach, 0.90 to 0.93 with the same
instructions, from where the code now lies.
