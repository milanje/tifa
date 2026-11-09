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
