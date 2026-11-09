# TIFA, unreleased

The changes on master since the last release, newest first. At the next
release, this file becomes that version's notes. [All releases](../news.md).

## 9 November 2026 &mdash; Every algorithm stops at the time limit

`rho()`, `squfof()`, `fermat()`, and the stages of `p_minus_1()` and
`p_plus_1()` above 256 bits never checked the time limit or the interrupt of
their options, nor did ECM while building the precomputations of large
bounds (the plans of the Montgomery curves, up to 0.43 s at B1 = 2^24 and
B2 = 2^28, which `tifa_factor()` builds when its search first reaches a
level above 60 bits, and the tables of stage 2 in GMP arithmetic). Under a
limit of 50 ms on an Apple M5 Pro, `rho()` with `max_niters` = 2^30 at 400
bits never stopped, `fermat()` ran for minutes at 120 bits, `squfof()` for
more than a minute at 122 bits, `p_plus_1()` at 1000 bits ended 0.62 s late
(0.05 s at 500 bits), `p_minus_1()` 0.19 s, `ecm()` with large bounds 0.14
to 0.22 s. They all stop within 3 ms of the limit now (below 256 bits, the
longest runs of p-1 and p+1 take 40 ms), at no measurable cost
(`bench/ab.sh`: `squfof()`, `p_minus_1()`, `p_plus_1()`, `ecm()` and
`tifa_factor()` 0.976 to 1.002 of the time), and `test_time_limit` checks
each.

The documentation now says what `fb_bound` of `tifa_cofactor()` is trusted
for (parts below its square are taken as prime without a test: a bound
below the true one is safe, one above it gives wrong answers), and that a
CFRAC factor base of very few primes makes the time grow exponentially
(two primes take a second at 60 bits). No ABI change.

## 9 November 2026 &mdash; ECM: Edwards stage 1 on 5 to 8 limbs

The curves of the table and those of the Z/12 family now run their stage 1
on the Edwards form on moduli of 5 to 8 limbs too (257 to 512 bits), long
curves included: `mont_ecm_edwards_wide.c` (last in `TIFA_SOURCES`)
instantiates it for those sizes, with the stop checks of the long curves.
On 1 to 4 limbs, the long curves of the levels from 90 bits keep the
ladder, and Suyama's curves keep it everywhere. In one thread on an Apple
M5 Pro, on numbers with a prime of 35 bits, a curve of the table at B1 =
3200 takes 0.905 of the time on 5 limbs, 0.84 on 6, 0.865 on 7 and 8. With
ECM alone, finding a factor of 45 and 55 bits (150 numbers, two seeds)
takes 0.87 and 0.84 of the time on 5 limbs, 0.89 and 0.80 on 6, 0.82 and
0.69 on 8; a search that finds nothing up to factors of 70 bits, on a
product of two primes, 0.88 at 300 bits and 0.84 at 400. Against the
previous master (`bench/ab.sh`, both orders), SIQS, `tifa_factor()` from 70
to 180 bits and the cofactor bench do not change (instructions within
1.2%). The library grows by 46 KB (745 KB of text on the M5); the point
operations are inlined on three to five limbs only, 3% faster on five, the
same from six.

## 9 November 2026 &mdash; ECM: Montgomery's curves with torsion Z/12

`z2z8_curves.h` holds 27 curves, and the search for more finds few: with
seven times as many candidates (`tools/z2z8`), 32. After them, the curves of
the searches for factors of 45 bits and more were Suyama's, on the ladder.
From B1 = 1600, on 1 to 4 limbs, they are now those of Montgomery's family
with torsion Z/12 (his thesis, 6.1, as in CADO-NFS): curve k from k.P on
v^2 = u^3 - 12u, computed modulo n, whose Montgomery form has (A + 2)/B =
(4uv/(v^2 + 4u^2))^2, a square, so that its stage 1 runs on an Edwards
form like the curves of the table. Against Suyama's curves on numbers with
a factor of 25 to 40 bits, they find 0.96 to 1.09 as many factors per
curve (±4%), in 0.83 to 0.96 of the time from B1 = 1600; below, on three
and four limbs, the setup (k.P and two inversions) makes them slower, and
Suyama's stay. With ECM alone, on numbers of 150 to 160 bits (a 100-bit
prime times the factor), finding a factor of 50 bits takes 0.84 of the time,
55 bits 0.81, 60 bits 0.93 (150 to 400 numbers, two seeds each). In one
thread (`bench/ab.sh`, both orders), `tifa_factor()` takes 0.85 of the time
on numbers of 160 bits with a factor of 50 bits (a third fewer
instructions); SIQS, `tifa_factor()` on balanced numbers, whose ECM search
stays within the table, and the cofactor bench do not change.

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
18 to 40 bits, 6 more found it, none fewer. In one thread on an Apple M5
Pro, on numbers without small factors (so that stage 2 runs, as it does in
most curves of a search), a curve takes 0.88 of the time on two limbs at B1
= 300 and 0.78 at B1 = 2000, 0.985 on three at B1 = 300, 0.90 on four at
B1 = 6400; below 2^64, `split_u64` on numbers of two primes of 25 to 32
bits 0.94. On the cofactors of CADO-NFS's
siever, `tifa_cofactor()` takes 0.94 to 0.96 of the time at c100 to c160
(two large primes, below 2^64), and 0.97 at c180 (three large primes,
two limbs, mostly SIQS), with the same relations.

Against 0.5.1, in one thread (`bench/ab.sh`, both orders): `tifa_factor()`
0.97 of the time at 70 bits, 0.99 at 100 to 180 bits and on unbalanced
numbers at 160 bits (1.9% fewer instructions there); the cofactor bench of
two large primes 0.87 (6% fewer instructions); SIQS 0.984 to 1.001, the
same instructions within 0.3%. For these curves, `run_curve` runs the
Edwards stage 1 of `mont_ecm_edwards.c` (last in `TIFA_SOURCES`) in place
of the ladder; its stage 2 stays the same, and inlined. The Edwards code
adds 45 KB to the library (its text, 654 KB with 0.5.1, is 699 KB on the
M5): its products are calls, not inlined, and its point operations are
inlined from three limbs only, which was also the fastest choice on two
and four limbs.

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
