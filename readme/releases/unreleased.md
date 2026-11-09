# TIFA, unreleased

The changes on master since the last release, newest first. At the next
release, this file becomes that version's notes. [All releases](../news.md).

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
