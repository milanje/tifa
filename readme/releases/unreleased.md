# TIFA, unreleased

The changes on master since the last release, newest first. At the next
release, this file becomes that version's notes. [All releases](../news.md).

## 8 November 2026 &mdash; ECM's curve counts counted again

The curve counts of the levels from 65 to 100 bits stay. Counted again on
155 to 340 factors per level (1951 in all), the means over each level's
range are 8 to 27% below them from 65 to 80 bits, 3% above at 85 and 6
to 12% below from 90 to 100; ECM alone with the means instead found 0.8%
fewer factors of 60 to 85 bits (on 400 numbers), in 3% less time.

## 8 November 2026 &mdash; ECM: the plans above 85 bits freed after use

The plans of the levels for factors of 90 to 100 bits, 57 MB, are freed
when the last search using them ends, rather than kept for the life of the
process; a search that reaches those levels builds them again (0.7 s in
all, against minutes of curves).

## 8 November 2026 &mdash; ECM in GMP arithmetic: a lean stage 2

Stage 2 in GMP arithmetic (above 512 bits, and without 128-bit integers)
kept a table of two bytes per unit of B2, 135 MB at 2^26, which capped B2
there: the levels for 95 and 100 bits ran with a B2 of 29 and 15 times
their B1, rather than 50, and found their factors less often. It now lists
the primes of each giant step in a byte each (the Montgomery plans take
two), 20 MB at 2^28, with the same products in the same order, and B2
goes up to 2^28 there too.

## 8 November 2026 &mdash; ECM up to factors of 100 bits

`ecm_max_bits` goes up to 100 (`TIFA_ECM_TOP_BITS`), rather than 85. ECM's
stage 1 took its primes from the table of `first_primes.h`, which ends at
821641, and so stopped there; it now goes on with a segmented sieve, as
stage 2 did. Three levels, for factors of 90, 95 and 100 bits, double B1
each time (up to 6553600 in the search of `tifa_factor()`, 4608000 in
`ecm()`), with the curve counts of Dickman's estimates, scaled to those
counted on 219 factors of 85 to 100 bits. The default stays 85 bits: each
5 bits more double the time or more. On an Apple M5 Pro, in 6 threads, a
400-bit RSA modulus gives `FACTORS_TOO_LARGE` after 9 s at 85 bits and
113 s at 100; a 400-bit number with a 95-bit factor, `FACTORS_TOO_LARGE`
after 9 s by default, and its factor in 14 s at 100 bits.

`ecm()` lowers the caller's `b1` to 2^24 rather than to 821641 (its stage
1 exponent takes 1.44 bits per unit of `b1`), and `b2` to 2^28 rather than
2^26. `p_minus_1()` and `p_plus_1()` keep their bounds. The curves
beyond the levels for 85-bit factors, which take up to seconds on any
number of limbs, check the time limit and the interrupt within themselves,
as those on five to eight limbs did. The plans build their stage 2 tables
directly in 16-bit entries, without a 64-bit copy.

Below `ecm_max_bits` of 85, the calls cost what they did: against 0.5.0,
in one thread on an Apple M5 Pro, 0.995 to 1.002 of the time for SIQS,
`tifa_factor()` and `tifa_cofactor()` (`bench/ab.sh`, both orders), with
at most 0.1% more instructions.
