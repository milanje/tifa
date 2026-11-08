# TIFA, unreleased

The changes on master since the last release, newest first. At the next
release, this file becomes that version's notes. [All releases](../news.md).

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
after 9 s by default, and its factor in 14 s at 100 bits. The plans of the
three levels, 57 MB, are built when a search first reaches them, and kept.

`ecm()` lowers the caller's `b1` to 2^24 rather than to 821641 (its stage
1 exponent takes 1.44 bits per unit of `b1`), and `b2` to 2^28 rather than
2^26, but in GMP arithmetic, whose stage 2 table takes two bytes per unit
of `b2`. `p_minus_1()` and `p_plus_1()` keep their bounds. The curves
beyond the levels for 85-bit factors, which take up to seconds on any
number of limbs, check the time limit and the interrupt within themselves,
as those on five to eight limbs did. The plans build their stage 2 tables
directly in 16-bit entries, without a 64-bit copy.

Below `ecm_max_bits` of 85, the calls cost what they did: against 0.5.0,
in one thread on an Apple M5 Pro, 0.995 to 1.002 of the time for SIQS,
`tifa_factor()` and `tifa_cofactor()` (`bench/ab.sh`, both orders), with
at most 0.1% more instructions.
