# TIFA, unreleased

The changes on master since the last release, newest first. At the next
release, this file becomes that version's notes. [All releases](../news.md).

## 8 November 2026 &mdash; A crash of 0.4.10 fixed

0.4.10 made `tifa_factor()` keep the root of a perfect power when stopped,
but read `multis` there in every mode: with `multis` NULL, as any mode but
`FIND_COMPLETE_FACTORIZATION` allows, it crashed on the powers of
composites without a prime factor below 1000 (such as (1000003 &middot;
1000033)^2). The Python package, which always asks for a complete
factorization, was not affected.
