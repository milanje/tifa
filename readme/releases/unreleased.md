# TIFA, unreleased

The changes on master since the last release, newest first. At the next
release, this file becomes that version's notes. [All releases](../news.md).

## 8 November 2026 &mdash; Options per call, outcomes by mode, calls that always end

The interface changes, and with it the ABI: the shared library becomes
`libtifa.so.0.5` at the release.

**Options per call.** The functions that can run long take a
`const tifa_options_t*` last (`tifa_options.h`), NULL for the defaults:
`tifa_factor()`, `tifa_prove_prime()`, `tifa_cofactor()`,
`tifa_batch_gcd()` and the algorithms (`siqs()`, `ecm()`, `cfrac()`,
`rho()`, `squfof()`, `fermat()`, `p_minus_1()`, `p_plus_1()`). All zero
is the defaults, so that a designated initializer gives the rest:
`tifa_options_t options = { .time_limit = 10, .threads = 4 };`. The
options hold for the call, the calls it makes included, and several
threads can share them; a call made from within another (from its
interrupt check, for one) has its own and gives the outer call back its
own. Reserved fields, which must be zero (the functions refuse others with
`FAILURE`), keep room for later options without breaking the ABI. The
thread-local setters, whose state a nested call could wipe, are gone.

**Outcomes by mode.** The outcome follows the mode: complete and partial
factorizations, and multiplicities, in `FIND_COMPLETE_FACTORIZATION`
alone. SIQS no longer answers `COMPLETE_FACTORIZATION_FOUND` in the other
modes when it splits a semiprime whole. `FIND_SOME_PRIME_FACTORS`, which
its own documentation doubted, and `SOME_PRIME_FACTORS_FOUND` are gone.

**Calls that always end.** Composites above the options' `siqs_max_bits`
(320 by default, about 40 minutes of SIQS on an Apple M5 Pro) never go to
SIQS: ECM alone searches them, for factors of up to `ecm_max_bits` bits
(85 by default; each 5 bits more about twice as long, up to 100), and a
composite it does not split is left, with the new outcome
`FACTORS_TOO_LARGE`. `tifa_factor()` and `tifa_prove_prime()` so always
end: a 400-bit balanced semiprime gives `FACTORS_TOO_LARGE` after 9 s in 6
threads, where SIQS ran without end. `ecm()` climbs up to `ecm_max_bits`
too, and the levels of both ECMs now go up to 100 bits (extrapolated).
Below `siqs_max_bits`, nothing changes.

Python keeps its interface, and gains `siqs_max_bits` and `ecm_max_bits`
in `factor()` and `prime_certificate()` (and `--siqs-max-bits`,
`--ecm-max-bits` in `python -m tifa` and `tifa-factor`);
`FactorizationError` says when composites are too large.

Migrating from 0.4:

| 0.4 | 0.5 |
|---|---|
| `tifa_factor(f, m, n, mode)` | `tifa_factor(f, m, n, mode, NULL)` |
| `siqs(f, m, n, &params, mode)` (and the other algorithms) | `siqs(f, m, n, &params, mode, NULL)` |
| `tifa_prove_prime(cert, n)`, `tifa_cofactor(f, m, n, s)`, `tifa_batch_gcd(g, ns)` | the same, with `NULL` (or options) last |
| `tifa_factor_with_hint(f, m, n, mode, hint)` | options with `.hint = hint` |
| `tifa_set_time_limit(seconds)` | options with `.time_limit = seconds` |
| `tifa_set_interrupt_check(check, arg)` | options with `.interrupt = check, .interrupt_arg = arg` |
| `tifa_set_threads(n)` | options with `.threads = n` |
| `tifa_time_limit_reached()`, `tifa_interrupted()` | gone: the caller knows its deadline and its interrupt |
| `time_limit.h` | `tifa_options.h` |
| `FIND_SOME_PRIME_FACTORS` | `FIND_SOME_COPRIME_FACTORS`, then a primality test |
| `COMPLETE_FACTORIZATION_FOUND` outside `FIND_COMPLETE_FACTORIZATION` | the mode's own outcome |

The values of the modes and of the exit codes after them change: compile
again against the new headers.

## 8 November 2026 &mdash; Two fixes for 0.4.10

- 0.4.10 made `tifa_factor()` keep the root of a perfect power when
  stopped, but read `multis` there in every mode: with `multis` NULL, as
  any mode but `FIND_COMPLETE_FACTORIZATION` allows, it crashed on the
  powers of composites without a prime factor below 1000 (such as
  (1000003 &middot; 1000033)^2). The Python package, which always asks for
  a complete factorization, was not affected.
- Stopped by a time limit (since 0.4.7) or an interrupt, an algorithm's
  run that had found nothing returned `NO_FACTOR_FOUND`, which
  `tifa_factor()` documents for primes: on composites of 65 to 95 bits
  under limits of a fraction of a millisecond, 20 times in 600. It returns
  `GIVING_UP` now. The Python package, which checks the time itself, was
  not affected.
