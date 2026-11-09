# TIFA news

The releases of TIFA, newest first: the main changes of each, and a link to
its release notes in [releases/](./releases/), which have all of them, with
their measurements. The changes on master since the last release are in
[releases/unreleased.md](./releases/unreleased.md), which becomes the next
version's notes. Times are on one thread of an Apple M5 Pro unless the notes
say otherwise.

## Unreleased

Nothing yet.

[Release notes](./releases/unreleased.md)

## 0.5.2 &mdash; 2026-11-09

- ECM's stage 1 on Edwards curves, from one word to eight limbs: the
  table of curves with torsion Z/2 x Z/8, then Montgomery's family with
  torsion Z/12 (from B1 = 1600) instead of Suyama's curves; 0.78 to 0.99
  of the time per curve, ECM alone 0.69 to 0.89 on 5 to 8 limbs.
- `tifa_cofactor()` on CADO-NFS's own cofactors: the same relations in
  0.50 to 0.69 of CADO's time.
- Without SIQS, three-large-prime cofactors no longer given up before ECM
  has spent its budget (0.35% of the relations lost at c180, not 8.4%).
- Every algorithm stops at the time limit and the interrupt: rho, SQUFOF,
  Fermat, p-1, p+1 and ECM's large precomputations too.
- Same ABI (`libtifa.so.0.5`).

[Release notes](./releases/0.5.2.md)

## 0.5.1 &mdash; 2026-11-09

- ECM up to factors of 100 bits: `ecm_max_bits` up to 100 (the default
  stays 85); stage 1 goes on past the table of primes.
- ECM in GMP arithmetic: a stage 2 of 20 MB at B2 = 2^28, rather than
  135 MB at 2^26.
- The plans for factors above 85 bits freed after use.
- `ecm()` caps `b1` at 2^24 and `b2` at 2^28.
- Same ABI (`libtifa.so.0.5`).

[Release notes](./releases/0.5.1.md)

## 0.5.0 &mdash; 2026-11-08

- Options per call (`tifa_options.h`) instead of thread-local setters:
  time limit, interrupt, threads, hint, sizes; a call made from within
  another keeps the outer one's.
- The outcome follows the factoring mode; `FIND_SOME_PRIME_FACTORS` gone.
- `tifa_factor()` always ends: composites above 320 bits go to ECM alone,
  and those it does not split give `FACTORS_TOO_LARGE`.
- The headers work from C++.
- Fixed: a crash of 0.4.10 (`multis` NULL, powers of composites), and
  `NO_FACTOR_FOUND` for composites stopped by a time limit.
- ABI changes: `libtifa.so.0.5`.

[Release notes](./releases/0.5.0.md)

## 0.4.10 &mdash; 2026-11-08

- SIQS, and the ECM search before it, in several threads
  (`tifa_set_threads()`, `--threads`, Python's `threads=`): in six
  threads, `tifa_factor()` takes 0.41 of the time at 140 bits, 0.22 at
  240.
- Fixes from a bug hunt over the whole code: block Lanczos without end on
  some CFRAC matrices, CFRAC's `FATAL_INTERNAL_ERROR` when `nprimes_tdiv`
  was below the base, Python calls from a signal handler, the programs'
  exit statuses, and more.
- ABI grows (`libtifa.so.0.4`).

[Release notes](./releases/0.4.10.md)

## 0.4.9 &mdash; 2026-11-07

- Block Lanczos (`BLOCK_LANCZOS`), for SIQS from about 230 bits: at 280
  bits, its linear algebra 0.16 of the time, the whole factorization 0.94,
  its peak memory 83 MB rather than 721.
- The default linear algebra of SIQS follows the factor base of a tuning
  file (`siqs_dflt_linalg_method()`).
- `tifa-cfrac` and `tifa-siqs` take block Lanczos and refuse an unknown
  linear algebra.
- ABI grows (`libtifa.so.0.4`).

[Release notes](./releases/0.4.9.md)

## 0.4.8 &mdash; 2026-11-06

- Montgomery ECM, p-1 and p+1 up to 512 bits (8 limbs); above 256 bits,
  `tifa_factor()` searches for factors of up to 85 bits before SIQS.
- Above 512 bits, ECM in GMP arithmetic before SIQS; the longest curves stop
  within a few ms of a time limit or an interrupt.
- `tifa_batch_gcd()`, `tifa.batch_gcd()` and the `tifa-batch-gcd` program:
  the integers (RSA moduli) sharing a prime with another.
- Python: `tifa.Cofactorizer`, `tifa_cofactor()` for one integer or a list;
  Ctrl-C interrupts `factor()` and `prime_certificate()`.
- `python -m tifa`, as GNU `factor`; the wheels install `tifa-factor` and
  `tifa-batch-gcd`; `--version` for every program.
- `tifa_factor()` 1% to 2% faster at 70 and 80 bits (0.4.6's speed back).
- The ECM search of `tifa_cofactor()` and `tifa_factor()` runs the
  Z/2 x Z/8 curves first: 1% to 2.4% less work on sets of cofactors.
- The programs take numbers of up to 4095 digits (255 before).
- ABI grows (`libtifa.so.0.4`).

[Release notes](./releases/0.4.8.md)

## 0.4.7 &mdash; 2026-11-05

- SIQS 1.1 to 4.1 times as fast from 232 to 280 bits: factor bases of up
  to 65536 primes, sparse exponents of the relations, new parameters.
- Primality certificates: `tifa_prove_prime()`, `tifa_check_prime_cert()`,
  as text both ways; `tifa-factor --prove`.
- `tifa_set_time_limit()`: a time limit per thread.
- `tifa-tune` for large sizes; tuning tables apply near their rows only.
- Python: proofs, `timeout=`, tuning files, `tifa-tune` in the wheels.
- Fixes for untrusted input: certificates, tuning files, factoring core.
- ABI grows (`libtifa.so.0.4`).

[Release notes](./releases/0.4.7.md)

## 0.4.6 &mdash; 2026-11-04

- Python package on PyPI, `pip install tifa-factor`: `tifa.factor(n)`,
  wheels for Linux and macOS.
- SIQS bucket fill with AVX-512: 0.84 to 0.87 of the time at 220 bits.
- ABI unchanged.

[Release notes](./releases/0.4.6.md)

## 0.4.5 &mdash; 2026-11-03

- SIQS 0.94 to 0.98 of the time from 140 to 220 bits (0.96 to 0.99 on
  x86-64): a cheaper candidate pre-check.
- `FILTERED_GAUSS_ELIM`: 0.29 to 0.46 of the time of `SMART_GAUSS_ELIM`
  from 180 to 220 bits.
- The table of primes from a Python script: no Perl left.

[Release notes](./releases/0.4.5.md)

## 0.4.4 &mdash; 2026-11-03

- SIQS 1.25 to 1.35 times as fast at 200 and 220 bits, M5 and x86-64:
  hopeless candidates dropped early, 40 small primes left out of the sieve,
  polynomial switches without divisions, SSE2 scan on x86-64.
- `tifa_factor()` 2% to 18% faster from 80 to 180 bits.

[Release notes](./releases/0.4.4.md)

## 0.4.3 &mdash; 2026-11-02

- SIQS 1.5 to 1.6 times as fast at 200 and 220 bits on x86-64, up to 1.3 on
  the M5: branch-free sieve kept by GCC, SSE2 sieve, Q2 polynomials,
  parameters retuned from 160 to 220 bits.
- `tifa_factor()` 5% to 30% faster from 80 to 180 bits.

[Release notes](./releases/0.4.3.md)

## 0.4.2 &mdash; 2026-11-02

- `tifa_factor()`: ECM before SIQS again from 74 bits.
- `ecm()` keeps the plans of its levels: 2% to 6% faster (4% to 10% on
  x86-64).

[Release notes](./releases/0.4.2.md)

## 0.4.1 &mdash; 2026-11-01

- SQUFOF 4 to 8 times as fast; p+1 12% to 32% faster.
- `tifa_factor()` 5% to 19% faster from 100 to 160 bits, 3 to 10 times on
  numbers with a factor of 26 to 40 bits.
- `tifa_factor_with_hint()`: what is known of the number.
- ABI grows.

[Release notes](./releases/0.4.1.md)

## 0.4.0 &mdash; 2026-10-30

- ECM 19% to 33% faster; SIQS 10% to 14% faster from 110 to 130 bits.
- `tifa_factor()` 7% to 11% faster from 128 to 160 bits.
- ABI changes (`libtifa.so.0.4`): the public API alone is exported and
  installed.

[Release notes](./releases/0.4.0.md)

## 0.3.4 &mdash; 2026-10-30

- SIQS 15% to 20% faster from 135 to 190 bits on x86-64 (32 KiB sieve
  chunks), 4% to 8% from 150 to 200 bits on the M5.

[Release notes](./releases/0.3.4.md)

## 0.3.3 &mdash; 2026-10-28

- SIQS 14% to 30% faster from 66 to 200 bits: smoothness tests eight at a
  time, a sieve without unpredictable branches.
- Benchmark programs in `bench/`.

[Release notes](./releases/0.3.3.md)

## 0.3.2 &mdash; 2026-10-28

- SIQS 11% to 32% faster from 66 to 110 bits: arrays kept between runs,
  fewer GMP calls.
- `tifa_factor()` 12% to 26% faster from 66 to 100 bits.

[Release notes](./releases/0.3.2.md)

## 0.3.1 &mdash; 2026-10-28

- SIQS 18% to 31% faster from 66 to 100 bits: less time outside the sieve,
  8 extra relations instead of 24.
- The programs take `-h` and `--help`.

[Release notes](./releases/0.3.1.md)

## 0.3.0 &mdash; 2026-10-27

- SIQS 7% to 20% faster from 66 to 110 bits.
- Fixes: perfect powers, hangs of `cfrac()` and `squfof()`, Fermat's square
  roots, ECM and p±1 on small composites.
- Programs named `tifa-<name>`, exit statuses 0 to 3.
- ABI changes (`libtifa.so.0.3`).

[Release notes](./releases/0.3.0.md)

## 0.2.0 &mdash; 2026-10-22

- Rework for speed, correctness and thread safety: `tifa_factor()` 6 to 100
  times as fast from 32 to 110 bits.
- `tifa_cofactor()`: the cofactorization of sieves.
- ABI changes (`libtifa.so.0.2`); GMP 6.0 or later.

[Release notes](./releases/0.2.0.md)

## 2011-06-10

TIFA has been registered in June 2011 at the French agency for software
protection ([APP](https://www.app.asso.fr/en)) with the Inter Deposit Digital
Number `IDDN.FR.001.220019.000.S.A.2011.000.31235`.
