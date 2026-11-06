# TIFA, unreleased

The changes on master since the last release, newest first. At the next
release, this file becomes that version's notes. [All releases](../news.md).

## 6 November 2026 &mdash; Block Lanczos

The linear algebra has a third method, `BLOCK_LANCZOS` (`methods.h`):
Montgomery's block Lanczos ("A block Lanczos algorithm for finding
dependencies over GF(2)", Eurocrypt 1995), after the presentation of
msieve's `lanczos.c`, on blocks of 64 vectors, on a sparse matrix. Columns
in no dependency go first (rows of their own), then the heaviest ones
beyond an excess of a hundred over the rows used; the dependencies found,
up to 64, are each checked. `find_dependencies()` takes it for any dense
matrix of relations. SIQS uses it from 12000 primes in the base (about 230
bits), with no dense matrix at all: the relations' exponents, kept sparse,
give the sparse matrix. On an Apple M5 Pro, one balanced semiprime per
size:

| bits, primes in the base | 220, 8520 | 240, 16143 | 256, 26456 | 264, 34700 | 272, 44250 | 280, 56300 |
|---|---:|---:|---:|---:|---:|---:|
| linear algebra, filtered Gauss | 0.034 s | 0.18 s | 0.72 s | 1.70 s | 5.00 s | 7.54 s |
| linear algebra, block Lanczos | 0.035 s | 0.11 s | 0.28 s | 0.45 s | 0.73 s | 1.10 s |
| peak memory, filtered Gauss | 98 MB | 323 MB | 187 MB | 305 MB | 473 MB | 731 MB |
| peak memory, block Lanczos | 19 MB | 28 MB | 43 MB | 59 MB | 70 MB | 124 MB |

(At 240 bits the 16143 primes still had dense decompositions, a byte per
relation and prime.) Below 12000 primes nothing changes: the trial division
of the relations has one more instance compiled, without the dense matrix,
and the code that runs is the same: against 0.4.8, the same instructions within
0.05% from 70 to 160 bits, and SIQS 0.98 to 1.01 of the time from 70 to 220
bits, `tifa_factor()` 1.00 from 70 to 100 bits (`bench/ab.sh`).
