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
give the sparse matrix. On an Apple M5 Pro, 19 balanced semiprimes of 240
to 280 bits, each factored with both methods (means; the whole
factorization is the ratio of its time with block Lanczos to its time with
filtered Gauss, with the lowest and highest ratios over the numbers):

| bits | 240 | 256 | 264 | 272 | 280 |
|---|---:|---:|---:|---:|---:|
| numbers | 5 | 4 | 4 | 3 | 3 |
| primes in the base | 16143&ndash;16600 | 26456&ndash;27200 | 33762&ndash;34700 | 43056&ndash;44250 | 54793&ndash;56300 |
| linear algebra, filtered Gauss | 0.19 s | 0.71 s | 1.80 s | 3.51 s | 7.15 s |
| linear algebra, block Lanczos | 0.11 s | 0.27 s | 0.45 s | 0.92 s | 1.12 s |
| whole factorization | 0.981 (0.962&ndash;0.991) | 0.977 (0.971&ndash;0.986) | 0.964 (0.961&ndash;0.966) | 0.955 (0.952&ndash;0.956) | 0.940 (0.938&ndash;0.943) |
| peak memory, filtered Gauss | 181 MB | 191 MB | 300 MB | 455 MB | 721 MB |
| peak memory, block Lanczos | 30 MB | 44 MB | 58 MB | 70 MB | 83 MB |

Every number split with both methods. At 240 bits the peak of filtered
Gauss depends on the decompositions: 323 MB where the base of 16143 primes
still kept them dense, a byte per relation and prime, 86 MB with 16600
primes. One 272-bit number took block Lanczos 1.40 s, the two others 0.65
and 0.69 s. At 220 bits, 8520 primes, below the threshold, both methods
took the same time on one number, 0.034 s and 0.035 s.

Below 12000 primes nothing changes: the trial division of the relations has
one more instance compiled, without the dense matrix, and the code that runs
is the same: against 0.4.8, the same instructions within 0.05% from 70 to
160 bits, and SIQS 0.98 to 1.01 of the time from 70 to 220
bits, `tifa_factor()` 1.00 from 70 to 100 bits (`bench/ab.sh`).
