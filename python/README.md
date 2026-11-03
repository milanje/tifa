# tifa

Fast factorization of integers of up to about 220 bits (66 digits), with
the [TIFA](https://github.com/milanje/tifa) C library: its quadratic sieve,
elliptic curves and smaller methods, chosen by the size of the number.

```python
>>> import tifa
>>> tifa.factor(360)
{2: 3, 3: 2, 5: 1}
>>> tifa.factor(2**64 + 1)
{274177: 1, 67280421310721: 1}
>>> tifa.factor(1)
{}
```

`factor(n)` returns the prime factorization of a positive integer as
`{prime: exponent}`, primes in increasing order, like SymPy's `factorint`.
`n` can be an `int` or any integer type with `__index__` (NumPy's, gmpy2's).

## Speed

One thread, balanced semiprimes (the hardest numbers of a size), on an AMD
EPYC 9V45:

| bits | 100 | 140 | 180 | 220 |
|---|---:|---:|---:|---:|
| `tifa.factor` | 0.7 ms | 9.6 ms | 0.18 s | 3.0 s |
| PARI/GP 2.15 `factorint` | 1.9 ms | 39 ms | 0.55 s | 9.0 s |

Numbers with small factors take much less. From 180 to 220 bits the time
doubles about every 10 bits, and keeps growing beyond.

## Hints

When something is known of `n`, a hint saves time:

```python
tifa.factor(n, hint="balanced")    # two primes of about the same size
tifa.factor(n, hint="unbalanced")  # a prime factor under a third of its size
```

A hint never changes the result, but a wrong one can cost much more time
than none.

## Threads

`factor` releases the GIL, so that threads factor in parallel:

```python
from concurrent.futures import ThreadPoolExecutor

with ThreadPoolExecutor() as pool:
    results = list(pool.map(tifa.factor, numbers))
```

## Errors

`TypeError` if `n` is not an integer, `ValueError` if it is not positive
or the hint is unknown. In the rare event that `n` can not be factored
completely, `tifa.FactorizationError` is raised, its `factors` attribute
holding what was found (factors that multiply to `n`, some of them not
prime).

## Licence

LGPL 2.1 or later, as TIFA.
