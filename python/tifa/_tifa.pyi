LIBRARY_VERSION: str

def _factor(
    n: int,
    hint: int,
    timeout: float,
    interruptible: bool,
    threads: int,
    siqs_max_bits: int,
    ecm_max_bits: int,
    /,
) -> tuple[bool, list[tuple[int, int]], bool, bool]: ...
def _prime_certificate(
    n: int,
    timeout: float,
    interruptible: bool,
    threads: int,
    siqs_max_bits: int,
    ecm_max_bits: int,
    /,
) -> tuple[int, list[tuple], bool]: ...
def _set_tuning_dir(path: bytes | None, /) -> None: ...
def _batch_gcd(
    moduli: list[int], timeout: float, interruptible: bool, /
) -> tuple[list[int], bool]: ...
def _cofactor_strategy(
    lpb: int, fb_bound: int, max_large_primes: int, max_siqs_bits: int, /
) -> object: ...
def _cofactor(
    strategy: object, ns: list[int], interruptible: bool, /
) -> list[dict[int, int] | None]: ...
