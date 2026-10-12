"""Exact-arithmetic primitives for Fermat-number research.

Everything here uses Python/gmpy2 arbitrary-precision integers; there is no
floating point anywhere in a decision path.

Conventions
-----------
F(n) = 2^(2^n) + 1.
A prime p divides F(n)  <=>  2^(2^n) == -1 (mod p)  <=>  ord_p(2) == 2^(n+1).
Writing p = k*2^m + 1 with k odd, a prime factor of F(n) (n >= 2) has
m >= n + 2 (Euler-Lucas); we call e = m - (n + 2) >= 0 the "excess".
"""
from __future__ import annotations

import gmpy2
from gmpy2 import mpz

__all__ = [
    "fermat", "split_k_m", "divides_fermat", "fermat_index_of_factor",
    "proth_certificate", "is_proven_prime", "pepin_residue", "pepin",
]


def fermat(n: int) -> mpz:
    """F_n = 2^(2^n) + 1 (exact)."""
    if n < 0:
        raise ValueError("n must be >= 0")
    return (mpz(1) << (1 << n)) + 1


def split_k_m(p: int) -> tuple[int, int]:
    """Return (k, m) with p - 1 = k * 2^m and k odd."""
    p = mpz(p)
    if p < 3 or p % 2 == 0:
        raise ValueError("p must be an odd integer >= 3")
    q = p - 1
    m = int(gmpy2.bit_scan1(q))
    return int(q >> m), m


def divides_fermat(p: int, n: int) -> bool:
    """Exact test: does p divide F_n?  Cost: n modular squarings mod p.

    Works for any integer p >= 2 (prime or not) because it tests
    2^(2^n) == -1 (mod p), which is equivalent to p | F_n.
    """
    p = mpz(p)
    if p < 2:
        raise ValueError("p must be >= 2")
    if p == 2:
        return False
    x = mpz(2) % p
    for _ in range(n):
        x = (x * x) % p
    return x == p - 1


def fermat_index_of_factor(p: int, n_max: int | None = None) -> int | None:
    """If p divides some F_n, return that n (unique, since Fermat numbers are
    pairwise coprime), else None.  Only n <= m-1 can occur where
    2^m || p-1, so the search is complete without n_max for odd p.
    """
    p = mpz(p)
    if p % 2 == 0:
        return None
    _, m = split_k_m(p)
    top = m - 1 if n_max is None else min(m - 1, n_max)
    x = mpz(2) % p
    for j in range(top + 1):
        if x == p - 1:
            return j
        x = (x * x) % p
    return None


def proth_certificate(p: int, max_base: int = 1000) -> int | None:
    """Proth's theorem: if p = k*2^m + 1 with k < 2^m and some a satisfies
    a^((p-1)/2) == -1 (mod p), then p is prime.  Returns the witness a, or
    None if the theorem does not apply / no witness was found (which does NOT
    prove compositeness unless a Fermat test fails).
    """
    p = mpz(p)
    k, m = split_k_m(p)
    if k >= (1 << m):
        return None
    half = (p - 1) // 2
    for a in range(3, max_base):
        if gmpy2.jacobi(a, p) != -1:
            continue
        r = gmpy2.powmod(a, half, p)
        if r == p - 1:
            return a
        # For a non-residue a, prime p would give -1; anything else => composite.
        return None
    return None


def is_proven_prime(p: int) -> tuple[bool, str]:
    """Rigorous primality with the method recorded.

    * Proth certificate when applicable (all Fermat-factor shapes with k < 2^m).
    * Otherwise PARI/GP isprime() (APRCL / ECPP-backed proof) via cypari2 if
      available, else the external `gp` binary.
    Returns (is_prime, method).  A False result means "proven composite" only
    when method says so.
    """
    p = mpz(p)
    if p < 2:
        return False, "trivial"
    if p < (1 << 64):
        # BPSW is deterministic (no counterexamples) below 2^64 [Feitsma/Galway table].
        return bool(gmpy2.is_bpsw_prp(p)), "BPSW (deterministic < 2^64)"
    if not gmpy2.is_fermat_prp(p, 3):
        return False, "Fermat base-3 test (proves compositeness)"
    a = proth_certificate(p)
    if a is not None:
        return True, f"Proth certificate a={a}"
    return _gp_isprime(p)


def _gp_isprime(p: int) -> tuple[bool, str]:
    import subprocess
    out = subprocess.run(["gp", "-q", "-f"], input=f"default(parisizemax, 4000000000);\nprint(isprime({int(p)}))\n",
                         capture_output=True, text=True, timeout=36000)
    res = out.stdout.strip()
    if res not in ("0", "1"):
        raise RuntimeError(f"gp failed: {out.stderr}")
    return res == "1", "PARI/GP isprime (APRCL/ECPP proof)"


def pepin_residue(n: int, base: int = 3) -> mpz:
    """base^((F_n - 1)/2) mod F_n, by 2^n - 1 squarings using the special
    reduction x mod (2^N + 1) = (x mod 2^N) - (x >> N).
    """
    N = 1 << n
    F = (mpz(1) << N) + 1
    mask = (mpz(1) << N) - 1
    x = mpz(base)
    for _ in range(N - 1):
        x = x * x
        x = (x & mask) - (x >> N)
        if x < 0:
            x += F
    return x


def pepin(n: int) -> bool:
    """Pépin's test (n >= 1): F_n prime <=> 3^((F_n-1)/2) == -1 mod F_n."""
    if n < 1:
        raise ValueError("Pépin's test with base 3 needs n >= 1")
    return pepin_residue(n) == fermat(n) - 1
