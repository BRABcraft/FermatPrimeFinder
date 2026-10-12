"""Stratified Monte-Carlo expectations under the random-order-of-2 model.

Region: for each m, all odd k <= kmax[m]; admissible n in [nmin, m-2].
Model: a prime p = k*2^m+1 divides F_n with probability 2^(n-m+1)/k.

We build a weighted sample of candidate primes in the region:
  * dyadic k-blocks [2^j, 2^(j+1)) intersected with [1, kmax];
  * blocks with few odd k are enumerated exactly; larger blocks are sampled
    uniformly (s odd k per block) with expansion factor N_odd_block / s;
  * primality indicator: exact (BPSW, 25 MR rounds) when p has <= EXACT_BITS bits
    or k < EXACT_K; otherwise the sieve-conditioned density
        [p has no prime factor <= 1000] * e^gamma * ln(1000) / ln p,
    the standard (Mertens) probability that a number of that size with no
    factor <= 1000 is prime.  (Validated against exact counts; see tests.)
Then E[#factors with feature f, index n] = sum fac * ind / k * 2^(n-m+1) * f.
"""
from __future__ import annotations

import math
import random
from dataclasses import dataclass

import gmpy2

EXACT_BITS = 600
EXACT_K = 64
_Q = 1000
_PRIMORIAL = 1
for _q in range(3, _Q + 1, 2):
    if gmpy2.is_prime(_q):
        _PRIMORIAL *= _q
_PRIMORIAL = gmpy2.mpz(_PRIMORIAL)
_SIEVE_BOOST = math.exp(0.5772156649015329) * math.log(_Q)


@dataclass
class WPrime:
    k: int
    m: int
    fac: float        # expansion factor (1 for exactly enumerated blocks)
    ind: float = 1.0  # primality indicator (1) or probability (density mode)
    exact: bool = True

    @property
    def w(self) -> float:  # expected-number weight before the n-split (model: 1/k)
        return self.fac * self.ind / self.k


def prime_indicator(k: int, m: int) -> tuple[float, bool]:
    p = (gmpy2.mpz(k) << m) + 1
    bits = m + k.bit_length()
    if bits <= EXACT_BITS or k < EXACT_K:
        return (1.0 if gmpy2.is_prime(p, 25) else 0.0), True
    if gmpy2.gcd(p, _PRIMORIAL) != 1:
        return 0.0, False
    return _SIEVE_BOOST / (bits * math.log(2)), False


def build_sample(m_kmax: dict[int, int], exact_max_odd: int = 512, s_small: int = 3000,
                 s_big: int = 48, big_bits: int = 64, seed: int = 1) -> list[WPrime]:
    rng = random.Random(seed)
    out: list[WPrime] = []
    for m, kmax in sorted(m_kmax.items()):
        j = 0
        while (1 << j) <= kmax:
            lo, hi = 1 << j, min((1 << (j + 1)) - 1, kmax)
            first = lo | 1
            n_odd = (hi - first) // 2 + 1 if hi >= first else 0
            if n_odd <= 0:
                j += 1
                continue
            if n_odd <= exact_max_odd:
                ks, fac = range(first, hi + 1, 2), 1.0
            else:
                s = s_small if (hi.bit_length() + m) <= big_bits else s_big
                fac = n_odd / s
                ks = [first + 2 * rng.randrange(n_odd) for _ in range(s)]
            for k in ks:
                ind, ex = prime_indicator(k, m)
                if ind > 0:
                    out.append(WPrime(k, m, fac, ind, ex))
            j += 1
    return out


class Arr:
    """Vectorised view of a sample (numpy).  k is kept as float64 for features
    (exact integer residues are computed from Python ints at build time)."""

    def __init__(self, sample: list[WPrime], nmin: int):
        import numpy as np
        self.np = np
        self.nmin = nmin
        self.kint = [w.k for w in sample]
        self.k = np.array([float(w.k) for w in sample])
        self.m = np.array([w.m for w in sample], dtype=np.int64)
        self.w = np.array([w.w for w in sample])
        self.fac = np.array([w.fac for w in sample])
        self.ind = np.array([w.ind for w in sample])
        self.exact = np.array([w.exact for w in sample])
        self.emax = self.m - 2 - nmin
        # total model mass over admissible n:  sum_{e=0}^{emax} 2^-(e+1) = 1 - 2^-(emax+1)
        self.mass = self.w * (1.0 - 2.0 ** (-(self.emax + 1).astype(float)))

    def by_key(self, keys) -> dict:
        """keys: array of per-entry class labels (depends on k, m only)."""
        np = self.np
        out = {}
        keys = np.asarray(keys)
        for g in np.unique(keys):
            out[g.item() if hasattr(g, "item") else g] = float(self.mass[keys == g].sum())
        return out

    def by_excess(self, cap: int = 5) -> dict:
        np = self.np
        out = {}
        for e in range(cap):
            out[e] = float((self.w * 2.0 ** (-(e + 1)) * (self.emax >= e)).sum())
        # e >= cap: 2^-(cap+1) + ... up to emax  = 2^-cap - 2^-(emax+1)
        tail = self.w * np.where(self.emax >= cap, 2.0 ** (-cap) - 2.0 ** (-(self.emax + 1).astype(float)), 0.0)
        out[cap] = float(tail.sum())
        return out

    def int_feature(self, f) -> "object":
        return self.np.array([f(k, m) for k, m in zip(self.kint, self.m.tolist())])


def expected_by(sample: list[WPrime], nmin: int, key) -> dict:
    """Sum of model expectations grouped by key(k, m, n) (key may return None to skip)."""
    acc: dict = {}
    for wp in sample:
        m, k = wp.m, wp.k
        base = wp.w
        # P(n) = 2^(n-m+1) = 2^-(e+1) with e = m-2-n; terms with e > 60 are < 2^-61, dropped
        for e in range(0, min(m - 2 - nmin, 60) + 1):
            n = m - 2 - e
            g = key(k, m, n)
            if g is None:
                continue
            acc[g] = acc.get(g, 0.0) + base * 2.0 ** (-(e + 1))
    return acc
