"""The 'random order of 2' model for Fermat-number factors, and its predictions.

Model (standard heuristic; Chebotarev-motivated):
  For a prime p = k*2^m + 1 (k odd, m >= 3), 2 is a quadratic residue, and we
  treat 2 as a uniformly random element of the subgroup of squares, which is
  cyclic of order k*2^(m-1).  Then
     P(ord_p(2) is a power of 2)            = 1/k,
     P(ord_p(2) = 2^(n+1) | power of 2)     = 2^(n-m+1),   0 <= n <= m-2,
  i.e. p | F_n with probability 2^(n-m+1)/k, and the excess e = m-n-2 is
  geometric: P(e = t) = 2^-(t+1).

Prime density among k*2^m+1 over odd k is ~2/ln p on average.  For small k the
actual prime indicator is used (exact), because the 1/k weight puts most of
the expected count on small k where densities fluctuate (Riesel/Sierpinski
effects).
"""
from __future__ import annotations

import math
from functools import lru_cache

import gmpy2

LN2 = math.log(2.0)


@lru_cache(maxsize=None)
def _is_prime_km(k: int, m: int) -> bool:
    return bool(gmpy2.is_prime((k << m) + 1, 25))


def weight_sum(m: int, kmax: int, kexact: int = 1 << 12, mexact: int = 600) -> float:
    """S(m) = sum over odd k <= kmax of [k*2^m+1 prime] / k.

    Exact prime indicator for k < kexact (when m <= mexact); beyond that the
    average density 2/ln(p) replaces the indicator and the sum over odd k is
    integrated:  sum_{odd k in (a,b]} 2/(k ln(k 2^m))
                ~ ∫_a^b dk / (k (ln k + m ln2))
                = ln((ln b + m ln2) / (ln a + m ln2)).
    """
    if kmax < 1:
        return 0.0
    s = 0.0
    lim = min(kmax, kexact - 1) if m <= mexact else 0
    for k in range(1, lim + 1, 2):
        if _is_prime_km(k, m):
            s += 1.0 / k
    a = max(lim, 1) if lim >= 1 else 0
    if kmax > a:
        if a == 0:
            # include k = 1 term by density, then integrate from 1
            s += 2.0 / (m * LN2)
            a = 1
        s += math.log((math.log(kmax) + m * LN2) / (math.log(a) + m * LN2))
    return s


def expected_factor_cells(m_kmax: dict[int, int], nmin: int, **kw) -> dict[tuple[int, int], float]:
    """Expected number of prime factors in each (n, m) cell of a region where,
    for each m, all odd k <= m_kmax[m] were searched and n in [nmin, m-2]."""
    cells = {}
    for m, kmax in m_kmax.items():
        S = weight_sum(m, kmax, **kw)
        for n in range(nmin, m - 1):
            cells[(n, m)] = S * 2.0 ** (n - m + 1)
    return cells


def prob_prime_given_no_factor(n: int, log2_bound: float, mertens_correction: bool = True) -> float:
    """P(F_n prime | F_n has no prime factor p < 2^log2_bound) under the model.

    Sum over admissible p = K*2^(n+2)+1 in (B, sqrt F_n] of 2^(n+1)/(p-1)
    ~ ln( ln sqrt(F_n) / ln B ), so P(no factor in that range) ~ ln B / ln sqrt(F_n)
    = log2_bound / 2^(n-1).  The optional factor e^gamma/2 ~ 0.8905 is the usual
    correction between 'no prime factor <= sqrt(x)' (Mertens product, ~2e^-gamma/ln x)
    and 'prime' (~1/ln x).
    """
    p = log2_bound / 2.0 ** (n - 1)
    if mertens_correction:
        p *= math.exp(0.5772156649015329) / 2.0
    return min(1.0, p)
