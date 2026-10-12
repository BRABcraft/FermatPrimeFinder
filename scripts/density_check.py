#!/usr/bin/env python3
"""Numerical check of Theorem D (density form of the excess law) and of the
cubic-residue density used in H-k3, over all primes p < X.
  * v2(ord_p 2) given v2(p-1)=m  -> predicted P(v2(ord) = m-1-t) = 2^-(t+1), t <= m-2
  * Hasse: density of p with ord_p(2) even = 17/24
  * among p with 3 | (p-1): P(2 is a cubic residue) = 1/3, also within each class of m = v2(p-1)
"""
import json, os, sys
import gmpy2
X = int(sys.argv[1]) if len(sys.argv) > 1 else 20_000_000
from collections import defaultdict
tab = defaultdict(lambda: defaultdict(int)); even = tot = 0
cub = defaultdict(lambda: [0, 0])
p = gmpy2.mpz(2)
while True:
    p = gmpy2.next_prime(p)
    if p >= X: break
    m = int(gmpy2.bit_scan1(p - 1)); q = (p - 1) >> m
    x = gmpy2.powmod(2, q, p)        # x has order 2^v2(ord_p(2))
    v = 0
    while x != 1:
        x = x * x % p; v += 1
    tab[m][v] += 1; tot += 1; even += (v > 0)
    if (p - 1) % 3 == 0:
        c = cub[min(m, 6)]; c[1] += 1; c[0] += (gmpy2.powmod(2, (p - 1) // 3, p) == 1)
res = {"X": X, "primes": tot, "frac_ord_even": even / tot, "hasse_17_24": 17 / 24, "excess_law": {}, "cubic": {}}
for m in range(3, 13):
    n = sum(tab[m].values())
    if n < 200: continue
    res["excess_law"][m] = {"count": n, "observed": {t: round(tab[m][m - 1 - t] / n, 4) for t in range(0, 5)},
                            "predicted": {t: (2 ** -(t + 1) if t <= m - 2 else 2 ** -(m - 1) if t == m - 1 else 0.0)
                                          for t in range(0, 5)}}
for m, (c, n) in sorted(cub.items()):
    res["cubic"][m] = {"count": n, "frac_2_cubic_residue": round(c / n, 4), "predicted": round(1 / 3, 4)}
os.makedirs("results", exist_ok=True)
json.dump(res, open("results/density_check.json", "w"), indent=1)
print(json.dumps(res, indent=1))
