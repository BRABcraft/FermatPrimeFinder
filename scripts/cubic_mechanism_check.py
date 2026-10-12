#!/usr/bin/env python3
"""Mechanism test for H-k3: among primes p = k*2^m + 1 with 3 | k (k odd), is 2 a cubic
residue with probability 1/3?  (The model's only k-mod-3-specific ingredient.)
Samples (k, m) uniformly with k odd, 3 | k, k < 2^20, 14 <= m <= 400; also reports the
small-k stratum k < 2^8 separately, since the 1/k weighting makes small k dominate."""
import json, random, sys, os
import gmpy2
rng = random.Random(7)
def run(kmax, target):
    hits = tot = 0
    while tot < target:
        k = 3 * (2 * rng.randrange(kmax // 6) + 1)
        m = rng.randrange(14, 401)
        p = (gmpy2.mpz(k) << m) + 1
        if not gmpy2.is_prime(p, 20): continue
        tot += 1
        hits += gmpy2.powmod(2, (p - 1) // 3, p) == 1
    return hits, tot
res = {}
for name, kmax, N in (("k<2^20", 1 << 20, 20000), ("k<2^8", 1 << 8, 5000)):
    h, t = run(kmax, N)
    se = (1 / 3 * 2 / 3 / t) ** 0.5
    res[name] = {"primes": t, "frac_2_cubic_residue": h / t, "z_vs_1/3": (h / t - 1 / 3) / se}
json.dump(res, open(os.path.join(os.path.dirname(__file__), "..", "results", "analysis", "cubic_mechanism.json"), "w"), indent=1)
print(json.dumps(res, indent=1))
