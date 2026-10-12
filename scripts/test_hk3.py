#!/usr/bin/env python3
"""Pre-registered test of H-k3 (research/hypotheses.md) on hold-out region R2c:
odd k in [65537, 1048575], 14 <= m <= 1000, n >= 12."""
import csv, json, math, os, sys
import numpy as np
from scipy import stats
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
from fermat.model_mc import Arr, build_sample
NMIN, KMIN, KMAX = 12, 65537, 1048575
region = {m: (KMIN, KMAX) for m in range(14, 1001)}
obs = []
with open(os.path.join(ROOT, "results", "search", "gmp_k65537-1048575_m14-1000.txt")) as f:
    for line in f:
        if line.startswith("FACTOR"):
            k, m, n = map(int, line.split()[1:4])
            import gmpy2
            if gmpy2.is_prime((k << m) + 1, 30):
                obs.append((k, m, n))
res = {"n_factors": len(obs)}
for seed in (11, 12):
    A = Arr(build_sample(region, s_small=1000, s_big=64, seed=seed), NMIN)
    exp = A.by_key(np.array([k % 3 for k in A.kint]))
    res[f"expected_seed{seed}"] = {int(c): round(v, 3) for c, v in exp.items()}
    res[f"expected_total_seed{seed}"] = round(float(A.mass.sum()), 3)
    if seed == 11:
        O = [sum(1 for k, m, n in obs if k % 3 == c) for c in range(3)]
        E = [exp.get(c, 0.0) for c in range(3)]
        Es = [e * sum(O) / sum(E) for e in E]
        G = 2 * sum(o * math.log(o / e) for o, e in zip(O, Es) if o > 0)
        lo = stats.poisson.cdf(O[0], E[0]); hi = stats.poisson.sf(O[0] - 1, E[0])
        res["observed_k_mod3"] = O
        res["i_poisson_p_k0"] = float(min(1, 2 * min(lo, hi)))
        res["ii_shape_G"] = G
        res["ii_shape_p"] = float(stats.chi2.sf(G, 2))
        res["direction_same_as_discovery"] = bool(O[0] < Es[0] and O[2] > Es[2])
        p = res["ii_shape_p"]
        res["decision"] = ("SUPPORTED" if (p < 0.01 and res["direction_same_as_discovery"]) else
                           "REJECTED as fluctuation" if p > 0.1 else "INCONCLUSIVE")
        k12 = sum(1 for k, m, n in obs if k % 12 == 5); k11 = sum(1 for k, m, n in obs if k % 11 == 5)
        e12 = A.by_key(np.array([k % 12 for k in A.kint])).get(5, 0.0)
        e11 = A.by_key(np.array([k % 11 for k in A.kint])).get(5, 0.0)
        res["secondary_k_mod12_eq5"] = [k12, round(e12, 3)]
        res["secondary_k_mod11_eq5"] = [k11, round(e11, 3)]
os.makedirs(os.path.join(ROOT, "results", "analysis"), exist_ok=True)
json.dump(res, open(os.path.join(ROOT, "results", "analysis", "hk3_holdout.json"), "w"), indent=1)
print(json.dumps(res, indent=1))
