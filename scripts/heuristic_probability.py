#!/usr/bin/env python3
"""Probability, under the random-order-of-2 model, that some F_n (n >= 33)
is prime, conditioned on (a) our own certified search regions and (b) a range of
hypothetical deeper literature search limits.  Writes results/heuristic_probability.{json,md}.

P(F_n prime | no factor in searched region R)
   = (e^gamma / 2) * exp( - Lambda_unsearched(n) ),
Lambda_unsearched(n) = sum_{m >= n+2} 2^(n+1-m) * [ S(m, K_all(m)) - S(m, kmax_R(m)) ],
K_all(m) = 2^(2^(n-1) - m)  (factors up to sqrt(F_n)),
S(m, K) = sum_{odd k <= K} [k 2^m + 1 prime] / k   (exact small k, density beyond).
"""
from __future__ import annotations

import csv
import json
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
from fermat.heuristics import weight_sum  # noqa: E402
from fermat.regions import load_regions, region_kmax  # noqa: E402

LN2 = math.log(2)
EG2 = math.exp(0.5772156649015329) / 2


def S_upto_sqrtF(n: int, m: int) -> float:
    """S(m, 2^(2^(n-1)-m)) without materialising K: exact part for k<2^12, then the
    density integral to ln K + m ln2 = 2^(n-1) ln2."""
    small = weight_sum(m, (1 << 12) - 1)
    a = (1 << 12) - 1
    return small + math.log((2 ** (n - 1) * LN2) / (math.log(a) + m * LN2))


def lam_unsearched(n: int, kmax_of_m) -> float:
    tot = 0.0
    for m in range(n + 2, n + 90):
        km = kmax_of_m(m)
        searched = weight_sum(m, km) if km >= 1 else 0.0
        tot += 2.0 ** (n + 1 - m) * max(0.0, S_upto_sqrtF(n, m) - searched)
    return tot


def main():
    regions = load_regions()
    known = set()
    with open(os.path.join(ROOT, "data", "fermat_factors.csv")) as f:
        for r in csv.DictReader(f):
            if r["role"] in ("prime factor", "prime cofactor"):
                known.add(int(r["n"]))
    N_MAX = 160
    scen = {"our certified regions": lambda m: max((region_kmax(r, m) for r in regions), default=0)}
    for b in (40, 50, 60):
        scen[f"hypothetical: all k < 2^{b} for every m"] = (lambda bb: (lambda m: (1 << bb) - 1))(b)
    out = {"n_excluded_known_factor": sorted(x for x in known if x >= 33), "scenarios": {}}
    for name, f in scen.items():
        rows = []
        total = 0.0
        p_none = 1.0
        for n in range(33, N_MAX + 1):
            if n in known:
                continue
            lam = lam_unsearched(n, f)
            p = EG2 * math.exp(-lam)
            total += p
            p_none *= (1 - p)
            if n <= 40:
                rows.append({"n": n, "Lambda_unsearched": lam, "P_prime": p})
        out["scenarios"][name] = {"per_n": rows, "expected_new_fermat_primes_n33_to_160": total,
                                  "P_at_least_one": 1 - p_none}
    # unconditional / calibration: n = 5..32 with no conditioning
    cal = []
    for n in range(5, 33):
        lam = lam_unsearched(n, lambda m: 0)
        cal.append({"n": n, "P_prime_unconditional": EG2 * math.exp(-lam)})
    out["calibration_unconditional_n5_32"] = cal
    out["expected_primes_n5_32_unconditional"] = sum(c["P_prime_unconditional"] for c in cal)
    p0 = 1.0
    for c in cal:
        p0 *= 1 - c["P_prime_unconditional"]
    out["P_no_prime_n5_32_unconditional"] = p0
    # naive Hardy-Wright: 2/ln F_n
    out["naive_HW_n33_160"] = sum(2 / ((2 ** n) * LN2) for n in range(33, N_MAX + 1) if n not in known)
    os.makedirs(os.path.join(ROOT, "results"), exist_ok=True)
    with open(os.path.join(ROOT, "results", "heuristic_probability.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps({k: (v if k != "scenarios" else {s: {kk: vv for kk, vv in d.items() if kk != 'per_n'} for s, d in v.items()})
                      for k, v in out.items() if k != "calibration_unconditional_n5_32"}, indent=1))
    for name, d in out["scenarios"].items():
        print(name)
        for r in d["per_n"][:4]:
            print(f"   n={r['n']}: Lambda={r['Lambda_unsearched']:.3f}  P={r['P_prime']:.3e}")


if __name__ == "__main__":
    main()
