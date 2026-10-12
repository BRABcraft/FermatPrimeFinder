#!/usr/bin/env python3
"""Phase 3/4: test the random-order-of-2 model on the EXHAUSTIVE (selection-
bias-free) part of the factor data, run an automated congruence-conjecture
search, fit a generalised model by maximum likelihood, and run ML baselines.

Region = union of completed regions in results/search/regions.json (n >= 12).
Observed = verified prime factors (data/fermat_factors.csv) inside the region.
Because the region was searched exhaustively, there is no selection bias:
every factor in it is observed.

Outputs results/analysis/summary.json.
"""
from __future__ import annotations

import csv
import json
import math
import os
import sys

import numpy as np
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
from fermat.heuristics import weight_sum  # noqa: E402
from fermat.model_mc import Arr, build_sample  # noqa: E402
from fermat.regions import load_regions, union_mkmax  # noqa: E402

OUT = os.path.join(ROOT, "results", "analysis")
NMIN = 12


def poisson_two_sided(obs: int, exp: float) -> float:
    if exp <= 0:
        return 1.0 if obs == 0 else 0.0
    return float(min(1.0, 2 * min(stats.poisson.cdf(obs, exp), stats.poisson.sf(obs - 1, exp))))


def omega(k: int) -> int:
    c, q = 0, 3
    while q * q <= k:
        if k % q == 0:
            c += 1
            while k % q == 0:
                k //= q
        q += 2
    return c + (1 if k > 1 else 0)


def load_observed(mk):
    obs = set()
    with open(os.path.join(ROOT, "data", "fermat_factors.csv")) as f:
        for r in csv.DictReader(f):
            if r["role"] != "prime factor":
                continue
            n, k, m = int(r["n"]), int(r["k"]), int(r["m"])
            if n >= NMIN and m in mk and k <= mk[m]:
                obs.add((k, m, n))
    return sorted(obs)


def table(exp: dict, cnt: dict, label: str) -> dict:
    rows = []
    for g in sorted(set(exp) | set(cnt), key=lambda x: (str(type(x)), x)):
        e, o = exp.get(g, 0.0), cnt.get(g, 0)
        rows.append({"class": g, "observed": o, "expected": round(e, 3), "p_two_sided": poisson_two_sided(o, e)})
    big = [r for r in rows if r["expected"] >= 2]
    so = sum(r["observed"] for r in rows if r["expected"] < 2)
    se = sum(r["expected"] for r in rows if r["expected"] < 2)
    O = [r["observed"] for r in big] + ([so] if se > 0 else [])
    E = [r["expected"] for r in big] + ([se] if se > 0 else [])
    Ot, Et = sum(O), sum(E)
    Es = [e * Ot / Et for e in E]
    G = 2 * sum(o * math.log(o / e) for o, e in zip(O, Es) if o > 0 and e > 0)
    dof = max(1, len(O) - 1)
    return {"label": label, "rows": rows, "G_shape": G, "dof": dof, "p_shape": float(stats.chi2.sf(G, dof))}


def counts(obs, f):
    c = {}
    for k, m, n in obs:
        g = f(k, m, n)
        c[g] = c.get(g, 0) + 1
    return c


def pmod_array(A: Arr, q: int):
    pw = np.array([pow(2, int(mm), q) for mm in range(int(A.m.max()) + 1)], dtype=np.int64)
    kq = np.array([k % q for k in A.kint], dtype=np.int64)
    return (kq * pw[A.m] + 1) % q


def congruence_search(A: Arr, obs, max_mod=24, n_sim=300, seed=5):
    """Automated conjecture generator: for every modulus M <= max_mod (skipping
    M == 2 mod 4, which for odd k duplicates M/2) and residue class r of k, and
    for p mod q (odd primes q <= 31), compare observed vs model-expected counts.
    The null distribution of the number of nominal hits is calibrated by running
    the identical search on n_sim synthetic datasets drawn from the model."""
    feats = []   # (name, per-entry labels, per-observation labels, modulus)
    for M in range(2, max_mod + 1):
        if M % 4 == 2 and M > 2:
            continue
        feats.append((f"k mod {M}", np.array([k % M for k in A.kint]), [k % M for k, m, n in obs], M))
    for q in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31):
        feats.append((f"p mod {q}", pmod_array(A, q), [((k << m) + 1) % q for k, m, n in obs], q))
    E = []
    for name, lab, olab, M in feats:
        E.append(np.array([A.mass[lab == r].sum() for r in range(M)]))

    def pvals(counts_list):
        ps = []
        for cnts, ex in zip(counts_list, E):
            for r in range(len(ex)):
                if ex[r] == 0 and cnts[r] == 0:
                    continue
                ps.append(poisson_two_sided(int(cnts[r]), float(ex[r])))
        return np.array(ps)
    obs_counts = [np.bincount(np.array(olab, dtype=int), minlength=M) for _, _, olab, M in feats]
    p_obs = pvals(obs_counts)
    tests = []
    i = 0
    for (name, lab, olab, M), cnts, ex in zip(feats, obs_counts, E):
        for r in range(M):
            if ex[r] == 0 and cnts[r] == 0:
                continue
            tests.append({"feature": f"{name} == {r}", "observed": int(cnts[r]), "expected": round(float(ex[r]), 3),
                          "p": float(p_obs[i])})
            i += 1
    nt = len(tests)
    for t in tests:
        t["p_bonferroni"] = min(1.0, t["p"] * nt)
    # simulation-calibrated null
    rng = np.random.default_rng(seed)
    prob = A.mass / A.mass.sum()
    tot = float(A.mass.sum())
    sim_hits, sim_minp = [], []
    for _ in range(n_sim):
        N = rng.poisson(tot)
        idx = rng.choice(len(prob), size=N, p=prob)
        cl = [np.bincount(lab[idx].astype(int), minlength=M) for _, lab, _, M in feats]
        ps = pvals(cl)
        sim_hits.append(int((ps < 0.01).sum()))
        sim_minp.append(float(ps.min()))
    nh = int((p_obs < 0.01).sum())
    tests.sort(key=lambda t: t["p"])
    return {"n_tests": nt, "n_nominal_p<0.01": nh,
            "null_sim_mean_hits_p<0.01": float(np.mean(sim_hits)),
            "null_sim_P(hits>=observed)": float(np.mean(np.array(sim_hits) >= nh)),
            "observed_min_p": float(p_obs.min()),
            "null_sim_P(min_p<=observed_min_p)": float(np.mean(np.array(sim_minp) <= p_obs.min())),
            "n_significant_bonferroni_0.05": sum(t["p_bonferroni"] < 0.05 for t in tests),
            "top10": tests[:10],
            "structural_zeros_respected": all(t["observed"] == 0 for t in tests if t["expected"] == 0)}


def mle_alpha_r(A: Arr, obs):
    """Generalised model: intensity C * k^-alpha * [p prime] * (1-r) r^e, e = m-n-2.
    (alpha, r) = (1, 1/2) is the random-order-of-2 model.  C is profiled out."""
    logk = np.log(A.k)
    base = A.fac * A.ind  # = w * k
    emax = A.emax.astype(float)
    ok = np.array([k for k, m, n in obs], dtype=float)
    oe = np.array([m - n - 2 for k, m, n in obs], dtype=float)
    olk = np.log(ok)
    N = len(ok)

    def negll(alpha, r):
        lam = np.sum(base * np.exp(-alpha * logk) * (1 - r ** (emax + 1)))
        ll = -alpha * olk.sum() + N * np.log(1 - r) + oe.sum() * np.log(r) - N * np.log(lam) + N * np.log(N) - N
        return -ll
    Agrid = np.linspace(0.70, 1.30, 121)
    Rgrid = np.linspace(0.25, 0.75, 101)
    grid = np.array([[negll(a, r) for r in Rgrid] for a in Agrid])
    i, j = np.unravel_index(np.argmin(grid), grid.shape)
    best = grid[i, j]
    pa, pr = grid.min(axis=1), grid.min(axis=0)
    lr = 2 * (negll(1.0, 0.5) - best)
    return {"n_factors": N, "alpha_hat": float(Agrid[i]),
            "alpha_95ci": [float(Agrid[pa <= best + 1.92].min()), float(Agrid[pa <= best + 1.92].max())],
            "r_hat": float(Rgrid[j]),
            "r_95ci": [float(Rgrid[pr <= best + 1.92].min()), float(Rgrid[pr <= best + 1.92].max())],
            "LR_vs_random_model": float(lr), "p_value_LR_2dof": float(stats.chi2.sf(lr, 2))}


def ml_experiments(A: Arr, obs, seed=0):
    """Case-control classification: positives = factors in the exactly-tested
    sub-region (p <= 600 bits); negatives = population-representative primes from
    the stratified sample (resampled with probability proportional to fac)."""
    from sklearn.ensemble import GradientBoostingClassifier, IsolationForest
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import average_precision_score, roc_auc_score
    rng = np.random.default_rng(seed)
    exact_idx = np.where(A.exact & (A.ind > 0))[0]
    pos = [(k, m) for k, m, n in obs if m + k.bit_length() <= 600]
    pos_set = set(pos)
    pr = A.fac[exact_idx] / A.fac[exact_idx].sum()
    pick = rng.choice(exact_idx, size=min(len(exact_idx), 60 * len(pos)), p=pr)
    neg = [(A.kint[i], int(A.m[i])) for i in pick if (A.kint[i], int(A.m[i])) not in pos_set]

    def feats(k, m):
        b = k.bit_length()
        p = (k << m) + 1
        return [math.log2(k), m, k % 3, k % 5, k % 7, k % 8, k % 16, bin(k).count("1") / b,
                omega(k), m % 2, m % 3, m % 4, p % 3, p % 5, p % 7]
    names = ["log2k", "m", "k%3", "k%5", "k%7", "k%8", "k%16", "popcount_frac", "omega(k)",
             "m%2", "m%3", "m%4", "p%3", "p%5", "p%7"]
    X = np.array([feats(k, m) for k, m in pos] + [feats(k, m) for k, m in neg], dtype=float)
    y = np.array([1] * len(pos) + [0] * len(neg))
    msc = np.array([-math.log(k) + math.log(1 - 2.0 ** (NMIN + 1 - m)) for k, m in pos + neg])
    groups = X[:, 1].astype(int) % 5
    preds = {nm: np.zeros(len(y)) for nm in ("logistic(model_score)", "logistic(model_score+features)",
                                             "gradient_boosting(features)")}
    for g in range(5):
        tr, te = groups != g, groups == g
        gbm = GradientBoostingClassifier(n_estimators=200, max_depth=3, learning_rate=0.05,
                                         subsample=0.8, random_state=seed).fit(X[tr], y[tr])
        preds["gradient_boosting(features)"][te] = gbm.predict_proba(X[te])[:, 1]
        lr = LogisticRegression(max_iter=2000).fit(msc[tr, None], y[tr])
        preds["logistic(model_score)"][te] = lr.predict_proba(msc[te, None])[:, 1]
        Z = np.c_[msc, X]
        mu, sd = Z[tr].mean(0), Z[tr].std(0) + 1e-9
        lra = LogisticRegression(max_iter=5000, C=1.0).fit((Z[tr] - mu) / sd, y[tr])
        preds["logistic(model_score+features)"][te] = lra.predict_proba((Z[te] - mu) / sd)[:, 1]
    preds["model_score_-log(k)_raw"] = msc
    preds["random"] = rng.random(len(y))

    def boot_ci(pred, B=1000):
        r = np.random.default_rng(seed + 1)
        vals = []
        for _ in range(B):
            ii = r.integers(0, len(y), len(y))
            if 0 < y[ii].sum() < len(ii):
                vals.append(roc_auc_score(y[ii], pred[ii]))
        return [float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))]
    res = {"n_pos": len(pos), "n_neg": len(neg), "split": "grouped 5-fold by m mod 5", "features": names}
    for nm, pd in preds.items():
        res[nm] = {"AUC": float(roc_auc_score(y, pd)), "AUC_95ci": boot_ci(pd),
                   "avg_precision": float(average_precision_score(y, pd))}
    gb = GradientBoostingClassifier(n_estimators=200, max_depth=3, learning_rate=0.05, subsample=0.8,
                                    random_state=seed).fit(X, y)
    res["gbm_feature_importance"] = dict(sorted(zip(names, map(float, gb.feature_importances_)), key=lambda t: -t[1]))
    # anomaly detection among all factors in region
    Xf = np.array([[math.log2(k), m, m - n - 2, n, k % 3, omega(k)] for k, m, n in obs], dtype=float)
    iso = IsolationForest(random_state=seed).fit(Xf)
    sc = iso.score_samples(Xf)
    res["isolation_forest_top_anomalies"] = [
        {"k": obs[i][0], "m": obs[i][1], "n": obs[i][2], "excess": obs[i][1] - obs[i][2] - 2, "score": float(sc[i])}
        for i in np.argsort(sc)[:6]]
    return res


def counting_function(Bmax=54):
    """N(2^B) = #{primes p < 2^B dividing some F_n, n >= 5}: exact from full factorizations of
    F5..F11 plus the exhaustive p < 2^54 search (n >= 12), vs model expectation."""
    ps = []
    with open(os.path.join(ROOT, "data", "fermat_factors.csv")) as f:
        for r in csv.DictReader(f):
            if int(r["n"]) >= 5 and r["role"] != "F_n is prime":
                ps.append(int(r["p"]))
    rows = []
    for B in range(10, Bmax + 1, 4):
        exp = sum(weight_sum(m, ((1 << B) - 2) >> m) * (1 - 2.0 ** (6 - m))
                  for m in range(7, B) if ((1 << B) - 2) >> m >= 1)
        o = sum(1 for p in ps if p < (1 << B))
        rows.append({"B": B, "observed_N": o, "model_expected": round(exp, 2), "p": poisson_two_sided(o, exp)})
    return rows


def main():
    regions = load_regions()
    m_hi = max(r.get("mmax", r.get("B", 0)) for r in regions)
    mk = union_mkmax(regions, NMIN, m_hi)
    obs = load_observed(mk)
    print(f"regions: {[r['name'] for r in regions]}; m in [{min(mk)}, {max(mk)}]; observed factors: {len(obs)}")
    out = {"regions": [r["name"] for r in regions], "n_observed": len(obs)}
    A1 = Arr(build_sample(mk, s_small=1000, seed=1), NMIN)
    A2 = Arr(build_sample(mk, s_small=1000, seed=2), NMIN)
    t1, t2 = float(A1.mass.sum()), float(A2.mass.sum())
    out["T1_total"] = {"observed": len(obs), "expected_seed1": t1, "expected_seed2": t2,
                       "p_two_sided": poisson_two_sided(len(obs), t1)}
    A = A1
    out["T2_excess"] = table(A.by_excess(5), counts(obs, lambda k, m, n: min(m - n - 2, 5)), "excess e=m-n-2 (5 means >=5)")
    out["T3_bitlen_k"] = table(A.by_key(np.array([min(k.bit_length() // 4, 9) for k in A.kint])),
                               counts(obs, lambda k, m, n: min(k.bit_length() // 4, 9)), "bitlength(k)//4 (1/k law)")
    for q in (3, 5, 7):
        out[f"T4_k_mod{q}"] = table(A.by_key(np.array([k % q for k in A.kint])),
                                    counts(obs, lambda k, m, n: k % q), f"k mod {q}")
    out["T4_omega"] = table(A.by_key(np.array([min(omega(k), 4) for k in A.kint])),
                            counts(obs, lambda k, m, n: min(omega(k), 4)), "omega(k) capped at 4")
    out["T5_m_bins"] = table(A.by_key(np.minimum(np.floor(np.log2(A.m)).astype(int), 11)),
                             counts(obs, lambda k, m, n: min(int(math.log2(m)), 11)), "floor(log2 m)")
    out["T6_congruence_search"] = congruence_search(A, obs)
    out["T7_mle"] = mle_alpha_r(A, obs)
    out["T7_mle_seed2"] = mle_alpha_r(A2, obs)
    out["T8_ml"] = ml_experiments(A, obs)
    out["T9_counting_function"] = counting_function()
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "summary.json"), "w") as f:
        json.dump(out, f, indent=1, default=str)
    print(json.dumps(out["T1_total"]))
    for key in [k for k in out if k[:2] in ("T2", "T3", "T4", "T5")]:
        d = out[key]
        print(f"{key}: p_shape={d['p_shape']:.3f}  " + " ".join(f"{r['class']}:{r['observed']}/{r['expected']:.1f}" for r in d["rows"]))
    c = out["T6_congruence_search"]
    print("T6:", {k: v for k, v in c.items() if k != "top10"}, c["top10"][:3])
    print("T7:", json.dumps(out["T7_mle"]))
    print("T7 seed2:", json.dumps(out["T7_mle_seed2"]))
    print("T8:", json.dumps({k: v for k, v in out["T8_ml"].items() if k not in ("features",)}, default=str))
    print("T9:", out["T9_counting_function"])


if __name__ == "__main__":
    main()
