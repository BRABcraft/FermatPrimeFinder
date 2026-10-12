#!/usr/bin/env python3
"""Build data/fermat_numbers.csv and data/fermat_factors.csv from
(1) literature factors (recalled from W. Keller's table / web snippets) that we
    RE-VERIFY here by exact arithmetic -- unverifiable claims go to
    data/unverified_reported_factors.csv, never into the main table;
(2) our own exhaustive searches (results/search/*.txt), which are complete
    inside their stated regions (results/search/regions.json);
(3) our own Pépin / Suyama runs (results/pepin/*.json).

Every prime factor row is checked:  p | F_n  (2^(2^n) == -1 mod p)  and  p prime
(BPSW below 2^64, Proth certificate or PARI/GP APRCL proof above).
"""
from __future__ import annotations

import csv
import glob
import json
import math
import os
import sys
from functools import reduce

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))

import gmpy2  # noqa: E402
from gmpy2 import mpz  # noqa: E402
from fermat.core import (divides_fermat, fermat, is_proven_prime,  # noqa: E402
                         split_k_m)

DATA = os.path.join(ROOT, "data")

# ---------------------------------------------------------------------------
# (1) Literature factors.  Values typed from memory of Keller's table or from
#     web-search snippets (see research/literature_review.md).  They are
#     *claims* until verified below.
# ---------------------------------------------------------------------------
FULLY_FACTORED_SMALL = {  # n: known prime factors except the largest cofactor
    5: [641],
    6: [274177],
    7: [59649589127497217],
    8: [1238926361552897],
    9: [2424833, 7455602825647884208337395736200454918783366342657],
    10: [45592577, 6487031809, 4659775785220018543264560743076778192897],
    11: [319489, 974849, 167988556341760475137, 3560841906445833920513],
}
LITERATURE_FACTORS = [  # (n, p, provenance)
    (12, 568630647535356955169033410940867804839360742060818433,
     "M. Vang 2010 (ECM) [snippet: caramel.loria.fr/f12.txt, Prime-Wiki]"),
    (13, 2663848877152141313, "R. Crandall 1991 [recalled, Keller table]"),
    (13, 3603109844542291969, "R. Crandall 1991 [recalled, Keller table]"),
    (13, 319546020820551643220672513, "R. Brent 1995 (ECM) [recalled, Keller table]"),
    (14, 116928085873074369829035993834596371340386703423373313,
     "2010 (ECM) [snippet: OEIS A093179]"),
    (15, 168768817029516972383024127016961, "Crandall et al. 1997 [recalled]"),
    (16, 188981757975021318420037633, "Crandall et al. 1996 [recalled]"),
    (17, 7751061099802522589358967058392886922693580423169, "2011 (ECM) [recalled]"),
    (18, 81274690703860512587777, "Crandall et al. 1999 [recalled]"),
    (19, 37590055514133754286524446080499713, "2009 (ECM) [recalled]"),
    (22, 64658705994591851009055774868504577,
     "2010-03-26 [snippet: mersenne.org F22 record]"),
    (25, 204393464266227713, "[recalled, Keller table]"),
    (25, 2170072644496392193, "[recalled, Keller table]"),
    (27, 231292694251438081, "[recalled, Keller table]"),
    (29, 1120048333814181889, "[recalled, Keller table] -- garbled recall; true factor 1120049*2^31+1 found by our search"),
    (28, 25463 * 2**30 + 1, "[recalled, Keller table] -- wrong recall"),
    (31, 46931635677864055013377, "A. Kruppa 2001-04-12 [snippet: aliquot list / Prime-Wiki]"),
    (251, 322072887044529 * 2**253 + 1, "G. Gostin 2025-05-24 [snippet: fermatsearch news]"),
    (3135, 36812439837 * 2**3144 + 1, "G. Gostin 2025-11 [snippet: fermatsearch news]"),
]
# Reported but too expensive to verify here (cost ~ n squarings of ~n-bit numbers).
REPORTED_UNVERIFIED = [
    (414677, "39701*2^414679+1", "Batalov & Propper 2025-01 [snippet: fermatsearch news]"),
    (3943261, "1085*2^3943263+1", "Dinkel / PrimeGrid 2025-12 [snippet: fermatsearch news]"),
    (5798447, "99*2^5798449+1", "V. Cavecchia 2025-01 [snippet: fermatsearch news]"),
]


def decimal_digits_of_fermat(n: int) -> int:
    """Number of decimal digits of F_n = digits of 2^(2^n) (never a power of 10)."""
    if n <= 22:
        return len(str(fermat(n)))
    import mpmath
    mpmath.mp.dps = int(n * 0.302) + 80  # enough for the fractional part of 2^n*log10(2)
    x = mpmath.mpf(2) ** n * mpmath.log10(2)
    fl = int(mpmath.floor(x))
    if min(x - fl, fl + 1 - x) < mpmath.mpf(10) ** -50:
        raise ArithmeticError("too close to an integer for safe floor")
    return fl + 1


def load_search_factors():
    rows = []
    for path in sorted(glob.glob(os.path.join(ROOT, "results", "search", "*.txt"))):
        with open(path) as f:
            for line in f:
                parts = line.split()
                if not parts or parts[0] != "FACTOR":
                    continue
                if os.path.basename(path).startswith("fsearch64"):
                    p, k, m, n = map(int, parts[1:5])
                else:
                    k, m, n = map(int, parts[1:4])
                    p = (k << m) + 1
                rows.append((n, p, f"exhaustive search (this work): {os.path.basename(path)}"))
    return rows


def load_pepin():
    res = {}
    for path in glob.glob(os.path.join(ROOT, "results", "pepin", "F*.json")):
        with open(path) as f:
            d = json.load(f)
        res[d["n"]] = d
    return res


def main():
    factors = {}      # (n, p) -> dict
    problems = []

    def add(n, p, prov, role="prime factor"):
        p = mpz(p)
        key = (n, int(p))
        if key in factors:
            factors[key]["provenance"] += " | " + prov
            return
        ok_div = divides_fermat(p, n)
        if not ok_div:
            problems.append(f"F{n}: claimed factor {p} does NOT divide F{n} ({prov})")
            return
        isp, method = is_proven_prime(p)
        if not isp:
            # Composite divisor (e.g. product of two factors found by a search).
            problems.append(f"F{n}: divisor {p} is composite ({method}); not a prime factor ({prov})")
            return
        k, m = split_k_m(p)
        factors[key] = dict(n=n, p=int(p), k=k, m=m, excess=m - (n + 2) if n >= 2 else None,
                            bits=int(p).bit_length(), digits=len(str(p)),
                            primality_proof=method, divisibility_verified=True,
                            provenance=prov, role=role)

    # Fully factored F5..F11: cofactor computed exactly and proven prime.
    for n, small in FULLY_FACTORED_SMALL.items():
        F = fermat(n)
        prod = reduce(lambda a, b: a * b, small, mpz(1))
        if F % prod != 0:
            problems.append(f"F{n}: product of listed small factors does not divide F{n}")
            continue
        cof = F // prod
        for p in small:
            add(n, p, "literature (Keller table) [recalled]; verified here")
        add(n, cof, "largest prime factor = F_n / (product of others), computed here",
            role="prime cofactor")
    for n, p, prov in LITERATURE_FACTORS:
        add(n, p, prov + "; verified here")
    for n, p, prov in load_search_factors():
        add(n, p, prov)

    # Fermat primes F0..F4 as their own (only) factor
    for n in range(5):
        add(n, fermat(n), "F_n itself (prime)", role="F_n is prime")

    fac_rows = sorted(factors.values(), key=lambda r: (r["n"], r["p"]))

    # ---------------------------------------------------------------- numbers
    pepin = load_pepin()
    by_n = {}
    for r in fac_rows:
        by_n.setdefault(r["n"], []).append(r)
    max_n = 40
    num_rows = []
    lit_composite_no_factor = {
        20: "Young & Buell 1988 (Pepin); no factor known",
        24: "Crandall, Mayer, Papadopoulos 2003 (Pepin); no factor known",
        28: "literature reports a known factor (Keller table) that we could not recall/verify; "
            "our exhaustive search proves F28 has no prime factor < 2^54",
    }
    for n in sorted(set(range(0, max_n + 1)) | set(by_n)):
        fl = [r for r in by_n.get(n, []) if r["role"] != "F_n is prime"]
        known = sorted(r["p"] for r in fl)
        full = any(r["role"] == "prime cofactor" for r in fl)
        if n <= 4:
            status, evidence = "prime", "Pepin test (this work, exact)" if n >= 1 else "trivial (3 is prime)"
        elif known:
            status, evidence = "composite", "known prime factor verified here"
        elif n in pepin and pepin[n]["verdict"] == "COMPOSITE":
            status, evidence = "composite", f"Pepin test (this work, {pepin[n]['impl']})"
        elif n in lit_composite_no_factor:
            status, evidence = "composite", "literature only: " + lit_composite_no_factor[n]
        else:
            status, evidence = "unknown", ""
        cof = ""
        if n in pepin and pepin[n].get("suyama"):
            cof = pepin[n]["suyama"]
        elif full:
            cof = "fully factored (cofactor proven prime here)"
        num_rows.append(dict(
            n=n, value=str(fermat(n)) if n <= 5 else "",
            binary_digits=(1 << n) + 1, decimal_digits=decimal_digits_of_fermat(n),
            status=status, status_evidence=evidence,
            fully_factored=full or n <= 4,
            num_known_prime_factors=len(known),
            known_prime_factors=" ".join(map(str, known)) if all(len(str(p)) < 200 for p in known)
            else " ".join(str(p) if len(str(p)) < 200 else f"<P{len(str(p))}>" for p in known),
            cofactor_status=cof,
            pepin_res64=pepin.get(n, {}).get("res64", ""),
        ))

    os.makedirs(DATA, exist_ok=True)
    with open(os.path.join(DATA, "fermat_factors.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(fac_rows[0].keys()))
        w.writeheader()
        w.writerows(fac_rows)
    with open(os.path.join(DATA, "fermat_numbers.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(num_rows[0].keys()))
        w.writeheader()
        w.writerows(num_rows)
    with open(os.path.join(DATA, "unverified_reported_factors.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["n", "p_expression", "provenance", "reason_not_verified"])
        for n, e, prov in REPORTED_UNVERIFIED:
            w.writerow([n, e, prov, f"verification needs ~{n} squarings of ~{n}-bit numbers; skipped"])
    with open(os.path.join(DATA, "build_problems.txt"), "w") as f:
        f.write("\n".join(problems) + ("\n" if problems else ""))
    print(f"{len(fac_rows)} verified prime-factor rows, {len(num_rows)} Fermat-number rows")
    print(f"{len(problems)} rejected claims/divisors (see data/build_problems.txt)")
    for pr in problems:
        print("  -", pr)


if __name__ == "__main__":
    main()
