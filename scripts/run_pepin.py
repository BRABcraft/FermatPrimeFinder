#!/usr/bin/env python3
"""Pépin test of F_n plus Suyama cofactor test, with two independent
implementations, writing results/pepin/F{n}.json.

  impl A: src/c/pepin_fft  (GMP mpn_mul_fft, products directly mod 2^N+1)
  impl B: fermat.core.pepin_residue (full gmpy2 squaring + shift/subtract reduction)

Usage: run_pepin.py n [--impl A|B|AB]
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
import tempfile
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
import gmpy2  # noqa: E402
from gmpy2 import mpz  # noqa: E402
from fermat.core import fermat, pepin_residue  # noqa: E402


def residue_A(n: int):
    with tempfile.NamedTemporaryFile(delete=False) as tf:
        path = tf.name
    t = time.time()
    out = subprocess.run([os.path.join(ROOT, "src", "c", "pepin_fft"), str(n), path],
                         capture_output=True, text=True, check=True)
    el = time.time() - t
    with open(path, "rb") as f:
        raw = f.read()
    os.unlink(path)
    # mpz_out_raw: 4-byte big-endian signed size, then big-endian magnitude
    size = int.from_bytes(raw[:4], "big", signed=True)
    R = mpz(int.from_bytes(raw[4:4 + abs(size)], "big"))
    return R, el, out.stdout.strip()


def residue_B(n: int):
    t = time.time()
    R = pepin_residue(n)
    return R, time.time() - t


def known_factors(n: int, roles=("prime factor",)):
    path = os.path.join(ROOT, "data", "fermat_factors.csv")
    ps = []
    if os.path.exists(path):
        for r in csv.DictReader(open(path)):
            if int(r["n"]) == n and r["role"] in roles:
                ps.append(mpz(r["p"]))
    return ps


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("n", type=int)
    ap.add_argument("--impl", default="AB")
    a = ap.parse_args()
    n = a.n
    F = fermat(n)
    rec = {"n": n}
    Rs = {}
    if "A" in a.impl and n >= 8:
        R, el, line = residue_A(n)
        Rs["A"] = R
        rec["time_A_s"] = round(el, 3)
        rec["impl_A_output"] = line
    if "B" in a.impl or n < 8:
        R, el = residue_B(n)
        Rs["B"] = R
        rec["time_B_s"] = round(el, 3)
    vals = list(Rs.values())
    if any(v != vals[0] for v in vals):
        rec["verdict"] = "MISMATCH"
        print(json.dumps(rec, indent=1))
        sys.exit(2)
    R = vals[0]
    rec["impl"] = "+".join(sorted(Rs)) + (" (agree)" if len(Rs) > 1 else "")
    rec["verdict"] = "PRIME" if R == F - 1 else "COMPOSITE"
    rec["res64"] = format(int(R % (1 << 64)), "016x")
    rec["SH_residues"] = [int(R % (1 << 36)), int(R % ((1 << 36) - 1)), int(R % ((1 << 35) - 1))]
    # Suyama: if C = F/P is prime then 3^(F-1) == 3^(P-1) (mod C).
    ps = known_factors(n)
    allp = known_factors(n, roles=("prime factor", "prime cofactor"))
    fully = len(allp) > len(ps)
    # Cheap independent check of the whole squaring chain: for a known prime p | F_n,
    # R mod p must equal 3^((F_n-1)/2 mod (p-1)) mod p, where (F_n-1)/2 = 2^(2^n - 1).
    if allp:
        checks = []
        for p in allp:
            e = gmpy2.powmod(2, (1 << n) - 1, p - 1)
            checks.append(bool(R % p == gmpy2.powmod(3, e, p)))
        rec["residue_mod_known_factors_consistent"] = all(checks)
        if not all(checks):
            rec["verdict"] = "INCONSISTENT"
            print(json.dumps(rec, indent=1))
            sys.exit(3)
    if fully:
        rec["suyama"] = "not needed: fully factored (all factors proven prime)"
    elif ps and rec["verdict"] == "COMPOSITE":
        P = mpz(1)
        for p in ps:
            P *= p
        assert F % P == 0
        C = F // P
        A = (R * R) % F % C
        B = gmpy2.powmod(3, P - 1, C)
        rec["known_factors"] = [str(p) for p in ps]
        rec["cofactor_digits"] = len(str(C))
        if A != B:
            rec["suyama"] = f"cofactor C{len(str(C))} proven composite (Suyama test, this work)"
        else:
            rec["suyama"] = f"cofactor C{len(str(C))} is a base-3 probable prime (Suyama) - needs proof"
    os.makedirs(os.path.join(ROOT, "results", "pepin"), exist_ok=True)
    with open(os.path.join(ROOT, "results", "pepin", f"F{n}.json"), "w") as f:
        json.dump(rec, f, indent=1)
    print(json.dumps(rec))


if __name__ == "__main__":
    main()
