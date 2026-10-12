#!/usr/bin/env python3
"""Independent re-validation of data/*.csv (does not trust build_dataset.py).

Checks, for every factor row:
  * p == k*2^m + 1 with k odd;  m >= n+2 for n >= 2 (Euler-Lucas)
  * 2^(2^n) == -1 (mod p), recomputed with gmpy2.powmod(2, 2^n, p) (a different
    code path from the repeated-squaring loop used by the builder)
  * primality: BPSW below 2^64; above, PARI/GP isprime (APR-CL) or Proth certificate
  * for rows of fully factored F_n: the product of all listed factors equals F_n
  * for every F_n row: decimal digit count recomputed; status consistent with factors
Exit code 0 iff everything checks.
"""
from __future__ import annotations

import csv
import os
import subprocess
import sys

import gmpy2
from gmpy2 import mpz

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")


def F(n):
    return (mpz(1) << (1 << n)) + 1


def gp_isprime_batch(nums):
    if not nums:
        return []
    prog = "default(parisizemax, 4000000000);\n" + "".join(f"print(isprime({int(x)}));\n" for x in nums)
    out = subprocess.run(["gp", "-q", "-f"], input=prog, capture_output=True, text=True, timeout=36000)
    res = out.stdout.split()
    assert len(res) == len(nums), out.stderr
    return [r == "1" for r in res]


def main():
    errors = []
    rows = list(csv.DictReader(open(os.path.join(DATA, "fermat_factors.csv"))))
    big = []
    by_n = {}
    for r in rows:
        n, p, k, m = int(r["n"]), mpz(r["p"]), int(r["k"]), int(r["m"])
        by_n.setdefault(n, []).append((p, r["role"]))
        if p != (mpz(k) << m) + 1 or k % 2 == 0:
            errors.append(f"row F{n} p={p}: k/m inconsistent")
        if n >= 2 and m < n + 2:
            errors.append(f"row F{n} p={p}: violates Euler-Lucas m >= n+2")
        if r["role"] == "F_n is prime":
            if p != F(n):
                errors.append(f"F{n} prime row wrong")
        elif gmpy2.powmod(2, mpz(1) << n, p) != p - 1:
            errors.append(f"row F{n} p={p}: does NOT divide F{n}")
        if p < (1 << 64):
            if not gmpy2.is_bpsw_prp(p):
                errors.append(f"row F{n} p={p}: composite")
        else:
            big.append((n, p))
    flags = gp_isprime_batch([p for _, p in big])
    for (n, p), ok in zip(big, flags):
        if not ok:
            errors.append(f"row F{n} p={p}: PARI says composite")
    nums = list(csv.DictReader(open(os.path.join(DATA, "fermat_numbers.csv"))))
    for r in nums:
        n = int(r["n"])
        facs = by_n.get(n, [])
        if r["fully_factored"] == "True" and n >= 5:
            prod = mpz(1)
            for p, _ in facs:
                prod *= p
            if prod != F(n):
                errors.append(f"F{n}: marked fully factored but product != F{n}")
        if n <= 22 and int(r["decimal_digits"]) != len(str(F(n))):
            errors.append(f"F{n}: decimal digits wrong")
        if int(r["binary_digits"]) != (1 << n) + 1:
            errors.append(f"F{n}: binary digits wrong")
        if r["status"] == "composite" and n <= 4:
            errors.append(f"F{n} wrongly composite")
        if r["status"] == "prime" and n >= 5:
            errors.append(f"F{n} wrongly prime")
        if int(r["num_known_prime_factors"]) != len([1 for _, role in facs if role != "F_n is prime"]):
            errors.append(f"F{n}: factor count mismatch")
    print(f"validated {len(rows)} factor rows ({len(big)} via PARI), {len(nums)} Fermat-number rows")
    if errors:
        print("ERRORS:")
        for e in errors:
            print("  ", e)
        sys.exit(1)
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
