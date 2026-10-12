"""Exact-arithmetic regression tests.  Run: python3 -m pytest -q tests/"""
import os
import subprocess
import sys

import gmpy2
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
from fermat.core import (divides_fermat, fermat, fermat_index_of_factor,  # noqa: E402
                         is_proven_prime, pepin, pepin_residue, proth_certificate,
                         split_k_m)


def test_fermat_values():
    assert [int(fermat(n)) for n in range(6)] == [3, 5, 17, 257, 65537, 4294967297]


def test_recurrences():
    prod = 1
    for n in range(12):
        assert fermat(n) == prod + 2          # F_n = F_0...F_{n-1} + 2
        prod *= fermat(n)
        if n >= 1:
            assert fermat(n) == (fermat(n - 1) - 1) ** 2 + 1
        if n >= 2:
            assert fermat(n) == fermat(n - 1) ** 2 - 2 * (fermat(n - 2) - 1) ** 2


def test_pairwise_coprime():
    for i in range(10):
        for j in range(i):
            assert gmpy2.gcd(fermat(i), fermat(j)) == 1


def test_euler_641():
    assert fermat(5) == 641 * 6700417
    assert divides_fermat(641, 5) and not divides_fermat(641, 4)
    assert fermat_index_of_factor(641) == 5
    assert split_k_m(641) == (5, 7)


def test_pepin_known_primes_and_composites():
    for n in range(1, 5):
        assert pepin(n)
    for n in range(5, 14):
        assert not pepin(n)


def test_pepin_matches_powmod():
    for n in range(1, 13):
        F = fermat(n)
        assert pepin_residue(n) == gmpy2.powmod(3, (F - 1) // 2, F)


def test_lucas_form_of_factors():
    # every prime factor of F_n (n>=2) is 1 mod 2^(n+2): check on full factorizations
    for n, ps in {5: [641, 6700417], 6: [274177, 67280421310721]}.items():
        prod = 1
        for p in ps:
            assert (p - 1) % (1 << (n + 2)) == 0
            prod *= p
        assert prod == fermat(n)


def test_proth():
    assert proth_certificate(641) is not None
    assert proth_certificate(5 * 2**39 + 1) is not None   # factor of F36
    assert divides_fermat(5 * 2**39 + 1, 36)


def test_is_proven_prime_small():
    assert is_proven_prime(6700417)[0]
    assert not is_proven_prime(114689 * 26017793)[0]


@pytest.mark.skipif(not os.path.exists(os.path.join(ROOT, "src/c/fsearch64")), reason="not built")
def test_fsearch64_matches_bruteforce():
    out = subprocess.run([os.path.join(ROOT, "src/c/fsearch64"), "32", "12"],
                         capture_output=True, text=True).stdout
    got = sorted(l for l in out.splitlines() if l.startswith("FACTOR"))
    exp = []
    for m in range(14, 32):
        for k in range(1, ((1 << 32) - 2 >> m) + 1, 2):
            p = (k << m) + 1
            n = fermat_index_of_factor(p)
            if n is not None and n >= 12:
                exp.append(f"FACTOR {p} {k} {m} {n}")
    # the C search also reports composite divisors; restrict to the prime ones
    got_prime = sorted(l for l in got if gmpy2.is_prime(int(l.split()[1])))
    exp_prime = sorted(l for l in exp if gmpy2.is_prime(int(l.split()[1])))
    assert got_prime == exp_prime


@pytest.mark.skipif(not os.path.exists(os.path.join(ROOT, "src/c/pepin_fft")), reason="not built")
def test_pepin_fft_matches_python():
    for n in (8, 10, 12):
        out = subprocess.run([os.path.join(ROOT, "src/c/pepin_fft"), str(n)],
                             capture_output=True, text=True).stdout
        res64 = out.split("res64=")[1].split()[0]
        assert int(res64, 16) == int(pepin_residue(n) % (1 << 64))


@pytest.mark.skipif(not os.path.exists(os.path.join(ROOT, "src/c/fsearch_gmp")), reason="not built")
def test_fsearch_gmp_inclusive_kmax():
    # regression for the off-by-one found in adversarial review: k = KMAX must be tested
    out = subprocess.run([os.path.join(ROOT, "src/c/fsearch_gmp"), "1", "7", "14", "14", "12", "1000"],
                         capture_output=True, text=True).stdout
    assert "FACTOR 7 14 12" in out
    out = subprocess.run([os.path.join(ROOT, "src/c/fsearch_gmp"), "3", "3", "41", "41"],
                         capture_output=True, text=True).stdout
    assert "FACTOR 3 41 38" in out
