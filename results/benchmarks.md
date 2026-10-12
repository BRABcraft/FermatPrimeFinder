# Benchmarks (Track C)

Hardware: 4-core Intel Xeon @ 2.8 GHz (cloud VM), GMP 6.3.0. Timings were taken while other jobs shared the CPU.
Treat them as ±30%. The ratios between methods measured in the same run are more reliable than absolute values.

## 1. One modular squaring mod F_n (`src/c/bench_sqr.c`, raw output in `results/bench_sqr_raw.txt`)

| n | bits of F_n | (A) `mpn_mul_fft` mod 2^N+1 | (B) `mpz_mul` + shift/subtract | (C) `mpz_mul` + `mpz_mod` | B/A | C/A |
|---|---|---|---|---|---|---|
| 20 | 1,048,576 | 1.25 ms | 2.78 ms | 15.3 ms | 2.2 | 12 |
| 22 | 4,194,304 | 7.97 ms | 18.3 ms | 93.4 ms | 2.3 | 12 |
| 24 | 16,777,216 | 36.5 ms | 80.4 ms | 436 ms | 2.2 | 12 |
| 25 | 33,554,432 | 116 ms | 233 ms | 956 ms | 2.0 | 8 |
| 26 | 67,108,864 | 245 ms | 468 ms | 2456 ms | 1.9 | 10 |

All three methods produce identical residues ("agree=1"). For n ≤ 16 the FFT path is *slower* than B, because GMP uses
Toom–Cook below its FFT threshold.

**Finding.** GMP already contains a negacyclic Schönhage–Strassen routine (`mpn_mul_fft`) that multiplies directly
modulo B^{pl}+1. For pl = 2^{n−6} that modulus *is* F_n. Calling it directly halves the transform length relative to
"square, then reduce" and gives a **≈2× speed-up for Pépin/Suyama arithmetic at n ≥ 18** with no new mathematics.
This is the classical observation behind Crandall–Fagin's weighted transforms (Mayer's and Crandall's F24 code), so it is
**not new**, but it is a practical, verified improvement over naive GMP use. The caveat: `mpn_mul_fft` is an internal GMP
symbol (exported in libgmp.so but not part of the documented API), so a future GMP release could change it.

## 2. Full Pépin tests (`scripts/run_pepin.py`; A = FFT path, B = gmpy2 square + special reduction)

| n | squarings | A (s) | B (s) | B/A | verdict | consistency checks |
|---|---|---|---|---|---|---|
| 16 | 65,535 | 3.7 | 4.8 | 1.3 | composite | A = B; R mod each known p correct |
| 17 | 131,071 | 16.9 | 26.5 | 1.6 | composite | ✓ |
| 18 | 262,143 | 69 | 122 | 1.8 | composite | ✓ |
| 19 | 524,287 | 353 | 619 | 1.8 | composite | ✓ |
| 20 | 1,048,575 | see results/pepin/F20.json | | | | |
| 21 | 2,097,151 | see results/pepin/F21.json | (A only) | | | |

## 3. Exhaustive factor search throughput

| Tool | Region | Candidates | Tested after sieve | Wall time (4 cores) | Rate |
|---|---|---|---|---|---|
| `fsearch64` (64-bit Montgomery) | p < 2^54, n ≥ 12 | 1.10·10^12 | 1.27·10^11 | 3458 s | 3.2·10^8 candidates/s |
| `fsearch_gmp` | k < 2^16, m ≤ 1000 | 3.3·10^7 | 2.6·10^6 | 131 s | — |
| `fsearch_gmp` | 2^16 < k < 2^20, m ≤ 1000 | 4.9·10^8 | 3.9·10^7 | 3872 s (shared CPU) | — |
| `fsearch_gmp` | k < 2^12, 1001 ≤ m ≤ 3000 | 4.1·10^6 | 3.3·10^5 | 355 s | — |

Validation: both tools reproduce an independent Python brute force exactly on small bounds. The adversarial reviewer
added an unsieved C brute force. After the off-by-one fix, a regression test covers the boundary k = KMAX.

## 4. What it would cost to decide F33 by Pépin's test

Per-squaring time grows by ≈ 2.1–2.3× per unit of n here (the size doubles, plus the log factor). Extrapolating
from n = 26 (0.245 s) with ×2.15 per step gives ≈ **52 s per squaring at n = 33**. Each operand is 8.6·10^9 bits = 1 GiB,
and the FFT working set is several GiB, which is at or beyond this VM's 15 GB RAM. Pépin needs 2^33 − 1 ≈ 8.6·10^9 squarings:

  8.6·10^9 × 52 s ≈ 4.5·10^11 s ≈ **14,000 core-years** on this hardware.

Perfect parallel scaling over 4 cores gives ≈ 3,500 years. A modern multi-GPU or large-cluster implementation
(e.g. a DWT-based code like Mayer's Mlucas, with Gerbicz error checking) might be 10^3–10^4× faster, which still means years of
dedicated hardware. **Conclusion: F33 is not decidable within this project by any amount of credit-funded
reasoning. The bottleneck is raw arithmetic, not cleverness, unless a genuinely new theorem appears (Track A found none).**
