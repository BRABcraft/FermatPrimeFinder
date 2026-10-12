# Experiment log (chronological)

Machine: 4-core Intel Xeon @ 2.8 GHz, 15 GB RAM, Python 3.13 + gmpy2 2.3.2, GMP 6.3.0,
PARI/GP 2.15.4. No paid external compute was used.

## E0 — Environment and source access (2026-10-12)
* Primary status pages (prothsearch, fermatsearch, OEIS, Wikipedia, arXiv) were blocked by the
  sandbox network policy. Only web-search snippets were available.
  **Decision:** every literature number is treated as a claim and re-verified by exact arithmetic.

## E1 — Exhaustive search tool `fsearch64` (p < 2^B, n ≥ 12)
* Design: for every m ≥ 14 and odd k with k·2^m+1 < 2^B, sieve by primes ≤ 16381, then test
  2^{2^j} ≡ −1 (mod p) for every j ∈ [12, m−2] (64-bit Montgomery arithmetic).
  The completeness proof is in hypotheses.md L1/L2.
* Validation: output identical to an independent Python brute force for B = 32, 36.
* Throughput: ≈3·10^8 candidates/s on 4 cores (B = 46 in 14 s).
* Observation: the search also reports *composite* divisors of F12 (114689·26017793 and
  114689·63766529). This is expected, because 2^{2^n} ≡ −1 (mod d) ⇔ d | F_n for any d.
  The dataset builder rejects them after BPSW (deterministic below 2^64).

## E2 — Exhaustive search tool `fsearch_gmp` (rectangle k < K, m ≤ M)
* Validated against Python brute force on k < 4096, 14 ≤ m ≤ 120: identical 22 factors.

## E3 — Dataset build, first pass
* All recalled factors of F5–F11 verified. Cofactors of F7 (P22), F8 (P62), F9 (P99), F10 (P252) and
  F11 (P564) were recomputed as F_n/∏(other factors) and proven prime by PARI/GP APR-CL.
  The products equal F_n exactly, so **F5–F11 are completely factored (verified)**.
* **One recalled claim rejected:** 1120048333814181889 was recalled as a factor of F29. It does not divide
  F29 and is in fact composite. The recall was wrong; the verifier caught it. (Lesson: recall is unreliable,
  so verification is mandatory.)
* Verified literature factors: F12 P54, F13 (3 factors), F14 P54, F15 P33, F16 P27, F17 P49, F18 P23, F19 P35,
  F22 P35, F25 (2), F27, F31 P23, F(251), F(3135).
* Not verifiable here (cost): factors of F(414677), F(3943261), F(5798447) → kept in
  `data/unverified_reported_factors.csv`, excluded from analysis.

## E4 — Pépin tests, two implementations
* Impl A: GMP `mpn_mul_fft` (products mod B^pl+1 = F_n directly). Impl B: gmpy2 full square
  + shift/subtract reduction. Residues agree for all n tested (n = 8…16 so far); see results/pepin/.
* F1–F4 prime; F5–F16 composite (reproduces known results).
* Suyama: cofactors of F12–F16 (relative to factors then in the dataset) are proven composite.
