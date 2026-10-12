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

## E5 — Exhaustive search R1: all prime factors p < 2^54 of F_n, n ≥ 12 (fsearch64, 58 min)
* 1.0995·10^12 candidates; 1.27·10^11 survived the sieve (primes ≤ 16381) and were tested.
* 29 divisors reported; 3 were composite products of F12 factors (rejected by BPSW). The other 26 are prime factors of
  F12, F13, F15–F19, F21, F23, F25–F27, F29, F30, F32, F36, F38, F39.
* **Recall corrections:** the true F29 factor is 1120049·2^31+1 = 2405286912458753. My recalled
  "1120048333814181889" was a garbled memory of it. **F28 has no prime factor < 2^54** (proved by this search),
  so my recalled "25463·2^30+1" (composite) was wrong, and F28's literature factor lies above 2^54.
* Regions certified: see `results/search/regions.json`.

## E6 — Exhaustive rectangles R2a (k < 2^16, m ≤ 1000: 131 s) and R2b (k < 2^12, 1001 ≤ m ≤ 3000: 355 s)
* 59 and 5 factors respectively. These include the classical small-k factors (5·2^39+1 | F36, 3·2^41+1 | F38,
  5·2^1947+1 | F1945, …). R2c (2^16 < k < 2^20, m ≤ 1000) was started as an independent hold-out set.

## E7 — Model tests on the exhaustive data (scripts/run_analysis.py)
* Totals, counting function, excess law, 1/k law, ω(k), m-distribution: all consistent with the model.
* MLE: α̂ = 0.955 [0.915, 1.005], r̂ = 0.53 [0.455, 0.605]; LR test vs (1, 1/2): p = 0.18.
* 369-test congruence search: 0 Bonferroni-significant hits.
* ML: gradient boosting (AUC 0.977) does not beat the model score −log k (AUC 0.995).
* **Anomaly:** k mod 3 table 28/18/30 vs 37.0/18.5/17.4 (p = 0.008, not significant after multiplicity).
  → Pre-registered as H-k3 with a decision rule, to be tested on hold-out R2c. Committed *before* R2c finished.
* Implementation lesson: the first version of the analysis looped in pure Python (≈10^9 iterations) and was
  killed. The vectorized version runs in ~12 min including Monte Carlo sample construction.

## E8 — Theorem D (density form of the excess law) and numerical checks
* Proof written in hypotheses.md (Chebotarev + Kummer; key point: 2^{1/4} ∉ ℚ(ζ_{2^m})).
* Over the 1,270,606 primes < 2·10^7: fraction with ord_p(2) even = 0.70827 vs Hasse's 17/24 = 0.70833; per-m
  excess distributions match 2^{−(t+1)}. P(2 is a cubic residue | 3 | p−1) = 1/3 in every m class (0.331–0.335).

## E9 — Heuristic probability of a new Fermat prime
* P(at least one prime F_n, 33 ≤ n ≤ 160 | our certified regions) = 2.1·10^{−8}. With hypothetical
  k < 2^60 clearance it is 3.7·10^{−8}. Calibration on n = 5..32: expected 0.89 primes, P(none) = 0.36.

## E10 — Pépin + Suyama, F17–F19 (both implementations), F20 (both, running), F21 (FFT only, running)
* F17: 16.9 s (FFT) vs 26.5 s (full square + reduce). F18: 69 s vs 122 s. F19: 353 s vs 619 s (CPU shared with searches).
* All residues agree between implementations and are consistent mod every known factor. Cofactors of F12–F19
  (relative to all factors in the dataset) are proven composite by Suyama's test.
