# Adversarial review (2026-10-12)

Scope: items 1–6 of the review brief. No existing file was modified. The scratch checks were run single-core in the session
scratchpad, outside the repo: `brute.c`, `edge.py`, `mont.c` and `fft.c`.

## 1. Completeness of the exhaustive searches

**R1 (`fsearch64`): OK.**
* m range `NMIN+2 … B-1` (fsearch64.c:61). `kmax = (2^B-2)>>m` (l.62) is the exact bound for k·2^m+1 < 2^B. `nodd=(kmax+1)/2` counts the odd k exactly. The total 2^40−1 matches the log.
* j range `NMIN … m-2` (l.100–109). The check runs before each squaring, so j = m−2 is tested. n ≥ 52 cannot have a factor below 2^54.
* Sieve index (l.75–84): `i0 = (k0−1)/2 mod q` and `start = (i0 − seg0) mod q` are correct. `k<<m < 2^B`, so nothing overflows.
* The `p != q` guard (l.88) is correct. With the default QMAX it is never needed.
* Montgomery: `mont_inv` reaches 96 bits of precision. t + m·p < 2^128 and r < 2p hold for p < 2^63. `target = p−R` is mont(−1).
* The OpenMP parallel loop has a correct reduction, per-thread buffers and a critical-section `printf`.
* *Tests:* I wrote an independent brute force with no sieve, `u128 %` arithmetic and deterministic MR. It iterates over every p ≡ 1 (mod 2^{NMIN+2}).
  * B=34, NMIN=5 with QMAX ∈ {3, 16381, 10^6}: the same 15 primes as fsearch64. This includes 641, 114689 and 274177, which exercise the p==q branch. With m=7 there are 16 segments, so segment boundaries are crossed.
  * B=40, NMIN=12: the same 10 primes.
  * A Montgomery fuzz test, 2·10^7 random products with p of 40–62 bits: 0 errors.
* The composite FACTOR lines (2983954661377 and two others) are filtered out by `build_dataset.py`. The observed set equals the prime-filtered search output exactly (76 = 76).

**R2a/R2b/R2c (`fsearch_gmp`): ISSUE (medium). Off-by-one at the upper k bound.**
`nodd = (KMAX−KMIN+1)/2` (fsearch_gmp.c:159) with odd KMIN and odd KMAX tests only k ≤ KMAX−2.
* Failing input: `fsearch_gmp 1 7 14 14 12 1000` does **not** find 7·2^14+1 = 114689 | F12. `fsearch_gmp 1 8 …` does find it.
* The header comment ("KMIN ≤ k < KMAX"), the README and E2/E6 say "k < K". But `regions.json` (R2a "k odd ≤ 65535", R2b "≤ 4095"), `regions.py:30` (`return r["kmax"]`, used as inclusive) and H-k3 claim inclusive bounds.
* Untested columns:
  * k=65535 for 39 ≤ m ≤ 1000 (R1 covers m ≤ 38).
  * k=4095 for 1001 ≤ m ≤ 3000.
  * k=1048575 in R2c. The "union gives all odd k ≤ 2^20−1" note is false twice over, since it misses 65535 and 1048575.
* *Patch check:* I tested k=65535 (m 39..1000) and k=4095 (m 1001..3000) directly. They contain 17 and 10 primes, none of which divides any F_n with n ≥ 12. So the stated R2a/R2b certificates are true after this check, and no data changes. Fix: either use `(KMAX−KMIN)/2+1`, or set kmax to 65533/4093 in `regions.json` and rerun the edge k.

Latent (low): fsearch_gmp.c:178–181 skips the p==q check for m ≥ 21. That is correct only when QMAX < 2^21+1. A user passing QMAX ≥ 2^21 could sieve out a prime p = k·2^m+1 ≤ QMAX. All runs used QMAX = 10^6.

## 2. Theorem D: OK with a proof gap (low)

* Degree, (b) disjointness, (c) Chebotarev, and the telescoping 2^{−t}−2^{−t−1}: all correct. The t = m−1 remainder is also correct.
* I recomputed the Hasse check: 1/4 + 1/4 + Σ_{m≥3}2^{−m}(1−2^{1−m}) = 1/2 + (1/4 − 1/24) = 17/24. Correct.
* **Gap in (a):** "2^{2^{j−2}} ∉ K^{*2^j} because 2^{1/4} ∉ K" does not follow. From 2^{2^{j−2}} = β^{2^j} one only gets β^4 = 2ζ with ζ ∈ μ_{2^{j−2}} ⊂ K. Roots of unity must be excluded, and (1+i)^4 = −4 shows they matter.
  * Clean fix: K/ℚ is Galois, so [K(2^{1/2^j}):K] = [ℚ(2^{1/2^j}) : ℚ(2^{1/2^j})∩K]. The intersection is a normal subfield of ℚ(2^{1/2^j}), hence ℚ or ℚ(√2). So the degree is 2^{j−1}.
* Correctly labelled not new. It is the classical computation behind Hasse (1966); Moree's papers on "d | ord_p(g)" contain equivalent densities. I cannot cite an exact theorem number offline.
* Minor: `density_check.py` prints "predicted 2^{−(t+1)}" even at t = m−1, where the theorem gives 2^{−(m−1)}. For example, m=3, t=2 is shown as observed 0.2501 vs "predicted" 0.125. The data agree with the theorem, not with the printed column.

## 3. Heuristic probability: OK (headline robust)

* P = (e^γ/2)·exp(−Λ_unsearched) is the standard Wagstaff-type conversion. P(N prime | no factor < y) ≈ e^γ ln y/ln N, against ln y/ln √N from the Poisson sum, gives a ratio of e^γ/2 ≈ 0.89. The direction is right.
* The odd-k integral is correct: dk/2 odd k, times density 2/ln p, gives ∫dk/(k ln p).
* The k < 4096 exact part and the integral start at a = 4095 consistently in both S_all and S_searched, so the difference is the pure integral.
* Sanity check: P(F33) ≈ e^γ·54/2^33 = 1.12e−8, which matches the table. The calibration also checks out: an exact prime sum gives P(F5) = 0.373 and P(F6) = 0.172, the same as the code.
* The m cutoff at n+89 and the extra k=65535 column are negligible. No bug changes the headline by ≥3×.
* ISSUE (low): `heuristic_probability.md` says the gain over the naive model is "≈(n+2)·e^γ/2". It is actually e^γ(n+2)ln2/2 ≈ 0.62(n+2).
* ISSUE (low–medium): Boklan–Conway comparison. The .md attributes the 20–40× gap to B&C not crediting the absence of factors below 2^{n+2}. But `literature_review.md` §4 says B&C *do* argue that the special form of factors makes F_n more likely prime. The two statements conflict, and the attribution is speculation. CANNOT VERIFY B&C's method offline; the attribution should be removed or marked as a guess.
* `literature_review.md` §4: the naive Σ_{n≥5} 2/ln F_n is 2/ln2·2^{−4} ≈ **0.18**, not 0.09.
* "Chebotarev … make this exact on average" overstates. Theorem D covers only the 2-adic (excess) factor. The 1/k factor has no density proof here.

## 4. Statistical analysis

* **Region identity: OK.** `load_observed` and `build_sample` use the same `union_mkmax(regions, 12, 3000)` with k ≤ mk[m]. Both inherit the inclusive-kmax error from item 1, which is empirically harmless.
* **Estimator: OK.** Dyadic blocks with uniform with-replacement sampling and fac = N_odd/s give an unbiased Horvitz–Thompson estimator. Seeds 1 and 2 give 73.0 vs 72.7.
* **Mass and excess: OK.**
  * mass = w(1−2^{−(emax+1)}) is correct.
  * The tail Σ_{e=cap}^{emax}2^{−(e+1)} = 2^{−cap}−2^{−(emax+1)} is correct.
  * The e ≤ 60 truncation drops less than 2^{−61}.
* **MLE: OK.** The normaliser is (1−r^{emax+1}), consistent with truncation at emax. C is correctly profiled: ll = N log N − N log Λ + …
  * Caveat: the α CI upper end of 1.005 sits at grid resolution, and the normaliser carries MC noise. "Consistent with 1" is borderline and should be shown to be seed-robust.
* **Congruence search: ISSUE (low).** k is always odd, so "k mod 2M" with M odd duplicates "k mod M" (M = 2, 6, 10, 14, 18, 22). That makes ~36 of the 369 tests exact copies (e.g. top-10 entries for mod 11/22 and mod 3/6), and many others are nested. The claim "11 nominal hits, roughly as expected" is unsupported: the naive expectation is ≤3.7, and the excess comes from correlated and duplicated tests. Bonferroni is over-conservative for the same reason. Dedupe, or use a permutation/max-statistic calibration.
* **ML: ISSUE (low).**
  * Negatives come from `exact` entries, i.e. ≤600 bits **or k<64**. Positives are restricted to ≤600 bits. The domains differ, but the k<64, m>600 strata have fac≈1 and are almost never drawn, so the practical effect is nil.
  * The real weakness is power. Negatives are dominated by k≈2^30–2^40 and positives by small k, so AUC≈0.99 just restates the 1/k law. A k-matched design would be needed to detect residual structure. The H5 "REJECTED" means "not detected with low power".
  * The H5 sentence "consistent with Theorem D" is a non sequitur.

## 5. Pépin/Suyama: OK

* `mpn_mul_fft` semantics: I checked 3000 squarings per n = 8…14 against `mpz` (0 errors). Carry = 1 occurred only for x² ≡ −1, always with the low limbs zero. So the `top` → "x = −1, next square = 1" shortcut (pepin_fft.c:40–45, 57) is valid. The iteration count 2^n−1 is correct.
* Python reduction (core.py:139–143): the range [0, 2^N] is kept. Proth (l.87–94) and Suyama (run_pepin.py:186–193; F−1 = P(C−1)+(P−1)) are correct.
* The residue check against known factors (e = 2^{2^n−1} mod (p−1)) is correct.
* ISSUE (low): impl A and impl B are not fully independent. gmpy2's big `mpz_mul` also dispatches to GMP's FFT. The mod-p residue checks partially compensate.
* `data/fermat_numbers.csv` is stale: the F6 and F8 rows say the cofactor "is a base-3 probable prime – needs proof", but both are fully factored with proven-prime cofactors. This happens because build_dataset.py:194 prefers the Pépin note over `full`.

## 6. Overclaims and other text

* `literature_review.md`: "Section 7 of FINAL_REPORT.md" points to a file that does not exist yet. `results/heuristic_probability.md` cites `results/benchmarks.md`, which is also missing.
* README: "none can be tested with present technology" is too broad. It should say F33 and above cannot be tested with Pépin today.
* README's reproduce line writes `gmp_k16_m1000.txt`, not the actual file name.
* H6.4: the addition-chain lower bound holds only for straight-line exponentiation. It is not a bound on "every test", although the status line does hedge.
* F28 (low; data improvement): `build_dataset.py` lists F28's factor as unrecalled. A candidate from my recall, 25709319373·2^36+1 = 1766730974551267606529, does divide F28 and is a PARI pseudoprime. I checked this with gp. It is consistent with R1, since it is about 2^70.6. I believe it is the Taura (1997) factor in Keller's table, but I could not verify that offline.
* "R2a = k < 2^16" (E6) vs "≤ 65535" (regions.json): inconsistent; see item 1.
* No probable prime is presented as proven: Suyama outcomes are labelled PRP or composite. Recalled values are re-verified.

## Ranked issues

1. **Medium:** the `fsearch_gmp` upper k bound is exclusive while regions.json/regions.py treat it as inclusive. k=65535, 4095 and 1048575 were untested. I closed the k=65535 and k=4095 gaps by direct check (no factors); k=1048575 is still open for R2c.
2. **Low–medium:** the Boklan–Conway discrepancy explanation is speculative and contradicts the literature review.
3. **Low:** proof gap in Theorem D (a), with the fix given above.
4. **Low:** the congruence search contains duplicate and nested tests, so "11 hits as expected" is unjustified.
5. **Low:** the ML case-control design has low power and a minor domain mismatch.
6. **Low:** impl A and impl B are not independent (shared GMP FFT).
7. **Low:** naive-sum arithmetic (0.09 → 0.18). The (n+2)e^γ/2 factor is missing ln2.
8. **Low:** the density_check "predicted" column is wrong at t=m−1. fermat_numbers.csv rows are stale. There are dangling file references. The fsearch_gmp p==q guard is latent for QMAX ≥ 2^21.
