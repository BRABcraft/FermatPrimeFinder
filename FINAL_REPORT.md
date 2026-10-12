# Final report — searching for a new Fermat prime

*Session of 2026-10-12. Everything here is reproducible from this repository (see README).*

## 0. Bottom line

* **No new Fermat prime was found, and none is claimed.** The primary mission, a verified prime F_n with n ≥ 5,
  was not achieved.
* The smallest Fermat number of unknown character is **F33** (2^33 + 1 ≈ 8.6·10^9 bits). Measured arithmetic costs
  put a Pépin test of F33 at ≈ **14,000 core-years** on this machine (`results/benchmarks.md` §4).
* Under a model whose predictions we **tested and confirmed** on 84 exhaustively found factors, the probability that
  *any* further Fermat prime exists (33 ≤ n ≤ 160; larger n are negligible) is ≈ **2·10^{−8}** (at most 4·10^{−8} under any
  plausible search depth). This is a heuristic, not a theorem.
* No theoretical shortcut (Track A) was found, and we explain why the obvious routes cannot work (H6).
* Secondary results: a fully verified, provenance-tracked dataset (127 prime factors); certified exhaustive factor
  searches; a confirmation, from data without selection bias, of the random-order-of-2 model (1/k law, geometric excess law,
  counting function); a proof of the density form of the excess law, which reproduces Hasse's 17/24; a ≈2× faster
  Pépin/Suyama implementation; Pépin and Suyama results reproduced for F1–F21. None of the factors is new.

## 1. What is new versus known

| Item | Status |
|---|---|
| Fermat factors found | **all previously known** (rediscovered; none new). This was expected: the regions we could afford were searched decades ago. |
| Exhaustive certificates (every prime factor p < 2^54 of F_n, n ≥ 12; every factor k·2^m+1 with k < 2^20, m ≤ 1000; k < 2^12, m ≤ 3000) | independent re-verification of known completeness, with our own code. Not new in scope (FermatSearch went far deeper) but independently reproduced. |
| Test of the random-order-of-2 model on **exhaustive, bias-free** data, incl. MLE "rediscovery" of α = 1 and r = 1/2 | to our knowledge a cleaner test than fits to the selection-biased literature table. **Modest; not claimed as new theory.** |
| Theorem D (density of v₂(ord_p 2) given v₂(p−1)) | **known in substance** (Hasse-type densities); our proof is self-contained and reproduces Hasse's 17/24. Not claimed as new. |
| Heuristic P(new Fermat prime) ≈ 2·10^{−8} conditioned on certified regions | a re-derivation. It is ≈ 20–40× larger than Boklan–Conway's headline < 10^{−9}. The discrepancy is unexplained because their text was not accessible. |
| `mpn_mul_fft` for arithmetic mod F_n (≈2× over square+reduce) | known idea (Crandall–Fagin weighted transforms); a practical verified engineering result. |

## 2. Proved results (complete proofs in `research/hypotheses.md`)
* L1–L8: the foundational lemmas the code relies on. These are p | F_n ⇔ ord_p 2 = 2^{n+1}; Euler–Lucas p ≡ 1 mod 2^{n+2};
  Pépin; Suyama; Proth; mod-(2^N+1) reduction; no Aurifeuillian/algebraic factorization of F_n for n ≥ 2;
  p² | F_n ⇒ p Wieferich. All are classical.
* **Theorem D.** Among primes with v₂(p−1) = m ≥ 3, the share with v₂(ord_p 2) = m−1−t is exactly 2^{−(t+1)}
  (0 ≤ t ≤ m−2). The proof uses Chebotarev, Kummer theory, and the fact that ℚ(2^{1/2^j}) ∩ ℚ(ζ_{2^m}) = ℚ(√2).
  Corollary: Hasse's density 17/24 of primes with ord_p 2 even. The proof was corrected once after adversarial review.
* Elementary lower bound: computing 3^{(F_n−1)/2} by any multiplication chain needs ≥ 2^n − 1 multiplications. This
  holds only for exponentiation-based tests.

## 3. Computational results (exact scope)
1. **Dataset** (`data/`): 127 prime factors of Fermat numbers. Each row is verified by exact divisibility (2^{2^n} ≡ −1 mod p)
   and a primality proof (BPSW below 2^64; Proth certificate or PARI/GP APR-CL above). Every row is re-checked by an
   independent validator. F5–F11 are verified to be completely factored (product = F_n, all cofactors up to P564 proven prime).
   Rejected claims are logged with reasons in `data/build_problems.txt`. Three were wrong recollections, including
   a garbled F29 factor and a wrong F28 factor; three were composite divisors. Literature factors too
   costly to verify (n = 414677, 3943261, 5798447) are kept separate in `data/unverified_reported_factors.csv`.
2. **Certified exhaustive searches** (`results/search/regions.json`):
   * R1: all p < 2^54 dividing any F_n, n ≥ 12 (1.1·10^12 candidates, 58 min on 4 cores).
     In particular **F28 has no prime factor < 2^54**. Its smallest known factor, 25709319373·2^36+1 (verified here), is ≈ 2^70.6.
   * R2a/R2c: all factors k·2^m+1 with odd k ≤ 2^20 − 1, 14 ≤ m ≤ 1000. R2b: odd k ≤ 4095, 1001 ≤ m ≤ 3000.
   * As a by-product, for example, F33, F34, F35, F37 and F40 have no prime factor below 2^54 and none of the form
     k·2^m+1 with odd k < 2^20 and m ≤ 1000 (the latter bound reaches p ≈ 2^1020). That is a small fraction of what
     FermatSearch has excluded, but independently certified here.
3. **Pépin + Suyama** (`results/pepin/`): F1–F4 prime; F5–F21 composite. Each residue was computed by two
   implementations for n ≤ 20 (FFT mod F_n, and full square + reduce). Every residue was checked mod every known prime factor:
   R ≡ 3^{2^{2^n−1} mod (p−1)} (mod p). Suyama's test proves the cofactors of F12–F19 composite (C1133, C2391, C4880, C9808,
   C19694, C39395, C78884, C157770, relative to the factors in our dataset). The F20 result reproduces Young–Buell (1988), which
   showed F20 is composite without exhibiting a factor.
4. **Model tests** (`results/analysis/summary.json`, 84 bias-free factors):
   * total count 84 vs 85.3 expected;
   * counting function N(2^B) tracks the model for B = 10…54 (36 vs 31.8 at 2^54);
   * MLE α̂ = 0.975 [0.93, 1.02] (model: 1) and r̂ = 0.545 [0.47, 0.615] (model: 1/2); LR p = 0.31;
   * Theorem D checked on all 1.27M primes < 2·10^7.
5. **Benchmarks** (`results/benchmarks.md`): FFT path ≈2× faster than square+reduce and ≈10× faster than generic mod for
   n ≥ 18. Search throughput is 3.2·10^8 candidates/s.

## 4. Conjectures and heuristics (not proved)
* **Random-order-of-2 model** for Fermat factors (H1–H3): consistent with every test we ran.
* **Counting function:** #{primes p ≤ x dividing some Fermat number} ~ c·log x (model). Our data fit it up to 2^54.
  The proven upper bounds known to us are far weaker (power-of-x type). We cite this from memory and could not check it.
* **No further Fermat primes:** heuristic probability ≈ 2·10^{−8} (`results/heuristic_probability.md`).

## 5. Negative results
* **H6 (theoretical shortcut):** every constraint we could derive or find concerns *candidate divisors*, not F_n itself.
  Recurrences give only periodic residues, there is no algebraic factorization, and exponentiation-based tests are bounded below by 2^n − 1
  multiplications. No route to deciding F33 was identified.
* **H5 (ML):** gradient boosting over 15 arithmetic features (AUC 0.985) does not beat the one-parameter model score −log k
  (AUC 0.995). Feature importance is 98% log₂k. With 72 positives this has low power. "Not detected" is the accurate
  phrasing. This matches theory: the remaining information is the power-residue symbol of 2 mod p, which is exactly what the
  direct test computes.
* **H-k3 (k mod 3 anomaly):** a p = 0.008 signal in the discovery data was pre-registered, then failed its hold-out test
  (p = 0.24). A direct mechanism test showed 2 is a cubic residue for 0.3285 of 20,000 primes k·2^m+1 with 3 | k, as
  predicted. It is classified as a fluctuation. The hold-out was underpowered, which is noted as a methodological lesson.
* **Automated congruence search** (334 tests, simulation-calibrated null): 0 Bonferroni-significant hits. Family-wise
  P ≈ 0.05–0.07, driven entirely by H-k3.

## 6. Errors made and caught (kept for transparency)
* Three recalled "literature" factors were wrong (F29, F28; composite or non-dividing). Exact verification caught all three.
* A real off-by-one in `fsearch_gmp` (k = KMAX never tested) was caught by the adversarial reviewer. It was fixed, the gaps were
  re-run (no factors), and a regression test was added.
* A proof gap in Theorem D, step (a), was corrected. Several numerical slips in the documents were corrected
  (naive sum 0.18, not 0.09; gain factor (n+2)e^γ ln2/2). An unsupported "as expected" claim about the congruence search
  was replaced by a simulation-calibrated null.

## 7. Limitations
* **Literature access was blocked** in this environment (prothsearch.com, fermatsearch.org, OEIS, arXiv, Wikipedia).
  Status facts come from search snippets and memory, and are labelled as such. The full ~375-factor table was not imported.
  Our F20/F21 residues were not compared with published Selfridge–Hurwitz residues.
* Both Pépin implementations rely on GMP multiplication. Independence comes from the residue checks against known factors,
  not from different arithmetic libraries.
* The model's 1/k factor is a heuristic for thin families. Only its 2-power part (Theorem D) is proved, and that only in density form.
* Benchmarks were taken on a shared, noisy VM (±30%).

## 8. Open problems
* Is there *any* property decidable faster than Pépin that is necessary for F_n to be prime and fails for F33?
* Are there infinitely many composite Fermat numbers? Is every F_n square-free? (Both open.)
* What explains the factor-of-20–40 gap between our conditional estimate and Boklan–Conway's? (Needs their text.)

## 9. Next steps
See `research/next_steps.md`. Highest value: restore literature access and diff our certified regions against Keller's
table, and compare our F20/F21 residues with published values. More CPU spent on factor hunting has low expected
yield and is dominated by existing GPU projects. More CPU spent on primality of F33 is infeasible.

## Sources
Status facts were obtained via web search (snippets) on 2026-10-12. Sources include fermatsearch.org (news and factors pages),
prothsearch.com (W. Keller, *Fermat factoring status*), OEIS A093179/A019434, mersenne.org (F22 record),
rieselprime.de Prime-Wiki, caramel.loria.fr/f12.txt, AMS Math. Comp. 72 (2003) 1555–1572 (Crandall–Mayer–Papadopoulos, F24),
and arXiv:1605.01371 (Boklan–Conway). Full list with provenance labels: `research/literature_review.md`.
