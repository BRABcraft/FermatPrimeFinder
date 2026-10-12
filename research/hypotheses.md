# Hypotheses, lemmas, proofs, and their status

Status vocabulary: **PROVED** (complete proof below, known result unless marked NEW),
**VERIFIED (computational, scope stated)**, **CONSISTENT (statistical test passed)**,
**REJECTED**, **OPEN**.

---

## Part I — Foundational lemmas used by the code (all classical; proofs included so the code's correctness arguments are self-contained)

### L1. p | F_n ⇔ ord_p(2) = 2^{n+1} (p odd prime). PROVED.
If p | F_n then 2^{2^n} ≡ −1, so 2^{2^{n+1}} ≡ 1. Hence ord_p(2) divides 2^{n+1}
but does not divide 2^n, so it equals 2^{n+1}. Conversely, if ord_p(2) = 2^{n+1}, then
x = 2^{2^n} satisfies x² ≡ 1 and x ≢ 1. Since p is prime, x ≡ −1. ∎
*Corollary.* A prime divides at most one F_n, so the Fermat numbers are pairwise coprime.
For a composite d, the test 2^{2^n} ≡ −1 (mod d) is still equivalent to d | F_n. This is why the
C searches also report composite divisors, which the dataset builder then filters out.

### L2 (Euler 1747, Lucas 1878). p | F_n, n ≥ 2 ⇒ p ≡ 1 (mod 2^{n+2}). PROVED.
By L1, 2^{n+1} | p−1, so p ≡ 1 (mod 8) because n ≥ 2. Then 2 is a quadratic residue mod p
(second supplement to quadratic reciprocity), so 2^{(p−1)/2} ≡ 1. Hence 2^{n+1} = ord_p(2)
divides (p−1)/2, i.e. 2^{n+2} | p−1. ∎
*Use in code.* Writing p = k·2^m + 1 with k odd, every prime factor of F_n has m ≥ n+2.
So for fixed p only n ≤ m−2 needs testing. This makes the exhaustive searches complete.

### L3 (Pépin 1877). For n ≥ 1, F_n prime ⇔ 3^{(F_n−1)/2} ≡ −1 (mod F_n). PROVED.
(⇐) Let q be any prime factor of F_n. Then 3^{(F_n−1)/2} ≡ −1 (mod q), so ord_q(3) = F_n − 1 = 2^{2^n}
(a power of 2, not dividing (F_n−1)/2). Hence F_n − 1 ≤ q − 1, so q = F_n.
(⇒) For n ≥ 1, F_n ≡ 2 (mod 3) and F_n ≡ 1 (mod 4). By reciprocity, (3/F_n) = (F_n/3) = (2/3) = −1.
Euler's criterion then gives 3^{(F_n−1)/2} ≡ −1. ∎

### L4 (Suyama). If F_n = P·C with gcd(P, C) = 1 and C prime, then 3^{F_n−1} ≡ 3^{P−1} (mod C). PROVED.
F_n − 1 = P(C−1) + (P−1). Since C ∤ 3, Fermat gives 3^{C−1} ≡ 1 (mod C). ∎
*Contrapositive (used):* inequality proves C composite. Cost: one powmod mod C, given
the Pépin residue R (3^{F_n−1} = R²).

### L5 (Proth 1878). If N = k·2^m + 1 with k < 2^m and a^{(N−1)/2} ≡ −1 (mod N) for some a, then N is prime. PROVED (standard; Pocklington with the fully factored part 2^m > √N).

### L6. Mod-(2^N+1) reduction. PROVED.
For x = h·2^N + l with 0 ≤ l < 2^N: 2^N ≡ −1, so x ≡ l − h. Used by the Python Pépin implementation.
GMP's `mpn_mul_fft` computes x·y mod (B^{pl}+1) with B = 2^64. For pl = 2^{n−6} this modulus is
F_n, so no reduction step is needed (C implementation).

### L7. No algebraic (cyclotomic/Aurifeuillian) factorization of F_n for n ≥ 2. PROVED (known).
F_n = Φ_{2^{n+1}}(2), and Φ_{2^{n+1}}(x) = x^{2^n}+1 is irreducible over ℚ (Eisenstein's criterion applied to Φ(x+1)).
Aurifeuillian identities for Φ_N(2) arise from 2 = 2·1², and give a factorization only when N ≡ 4 (mod 8).
For example, Φ_4(2) = 5 = (2+2+1)(2−2+1) is the trivial instance. Since 2^{n+1} ≡ 0 (mod 8) for n ≥ 2, none applies.
Directly: the Sophie Germain identity a⁴ + 4b⁴ = (a²+2b²+2ab)(a²+2b²−2ab) would need
2^{2^n} = 4b⁴, i.e. 2^n − 2 ≡ 0 (mod 4). That is false for n ≥ 2. ∎
*Consequence.* There is no "free" structural factorization. Every factor has to be found arithmetically.

### L8. p² | F_n ⇒ 2^{p−1} ≡ 1 (mod p²) (p is a Wieferich prime). PROVED.
If p² | F_n then ord_{p²}(2) = 2^{n+1}, which divides p−1 by L2. ∎

---

## Part II — Research hypotheses (one record each)

Data for all statistical records: the **exhaustive regions** in `results/search/regions.json`.
Every prime factor of F_n (n ≥ 12) in these regions is known, so there is no selection bias.
Statistics are in `results/analysis/summary.json`.

### H1. Random-order-of-2 model for factor *counts*. CONSISTENT.
*Statement.* A prime p = k·2^m+1 (k odd) divides some F_n, n ∈ [12, m−2], with probability (1 − 2^{13−m})/k.
*Motivation.* 2 is a "random" quadratic residue; Chebotarev for fixed-size families.
*Evidence.* Total: 84 observed vs 85.3 expected in R1 ∪ R2a ∪ R2b ∪ R2c (final). 76 vs 73.0 before R2c was added. Counting function
N(2^B) = #{p < 2^B, p | some F_n, n ≥ 5}: observed 1, 1, 2, 6, 11, 13, 15, 17, 20, 27, 31, 36 vs expected
1.1, 2.1, 5.0, 7.2, 9.6, 12.4, 15.5, 18.6, 21.6, 25.0, 28.3, 31.8 for B = 10, 14, …, 54. All p > 0.24.
Distribution over ⌊log₂ m⌋: p = 0.67.
*Proof status.* A heuristic, so it cannot be proved for individual p. The density analogue is Theorem D below.
*Consequence.* N(x) ≈ c·log x, so Fermat-divisor primes are exceedingly sparse. The model predicts N(x) ~ log₂ x.

### H2. Geometric excess law P(e = t) = 2^{−(t+1)}, e = m − n − 2. CONSISTENT; density version PROVED (Theorem D).
*Evidence (final, 84 factors).* e = 0, 1, 2, 3, 4, ≥5: observed 39/24/11/1/3/6 vs expected 43.3/21.3/10.5/5.1/2.5/2.5
(shape p = 0.09). MLE of the ratio: r̂ = 0.545, 95% CI [0.47, 0.615] ∋ 1/2.
**Theorem D (density form).** Fix m ≥ 3. Among primes p with v₂(p−1) = m, the relative natural density of
those with v₂(ord_p 2) = m−1−t is 2^{−(t+1)} for 0 ≤ t ≤ m−2, and 2^{−(m−1)} for t = m−1 (odd order).
*Proof.* Let K = ℚ(ζ_{2^m}) and K' = ℚ(ζ_{2^{m+1}}). The primes with v₂(p−1) = m are those split in K
and not split in K' (density 2^{−m} by Dirichlet). For such p and j ≤ m, 2 is a 2^j-th power mod p iff p splits completely
in L_j = K(2^{1/2^j}). This holds because F_p^* ⊃ μ_{2^m} and the 2^j-th power map is well defined; Kummer theory applies since μ_{2^j} ⊂ K.
(a) Degree. Let α = 2^{1/2^j} > 0. Then x^{2^j} − 2 is Eisenstein at 2, so [ℚ(α):ℚ] = 2^j. Since K/ℚ is
Galois, [K(α):K] = [ℚ(α) : ℚ(α) ∩ K]. ℚ(α) ∩ K is a subfield of the abelian field K, so it is normal over ℚ. The subfields of
ℚ(α) are the ℚ(α^{2^i}), and among them only ℚ and ℚ(√2) are normal over ℚ: ℚ(2^{1/4}) and the larger ones are real
fields missing their complex conjugate roots. Since √2 = ζ_8 + ζ_8^{−1} ∈ K (m ≥ 3), ℚ(α) ∩ K = ℚ(√2) for j ≥ 1, and so
[L_j : K] = 2^j / 2 = 2^{j−1}. *(Corrected after adversarial review. The earlier wording argued only that 2^{1/4} ∉ K,
which by itself ignores roots of unity.)* L_j = K(α) is the splitting field of x^{2^j} − 2 over K (μ_{2^j} ⊂ K), hence Galois.
(b) Disjointness: L_j/K is cyclic, so its unique quadratic subextension is K(2^{1/4}) (for j ≥ 2). If K(2^{1/4}) = K',
then ℚ(2^{1/4}) ⊂ K' would be abelian. That is impossible as in (a). Hence L_j ∩ K' = K and [L_j K' : K] = 2^{j}.
(c) Chebotarev in the compositum L_j K': the density of primes that split completely in L_j but not in K',
relative to primes split in K and not in K', is [L_j : K]^{−1} = 2^{−(j−1)}. So, given v₂(p−1) = m,
P(2 is a 2^j-th power) = 2^{1−j} (1 ≤ j ≤ m). Since v₂(ord_p 2) = m − (largest such j), we get
P(v₂(ord) ≤ m−j) = 2^{1−j}, hence P(v₂(ord) = m−1−t) = 2^{−t} − 2^{−t−1} = 2^{−(t+1)}. The remainder at
t = m−1 is 2^{−(m−1)}. ∎
*Check against the literature.* Summing over m with weights 2^{−m}, and adding the m = 1, 2 cases (where 2 may be a non-residue),
gives the density of primes with ord_p(2) even: 1/4 + 1/4 + Σ_{m≥3} 2^{−m}(1 − 2^{1−m}) = 17/24.
This is **Hasse's theorem (1966)**. Numerically, over the 1,270,606 primes < 2·10^7 the fraction is 0.70827 (17/24 = 0.70833),
and the per-m distributions match 2^{−(t+1)} to within sampling error (`results/density_check.json`).
*Status.* PROVED (standard Chebotarev/Kummer argument; almost certainly known in the Hasse–Wiertelak–Moree
literature on divisibility of orders, so **not claimed as new**). The Fermat-factor version (thin families with
fixed k) remains heuristic but is statistically confirmed (H2 evidence above).

### H3. The 1/k law: factor intensity ∝ k^{−α} with α = 1. CONSISTENT (rediscovered from data).
*Evidence (final, 84 factors).* Maximum likelihood with C profiled out: α̂ = 0.975, 95% CI [0.93, 1.02]
(identical under an independent Monte Carlo seed). Joint likelihood-ratio test of (α, r) = (1, 1/2): LR = 2.32, p = 0.31.
Bit-length(k) distribution: shape p = 0.77. (First pass on 76 factors: α̂ = 0.955 [0.915, 1.005].)
*Interpretation.* The exponent 1 comes from "2 is a k-th power residue mod p" having probability 1/k. The data
recover it to within ±5% without being told.

### H4. Small-modulus structure of k beyond the model. NO EVIDENCE (one pre-registered follow-up: H-k3).
*Experiment.* Automated conjecture generator: 369 tests (k mod M for 2 ≤ M ≤ 24, all residues; p mod q for
odd primes q ≤ 31). Observed vs expected under the model, with exact Poisson tails.
*Result (final, after review fixes).* Moduli M ≡ 2 (mod 4) were removed (for odd k they duplicate M/2), leaving
334 tests. The null was **calibrated by simulation**: the identical search on 300 synthetic datasets drawn from the model.
Observed: 5 tests with nominal p < 0.01, against a null mean of 1.46 (P(≥ 5 | null) = 0.053). Smallest p = 8.9·10^{−4}
(k ≡ 5 mod 12: 24 vs 11.0), with family-wise calibrated P(min p ≤ observed) = 0.073. 0 Bonferroni-significant tests.
All structural zeros were respected. All hits are facets of the k mod 3 imbalance (k ≡ 5 mod 12 ⊂ k ≡ 2 mod 3).
*Correction:* the first pass said "11 nominal hits, roughly as expected". The reviewer showed that claim was unsupported
(duplicated tests, no calibrated null). It is replaced by the calibrated statement above.
**Status: borderline (P ≈ 0.05–0.07 family-wise), driven entirely by H-k3, whose mechanism test is negative.**

### H5. Machine-learned features improve on the model score for predicting which primes divide Fermat numbers. NOT DETECTED (low power).
*Design.* Case-control: 66 positives (all factors with p ≤ 600 bits in the exhaustive regions) and 3,960
population-representative primes from the same region. Grouped 5-fold CV by m mod 5. 15 features
(log₂k, m, k mod 3/5/7/8/16, popcount, ω(k), m mod 2/3/4, p mod 3/5/7).
*Result (AUC, 95% bootstrap CI).* Model score −log k: 0.995 [0.985, 1.0]. Logistic on model score: 0.992.
Logistic on model score + all features: 0.991. Gradient boosting on all features: 0.977 [0.939, 1.0]. Random: 0.48.
Gradient boosting feature importance is 99.95% log₂k. **No learnable structure beyond the 1/k law was found.**
*Power caveat (from adversarial review):* with 66 positives and AUC ≈ 0.99 already achieved by −log k, this design can
only detect large effects. "Not detected" here means a weak upper bound, not proof of absence. This
is consistent with Theorem D: the remaining information (power-residue symbols of 2 mod p) is invisible to
low-complexity features of (k, m), and computing it is exactly the cost of the direct test.
*Anomaly detection* (isolation forest) flags only the factors with extreme m (e.g. 431·2^{2099}+1, excess 8).
Their excess values are within model expectations (P(e ≥ 8) = 1/256 per factor; with 76 factors, P(at least one) ≈ 0.26).

### H6. A "theoretical shortcut" to the primality of F33. NOT FOUND (negative result with reasons).
Approaches examined, and why each fails:
1. *Congruence obstructions:* all known constraints (L2, power-residue conditions, Theorem D) concern divisors p.
   They can rule candidate factors in or out, but say nothing about whether F_n itself has a factor at all.
2. *Algebraic factorization:* none exists (L7).
3. *Recurrences* (F_n = ∏F_i + 2 etc.): these determine F_n mod small numbers, which is uninformative about primality,
   since F_n mod q is periodic in n with no prime-specific signal.
4. *Faster Pépin:* Pépin computes 3^{2^{2^n−1}}. Any straight-line exponentiation of a fixed base to the exponent e = 2^{N−1}
   needs at least N−1 multiplications, because an addition chain for e has length ≥ log₂ e.
   So *this exponentiation*, and any test that consists of computing it (Pépin, and Proth/Lucas tests with
   exponent ≈ F_n), costs ≥ 2^n − 1 modular multiplications of 2^n-bit numbers. This is elementary and known.
   It says nothing about tests of an entirely different nature, which is exactly the open question.
5. *Generic primality proofs* (ECPP, APR-CL, AKS) are asymptotically worse than Pépin at this size.
6. *Heuristic outlook:* P(any new Fermat prime) ≈ 2·10^{−8} (`results/heuristic_probability.md`).
*Status.* No route to a shortcut identified. This is not a proof that none exists.

### H-k3 (PRE-REGISTERED 2026-10-12, before inspecting region R2c output)
**Statement under test.** Among prime factors p = k·2^m+1 of F_n (n ≥ 12), the residue k mod 3
deviates from the random-order-of-2 model. The model predicts that the expected share of k ≡ 0 (mod 3) is ≈ 1/2.
**Motivation.** In regions R1 ∪ R2a ∪ R2b (76 factors) the counts for k ≡ 0, 1, 2 (mod 3) were 28/18/30, against
model expectations 37.0/18.5/17.4 (shape G-test p = 0.008). This is the most significant of ~8 tabulated features.
The related congruence-search hits were "k ≡ 5 mod 12" (23 vs 9.5) and "k ≡ 5 mod 11" (17 vs 7.1),
neither Bonferroni-significant (adjusted p = 0.11 and 0.81).
**Pre-registered test.** Use region R2c only (odd k ∈ [65537, 1048575], 14 ≤ m ≤ 1000, n ≥ 12). Its data are
independent of the discovery sample. Compute model expectations for k mod 3 by the same Monte Carlo, then report
(i) the two-sided Poisson p-value for the count with k ≡ 0 (mod 3) and (ii) the shape G-test over the 3 classes.
Decision rule: H-k3 is *supported* if (ii) has p < 0.01 with the deviation in the same direction (k ≡ 0 deficient,
k ≡ 2 in excess); *rejected as a fluctuation* if p > 0.1. Otherwise it stays *inconclusive*.
Secondary (exploratory, no decision rule): the shares of k ≡ 5 (mod 12) and k ≡ 5 (mod 11).
**Prior.** Low. Chebotarev density arguments (Kummer fields ℚ(ζ_{3·2^m}, 2^{1/3}) are linearly disjoint from
ℚ(ζ_{2^{m+1}})) support the model's prediction that 2 is a cubic residue with probability 1/3.
**RESULT (hold-out R2c, run after pre-registration; `results/analysis/hk3_holdout.json`).** 14 factors (model total 13.97).
k ≡ 0/1/2 (mod 3): observed 4/5/5 vs expected 7.10/3.45/3.42. (i) Poisson p for k ≡ 0: 0.33; (ii) shape G-test p = 0.24.
**Decision by the pre-registered rule: REJECTED as a fluctuation** (p > 0.1). Secondary: k ≡ 5 (mod 12) 4 vs 1.67;
k ≡ 5 (mod 11) 0 vs 1.25, so the mod-11 signal did not replicate at all.
**Honest caveat (power analysis done after the fact, which is a lesson):** with only 14 hold-out factors, the
test reaches p < 0.01 in only 8% of simulations even if the discovery-size effect were real, and lands in the "reject"
zone 69% of the time. The hold-out likelihood ratio is ≈ 2.4 in favour of the discovery-size effect, which is weak evidence.
So the data do not settle H-k3. Its status rests on the strong prior: the density of 2 being a cubic residue is 1/3
by Chebotarev and is verified numerically to 0.3% in every m class (E8). **Mechanism test (more powerful than factor counts):** the model's only k-mod-3-specific ingredient is
"2 is a cubic residue mod p" for primes p = k·2^m+1 with 3 | k. Over 20,000 such primes (k < 2^20, 14 ≤ m ≤ 400)
the fraction is 0.3285 (z = −1.45 vs 1/3). Over 5,000 with k < 2^8 it is 0.3366 (z = +0.49)
(`results/analysis/cubic_mechanism.json`). So the mechanism behaves exactly as modelled.
Final full-data table (all 84 exhaustive factors, including the discovery data): 30/22/32 vs 43.1/21.8/20.4. This is not an
independent test. **Overall status: classified as a statistical fluctuation; not confirmed.** Settling it would need ≈ 10× more exhaustive hold-out factors. That costs ~10^2 core-hours with the current CPU code,
and less with GPU sieves.
