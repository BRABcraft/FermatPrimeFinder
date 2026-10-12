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

Records H1–H8 are filled in as experiments complete; see `research/experiment_log.md`.

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
