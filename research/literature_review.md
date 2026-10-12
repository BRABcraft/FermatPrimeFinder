# Literature review: the research frontier for Fermat numbers

*Compiled 2026-10-12. Every numerical fact that enters the dataset is recomputed
locally (see `data/` and `scripts/validate_dataset.py`), so transcription errors here
cannot turn into false claims in the results.*

## 0. Source-access caveat

This session's network policy blocked direct access to the primary status pages
(`prothsearch.com/fermat.html` (W. Keller's *Fermat factoring status*), `fermatsearch.org`,
`oeis.org`, `wikipedia.org`, `mathworld.wolfram.com`, `arxiv.org`). The connections were
refused by the proxy (HTTP 403 / DNS failure). Web-search snippets were reachable. Provenance labels used below:

* **[S]** stated in a web-search snippet of the named source (retrieved 2026-10-12);
* **[K]** standard textbook/literature knowledge (e.g. Křížek–Luca–Somer, *17 Lectures on
  Fermat Numbers*, Springer 2001; Crandall–Pomerance, *Prime Numbers: A Computational
  Perspective*, 2nd ed. 2005), not re-fetched in this session;
* **[V]** verified by our own exact computation in this repository.

## 1. Current status (as far as could be established)

| Item | Status | Provenance |
|---|---|---|
| Known Fermat primes | F0..F4 (3, 5, 17, 257, 65537) | [K], [V] |
| F5 … F32 | all composite | [S: J. D. Cook blog; Prime-Wiki; Wikipedia snippet "as of 2017, F_n composite for 5 ≤ n ≤ 32"] |
| Completely factored | F5 … F11 (7 numbers) | [S: fermatsearch.org factors page "7 complete factorizations"], [V] |
| Smallest unfactored | F12 (6 known prime factors, composite C1133 cofactor) | [S: caramel.loria.fr/f12.txt, Prime-Wiki], [V] |
| Composite with no known factor | F20, F24 | [S: fermatsearch.org "2 composites with unknown factors"] |
| F14, F22 | first factors found 2010 (F14: P54; F22: P35, 2010-03-26) | [S: OEIS A093179, mersenne.org F22 record, Yamasaki page], [V] |
| Smallest F_n of unknown character | **F33** (2,585,827,973 decimal digits) | [S: Prime-Wiki; Kruppa found 46931635677864055013377 \| F31 on 2001-04-12] |
| Number of known prime factors | 374 (fermatsearch.org factors page, "last updated 5 Feb 2026") vs 375 (fermatsearch news, Dec 2025, factor 1085·2^3943263+1 of F(3943261)); the two snippets disagree | [S] |
| Recent finds (2025) | 39701·2^414679+1 \| F(414677) (Batalov, Propper, Jan 2025); 99·2^5798449+1 \| F(5798447) (Cavecchia, Jan 2025); 322072887044529·2^253+1 \| F(251) (Gostin, May 2025); 36812439837·2^3144+1 \| F(3135) (Gostin, Nov 2025); 1085·2^3943263+1 \| F(3943261) (Dinkel / PrimeGrid, Dec 2025) | [S: fermatsearch.org news] |
| Pépin tests of F_n of previously unknown character | F14 (Selfridge–Hurwitz 1963), F20 (Young–Buell 1988), F22 (Crandall–Doenias–Norrie–Young 1995), F24 (Crandall–Mayer–Papadopoulos, Math. Comp. 72 (2003) 1555–1572, two independent FFT/DWT runs with agreeing residues) | [S: AMS abstract], [K] |
| F30 | Pépin test (cofactor compositeness) by E. W. Mayer, 2022 | [S: Prime-Wiki snippet] |

**Consequence for the mission.** Any new Fermat prime must be F_n with n ≥ 33 and n not
in the set of indices with a known factor. F33 has 2^33 ≈ 8.6·10^9 bits. A Pépin test
needs 2^33 − 1 modular squarings of 8.6-gigabit numbers. Section 7 of `FINAL_REPORT.md`
gives a measured cost extrapolation.

## 2. Established theorems (all [K]; short proofs reproduced in `research/hypotheses.md`)

1. **Recurrences.** F_n = F_0 F_1 ⋯ F_{n−1} + 2; F_n = (F_{n−1} − 1)^2 + 1;
   F_n = F_{n−1}^2 − 2(F_{n−2} − 1)^2. Hence the F_n are pairwise coprime (Goldbach 1730).
   This gives a proof that there are infinitely many primes, but says nothing about the
   primality of an individual F_n.
2. **Euler (1747) and Lucas (1878).** If a prime p divides F_n and n ≥ 2, then
   p ≡ 1 (mod 2^{n+2}). Euler's weaker form is p ≡ 1 (mod 2^{n+1}). In fact ord_p(2) = 2^{n+1}. Lucas's extra factor of 2
   holds because 2 is a quadratic residue mod p ≡ 1 (mod 8). Every *divisor* of F_n has the same form.
3. **Pépin (1877).** For n ≥ 1: F_n is prime ⇔ 3^{(F_n−1)/2} ≡ −1 (mod F_n). Bases 5 and 10
   also work for n ≥ 2, and in general any a with Jacobi symbol (a/F_n) = −1 works. The test is a special case of Proth's
   theorem (1878).
4. **Every F_n is a base-2 Fermat pseudoprime or prime**: 2^{F_n−1} ≡ 1 (mod F_n),
   because 2^{n+1} | 2^{2^n}. So base-2 tests are useless here.
5. **Gauss–Wantzel.** A regular N-gon is constructible ⇔ N = 2^a · (product of distinct Fermat primes).
6. **Cyclotomic view.** F_n = Φ_{2^{n+1}}(2), and x^{2^n} + 1 is irreducible over ℚ.
   Aurifeuillian factorizations of Φ_N(2) exist only for N ≡ 4 (mod 8). For N = 2^{n+1} this means
   n = 1, so **F_n has no algebraic factorization for n ≥ 2**.
7. **Square-freeness (open).** p^2 | F_n implies 2^{p−1} ≡ 1 (mod p^2), i.e. p is a
   Wieferich prime. The only Wieferich primes found in searches to very large bounds (Dorais–Klyve 2011:
   below 6.7·10^15; later distributed searches reportedly to 2^64) are 1093 and 3511.
   Neither divides a Fermat number. So no known F_n is divisible by the square of a known prime factor.
8. **Further congruence information.** If p = k·2^m + 1 (k odd) divides F_n, then 2 is a
   (k·2^{m−n−1})-th power residue mod p. For each odd prime q | k, 2 must be a q-th power residue mod p.
   For example, q = 3 requires p = x^2 + 27y^2 (Gauss). If m ≥ n+3, then 2 must be a quartic residue, which requires
   p = x^2 + 64y^2 (Gauss / Dirichlet). These are necessary conditions on a candidate factor and give no information about
   F_n itself.
9. **Suyama's test (1984).** Let F_n = P·C, where P is a product of known primes. If C is prime, then
   3^{F_n−1} ≡ 3^{P−1} (mod C). This is used to prove cofactors composite cheaply, given the Pépin residue.
10. **Known results on prime factors of composite F_n** (Křížek–Luca–Somer): F_n is never a
    perfect power and never a sum of two primes for n ≥ 2. Almost nothing is known
    unconditionally about how many prime factors a composite F_n has.

## 3. Algorithms and their limits [K]

| Method | What it gives | Cost for F_n | Notable use |
|---|---|---|---|
| Sieved trial division over k·2^{n+2}+1 (Proth.exe, FermatSearch, mmff on GPUs) | factors p with small k | ~n squarings mod p per candidate (cheap). The yield per doubling of k is only ~1/(n) | the source of most of the ~375 known factors |
| Pollard rho | small factors | √p | F8 (Brent–Pollard 1980) |
| Pollard p−1 | p with smooth p−1 | — | Fermat numbers' factors have p−1 = k·2^m, so p−1 is smooth if k is smooth |
| ECM (GMP-ECM with Kruppa's Fermat code) | factors up to ~60 digits | needs arithmetic mod F_n | F10, F11 (Brent), F12 P54 (Vang 2010), F14 P54 (2010) |
| NFS (special) | complete factorization | sub-exponential in size of F_n | F9 (Lenstra, Lenstra, Manasse, Pollard 1990) |
| Pépin test | proves primality or compositeness | 2^n − 1 squarings mod F_n: Θ(2^n · M(2^n)) bit ops | F14, F20, F22, F24 |
| Suyama test | cofactor compositeness | one extra exponentiation mod C | F12–F30 cofactors |

The core obstruction: **no known primality test for F_n is cheaper than about 2^n
multiplications of 2^n-bit numbers.** Pépin's test is already optimal among Lucas-type
tests, since F_n − 1 = 2^{2^n} is fully factored. General-purpose methods (ECPP, APR-CL, AKS) are
asymptotically far worse at this size. For multiplication modulo 2^N + 1, the fastest
implementations use the weighted (negacyclic) FFT/DWT of Crandall–Fagin, which avoids a separate reduction step.

## 4. Heuristics on the number of Fermat primes

* **Hardy–Wright style.** Treat F_n as a random odd integer of its size. Then
  P(prime) ≈ 2/ln F_n ≈ 2.885/2^n. Summing over n ≥ 5 gives ≈ 0.18 expected new primes. Restricting to
  n ≥ 33 gives ≈ 6.7·10^{-10}. *(Corrected from 0.09 after review.)*
* **Factor-shape-aware heuristic.** Here p | F_n is only possible for p ≡ 1 (mod 2^{n+2}), and
  Prob(p | F_n) ≈ 2^{n+1}/(p−1) for such p. A Mertens-type product then gives
  P(F_n prime | no factor ≤ B) ≈ ln B / ln √F_n = 2·log₂B / 2^n. The unconditional version
  with B ≈ 2^{n+2} gives P ≈ (e^γ/2)(n+2)/2^{n−1}. That is a factor (n+2)·e^γ·ln2/2 ≈ 0.62(n+2) larger than the
  naive 2.885/2^n, because Fermat numbers structurally have no prime factors below 2^{n+2}.
* **Boklan & Conway (2016/2017),** *"Expect at most one billionth of a new Fermat prime!"*, Math.
  Intelligencer 39(1), arXiv:1605.01371 [S]. The search snippets indicate they discuss the special form of factors
  (Euler–Lucas) as a reason the naive Hardy–Wright argument needs correcting. How exactly they use it could not be checked. Their headline estimate is that, given the known compositeness
  data, the probability of any further Fermat prime is < 10^{-9}. We could not read the
  full text, so in `results/heuristic_probability.md` we re-derive a version from first
  principles using current search limits.
* **Historical claims.** Fermat (1640) conjectured that all F_n are prime; Euler refuted this in 1732 (641 | F5).
  Eisenstein (1844) proposed proving there are infinitely many Fermat primes. Several authors have published
  "proofs" that only finitely many Fermat numbers are prime, or that all F_n with n ≥ 5 are
  composite. None is accepted. Both questions — infinitely many Fermat primes? infinitely many composite Fermat numbers? — are
  open. Even the second, which is "obviously" true heuristically, is unproved.

## 5. Machine learning and AI-assisted mathematical discovery [K]/[S]

* Davies et al., *Advancing mathematics by guiding human intuition with AI*, Nature 600
  (2021): supervised learning plus attribution to find relationships between knot invariants and
  in representation theory. It needs large, unbiased datasets of objects.
* Raayoni et al., *The Ramanujan Machine* (Nature 2021): conjecture generation for continued fractions.
* FunSearch (Romera-Paredes et al., Nature 2024) and AlphaEvolve (2025) [S]: LLM-guided program
  evolution with an exact evaluator. They succeed where a cheap verifier scores candidates
  (cap sets, packing, Kakeya-type bounds).
* Attempts to learn primality or factorization from integer representations consistently
  reduce to small-modulus residues. This is expected: primality of a random integer is
  pseudo-random with respect to all low-complexity features, apart from divisibility by small primes.

**Implication.** For Fermat primality itself there are 5 positive and 28 negative examples,
and the negatives are determined by deterministic factor searches. That is far too little
data for a credible predictor. ML is appropriate only for the ~375 known *factors* and,
better, for an **unbiased exhaustive sample** of factors in a fixed region. We generate such a sample ourselves.

## 6. Open questions and research opportunities identified

1. Is there any property P(n), computable much faster than a Pépin test, with
   "F_n prime ⇒ P(n)" and P(33) false? (Track A.) Known necessary conditions all concern
   *divisors*, not n itself.
2. Do the known factors deviate from the "random order of 2" model? (Tracks B and E.) Examples:
   the distribution of the excess e = m − n − 2, of k mod small primes, and the counts per region
   compared with the ln ln prediction. A detectable deviation would be a genuine structural signal.
3. Re-derive the conditional probability of a new Fermat prime from current search limits.
4. Engineering (Track C): Pépin/Suyama arithmetic directly modulo 2^N + 1 with a negacyclic
   transform. GMP contains such a routine internally (`mpn_mul_fft`); is it usable, and how
   much faster is it than square-then-reduce?
5. Independent completeness checks of the factor tables in bounded regions.
