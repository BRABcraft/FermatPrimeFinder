# How likely is another Fermat prime? (heuristic, re-derived)

Source: `scripts/heuristic_probability.py` → `results/heuristic_probability.json`.
**This is a heuristic probability under an explicit model. It is not a theorem.**

## Model
For a prime p = k·2^m + 1 (k odd), treat 2 as a uniformly random quadratic residue. Then
P(p | F_n) = 2^{n+1−m}/k for n ≤ m−2. Chebotarev/Kummer arguments prove the *2-power part* of this (the excess
law, Theorem D in research/hypotheses.md) as a density statement. The 1/k factor (2 is a k-th power residue for the
specific prime p = k·2^m+1) is a heuristic for thin families and is only *tested empirically* (H3). The prime factors of F_n are taken to arise independently, which is the standard heuristic step.
Then

  P(F_n prime | no factor in a searched region R) = (e^γ/2) · exp(−Λ_unsearched(n)),

where Λ_unsearched(n) is the expected number of prime factors of F_n in [2^{n+2}, √F_n] outside R. It is computed with
exact primality for k < 4096 (where the 1/k weight concentrates) and the 2/ln p density beyond.
The factor e^γ/2 converts "no prime factor ≤ √N" into "prime" (Mertens).

## Calibration on known data
* Unconditional model probabilities for n = 5…32: P(F5) = 0.37, P(F6) = 0.17, P(F7) = 0.17, P(F8) = 0.085, …
  The expected number of primes among F5…F32 is **0.89**, and P(none) = **0.36**. Observing none is unremarkable.
* The same model predicts the *factor* data well (see `results/analysis/`). It predicts the number of Fermat-divisor primes
  below 2^B for B = 10…54 (e.g. 31.8 expected vs 36 observed below 2^54), the 1/k law, and the
  geometric excess law. So the model is not merely assumed; it is tested on the part of the data that can test it.

## Result (n = 33…160, excluding n with a factor verified in this repository)

| Conditioning | P(F33 prime) | P(F34) | P(F35) | P(at least one new Fermat prime, 33 ≤ n ≤ 160) |
|---|---|---|---|---|
| our certified regions (p < 2^54; k < 2^16 for m ≤ 1000; k < 2^12 for m ≤ 3000) | 1.12e−8 | 5.6e−9 | 2.8e−9 | **2.1e−8** |
| hypothetical: all k < 2^40 cleared | 1.58e−8 | 8.0e−9 | 4.0e−9 | 2.9e−8 |
| hypothetical: all k < 2^50 cleared | 1.78e−8 | 9.0e−9 | 4.6e−9 | 3.3e−8 |
| hypothetical: all k < 2^60 cleared | 1.99e−8 | 1.0e−8 | 5.1e−9 | 3.7e−8 |
| naive "random odd number" (2/ln F_n) | — | — | — | 6.1e−10 |

Notes:
* Deeper clean searches *raise* the conditional probability, but only logarithmically: P ∝ log(search bound).
  Even an absurd bound (k < 2^60 for every n) leaves the total below 4·10^{−8}.
* n > 160 contributes < 10^{−40}. Terms fall by a factor of 2 per index.
* n with a factor known in the literature but not verified here are *kept* as candidates (e.g. F28-type cases
  above n = 33 that we could not recall). This biases the total slightly *upward*, i.e. conservatively.
* **Comparison with Boklan & Conway (2016), "< 10^{−9}".** Our estimate is ≈ 20–40× larger than their headline figure.
  Relative to the *naive* random-odd-number estimate, our model gains a factor ≈ (n+2)·e^γ·ln2/2 ≈ 22 at n = 33,
  because F_n has no prime factors below 2^{n+2}. Snippets suggest Boklan–Conway also discuss this point, so we do
  **not** know the source of the discrepancy. It may lie in different conditioning, different bounds, or a different
  treatment of n ≥ 33. We could not read their text in this session, so no error is claimed on either side.
  Both estimates make a new Fermat prime heuristically negligible, and neither is a proof.

## Interpretation for the mission
A heuristic probability of ~2·10^{−8} that *any* further Fermat prime exists is the single most important
quantitative fact for this project. Together with the Pépin cost for F33 (see `results/benchmarks.md`), it means
a brute-force or "AI-guided" search for a new Fermat prime has an expected yield indistinguishable from zero.
A positive result would require a genuinely new theorem, and no evidence of a route to one was found (Track A).
