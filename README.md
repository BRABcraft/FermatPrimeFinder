# FermatPrimeFinder — an AI-assisted research attempt on Fermat primes

Goal: find a Fermat prime F_n = 2^(2^n)+1 with n ≥ 5, or failing that, produce
verified mathematical and computational results about Fermat numbers.
**Result in one line:** no new Fermat prime was found. The smallest candidate, F33 (8.6·10^9 bits), is far beyond
what can be Pépin-tested with the resources here (see the cost estimate in the report).
See [`FINAL_REPORT.md`](FINAL_REPORT.md) for what *was* established.

## Layout
| Path | Content |
|---|---|
| `research/literature_review.md` | state of the art, with source provenance |
| `research/hypotheses.md` | lemmas with proofs; hypothesis records with status |
| `research/experiment_log.md` | chronological log of experiments |
| `research/next_steps.md` | prioritised follow-ups |
| `src/fermat/` | exact-arithmetic library (`core.py`), heuristic model (`heuristics.py`, `model_mc.py`), region bookkeeping |
| `src/c/fsearch64.c` | exhaustive search for all prime factors p < 2^B (B ≤ 62) of F_n, n ≥ 12 |
| `src/c/fsearch_gmp.c` | exhaustive search over rectangles k < K, m ≤ M |
| `src/c/pepin_fft.c` | Pépin test using GMP's negacyclic FFT (`mpn_mul_fft`), working directly mod F_n |
| `scripts/` | dataset build/validate, Pépin+Suyama runner, analysis, heuristic probability, benchmarks |
| `data/` | `fermat_numbers.csv`, `fermat_factors.csv` (every row verified), unverified claims kept separate |
| `results/` | search outputs + certified regions, Pépin residues, analysis, benchmarks |
| `tests/` | exact-arithmetic regression tests |

## Reproduce
```bash
pip install gmpy2 sympy numpy scipy scikit-learn mpmath pytest   # plus PARI/GP (gp) and libgmp-dev
make            # builds the C tools
make test
./src/c/fsearch64 54 12 > results/search/fsearch64_B54.txt      # ~1 h on 4 cores
./src/c/fsearch_gmp 1 65535 14 1000 > results/search/gmp_k16_m1000.txt
python3 scripts/build_dataset.py && python3 scripts/validate_dataset.py
for n in $(seq 1 20); do python3 scripts/run_pepin.py $n; done
python3 scripts/run_analysis.py
python3 scripts/heuristic_probability.py
```
