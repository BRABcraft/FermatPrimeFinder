# Next steps (prioritised by expected information gain per unit cost)

1. **Restore literature access and reconcile with Keller's table** (cheap, high value for data quality).
   Allow `prothsearch.com`, `fermatsearch.org`, `oeis.org`, `arxiv.org` in the environment's network settings. Then:
   import the full table of ~375 factors, run `scripts/validate_dataset.py` on all of it, and diff it against our exhaustive
   regions. Any factor in our regions that is missing from the table, or vice versa, would be a genuine discrepancy.
   Read Boklan–Conway's paper to reconcile their < 10^{−9} with our ≈ 2·10^{−8}.
2. **Grow the unbiased factor sample** to tighten the model tests (H1–H3) and to settle H-k3 with more power:
   * R1 extension to p < 2^58 costs ≈ 16× the 2^54 run (≈ 15 h on 4 cores) and adds ≈ 3 expected factors. That is poor value.
   * Rectangles k < 2^24 for m ≤ 1000 cost ≈ 16× R2c and add ≈ 15–20 factors. This is better value, but most of these factors are already
     in Keller's table (the searches are classical). Their value here is *exhaustiveness*, not novelty.
3. **Pépin residues for F20–F22 published comparison.** Our F20 (and F21) residues should be compared with the
   published Selfridge–Hurwitz residues (Young–Buell 1988 for F20; Crandall et al. 1995 for F22). This needs the
   literature (step 1). A Pépin test of F22 with the FFT code would take ≈ 8 h on one core.
4. **Gerbicz-style error checking** for long Pépin runs (needed beyond F22). It is already standard in the literature;
   implement it only if long runs are attempted.
5. **New-factor hunting** (Track B) has an expected yield of ≈ (bits of k searched) × Σ 1/m new factors *outside*
   ranges already covered by FermatSearch/PrimeGrid. Without their coverage map (blocked) we cannot target gaps,
   and blind rediscovery is wasteful. With access, pick n-ranges whose k-limits lag (typically very large n) and run
   `fsearch_gmp` there. GPU tools (mmff, GFN-sieve) are ~100× faster and are the right tools for this.
6. **Theory (Track A).** The only genuinely new direction would be an argument about F_n *itself* rather than about
   its possible divisors (H6). No candidate idea survived examination in this session. The most honest next step
   is a literature search on known conditional results (e.g. under GRH, results on the number of prime factors of F_n
   and on the counting function of Fermat-divisor primes) rather than further speculation.
