/*
 * pepin_fft: Pépin test for F_n = 2^(2^n)+1 using GMP's internal negacyclic
 * Schönhage–Strassen multiplication mpn_mul_fft, which returns products
 * directly modulo B^pl + 1 (B = 2^64).  With pl = 2^(n-6) limbs, B^pl + 1 = F_n,
 * so each Pépin step is ONE half-size transform and no separate reduction.
 *
 * Usage: pepin_fft n [outfile]
 *   prints: n, verdict, Selfridge–Hurwitz residues (R mod 2^36, 2^36-1, 2^35-1),
 *   and a 64-bit residue; optionally writes the full residue 3^((F_n-1)/2) mod F_n
 *   in GMP raw format (mpz_out_raw) to outfile for Suyama tests.
 * Requires n >= 8 (mpn_mul_fft needs pl >= 4 limbs).
 */
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* GMP internals (exported in libgmp.so, declared in gmp-impl.h) */
extern mp_limb_t __gmpn_mul_fft(mp_ptr op, mp_size_t pl, mp_srcptr n, mp_size_t nl,
                                mp_srcptr m, mp_size_t ml, int k);
extern int __gmpn_fft_best_k(mp_size_t n, int sqr);

static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s n [outfile]\n", argv[0]); return 1; }
    int n = atoi(argv[1]);
    if (n < 8 || n > 40) { fprintf(stderr, "need 8 <= n <= 40\n"); return 1; }
    mp_size_t pl = (mp_size_t)1 << (n - 6);
    int k = __gmpn_fft_best_k(pl, 1);
    while (k > 0 && (pl % ((mp_size_t)1 << k)) != 0) k--;
    mp_limb_t *x = calloc(pl + 1, sizeof(mp_limb_t));
    mp_limb_t *y = calloc(pl + 1, sizeof(mp_limb_t));
    x[0] = 3;
    int top = 0;            /* value = {x,pl} + top*B^pl, top in {0,1}; B^pl == -1 */
    unsigned long long iters = (1ULL << n) - 1;  /* (F_n - 1)/2 = 2^(2^n - 1) */
    double t0 = now();
    for (unsigned long long i = 0; i < iters; i++) {
        if (top) {          /* x == -1 => x^2 == 1 */
            memset(x, 0, pl * sizeof(mp_limb_t)); x[0] = 1; top = 0; continue;
        }
        mp_limb_t c = __gmpn_mul_fft(y, pl, x, pl, x, pl, k);
        mp_limb_t *t = x; x = y; y = t;
        top = (int)c;
        if ((i & 0xFFFF) == 0xFFFF) {
            double el = now() - t0;
            fprintf(stderr, "n=%d iter %llu/%llu  %.1fs  eta %.0fs\n", n, i + 1, iters, el,
                    el / (i + 1) * (iters - i - 1));
        }
    }
    double el = now() - t0;
    /* convert to mpz and reduce mod F_n */
    mpz_t R, F, t;
    mpz_inits(R, F, t, NULL);
    mpz_import(R, pl, -1, sizeof(mp_limb_t), 0, 0, x);
    if (top) { mpz_set_ui(t, 1); mpz_mul_2exp(t, t, (mp_bitcnt_t)1 << n); mpz_add(R, R, t); }
    mpz_set_ui(F, 1); mpz_mul_2exp(F, F, (mp_bitcnt_t)1 << n); mpz_add_ui(F, F, 1);
    mpz_mod(R, R, F);
    mpz_sub_ui(t, F, 1);
    int prime = mpz_cmp(R, t) == 0;
    mpz_t m1, m2, m3; mpz_inits(m1, m2, m3, NULL);
    mpz_set_ui(t, 1); mpz_mul_2exp(t, t, 36); mpz_mod(m1, R, t);
    mpz_sub_ui(t, t, 1); mpz_mod(m2, R, t);
    mpz_set_ui(t, 1); mpz_mul_2exp(t, t, 35); mpz_sub_ui(t, t, 1); mpz_mod(m3, R, t);
    mpz_set_ui(t, 1); mpz_mul_2exp(t, t, 64);
    mpz_t r64; mpz_init(r64); mpz_mod(r64, R, t);
    gmp_printf("n=%d verdict=%s SH_residues=(%Zd, %Zd, %Zd) res64=%016Zx k_fft=%d squarings=%llu time=%.3fs per_sq=%.3gms\n",
               n, prime ? "PRIME" : "COMPOSITE", m1, m2, m3, r64, k, iters, el, 1000.0 * el / iters);
    if (argc > 2) {
        FILE *f = fopen(argv[2], "wb"); mpz_out_raw(f, R); fclose(f);
    }
    return 0;
}
