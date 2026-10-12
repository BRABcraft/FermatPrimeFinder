/*
 * bench_sqr: time one modular squaring modulo F_n = 2^(2^n)+1 for
 *   (A) mpn_mul_fft directly mod B^pl+1   (negacyclic Schönhage–Strassen)
 *   (B) mpz_mul (full 2N-bit product) + shift/subtract reduction
 *   (C) mpz_mul + generic mpz_mod (division)
 * Usage: bench_sqr n reps
 * Prints: n, N bits, seconds per squaring for A, B, C, and checks A == B == C.
 */
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
extern mp_limb_t __gmpn_mul_fft(mp_ptr, mp_size_t, mp_srcptr, mp_size_t, mp_srcptr, mp_size_t, int);
extern int __gmpn_fft_best_k(mp_size_t, int);
static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }

int main(int argc, char **argv) {
    int n = atoi(argv[1]); int reps = atoi(argv[2]);
    mp_bitcnt_t N = (mp_bitcnt_t)1 << n;
    mp_size_t pl = (mp_size_t)1 << (n - 6);
    int k = __gmpn_fft_best_k(pl, 1);
    while (k > 0 && (pl % ((mp_size_t)1 << k))) k--;
    gmp_randstate_t rs; gmp_randinit_default(rs); gmp_randseed_ui(rs, 12345);
    mpz_t F, x0, xb, xc, t, lo; mpz_inits(F, x0, xb, xc, t, lo, NULL);
    mpz_set_ui(F, 1); mpz_mul_2exp(F, F, N); mpz_add_ui(F, F, 1);
    mpz_urandomb(x0, rs, N);
    /* A */
    mp_limb_t *xa = calloc(pl + 1, 8), *ya = calloc(pl + 1, 8);
    size_t cnt; mpz_export(xa, &cnt, -1, 8, 0, 0, x0);
    double t0 = now();
    for (int i = 0; i < reps; i++) {
        mp_limb_t c = __gmpn_mul_fft(ya, pl, xa, pl, xa, pl, k);
        if (c) { fprintf(stderr, "hit x == -1 (prob ~2^-N); rerun with another seed\n"); return 2; }
        mp_limb_t *s = xa; xa = ya; ya = s;
    }
    double tA = (now() - t0) / reps;
    /* B */
    mpz_set(xb, x0);
    t0 = now();
    for (int i = 0; i < reps; i++) {
        mpz_mul(xb, xb, xb);
        mpz_tdiv_r_2exp(lo, xb, N); mpz_tdiv_q_2exp(xb, xb, N);
        mpz_sub(xb, lo, xb); if (mpz_sgn(xb) < 0) mpz_add(xb, xb, F);
    }
    double tB = (now() - t0) / reps;
    /* C */
    mpz_set(xc, x0);
    t0 = now();
    for (int i = 0; i < reps; i++) { mpz_mul(xc, xc, xc); mpz_mod(xc, xc, F); }
    double tC = (now() - t0) / reps;
    mpz_import(t, pl, -1, 8, 0, 0, xa); mpz_mod(t, t, F);
    int agree = (mpz_cmp(t, xb) == 0) && (mpz_cmp(xb, xc) == 0);
    printf("n=%d bits=%lu fft_k=%d A_mulfft=%.4e B_mul+special=%.4e C_mul+mod=%.4e agree=%d\n",
           n, (unsigned long)N, k, tA, tB, tC, agree);
    return 0;
}
