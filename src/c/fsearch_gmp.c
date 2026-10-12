/*
 * fsearch_gmp: exhaustive search over the rectangle
 *     k odd, KMIN <= k < KMAX,   MMIN <= m <= MMAX,   p = k*2^m + 1,
 * reporting every p that divides some F_n with n in [NMIN, m-2].
 *
 * Completeness: same argument as fsearch64.c.  A candidate is discarded only
 * if a prime q <= QMAX, q < p, divides it (then p is not prime).  Every
 * admissible n (NMIN <= n <= m-2) is checked by repeated squaring.
 *
 * Usage: fsearch_gmp KMIN KMAX MMIN MMAX [NMIN] [QMAX]
 * Output: FACTOR k m n
 */
#include <gmp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "usage\n"); return 1; }
    uint64_t KMIN = strtoull(argv[1], 0, 10), KMAX = strtoull(argv[2], 0, 10);
    int MMIN = atoi(argv[3]), MMAX = atoi(argv[4]);
    int NMIN = argc > 5 ? atoi(argv[5]) : 12;
    int QMAX = argc > 6 ? atoi(argv[6]) : 1000000;
    if (MMIN < NMIN + 2) MMIN = NMIN + 2;
    if (KMIN % 2 == 0) KMIN++;

    char *comp = calloc(QMAX + 1, 1);
    int *qs = malloc(sizeof(int) * (QMAX + 1)); int nq = 0;
    for (int i = 2; i <= QMAX; i++) if (!comp[i]) {
        if (i > 2) qs[nq++] = i;
        for (long j = (long)i * i; j <= QMAX; j += i) comp[j] = 1;
    }
    uint64_t nodd = (KMAX > KMIN) ? (KMAX - KMIN + 1) / 2 : 0;  /* k = KMIN + 2i */
    time_t t0 = time(NULL);
    uint64_t tested_total = 0;

    #pragma omp parallel for schedule(dynamic, 1) reduction(+:tested_total)
    for (int m = MMAX; m >= MMIN; m--) {   /* big m first for load balance */
        unsigned char *sv = malloc(nodd ? nodd : 1);
        memset(sv, 1, nodd);
        for (int t = 0; t < nq; t++) {
            uint64_t q = qs[t];
            uint64_t tm = 1, b = 2, e = m;   /* 2^m mod q */
            while (e) { if (e & 1) tm = tm * b % q; b = b * b % q; e >>= 1; }
            uint64_t inv = 1; b = tm; e = q - 2;
            while (e) { if (e & 1) inv = inv * b % q; b = b * b % q; e >>= 1; }
            uint64_t k0 = (q - inv) % q;            /* k == k0 mod q => q | p */
            /* k = KMIN + 2i == k0  =>  i == (k0 - KMIN) * inv2 mod q */
            uint64_t inv2 = (q + 1) / 2;
            uint64_t i0 = ((k0 + q - KMIN % q) % q) * inv2 % q;
            for (uint64_t i = i0; i < nodd; i += q) {
                /* p == q only possible if p small; p >= 2^m+1 > q when m >= 20 */
                if (m < 21) {
                    uint64_t k = KMIN + 2 * i;
                    if (k < (1ULL << 40) && ((k << m) + 1) == q) continue;
                }
                sv[i] = 0;
            }
        }
        mpz_t p, x, pm1;
        mpz_inits(p, x, pm1, NULL);
        uint64_t tested = 0;
        for (uint64_t i = 0; i < nodd; i++) {
            if (!sv[i]) continue;
            uint64_t k = KMIN + 2 * i;
            mpz_set_ui(p, k); mpz_mul_2exp(p, p, m); mpz_add_ui(p, p, 1);
            mpz_sub_ui(pm1, p, 1);
            mpz_set_ui(x, 2);
            tested++;
            for (int j = 0; j <= m - 2; j++) {
                if (j >= NMIN && mpz_cmp(x, pm1) == 0) {
                    #pragma omp critical
                    { printf("FACTOR %llu %d %d\n", (unsigned long long)k, m, j); fflush(stdout); }
                    break;
                }
                mpz_mul(x, x, x);
                mpz_tdiv_r(x, x, p);
            }
        }
        tested_total += tested;
        mpz_clears(p, x, pm1, NULL);
        free(sv);
    }
    printf("DONE KMIN=%llu KMAX=%llu MMIN=%d MMAX=%d NMIN=%d QMAX=%d tested=%llu elapsed=%lds\n",
           (unsigned long long)KMIN, (unsigned long long)KMAX, MMIN, MMAX, NMIN, QMAX,
           (unsigned long long)tested_total, (long)(time(NULL) - t0));
    return 0;
}
