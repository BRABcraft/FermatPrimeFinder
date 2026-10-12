/*
 * fsearch64: exhaustive search for prime factors p < 2^B (B <= 62) of Fermat
 * numbers F_n with n >= NMIN (default 12).
 *
 * Completeness argument:
 *   p | F_n (n >= 2)  =>  p = k*2^m + 1, k odd, m >= n+2.
 *   For each m in [NMIN+2, B-1] and every odd k with k*2^m+1 < 2^B we test
 *   2^(2^j) == -1 (mod p) for every j in [NMIN, m-2], i.e. every admissible n.
 *   Candidates are skipped only if they have a prime divisor q <= QMAX with
 *   q < p (then p is composite, and a composite cannot be a *prime* factor).
 *   Hence every prime factor p < 2^B of every F_n, n >= NMIN, is reported.
 *
 * Arithmetic: 64-bit Montgomery multiplication with unsigned __int128.
 * Requires p < 2^63.
 *
 * Usage: fsearch64 B [NMIN] [QMAX]
 * Output lines: FACTOR p k m n
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef unsigned __int128 u128;
typedef uint64_t u64;

static inline u64 mont_inv(u64 p) { /* -p^{-1} mod 2^64 */
    u64 x = p;                       /* p*p == 1 mod 8 */
    for (int i = 0; i < 5; i++) x *= 2 - p * x;
    return (u64)0 - x;
}
static inline u64 mont_mul(u64 a, u64 b, u64 p, u64 pinv) {
    u128 t = (u128)a * b;
    u64 m = (u64)t * pinv;
    u128 u = (t + (u128)m * p) >> 64;
    u64 r = (u64)u;
    return r >= p ? r - p : r;
}

#define SEG (1u << 22)  /* odd k per segment */

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s B [NMIN] [QMAX]\n", argv[0]); return 1; }
    int B = atoi(argv[1]);
    int NMIN = argc > 2 ? atoi(argv[2]) : 12;
    int QMAX = argc > 3 ? atoi(argv[3]) : 16381;
    if (B > 62 || B < NMIN + 3) { fprintf(stderr, "bad B\n"); return 1; }

    /* small odd primes */
    char *comp = calloc(QMAX + 1, 1);
    int *qs = malloc(sizeof(int) * (QMAX + 1)); int nq = 0;
    for (int i = 2; i <= QMAX; i++) if (!comp[i]) {
        if (i > 2) qs[nq++] = i;
        for (long j = (long)i * i; j <= QMAX; j += i) comp[j] = 1;
    }

    u64 tested_total = 0, cand_total = 0;
    time_t t0 = time(NULL);

    for (int m = NMIN + 2; m <= B - 1; m++) {
        u64 kmax = ((((u64)1) << B) - 2) >> m;       /* k*2^m+1 < 2^B  <=> k <= (2^B-2)/2^m */
        u64 nodd = (kmax + 1) / 2;                    /* odd k = 2i+1, i < nodd */
        int jlo = NMIN, jhi = m - 2;
        u64 found_m = 0, tested_m = 0;
        #pragma omp parallel for schedule(dynamic, 1) reduction(+:tested_m)
        for (u64 seg0 = 0; seg0 < nodd; seg0 += SEG) {
            u64 len = nodd - seg0 < SEG ? nodd - seg0 : SEG;
            unsigned char *sv = malloc(len);
            memset(sv, 1, len);
            for (int t = 0; t < nq; t++) {
                u64 q = qs[t];
                /* want k*2^m+1 == 0 mod q, k = 2i+1  =>  i == (k0-1)/2 mod q,
                   k0 = -(2^m)^{-1} mod q */
                u64 tm = 1; for (int s = 0; s < m; s++) tm = (tm * 2) % q;
                /* inverse via Fermat little theorem */
                u64 inv = 1, b = tm, e = q - 2;
                while (e) { if (e & 1) inv = inv * b % q; b = b * b % q; e >>= 1; }
                u64 k0 = (q - inv) % q;
                u64 inv2 = (q + 1) / 2;
                u64 i0 = ((k0 + q - 1) % q) * inv2 % q;
                /* first index >= seg0 */
                u64 r = seg0 % q;
                u64 start = (i0 + q - r) % q;
                for (u64 i = start; i < len; i += q) {
                    u64 k = 2 * (seg0 + i) + 1;
                    u64 p = (k << m) + 1;
                    if (p != q) sv[i] = 0;
                }
            }
            for (u64 i = 0; i < len; i++) {
                if (!sv[i]) continue;
                u64 k = 2 * (seg0 + i) + 1;
                u64 p = (k << m) + 1;
                tested_m++;
                u64 pinv = mont_inv(p);
                u64 R = (u64)(((u128)1 << 64) % p);
                u64 x = (u64)(((u128)R * 2) % p);             /* mont(2)  */
                u64 target = p - R;                            /* mont(p-1) */
                for (int j = 0; j < jlo; j++) x = mont_mul(x, x, p, pinv); /* mont(2^(2^jlo)) */
                for (int j = jlo; j <= jhi; j++) {
                    if (x == target) {
                        #pragma omp critical
                        { printf("FACTOR %llu %llu %d %d\n", (unsigned long long)p,
                                 (unsigned long long)k, m, j); fflush(stdout); }
                        break;  /* order is unique */
                    }
                    x = mont_mul(x, x, p, pinv);
                }
            }
            free(sv);
        }
        tested_total += tested_m; cand_total += nodd;
        fprintf(stderr, "m=%d candidates=%llu tested=%llu elapsed=%lds\n", m,
                (unsigned long long)nodd, (unsigned long long)tested_m, (long)(time(NULL) - t0));
        (void)found_m;
    }
    fprintf(stderr, "DONE B=%d NMIN=%d QMAX=%d candidates=%llu tested=%llu elapsed=%lds\n",
            B, NMIN, QMAX, (unsigned long long)cand_total, (unsigned long long)tested_total,
            (long)(time(NULL) - t0));
    printf("DONE B=%d NMIN=%d QMAX=%d candidates=%llu tested=%llu\n", B, NMIN, QMAX,
           (unsigned long long)cand_total, (unsigned long long)tested_total);
    return 0;
}
