/*
 * verify_zp.c -- independent check of the block-colouring bounds for Erdos #1186.
 *
 * This program shares no code and no method with src/exact.py.  It does not know
 * about carry regions, cell areas, or weight tables.  It takes a prime p and a
 * binary word c of length m, builds the colouring of Z_p directly as
 *
 *     colour(x) = c[ floor(m*x/p) ],      x = 0,...,p-1,
 *
 * and counts monochromatic k-term APs in Z_p by brute force over all p^2 pairs
 * (x,d).  It then reports the measured density, which must converge to the
 * exactly computed Phi_k(c)/2 at rate O(1/p).
 *
 * Counting convention, matching the problem statement: a k-AP is an unordered
 * non-degenerate progression.  The pairs (x,d) and (x+(k-1)d, -d) give the same
 * AP, and d=0 is excluded, so
 *
 *     #mono k-APs = (ordered monochromatic (x,d) pairs - p) / 2.
 *
 * Build:  cc -O3 -march=native -pthread -o src/verify_zp src/verify_zp.c
 * Usage:  ./src/verify_zp <p> <word> <k> [nthreads]
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <pthread.h>

static long long P;
static int K, M;
static unsigned char *COL;

typedef struct { long long d0, d1; long long count; } job;

static void *worker(void *vp) {
    job *J = (job *)vp;
    long long acc = 0;
    const long long p = P;
    const int k = K;
    for (long long d = J->d0; d < J->d1; d++) {
        for (long long x = 0; x < p; x++) {
            unsigned char c0 = COL[x];
            long long y = x;
            int ok = 1;
            for (int j = 1; j < k; j++) {
                y += d;
                if (y >= p) y -= p;
                if (COL[y] != c0) { ok = 0; break; }
            }
            acc += ok;
        }
    }
    J->count = acc;
    return NULL;
}

static int is_prime(long long n) {
    if (n < 2) return 0;
    for (long long f = 2; f * f <= n; f++) if (n % f == 0) return 0;
    return 1;
}

int main(int argc, char **argv) {
    if (argc < 4) {
        fprintf(stderr, "usage: %s <p> <word> <k> [nthreads]\n", argv[0]);
        return 2;
    }
    P = atoll(argv[1]);
    const char *word = argv[2];
    M = (int)strlen(word);
    K = atoi(argv[3]);
    int nt = (argc > 4) ? atoi(argv[4]) : 8;
    if (!is_prime(P)) { fprintf(stderr, "p must be prime\n"); return 2; }

    COL = malloc(P);
    for (long long x = 0; x < P; x++) {
        long long blk = ((long long)M * x) / P;          /* floor(m*x/p) */
        COL[x] = (unsigned char)(word[blk] - '0');
    }

    pthread_t *th = malloc(sizeof(pthread_t) * nt);
    job *jobs = malloc(sizeof(job) * nt);
    for (int i = 0; i < nt; i++) {
        jobs[i].d0 = P * i / nt;
        jobs[i].d1 = P * (i + 1) / nt;
        jobs[i].count = 0;
        pthread_create(&th[i], NULL, worker, &jobs[i]);
    }
    long long ordered = 0;
    for (int i = 0; i < nt; i++) { pthread_join(th[i], NULL); ordered += jobs[i].count; }

    long long nondeg = (ordered - P) / 2;               /* d=0 contributes exactly p */
    double pp = (double)P * (double)P;
    printf("{\"p\": %lld, \"m\": %d, \"k\": %d, \"ordered_mono\": %lld, "
           "\"mono_aps\": %lld, \"phi_measured\": %.12f, \"delta_measured\": %.12f}\n",
           P, M, K, ordered, nondeg, (double)ordered / pp, (double)nondeg / pp);
    free(COL); free(th); free(jobs);
    return 0;
}
