/*
 * verify_n.c -- independent check of the {1,...,n} block-colouring bounds (Erdos #1186).
 *
 * Shares no code and no method with src/exact.py: no carry regions, no cell
 * areas, no weight tables.  It takes n and a binary word c of length m, builds
 *
 *     colour(i) = c[ floor(m*i/n) ],      i = 0,...,n-1,
 *
 * and counts monochromatic k-term APs by brute force over every progression
 * (a, a+d, ..., a+(k-1)d) with d >= 1 and a + (k-1)d <= n-1.  That is the count
 * the problem normalises by n^2, so the measured density must converge to the
 * exactly computed Psi_k(c) at rate O(1/n).
 *
 * The total number of k-APs is also reported, as a check on the convention: it
 * must be ~ n^2/(2(k-1)).
 *
 * Build:  cc -O3 -march=native -pthread -o src/verify_n src/verify_n.c
 * Usage:  ./src/verify_n <n> <word> <k> [nthreads]
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <pthread.h>

static long long N;
static int K, M;
static unsigned char *COL;

typedef struct { long long d0, d1, mono, total; } job;

static void *worker(void *vp) {
    job *J = (job *)vp;
    long long mono = 0, total = 0;
    const int k = K;
    for (long long d = J->d0; d < J->d1; d++) {
        long long amax = N - 1 - (long long)(k - 1) * d;
        if (amax < 0) break;
        total += amax + 1;
        for (long long a = 0; a <= amax; a++) {
            unsigned char c0 = COL[a];
            int ok = 1;
            long long y = a;
            for (int j = 1; j < k; j++) {
                y += d;
                if (COL[y] != c0) { ok = 0; break; }
            }
            mono += ok;
        }
    }
    J->mono = mono; J->total = total;
    return NULL;
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: %s <n> <word> <k> [nthreads]\n", argv[0]); return 2; }
    N = atoll(argv[1]);
    const char *word = argv[2];
    M = (int)strlen(word);
    K = atoi(argv[3]);
    int nt = (argc > 4) ? atoi(argv[4]) : 8;

    COL = malloc(N);
    for (long long i = 0; i < N; i++) COL[i] = (unsigned char)(word[(M * i) / N] - '0');

    long long dmax = (N - 1) / (K - 1);          /* largest usable common difference */
    pthread_t *th = malloc(sizeof(pthread_t) * nt);
    job *jobs = malloc(sizeof(job) * nt);
    for (int i = 0; i < nt; i++) {
        /* split on d with equal work: the work in d is proportional to n-(k-1)d,
           so cut the triangle into equal-area strips */
        double f0 = (double)i / nt, f1 = (double)(i + 1) / nt;
        jobs[i].d0 = 1 + (long long)(dmax * (1.0 - __builtin_sqrt(1.0 - f0)));
        jobs[i].d1 = 1 + (long long)(dmax * (1.0 - __builtin_sqrt(1.0 - f1)));
        if (i == nt - 1) jobs[i].d1 = dmax + 1;
        jobs[i].mono = jobs[i].total = 0;
        pthread_create(&th[i], NULL, worker, &jobs[i]);
    }
    long long mono = 0, total = 0;
    for (int i = 0; i < nt; i++) { pthread_join(th[i], NULL); mono += jobs[i].mono; total += jobs[i].total; }

    double nn = (double)N * (double)N;
    printf("{\"n\": %lld, \"m\": %d, \"k\": %d, \"mono_aps\": %lld, \"total_aps\": %lld, "
           "\"psi_measured\": %.12f, \"total_density\": %.12f}\n",
           N, M, K, mono, total, (double)mono / nn, (double)total / nn);
    free(COL); free(th); free(jobs);
    return 0;
}
