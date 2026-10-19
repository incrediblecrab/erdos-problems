/*
 * brute_count: histogram of doubled triangle areas in the a x b grid by enumerating every triple, O(n^3).
 *
 *   usage: brute_count A B [THREADS]
 *
 * Output is byte-identical in format to `grid_count A B 0`, so `cmp` compares the two.  Shares no code or method with
 * grid_count.c, which sums over triangle shapes instead of visiting triangles.
 */
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>

typedef unsigned long long u64;

static long NPTS, DMAX;
static long *X, *Y;
static long next_i = 0;
static pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;

static void *work(void *arg) {
    u64 *h = (u64 *)arg;
    for (;;) {
        long i;
        pthread_mutex_lock(&lock);
        i = next_i++;
        pthread_mutex_unlock(&lock);
        if (i >= NPTS - 2) break;
        for (long j = i + 1; j < NPTS; j++) {
            long ux = X[j] - X[i], uy = Y[j] - Y[i];
            for (long k = j + 1; k < NPTS; k++) {
                long d = ux * (Y[k] - Y[i]) - uy * (X[k] - X[i]);
                h[d < 0 ? -d : d]++;
            }
        }
    }
    return NULL;
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s A B [THREADS]\n", argv[0]); return 2; }
    long A = atol(argv[1]), B = atol(argv[2]);
    if (A < 1 || B < 1) { fprintf(stderr, "A, B must be >= 1\n"); return 2; }
    int T = argc > 3 ? atoi(argv[3]) : 8;
    if (T < 1) T = 1;
    NPTS = A * B;
    DMAX = (A - 1) * (B - 1);
    X = malloc((size_t)NPTS * sizeof *X);
    Y = malloc((size_t)NPTS * sizeof *Y);
    if (!X || !Y) { fprintf(stderr, "out of memory\n"); return 1; }
    for (long p = 0; p < NPTS; p++) { X[p] = p / B; Y[p] = p % B; }

    u64 **hs = calloc((size_t)T, sizeof(u64 *));
    pthread_t *th = calloc((size_t)T, sizeof(pthread_t));
    for (int t = 0; t < T; t++) {
        hs[t] = calloc((size_t)DMAX + 1, sizeof(u64));
        if (!hs[t]) { fprintf(stderr, "out of memory\n"); return 1; }
        pthread_create(&th[t], NULL, work, hs[t]);
    }
    for (int t = 0; t < T; t++) pthread_join(th[t], NULL);

    printf("# a=%ld b=%ld n=%ld dmax=%ld\n", A, B, NPTS, DMAX);
    for (long d = 1; d <= DMAX; d++) {
        u64 s = 0;
        for (int t = 0; t < T; t++) s += hs[t][d];
        if (s) printf("%ld %llu\n", d, s);
    }
    return 0;
}
