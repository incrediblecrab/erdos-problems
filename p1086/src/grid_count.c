/*
 * grid_count: exact histogram of doubled triangle areas in the a x b integer grid {0..a-1} x {0..b-1}.
 *
 *   usage: grid_count A B [DMAX [THREADS]]      (DMAX <= 0 or omitted: every D up to (A-1)(B-1); THREADS default 8)
 *   output: a "# a=A b=B n=AB dmax=DMAX" header, then one "D count" line for each D = 1..DMAX with count > 0,
 *           where count is the number of unordered triangles {p,q,r} of the grid with |det(q-p, r-p)| = D.
 *
 * Method.  An ordered triple (p,q,r) is a placement of the shape {0,u,v}, u = q-p, v = r-p.  The shape fits
 * (a - Wx)(b - Wy) times, where Wx = max(0,ux,vx) - min(0,ux,vx) and Wy likewise, when Wx < a and Wy < b.
 * Summing over all (u,v) with det != 0 counts every triangle 6 times.  (u,v) and (-u,-v) contribute the
 * same, so u runs over the half-plane ux > 0 or (ux = 0, uy > 0) and the sum is doubled.
 */
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

typedef unsigned long long u64;

static long A, B, DMAX;
static long NU;            /* number of u vectors in the half-plane */
static long next_u = 0;    /* work counter, guarded by lock */
static pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;

static void u_of(long idx, long *ux, long *uy) {
    /* idx 0..B-2 -> (0, 1..B-1); then ux = 1..A-1 with uy = -(B-1)..B-1 */
    if (idx < B - 1) { *ux = 0; *uy = idx + 1; return; }
    idx -= B - 1;
    *ux = 1 + idx / (2 * B - 1);
    *uy = idx % (2 * B - 1) - (B - 1);
}

static void *work(void *arg) {
    u64 *h = (u64 *)arg;
    for (;;) {
        long idx;
        pthread_mutex_lock(&lock);
        idx = next_u++;
        pthread_mutex_unlock(&lock);
        if (idx >= NU) break;
        long ux, uy;
        u_of(idx, &ux, &uy);
        for (long vx = -(A - 1); vx <= A - 1; vx++) {
            long lo = 0, hi = 0;
            if (ux < lo) lo = ux; if (vx < lo) lo = vx;
            if (ux > hi) hi = ux; if (vx > hi) hi = vx;
            long wx = hi - lo;
            if (wx >= A) continue;
            u64 fx = (u64)(A - wx);
            for (long vy = -(B - 1); vy <= B - 1; vy++) {
                long d = ux * vy - uy * vx;
                if (d < 0) d = -d;
                if (d == 0 || d > DMAX) continue;
                long ylo = 0, yhi = 0;
                if (uy < ylo) ylo = uy; if (vy < ylo) ylo = vy;
                if (uy > yhi) yhi = uy; if (vy > yhi) yhi = vy;
                long wy = yhi - ylo;
                if (wy >= B) continue;
                h[d] += fx * (u64)(B - wy);
            }
        }
    }
    return NULL;
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s A B [DMAX [THREADS]]\n", argv[0]); return 2; }
    A = atol(argv[1]);
    B = atol(argv[2]);
    if (A < 1 || B < 1) { fprintf(stderr, "A, B must be >= 1\n"); return 2; }
    long full = (A - 1) * (B - 1);
    DMAX = argc > 3 ? atol(argv[3]) : full;
    if (DMAX <= 0 || DMAX > full) DMAX = full;
    int T = argc > 4 ? atoi(argv[4]) : 8;
    if (T < 1) T = 1;
    NU = (B - 1) + (A - 1) * (2 * B - 1);

    u64 **hs = calloc((size_t)T, sizeof(u64 *));
    pthread_t *th = calloc((size_t)T, sizeof(pthread_t));
    for (int t = 0; t < T; t++) {
        hs[t] = calloc((size_t)DMAX + 1, sizeof(u64));
        if (!hs[t]) { fprintf(stderr, "out of memory\n"); return 1; }
        pthread_create(&th[t], NULL, work, hs[t]);
    }
    for (int t = 0; t < T; t++) pthread_join(th[t], NULL);

    printf("# a=%ld b=%ld n=%ld dmax=%ld\n", A, B, A * B, DMAX);
    for (long d = 1; d <= DMAX; d++) {
        u64 s = 0;
        for (int t = 0; t < T; t++) s += hs[t][d];
        s *= 2;                      /* the (-u,-v) half */
        if (s % 6 != 0) { fprintf(stderr, "count for D=%ld not divisible by 6: %llu\n", d, s); return 1; }
        if (s) printf("%ld %llu\n", d, s / 6);
    }
    return 0;
}
