/*
 * verify_zpx.c -- exact, independent value of Phi_k for a block word, by
 * brute-force counting and Ehrhart interpolation.
 *
 * This program shares no code and no method with src/exact.py or
 * src/circle.py.  It knows nothing about carry regions, polygons or areas.
 * Given a binary word c of length m and k, it colours Z_N, N = m*t, by
 *
 *     colour(x) = c[ floor(x / t) ],       x = 0,...,N-1,
 *
 * (the equal-block colouring with blocks of size t) and counts the ordered
 * pairs (x,d) in Z_N^2, d = 0 included, with x, x+d, ..., x+(k-1)d all one
 * colour.  Call this count(t).
 *
 * Why this is exact.  For fixed m and word, the pairs counted at block size t
 * are the lattice points of t*P, where P is a finite disjoint union of
 * half-open polygons in [0,m)^2 cut out by the lines x + j d = integer
 * (j < k).  Each vertex of P is the meet of two such lines, so its
 * coordinates have denominators dividing some j - j' < k.  By Ehrhart's
 * theorem count(t) is therefore a quasi-polynomial of degree 2 in t whose
 * period divides L = lcm(1,...,k-1), with a constant leading coefficient
 * equal to the area of P.  Hence, for any t0 >= 1,
 *
 *     area(P) = [count(t0+2L) - 2 count(t0+L) + count(t0)] / (2 L^2),
 *     Phi_k   = area(P) / m^2   exactly,
 *
 * which is the limiting density of monochromatic pairs (x,d) of the block
 * colouring x -> c[floor(m x/p)] of F_p as p -> infinity, the quantity whose
 * half bounds delta~_k.  The program checks the quasi-polynomial structure
 * instead of assuming it: for every t0 in 1..L the second difference must be
 * the same, and the third difference count(t0+3L) - 3 count(t0+2L)
 * + 3 count(t0+L) - count(t0) must vanish.  Any failure exits 1.
 *
 * Build:  cc -O3 -march=native -pthread -o src/verify_zpx src/verify_zpx.c
 * Usage:  ./src/verify_zpx <word> <k> [nthreads] [residues]
 *         residues (default L) limits the t0 checked to 1..residues; the
 *         value is exact with one residue, the others test the structure.
 * Output: one JSON line with the exact Phi as "num/den".
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>

static int K, M, T;
static long long N;
static unsigned char *COL;

typedef struct { long long d0, d1, count; } job;

static void *worker(void *vp) {
    job *J = (job *)vp;
    long long acc = 0;
    for (long long d = J->d0; d < J->d1; d++) {
        for (long long x = 0; x < N; x++) {
            unsigned char c0 = COL[x];
            long long y = x;
            int ok = 1;
            for (int j = 1; j < K; j++) {
                y += d;
                if (y >= N) y -= N;
                if (COL[y] != c0) { ok = 0; break; }
            }
            acc += ok;
        }
    }
    J->count = acc;
    return NULL;
}

static long long count_at(const char *word, int t, int nt) {
    T = t;
    N = (long long)M * t;
    COL = malloc(N);
    for (long long x = 0; x < N; x++) COL[x] = (unsigned char)(word[x / t] - '0');
    pthread_t *th = malloc(sizeof(pthread_t) * nt);
    job *jobs = malloc(sizeof(job) * nt);
    for (int i = 0; i < nt; i++) {
        jobs[i].d0 = N * i / nt;
        jobs[i].d1 = N * (i + 1) / nt;
        jobs[i].count = 0;
        pthread_create(&th[i], NULL, worker, &jobs[i]);
    }
    long long tot = 0;
    for (int i = 0; i < nt; i++) { pthread_join(th[i], NULL); tot += jobs[i].count; }
    free(COL); free(th); free(jobs);
    return tot;
}

static long long gcdll(long long a, long long b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b) { long long r = a % b; a = b; b = r; }
    return a;
}

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "usage: %s <word> <k> [nthreads] [residues]\n", argv[0]);
        return 2;
    }
    const char *word = argv[1];
    M = (int)strlen(word);
    K = atoi(argv[2]);
    int nt = (argc > 3) ? atoi(argv[3]) : 8;
    for (int i = 0; i < M; i++)
        if (word[i] != '0' && word[i] != '1') { fprintf(stderr, "word must be 0/1\n"); return 2; }
    if (K < 2 || M < 1) { fprintf(stderr, "bad k or word\n"); return 2; }

    long long L = 1;
    for (int j = 2; j < K; j++) L = L / gcdll(L, j) * j;

    long long R = (argc > 4) ? atoll(argv[4]) : L;
    if (R < 1 || R > L) R = L;
    long long d2_ref = 0;
    int bad = 0;
    for (long long t0 = 1; t0 <= R; t0++) {
        long long c[4];
        for (int s = 0; s < 4; s++) c[s] = count_at(word, (int)(t0 + s * L), nt);
        long long d2 = c[2] - 2 * c[1] + c[0];
        long long d3 = c[3] - 3 * c[2] + 3 * c[1] - c[0];
        if (t0 == 1) d2_ref = d2;
        if (d2 != d2_ref || d3 != 0) {
            fprintf(stderr, "quasi-polynomial check failed at t0=%lld: d2=%lld (ref %lld) d3=%lld\n",
                    t0, d2, d2_ref, d3);
            bad = 1;
        }
    }
    /* Phi = d2 / (2 L^2 m^2) */
    long long num = d2_ref, den = 2 * L * L * (long long)M * (long long)M;
    long long g = gcdll(num, den);
    num /= g; den /= g;
    printf("{\"m\": %d, \"k\": %d, \"period\": %lld, \"residues_checked\": %lld, "
           "\"phi\": \"%lld/%lld\", \"phi_float\": %.12f, \"structure_ok\": %s}\n",
           M, K, L, R, num, den, (double)num / (double)den, bad ? "false" : "true");
    return bad ? 1 : 0;
}
