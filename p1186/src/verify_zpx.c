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
 * Counting.  count(t) is computed with bitsets: for each d, the colour
 * sequence rotated by j*d is extracted 64 positions at a time from a doubled
 * copy, and the positions where all k rotations agree are popcounted.  The
 * reflection (x,d) -> (x+(k-1)d, -d) maps monochromatic pairs to
 * monochromatic pairs, so count(d) = count(N-d) and only d <= N/2 is scanned.
 * Neither shortcut is assumed: at t = 1 and t = 2 the count is also taken by
 * the plain double loop over every (x,d), d < N, and any disagreement exits 1.
 * (The plain loop was the whole program before November 2, 2026.)
 *
 * Build:  cc -O3 -pthread -o src/verify_zpx src/verify_zpx.c
 * Usage:  ./src/verify_zpx <word> <k> [nthreads] [residues]
 *         residues (default L) limits the t0 checked to 1..residues; the
 *         value is exact with one residue, the others test the structure.
 * Output: one JSON line with the exact Phi as "num/den".
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <pthread.h>

static int K, M, NT;
static long long N, W;              /* N = M*t positions, W = ceil(N/64) words */
static unsigned char *COL;          /* colour of each position (plain loop) */
static uint64_t *B2;                /* colour bits of positions 0..2N+127, mod N */
static uint64_t LASTMASK;

/* ---------------------------------------------------------------- plain loop */

typedef struct { long long d0, d1, count; int id; } job;

static void *plain_worker(void *vp) {
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

/* ------------------------------------------------------------------- bitsets */

static inline uint64_t rot_word(long long s, long long w) {
    long long bit = s + 64 * w;
    long long q = bit >> 6;
    int r = (int)(bit & 63);
    if (r == 0) return B2[q];
    return (B2[q] >> r) | (B2[q + 1] << (64 - r));
}

static long long count_d(long long d) {
    long long off[16];
    for (int j = 1; j < K; j++) off[j] = (long long)((__int128)j * d % N);
    long long acc = 0;
    for (long long w = 0; w < W; w++) {
        uint64_t a = B2[w];
        uint64_t one = a, zero = ~a;
        for (int j = 1; j < K; j++) {
            uint64_t v = rot_word(off[j], w);
            one &= v;
            zero &= ~v;
        }
        uint64_t m = (w == W - 1) ? LASTMASK : ~(uint64_t)0;
        acc += __builtin_popcountll(one & m) + __builtin_popcountll(zero & m);
    }
    return acc;
}

static void *bit_worker(void *vp) {
    job *J = (job *)vp;
    long long acc = 0;
    long long half = N / 2;           /* d = 1 .. ceil(N/2)-1 doubled, N/2 once if N even */
    for (long long d = 1 + J->id; d <= half; d += NT) {
        long long c = count_d(d);
        acc += (2 * d == N) ? c : 2 * c;
    }
    J->count = acc;
    return NULL;
}

static void setup(const char *word, int t) {
    N = (long long)M * t;
    COL = malloc(N);
    for (long long x = 0; x < N; x++) COL[x] = (unsigned char)(word[x / t] - '0');
    W = (N + 63) / 64;
    long long nw = 2 * W + 4;
    B2 = calloc(nw, sizeof(uint64_t));
    for (long long i = 0; i < 64 * nw; i++)
        if (COL[i % N]) B2[i >> 6] |= (uint64_t)1 << (i & 63);
    LASTMASK = (N % 64) ? (((uint64_t)1 << (N % 64)) - 1) : ~(uint64_t)0;
}

static long long run(void *(*fn)(void *), long long dmax) {
    pthread_t *th = malloc(sizeof(pthread_t) * NT);
    job *jobs = malloc(sizeof(job) * NT);
    for (int i = 0; i < NT; i++) {
        jobs[i].d0 = dmax * i / NT;
        jobs[i].d1 = dmax * (i + 1) / NT;
        jobs[i].id = i;
        jobs[i].count = 0;
        pthread_create(&th[i], NULL, fn, &jobs[i]);
    }
    long long tot = 0;
    for (int i = 0; i < NT; i++) { pthread_join(th[i], NULL); tot += jobs[i].count; }
    free(th); free(jobs);
    return tot;
}

static int plain_mismatch = 0;

static long long count_at(const char *word, int t) {
    setup(word, t);
    long long tot = count_d(0) + run(bit_worker, 0);
    if (t <= 2) {
        long long plain = run(plain_worker, N);
        if (plain != tot) {
            fprintf(stderr, "bitset count %lld != plain count %lld at t=%d\n", tot, plain, t);
            plain_mismatch = 1;
        }
    }
    free(COL); free(B2);
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
    NT = (argc > 3) ? atoi(argv[3]) : 8;
    if (NT < 1) NT = 1;
    for (int i = 0; i < M; i++)
        if (word[i] != '0' && word[i] != '1') { fprintf(stderr, "word must be 0/1\n"); return 2; }
    if (K < 2 || K > 16 || M < 1) { fprintf(stderr, "bad k or word\n"); return 2; }

    long long L = 1;
    for (int j = 2; j < K; j++) L = L / gcdll(L, j) * j;

    long long R = (argc > 4) ? atoll(argv[4]) : L;
    if (R < 1 || R > L) R = L;
    long long d2_ref = 0;
    int bad = 0;
    for (long long t0 = 1; t0 <= R; t0++) {
        long long c[4];
        for (int s = 0; s < 4; s++) c[s] = count_at(word, (int)(t0 + s * L));
        long long d2 = c[2] - 2 * c[1] + c[0];
        long long d3 = c[3] - 3 * c[2] + 3 * c[1] - c[0];
        if (t0 == 1) d2_ref = d2;
        if (d2 != d2_ref || d3 != 0) {
            fprintf(stderr, "quasi-polynomial check failed at t0=%lld: d2=%lld (ref %lld) d3=%lld\n",
                    t0, d2, d2_ref, d3);
            bad = 1;
        }
    }
    if (plain_mismatch) bad = 1;
    /* Phi = d2 / (2 L^2 m^2) */
    long long num = d2_ref, den = 2 * L * L * (long long)M * (long long)M;
    long long g = gcdll(num, den);
    num /= g; den /= g;
    printf("{\"m\": %d, \"k\": %d, \"period\": %lld, \"residues_checked\": %lld, "
           "\"phi\": \"%lld/%lld\", \"phi_float\": %.12f, \"plain_loop_agrees\": %s, "
           "\"structure_ok\": %s}\n",
           M, K, L, R, num, den, (double)num / (double)den,
           plain_mismatch ? "false" : "true", bad ? "false" : "true");
    return bad ? 1 : 0;
}
