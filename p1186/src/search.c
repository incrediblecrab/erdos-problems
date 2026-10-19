/*
 * search.c -- optimiser over block colourings for Erdos #1186.
 *
 * This program does no mathematics of its own.  It reads an exact integer
 * weight table produced by src/export_weights.py, where
 *
 *     Phi_k(c) = ( sum over sets U on which c is constant of w_U ) / DEN
 *
 * with all w_U and DEN integers, and searches binary words c in {0,1}^m for one
 * minimising that sum.  Every colouring it reports is re-evaluated from scratch
 * in exact rational arithmetic by src/certify.py, which shares no code with this
 * file.  Nothing printed here is a certificate.
 *
 * Sets are stored as explicit index lists, not bit masks, so m is not capped at
 * 64.  Exhaustive mode needs a Gray-code counter and is limited to m <= 62.
 *
 * Modes:
 *   exhaustive [NT]             Gray-code sweep of all 2^m words over NT threads.
 *                               Uses Phi(c) = Phi(~c) to halve the space.
 *   anneal R S SEED [NT]        restarted simulated annealing with single-bit
 *                               flips; every restart ends with a greedy 1- and
 *                               2-flip descent to a local minimum.
 *   seed WORD R S SEED [NT]     as anneal, but thread 0 restart 0 starts from
 *                               WORD.  Used to warm-start m' = r*m from the best
 *                               word at m, whose blow-up has exactly the same Phi.
 *
 * Build:  cc -O3 -march=native -pthread -o src/search src/search.c -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <pthread.h>
#include <math.h>

static int M, K, NSETS;
static long long DEN;
static int *SOFF, *SIDX, *SSIZE;   /* set i occupies SIDX[SOFF[i] .. +SSIZE[i]) */
static long long *SW;
static int **BITSETS, *BITCNT;     /* BITSETS[b] lists the sets containing block b */

static pthread_mutex_t LOCK = PTHREAD_MUTEX_INITIALIZER;
static long long BEST = 0x7fffffffffffffffLL;
static char *BESTC;

static void offer(long long val, const char *c) {
    pthread_mutex_lock(&LOCK);
    if (val < BEST) { BEST = val; memcpy(BESTC, c, M); }
    pthread_mutex_unlock(&LOCK);
}

static void load(const char *path) {
    FILE *f = fopen(path, "r");
    if (!f) { perror("open weights"); exit(1); }
    if (fscanf(f, "%d %d %d %lld", &M, &K, &NSETS, &DEN) != 4) {
        fprintf(stderr, "bad header\n"); exit(1);
    }
    SSIZE = malloc(sizeof(int) * NSETS);
    SOFF  = malloc(sizeof(int) * NSETS);
    SW    = malloc(sizeof(long long) * NSETS);
    SIDX  = malloc(sizeof(int) * (size_t)NSETS * (size_t)K);
    BITCNT = calloc(M, sizeof(int));
    int off = 0;
    for (int i = 0; i < NSETS; i++) {
        int sz; long long w;
        if (fscanf(f, "%d %lld", &sz, &w) != 2) { fprintf(stderr, "bad row %d\n", i); exit(1); }
        SSIZE[i] = sz; SW[i] = w; SOFF[i] = off;
        for (int t = 0; t < sz; t++) {
            int b;
            if (fscanf(f, "%d", &b) != 1) { fprintf(stderr, "bad index\n"); exit(1); }
            SIDX[off++] = b; BITCNT[b]++;
        }
    }
    fclose(f);
    BITSETS = malloc(sizeof(int *) * M);
    int *fill = calloc(M, sizeof(int));
    for (int b = 0; b < M; b++) BITSETS[b] = malloc(sizeof(int) * BITCNT[b]);
    for (int i = 0; i < NSETS; i++)
        for (int t = 0; t < SSIZE[i]; t++) {
            int b = SIDX[SOFF[i] + t];
            BITSETS[b][fill[b]++] = i;
        }
    free(fill);
    BESTC = calloc(M, 1);
}

static long long eval_full(const char *c, int *cnt1) {
    long long tot = 0;
    for (int i = 0; i < NSETS; i++) {
        int n1 = 0;
        for (int t = 0; t < SSIZE[i]; t++) n1 += c[SIDX[SOFF[i] + t]];
        cnt1[i] = n1;
        if (n1 == 0 || n1 == SSIZE[i]) tot += SW[i];
    }
    return tot;
}

/* flip block b of c (c is updated in place); returns the change in the objective */
static inline long long flip(int b, char *c, int *cnt1) {
    long long d = 0;
    int up = c[b] ? -1 : +1;
    const int *lst = BITSETS[b];
    const int n = BITCNT[b];
    for (int t = 0; t < n; t++) {
        int i = lst[t];
        int old = cnt1[i], sz = SSIZE[i], nw = old + up;
        cnt1[i] = nw;
        int was = (old == 0 || old == sz), is = (nw == 0 || nw == sz);
        if (was != is) d += is ? SW[i] : -SW[i];
    }
    c[b] ^= 1;
    return d;
}

/* greedy descent over single and pair flips until no improving move remains */
static long long polish(char *c, int *cnt1, long long cur) {
    int improved = 1;
    while (improved) {
        improved = 0;
        for (int b = 0; b < M; b++) {
            long long d = flip(b, c, cnt1);
            if (d < 0) { cur += d; improved = 1; } else flip(b, c, cnt1);
        }
        if (improved) continue;
        for (int b = 0; b < M && !improved; b++) {
            long long d1 = flip(b, c, cnt1);
            for (int e = b + 1; e < M; e++) {
                long long d2 = flip(e, c, cnt1);
                if (d1 + d2 < 0) { cur += d1 + d2; improved = 1; break; }
                flip(e, c, cnt1);
            }
            if (!improved) flip(b, c, cnt1);
        }
    }
    return cur;
}

/* ---------------- exhaustive ---------------- */
typedef struct { int hi_bits; uint64_t hi_val; } exh_arg;

static void *exh_worker(void *vp) {
    exh_arg *A = (exh_arg *)vp;
    int lo = M - A->hi_bits;
    char *c = calloc(M, 1);
    for (int b = 0; b < A->hi_bits; b++) c[lo + b] = (char)((A->hi_val >> b) & 1ULL);
    int *cnt1 = malloc(sizeof(int) * NSETS);
    long long cur = eval_full(c, cnt1);
    long long best = cur;
    char *bc = malloc(M);
    memcpy(bc, c, M);
    uint64_t total = 1ULL << lo;
    for (uint64_t g = 1; g < total; g++) {
        cur += flip(__builtin_ctzll(g), c, cnt1);
        if (cur < best) { best = cur; memcpy(bc, c, M); }
    }
    offer(best, bc);
    free(c); free(bc); free(cnt1);
    return NULL;
}

static void run_exhaustive(int nthreads) {
    if (M > 62) { fprintf(stderr, "exhaustive requires m <= 62\n"); exit(1); }
    int hb = 0;
    while ((1 << hb) < nthreads) hb++;
    if (hb > M - 4) hb = (M > 4) ? M - 4 : 0;
    int chunks = 1 << hb;
    pthread_t *th = malloc(sizeof(pthread_t) * chunks);
    exh_arg *args = malloc(sizeof(exh_arg) * chunks);
    int n = 0;
    for (int h = 0; h < chunks; h++) {
        if (hb > 0 && ((h >> (hb - 1)) & 1)) continue;   /* Phi(c) = Phi(~c) */
        args[n].hi_bits = hb; args[n].hi_val = (uint64_t)h;
        pthread_create(&th[n], NULL, exh_worker, &args[n]);
        n++;
    }
    for (int i = 0; i < n; i++) pthread_join(th[i], NULL);
    free(th); free(args);
}

/* ---------------- annealing ---------------- */
typedef struct { int restarts, sweeps; uint64_t seed; const char *start; } ann_arg;

static inline uint64_t xs(uint64_t *s) {
    uint64_t x = *s; x ^= x << 13; x ^= x >> 7; x ^= x << 17; return *s = x;
}

static void *ann_worker(void *vp) {
    ann_arg *A = (ann_arg *)vp;
    uint64_t rs = A->seed ? A->seed : 0x9E3779B97F4A7C15ULL;
    int *cnt1 = malloc(sizeof(int) * NSETS);
    char *c = malloc(M), *bc = malloc(M), *lb = malloc(M);
    long long best = 0x7fffffffffffffffLL;
    double scale = (double)DEN;
    memset(bc, 0, M);
    for (int r = 0; r < A->restarts; r++) {
        int warm = (r == 0 && A->start);
        if (warm) for (int b = 0; b < M; b++) c[b] = (char)(A->start[b] - '0');
        else      for (int b = 0; b < M; b++) c[b] = (char)(xs(&rs) & 1ULL);
        long long cur = eval_full(c, cnt1);
        long long lbest = cur;
        memcpy(lb, c, M);
        double T0 = warm ? 0.012 : 0.06, T1 = 0.0002;
        for (int s = 0; s < A->sweeps; s++) {
            double T = T0 * pow(T1 / T0, (double)s / (double)A->sweeps);
            for (int q = 0; q < M; q++) {
                int b = (int)(xs(&rs) % (uint64_t)M);
                long long d = flip(b, c, cnt1);
                int accept;
                if (d <= 0) accept = 1;
                else {
                    double u = (double)(xs(&rs) >> 11) / 9007199254740992.0;
                    accept = (u < exp(-((double)d / scale) / T));
                }
                if (accept) { cur += d; if (cur < lbest) { lbest = cur; memcpy(lb, c, M); } }
                else flip(b, c, cnt1);           /* undo; cnt1 is restored exactly */
            }
        }
        memcpy(c, lb, M);
        cur = eval_full(c, cnt1);
        cur = polish(c, cnt1, cur);
        if (cur < best) { best = cur; memcpy(bc, c, M); }
    }
    offer(best, bc);
    free(cnt1); free(c); free(bc); free(lb);
    return NULL;
}

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr,
            "usage: %s <weights> exhaustive [nthreads]\n"
            "       %s <weights> anneal <restarts> <sweeps> <seed> [nthreads]\n"
            "       %s <weights> seed <word> <restarts> <sweeps> <seed> [nthreads]\n",
            argv[0], argv[0], argv[0]);
        return 2;
    }
    load(argv[1]);
    if (!strcmp(argv[2], "exhaustive")) {
        run_exhaustive(argc > 3 ? atoi(argv[3]) : 8);
    } else {
        const char *start = NULL;
        int ai = 3;
        if (!strcmp(argv[2], "seed")) {
            start = argv[3]; ai = 4;
            if ((int)strlen(start) != M) { fprintf(stderr, "seed word length != m\n"); return 2; }
        } else if (strcmp(argv[2], "anneal")) {
            fprintf(stderr, "unknown mode\n"); return 2;
        }
        int restarts = atoi(argv[ai]), sweeps = atoi(argv[ai + 1]);
        uint64_t seed = strtoull(argv[ai + 2], NULL, 10);
        int nt = (argc > ai + 3) ? atoi(argv[ai + 3]) : 8;
        pthread_t *th = malloc(sizeof(pthread_t) * nt);
        ann_arg *args = malloc(sizeof(ann_arg) * nt);
        for (int i = 0; i < nt; i++) {
            args[i].restarts = restarts; args[i].sweeps = sweeps;
            args[i].seed = seed + 0x1234567ULL * (uint64_t)(i + 1);
            args[i].start = (i == 0) ? start : NULL;
            pthread_create(&th[i], NULL, ann_worker, &args[i]);
        }
        for (int i = 0; i < nt; i++) pthread_join(th[i], NULL);
        free(th); free(args);
    }

    int *cnt1 = malloc(sizeof(int) * NSETS);
    long long check = eval_full(BESTC, cnt1);
    free(cnt1);
    printf("{\"m\": %d, \"k\": %d, \"den\": %lld, \"num\": %lld, \"recheck\": %lld, \"word\": \"",
           M, K, DEN, BEST, check);
    for (int b = 0; b < M; b++) putchar(BESTC[b] ? '1' : '0');
    printf("\", \"phi\": %.12f}\n", (double)BEST / (double)DEN);
    return (check == BEST) ? 0 : 1;
}
