// Lower bounds for g(n): configurations of n points of the k x k grid with many triangles of one doubled area.
//
//   small_search exhaustive N K [THREADS]        best largest equal-area class over all N-subsets of the K x K grid
//   small_search anneal N K D ITERS SEED         simulated annealing on the number of triangles of doubled area D
//
// Both print "best COUNT D  x,y x,y ..." for the best configuration found. The counts are rechecked independently by
// recount.py; nothing here is trusted on its own.
#include <math.h>
#include <pthread.h>
#include <stdatomic.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXP 4096
#define MAXN 64

static int N, K, NP;
static int PX[MAXP], PY[MAXP];

static inline int adet(int a, int b, int c) {
    int d = (PX[b] - PX[a]) * (PY[c] - PY[a]) - (PY[b] - PY[a]) * (PX[c] - PX[a]);
    return d < 0 ? -d : d;
}

// ---- exhaustive ----
static atomic_int next_task;
static int ntask, task_a[MAXP * 64], task_b[MAXP * 64];
static pthread_mutex_t mu = PTHREAD_MUTEX_INITIALIZER;
static int g_best = 0, g_d = 0, g_cfg[MAXN];

typedef struct { int sel[MAXN]; int hist[2 * MAXP]; int best, bd, cfg[MAXN]; } work_t;

// Ties are broken by smaller D, then by the lexicographically smaller point list, so the output is deterministic.
static int better(int c, int d, const int *sel, int best, int bd, const int *cfg) {
    if (c != best) return c > best;
    if (d != bd) return d < bd;
    for (int i = 0; i < N; i++)
        if (sel[i] != cfg[i]) return sel[i] < cfg[i];
    return 0;
}

static void dfs(work_t *w, int depth, int start) {
    if (depth == N) {
        int minx = K;
        for (int i = 0; i < N; i++) if (PX[w->sel[i]] < minx) minx = PX[w->sel[i]];
        if (minx != 0) return;  // translates: keep only sets touching column 0 (row 0 holds the first point)
        for (int d = 1; d <= (K - 1) * (K - 1); d++)
            if (w->hist[d] >= w->best && better(w->hist[d], d, w->sel, w->best, w->bd, w->cfg)) {
                w->best = w->hist[d];
                w->bd = d;
                memcpy(w->cfg, w->sel, sizeof(int) * N);
            }
        return;
    }
    for (int p = start; p <= NP - (N - depth); p++) {
        for (int i = 1; i < depth; i++)
            for (int j = 0; j < i; j++) w->hist[adet(w->sel[j], w->sel[i], p)]++;
        w->sel[depth] = p;
        dfs(w, depth + 1, p + 1);
        for (int i = 1; i < depth; i++)
            for (int j = 0; j < i; j++) w->hist[adet(w->sel[j], w->sel[i], p)]--;
    }
}

static void *exh_worker(void *arg) {
    (void)arg;
    work_t *w = calloc(1, sizeof(work_t));
    for (;;) {
        int t = atomic_fetch_add(&next_task, 1);
        if (t >= ntask) break;
        memset(w->hist, 0, sizeof(w->hist));
        w->sel[0] = task_a[t];
        w->sel[1] = task_b[t];
        dfs(w, 2, task_b[t] + 1);
    }
    pthread_mutex_lock(&mu);
    if (better(w->best, w->bd, w->cfg, g_best, g_d, g_cfg)) {
        g_best = w->best;
        g_d = w->bd;
        memcpy(g_cfg, w->cfg, sizeof(int) * N);
    }
    pthread_mutex_unlock(&mu);
    free(w);
    return NULL;
}

// ---- annealing ----
static unsigned long long rs;
static inline unsigned long long rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
static inline double urand(void) { return (rnd() >> 11) * (1.0 / 9007199254740992.0); }

static int contrib(const int *pos, int skip, int q, int D) {
    int c = 0;
    for (int i = 0; i < N; i++) {
        if (i == skip) continue;
        for (int j = i + 1; j < N; j++) {
            if (j == skip) continue;
            if (adet(pos[i], pos[j], q) == D) c++;
        }
    }
    return c;
}

static int count_all(const int *pos, int D) {
    int c = 0;
    for (int i = 0; i < N; i++)
        for (int j = i + 1; j < N; j++)
            for (int l = j + 1; l < N; l++) c += adet(pos[i], pos[j], pos[l]) == D;
    return c;
}

static void anneal(int D, long iters, int restarts) {
    int pos[MAXN], best_pos[MAXN], occ[MAXP];
    int best = -1;
    for (int r = 0; r < restarts; r++) {
        memset(occ, 0, sizeof(int) * NP);
        for (int i = 0; i < N; i++) {
            int c;
            do c = (int)(rnd() % NP); while (occ[c]);
            occ[c] = 1;
            pos[i] = c;
        }
        int cur = count_all(pos, D);
        if (cur > best) { best = cur; memcpy(best_pos, pos, sizeof(int) * N); }
        double T0 = 2.0, T1 = 0.05;
        for (long it = 0; it < iters; it++) {
            double T = T0 * pow(T1 / T0, (double)it / iters);
            int i = (int)(rnd() % N), c;
            if (urand() < 0.5) {
                c = (int)(rnd() % NP);
            } else {
                int x = PX[pos[i]] + (int)(rnd() % 5) - 2, y = PY[pos[i]] + (int)(rnd() % 5) - 2;
                if (x < 0 || y < 0 || x >= K || y >= K) continue;
                c = y * K + x;
            }
            if (occ[c]) continue;
            int delta = contrib(pos, i, c, D) - contrib(pos, i, pos[i], D);
            if (delta >= 0 || urand() < exp(delta / T)) {
                occ[pos[i]] = 0;
                occ[c] = 1;
                pos[i] = c;
                cur += delta;
                if (cur > best) { best = cur; memcpy(best_pos, pos, sizeof(int) * N); }
            }
        }
    }
    if (count_all(best_pos, D) != best) { fprintf(stderr, "internal count mismatch\n"); exit(2); }
    printf("best %d %d ", best, D);
    for (int i = 0; i < N; i++) printf(" %d,%d", PX[best_pos[i]], PY[best_pos[i]]);
    printf("\n");
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: see header\n"); return 1; }
    N = atoi(argv[2]);
    K = atoi(argv[3]);
    NP = K * K;
    if (N < 3 || N > MAXN || K < 2 || NP > MAXP || N > NP) { fprintf(stderr, "bad N or K\n"); return 1; }
    for (int p = 0; p < NP; p++) { PX[p] = p % K; PY[p] = p / K; }
    if (!strcmp(argv[1], "exhaustive")) {
        int threads = argc > 4 ? atoi(argv[4]) : 8;
        ntask = 0;
        for (int a = 0; a < K; a++)  // the smallest index lies in row 0
            for (int b = a + 1; b < NP; b++) { task_a[ntask] = a; task_b[ntask] = b; ntask++; }
        pthread_t th[64];
        for (int t = 0; t < threads; t++) pthread_create(&th[t], NULL, exh_worker, NULL);
        for (int t = 0; t < threads; t++) pthread_join(th[t], NULL);
        printf("best %d %d ", g_best, g_d);
        for (int i = 0; i < N; i++) printf(" %d,%d", PX[g_cfg[i]], PY[g_cfg[i]]);
        printf("\n");
        return 0;
    }
    if (!strcmp(argv[1], "anneal") && argc >= 7) {
        int D = atoi(argv[4]);
        long iters = atol(argv[5]);
        rs = 0x9E3779B97F4A7C15ULL ^ (unsigned long long)atoll(argv[6]) * 0xD1B54A32D192ED03ULL;
        if (!rs) rs = 1;
        anneal(D, iters, argc > 7 ? atoi(argv[7]) : 4);
        return 0;
    }
    fprintf(stderr, "usage: see header\n");
    return 1;
}
