// Exhaustive scan of the K-subsets of an A x B lattice box for large classes of one doubled area.
//
//   box_scan A B K T [THREADS]
//
// Prints one line "D COUNT x,y x,y ..." for every K-subset (up to translation inside the box: only subsets that
// meet column 0 and row 0) and every doubled area D whose class has at least T triangles, then "done SUBSETS".
// Written for src/g9_box_check.py; it shares no code with small_search.c, g7_exact.py or g9_exact.py.
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXK 16
#define MAXPTS 400

static int A, B, K, T, NPTS;
static int X[MAXPTS], Y[MAXPTS];
static int maxd;
static pthread_mutex_t out_mu = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t task_mu = PTHREAD_MUTEX_INITIALIZER;
static int next_first = 0;
static long long total_subsets = 0;

typedef struct {
    int chosen[MAXK];
    int *hist;
    long long subsets;
} state_t;

static int area2(int i, int j, int l) {
    int v = (X[j] - X[i]) * (Y[l] - Y[i]) - (Y[j] - Y[i]) * (X[l] - X[i]);
    return v < 0 ? -v : v;
}

static void leaf(state_t *s) {
    int has_col0 = 0, has_row0 = 0;
    for (int i = 0; i < K; i++) {
        if (X[s->chosen[i]] == 0) has_col0 = 1;
        if (Y[s->chosen[i]] == 0) has_row0 = 1;
    }
    if (!has_col0 || !has_row0) return;
    s->subsets++;
    for (int d = 1; d <= maxd; d++) {
        if (s->hist[d] < T) continue;
        pthread_mutex_lock(&out_mu);
        printf("%d %d", d, s->hist[d]);
        for (int i = 0; i < K; i++) printf(" %d,%d", X[s->chosen[i]], Y[s->chosen[i]]);
        printf("\n");
        pthread_mutex_unlock(&out_mu);
    }
}

static void grow(state_t *s, int size, int from) {
    if (size == K) {
        leaf(s);
        return;
    }
    for (int p = from; p <= NPTS - (K - size); p++) {
        for (int i = 0; i < size; i++)
            for (int j = i + 1; j < size; j++) s->hist[area2(s->chosen[i], s->chosen[j], p)]++;
        s->chosen[size] = p;
        grow(s, size + 1, p + 1);
        for (int i = 0; i < size; i++)
            for (int j = i + 1; j < size; j++) s->hist[area2(s->chosen[i], s->chosen[j], p)]--;
    }
}

static void *worker(void *unused) {
    (void)unused;
    state_t s;
    s.hist = calloc((size_t)maxd + 1, sizeof(int));
    s.subsets = 0;
    for (;;) {
        pthread_mutex_lock(&task_mu);
        int first = next_first++;
        pthread_mutex_unlock(&task_mu);
        if (first > NPTS - K) break;
        s.chosen[0] = first;
        grow(&s, 1, first + 1);
    }
    pthread_mutex_lock(&task_mu);
    total_subsets += s.subsets;
    pthread_mutex_unlock(&task_mu);
    free(s.hist);
    return NULL;
}

int main(int argc, char **argv) {
    if (argc < 5) {
        fprintf(stderr, "usage: box_scan A B K T [THREADS]\n");
        return 1;
    }
    A = atoi(argv[1]);
    B = atoi(argv[2]);
    K = atoi(argv[3]);
    T = atoi(argv[4]);
    int threads = argc > 5 ? atoi(argv[5]) : 1;
    NPTS = A * B;
    if (A < 1 || B < 1 || NPTS > MAXPTS || K < 3 || K > MAXK || K > NPTS || T < 1 || threads < 1 || threads > 64) {
        fprintf(stderr, "bad arguments\n");
        return 1;
    }
    // points ordered row by row, so every subset meeting row 0 has its first point in row 0
    for (int y = 0, i = 0; y < B; y++)
        for (int x = 0; x < A; x++, i++) {
            X[i] = x;
            Y[i] = y;
        }
    maxd = (A - 1) * (B - 1);
    pthread_t th[64];
    for (int t = 0; t < threads; t++) pthread_create(&th[t], NULL, worker, NULL);
    for (int t = 0; t < threads; t++) pthread_join(th[t], NULL);
    printf("done %lld\n", total_subsets);
    return 0;
}
