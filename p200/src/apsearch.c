/* apsearch.c -- exhaustive enumeration of k-term arithmetic progressions of
 * primes whose last term is <= B.
 *
 * Structure theorem used (proved in NOTES.md, and re-derived by the brute-force
 * reference implementation in ref_bruteforce.py):
 *
 *   Let a, a+d, ..., a+(k-1)d be primes, d >= 1, k >= 3, and let p <= k be prime.
 *   If p does not divide d then the k terms cover every residue class mod p, so
 *   some term is divisible by p; being prime it equals p.  Hence
 *        p | d   or   p is a term of the progression.
 *   If p is a term then a <= p <= k.  Writing a = p:  p | d would make a+d a
 *   proper multiple of p, so p does not divide d, and then a+pd is a multiple of
 *   p exceeding p unless k <= p.  So p = k.  Therefore exactly two cases:
 *
 *   (A)  a > k  and  W_k := prod_{p<=k} p  divides d.  Every term is then prime
 *        and coprime to W_k, so all k terms lie in one residue class r mod W_k
 *        with gcd(r, W_k) = 1.
 *   (B)  k is prime, a = k, and (W_k / k) | d.
 *
 * Case A is enumerated residue class by residue class: for each r coprime to
 * W_k we sieve the progression r, r+W, r+2W, ... <= B into a bitset and look for
 * k bits in arithmetic progression.  The working set per class is tiny
 * (B/(8 W_k) bytes), so the search runs out of L1/L2 cache.
 *
 * Case B is a short separate loop using deterministic Miller-Rabin.
 *
 * Output: the minimal last term over all k-APs found, the total number of
 * k-APs with last term <= B, and a cumulative histogram of last terms in
 * buckets of 2^(1/16) so that counts up to intermediate thresholds can be
 * recovered.
 *
 * usage: ./apsearch k B [nthreads]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <pthread.h>
#include <time.h>

typedef unsigned __int128 u128;

#define NBUCK 1200 /* covers last terms < 2^75 */

static uint64_t K, B, W;
static uint64_t *sp;   /* small primes q <= sqrt(B) with q not dividing W */
static uint64_t *spinv;/* W^{-1} mod q */
static uint64_t *spjq; /* if q == r + j*W for the current class, that j; else UINT64_MAX */
static int nsp;
static uint64_t Mmax;  /* max positions in a class */
static int NT = 1;

/* ---------------- bit helpers ---------------- */
static inline void clearbit(uint64_t *b, uint64_t i) { b[i >> 6] &= ~(1ULL << (i & 63)); }
static inline uint64_t loadw(const uint64_t *b, uint64_t p) {
    uint64_t i = p >> 6, s = p & 63;
    uint64_t lo = b[i] >> s;
    return s ? (lo | (b[i + 1] << (64 - s))) : lo;
}
static inline int bucket_of(uint64_t x) {
    int e = 63 - __builtin_clzll(x);
    uint64_t m = (e >= 4) ? ((x >> (e - 4)) & 15u) : ((x << (4 - e)) & 15u);
    int idx = e * 16 + (int)m;
    if (idx >= NBUCK) idx = NBUCK - 1;
    return idx;
}

/* ---------------- deterministic Miller-Rabin for n < 3.3e24 ---------------- */
static uint64_t mulmod(uint64_t a, uint64_t b, uint64_t m) { return (uint64_t)((u128)a * b % m); }
static uint64_t powmod(uint64_t a, uint64_t e, uint64_t m) {
    uint64_t r = 1; a %= m;
    while (e) { if (e & 1) r = mulmod(r, a, m); a = mulmod(a, a, m); e >>= 1; }
    return r;
}
static int is_prime_u64(uint64_t n) {
    if (n < 2) return 0;
    static const uint64_t sm[] = {2,3,5,7,11,13,17,19,23,29,31,37};
    for (int i = 0; i < 12; i++) { if (n % sm[i] == 0) return n == sm[i]; }
    uint64_t d = n - 1; int s = 0;
    while (!(d & 1)) { d >>= 1; s++; }
    for (int i = 0; i < 12; i++) {
        uint64_t x = powmod(sm[i], d, n);
        if (x == 1 || x == n - 1) continue;
        int ok = 0;
        for (int j = 1; j < s; j++) { x = mulmod(x, x, n); if (x == n - 1) { ok = 1; break; } }
        if (!ok) return 0;
    }
    return 1;
}

/* ---------------- per-thread state ---------------- */
typedef struct {
    int id;
    uint64_t *bits;
    uint64_t *buck;
    uint64_t minlast, mina, mind;
    uint64_t total;
    uint64_t *jq;
} thr_t;

static uint64_t *rlist; static uint64_t nr; /* residues coprime to W */

static void *worker(void *arg) {
    thr_t *T = (thr_t *)arg;
    uint64_t *bits = T->bits, *jq = T->jq;
    uint64_t nwmax = Mmax / 64 + 4;
    for (uint64_t ri = T->id; ri < nr; ri += NT) {
        uint64_t r = rlist[ri];
        if (r > B) continue;
        uint64_t M = (B - r) / W + 1;              /* j = 0..M-1, n = r + jW <= B */
        if (M < K) continue;
        uint64_t nw = M / 64 + 4;
        if (nw > nwmax) nw = nwmax;
        memset(bits, 0xFF, nw * 8);
        for (uint64_t i = M; i < nw * 64; i++) clearbit(bits, i);
        if (r == 1) clearbit(bits, 0);             /* n = 1 is not prime */
        for (int t = 0; t < nsp; t++) {
            uint64_t q = sp[t];
            uint64_t rq = r % q;
            uint64_t j0 = rq ? ((q - rq) * spinv[t]) % q : 0;
            if (j0 >= M) { jq[t] = UINT64_MAX; continue; }
            uint64_t skip = UINT64_MAX;
            if (q >= r && (q - r) % W == 0) skip = (q - r) / W;   /* n == q itself */
            for (uint64_t j = j0; j < M; j += q) if (j != skip) clearbit(bits, j);
            jq[t] = skip;
        }
        /* search for K bits in AP */
        const uint64_t *__restrict b = bits;
        uint64_t maxm = (M - 1) / (K - 1);
        for (uint64_t m = 1; m <= maxm; m++) {
            uint64_t lim = M - (K - 1) * m;        /* j in [0, lim) */
            for (uint64_t p0 = 0; p0 < lim; p0 += 64) {
                uint64_t x = b[p0 >> 6];
                if (p0 + 64 > lim) x &= (~0ULL) >> (64 - (lim - p0));
                if (!x) continue;
                x &= loadw(b, p0 + m);
                if (!x) continue;
                uint64_t off = p0 + 2 * m;
                for (uint64_t i = 2; i < K; i++, off += m) {
                    x &= loadw(b, off);
                    if (!x) break;
                }
                while (x) {
                    uint64_t j = p0 + (uint64_t)__builtin_ctzll(x);
                    x &= x - 1;
                    uint64_t a = r + j * W, d = m * W;
                    uint64_t last = a + (K - 1) * d;
                    T->total++;
                    T->buck[bucket_of(last)]++;
                    if (last < T->minlast) { T->minlast = last; T->mina = a; T->mind = d; }
                }
            }
        }
    }
    return NULL;
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s k B [nthreads]\n", argv[0]); return 1; }
    K = strtoull(argv[1], 0, 10);
    B = strtoull(argv[2], 0, 10);
    NT = (argc > 3) ? atoi(argv[3]) : 1;
    if (K < 3) { fprintf(stderr, "need k >= 3\n"); return 1; }

    /* W = product of primes <= K */
    W = 1;
    for (uint64_t p = 2; p <= K; p++) { int ip = 1; for (uint64_t t = 2; t * t <= p; t++) if (p % t == 0) { ip = 0; break; } if (ip) W *= p; }

    /* small primes up to sqrt(B) */
    uint64_t S = (uint64_t)sqrtl((long double)B) + 2;
    while (S * S > B) S--;
    S += 1;
    char *cs = calloc(S + 1, 1);
    for (uint64_t i = 2; i * i <= S; i++) if (!cs[i]) for (uint64_t j = i * i; j <= S; j += i) cs[j] = 1;
    nsp = 0;
    for (uint64_t i = 2; i <= S; i++) if (!cs[i] && W % i) nsp++;
    sp = malloc(sizeof(uint64_t) * nsp); spinv = malloc(sizeof(uint64_t) * nsp);
    { int t = 0; for (uint64_t i = 2; i <= S; i++) if (!cs[i] && W % i) sp[t++] = i; }
    for (int t = 0; t < nsp; t++) {           /* W^{-1} mod q by Fermat */
        uint64_t q = sp[t], wm = W % q, inv = 1, b = wm, e = q - 2;
        while (e) { if (e & 1) inv = inv * b % q; b = b * b % q; e >>= 1; }
        spinv[t] = inv;
    }
    free(cs);

    /* residues coprime to W */
    rlist = malloc(sizeof(uint64_t) * (W / 2 + 8)); nr = 0;
    for (uint64_t r = 1; r < W; r++) {
        uint64_t a = r, b = W; while (b) { uint64_t t = a % b; a = b; b = t; }
        if (a == 1) rlist[nr++] = r;
    }
    Mmax = B / W + 2;

    fprintf(stderr, "k=%llu B=%llu W=%llu phi(W)=%llu Mmax=%llu nsp=%d threads=%d\n",
            (unsigned long long)K, (unsigned long long)B, (unsigned long long)W,
            (unsigned long long)nr, (unsigned long long)Mmax, nsp, NT);

    struct timespec t0, t1; clock_gettime(CLOCK_MONOTONIC, &t0);

    pthread_t *th = malloc(sizeof(pthread_t) * NT);
    thr_t *TS = calloc(NT, sizeof(thr_t));
    for (int i = 0; i < NT; i++) {
        TS[i].id = i;
        TS[i].bits = malloc((Mmax / 64 + 4) * 8);
        TS[i].buck = calloc(NBUCK, sizeof(uint64_t));
        TS[i].jq = malloc(sizeof(uint64_t) * nsp);
        TS[i].minlast = UINT64_MAX;
        pthread_create(&th[i], NULL, worker, &TS[i]);
    }
    for (int i = 0; i < NT; i++) pthread_join(th[i], NULL);

    uint64_t minlast = UINT64_MAX, mina = 0, mind = 0, total = 0;
    uint64_t *buck = calloc(NBUCK, sizeof(uint64_t));
    for (int i = 0; i < NT; i++) {
        total += TS[i].total;
        for (int b = 0; b < NBUCK; b++) buck[b] += TS[i].buck[b];
        if (TS[i].minlast < minlast) { minlast = TS[i].minlast; mina = TS[i].mina; mind = TS[i].mind; }
    }

    /* Case B: k prime, a = k, (W/k) | d */
    int kprime = (K >= 2); for (uint64_t t = 2; t * t <= K; t++) if (K % t == 0) { kprime = 0; break; }
    uint64_t caseB = 0;
    if (kprime) {
        uint64_t Wk = W / K;
        for (uint64_t d = Wk; d <= (B - K) / (K - 1); d += Wk) {
            int ok = 1;
            for (uint64_t i = 1; i < K; i++) if (!is_prime_u64(K + i * d)) { ok = 0; break; }
            if (ok) {
                uint64_t last = K + (K - 1) * d;
                caseB++; total++; buck[bucket_of(last)]++;
                if (last < minlast) { minlast = last; mina = K; mind = d; }
            }
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &t1);
    double el = (t1.tv_sec - t0.tv_sec) + 1e-9 * (t1.tv_nsec - t0.tv_nsec);

    printf("{\n  \"k\": %llu,\n  \"B\": %llu,\n  \"W\": %llu,\n  \"caseB_count\": %llu,\n",
           (unsigned long long)K, (unsigned long long)B, (unsigned long long)W, (unsigned long long)caseB);
    if (minlast == UINT64_MAX) printf("  \"min_last\": null,\n");
    else printf("  \"min_last\": %llu,\n  \"min_first\": %llu,\n  \"min_diff\": %llu,\n",
                (unsigned long long)minlast, (unsigned long long)mina, (unsigned long long)mind);
    printf("  \"total\": %llu,\n  \"seconds\": %.2f,\n  \"hist\": {", (unsigned long long)total, el);
    int first = 1;
    for (int b = 0; b < NBUCK; b++) if (buck[b]) {
        double lo = ldexp(1.0 + (b % 16) / 16.0, b / 16);
        printf("%s\n    \"%.10g\": %llu", first ? "" : ",", lo, (unsigned long long)buck[b]);
        first = 0;
    }
    printf("\n  }\n}\n");
    return 0;
}
