/* Exhaustive enumeration of the full subgroup lattice of S_n.
 *
 * Emits the exact multiset of subgroup orders -- the object Erdos asked about
 * in problem #1162 ("is there a statistical theorem on their order?") -- plus,
 * as a by-product, the exact number of elementary abelian 2-subgroups of each
 * rank.  The latter is an independent check on src/elemab.py, which computes
 * the same numbers from a generating function and shares no code with this.
 *
 * ---------------------------------------------------------------- algorithm
 *
 * Work class by class, not subgroup by subgroup.  S_8 has 151221 subgroups but
 * only 296 conjugacy classes of them, so expanding one representative per
 * class is ~500x less work.
 *
 *   expand(H):  find every subgroup of the form <H, g>.
 *               Since <H, g> = <H, hg> for h in H, it is enough to let g run
 *               over a transversal of the right cosets Hg.  We therefore scan
 *               g over G, skip any g already in a marked coset, and after each
 *               closure mark the whole coset Hg.  This is sound -- it never
 *               discards a value of <H, g> -- and costs |G| products in total
 *               per representative, with |G|/|H| closures.
 *
 *   register(K): if K is new, walk its entire conjugacy class by BFS under
 *               conjugation by two generators of S_n, recording every member
 *               in the global table.  Order, abelian-ness and elementary
 *               abelian rank are class invariants, computed once per class.
 *
 * Every subgroup L lies at the top of a maximal chain 1 = L_0 < ... < L_m = L
 * in which each L_{i+1} is a minimal overgroup of L_i, hence of the form
 * <L_i, g>; and minimal overgroups of L_i are conjugate to minimal overgroups
 * of the representative of its class.  So expanding representatives reaches
 * every class, and the orbit walk then reaches every subgroup.
 *
 * Subgroups are identified by a 128-bit order-independent hash of their
 * element set (two independent random 64-bit tags summed over the elements).
 * At under 2*10^6 subgroups a collision has probability below 10^-26.  The
 * results are in any case checked against two independently published
 * sequences: OEIS A005432 (subgroups) and A000638 (classes).
 *
 * Build:  cc -O3 -march=native -o subgroups subgroups.c
 * Usage:  ./subgroups <n> [out.json]
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

#if defined(__ARM_NEON)
#include <arm_neon.h>
#endif

#define MAXGENS 24

typedef struct { uint64_t lo, hi; } perm_t; /* 16 packed bytes, p(i) in byte i */

static int N;
static uint32_t NF;
static perm_t *PERM;
static uint64_t *R1, *R2;

/* ------------------------------------------------------------------ perms */

static inline perm_t pmul(perm_t p, perm_t q) /* (p*q)(i) = p(q(i)) */
{
#if defined(__ARM_NEON)
    uint8x16_t tbl = vreinterpretq_u8_u64(vcombine_u64(vcreate_u64(p.lo), vcreate_u64(p.hi)));
    uint8x16_t idx = vreinterpretq_u8_u64(vcombine_u64(vcreate_u64(q.lo), vcreate_u64(q.hi)));
    uint8x16_t res = vqtbl1q_u8(tbl, idx);
    perm_t r;
    r.lo = vgetq_lane_u64(vreinterpretq_u64_u8(res), 0);
    r.hi = vgetq_lane_u64(vreinterpretq_u64_u8(res), 1);
    return r;
#else
    perm_t r; const uint8_t *pb = (const uint8_t *)&p, *qb = (const uint8_t *)&q;
    uint8_t *rb = (uint8_t *)&r;
    for (int i = 0; i < 16; i++) rb[i] = pb[qb[i]];
    return r;
#endif
}

static inline int peq(perm_t a, perm_t b) { return a.lo == b.lo && a.hi == b.hi; }

static inline perm_t pinv(perm_t p)
{
    perm_t r; const uint8_t *pb = (const uint8_t *)&p; uint8_t *rb = (uint8_t *)&r;
    for (int i = 0; i < 16; i++) rb[pb[i]] = (uint8_t)i;
    return r;
}

static inline uint64_t mix64(uint64_t x)
{
    x ^= x >> 33; x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL;
    x ^= x >> 33; return x;
}

/* IDXTAB stores index+1, so 0 means "empty slot" */
static uint32_t *IDXTAB; static uint32_t IDXMASK;

static inline uint32_t perm_index(perm_t p)
{
    uint64_t h = mix64(p.lo ^ mix64(p.hi));
    uint32_t s = (uint32_t)(h & IDXMASK);
    for (;;) {
        uint32_t v = IDXTAB[s];
        if (v && peq(PERM[v - 1], p)) return v - 1;
        s = (s + 1) & IDXMASK;
    }
}

/* ------------------------------------------------- global table of subgroups */

typedef struct { uint64_t h1, h2; } key_t;
static key_t *SEEN; static uint64_t SEEN_MASK, SEEN_COUNT;

static void seen_init(uint64_t bits)
{
    SEEN_MASK = (1ULL << bits) - 1;
    SEEN = calloc(SEEN_MASK + 1, sizeof(key_t));
    SEEN_COUNT = 0;
}

static void seen_grow(void);

/* returns 1 if newly inserted, 0 if already present */
static int seen_add(uint64_t h1, uint64_t h2)
{
    uint64_t s = mix64(h1) & SEEN_MASK;
    for (;;) {
        key_t *e = &SEEN[s];
        if (e->h1 == 0 && e->h2 == 0) {
            e->h1 = h1; e->h2 = h2; SEEN_COUNT++;
            if (SEEN_COUNT * 2 > SEEN_MASK) seen_grow();
            return 1;
        }
        if (e->h1 == h1 && e->h2 == h2) return 0;
        s = (s + 1) & SEEN_MASK;
    }
}

static void seen_grow(void)
{
    key_t *old = SEEN; uint64_t oldmask = SEEN_MASK;
    SEEN_MASK = SEEN_MASK * 2 + 1;
    SEEN = calloc(SEEN_MASK + 1, sizeof(key_t));
    for (uint64_t i = 0; i <= oldmask; i++) {
        if (old[i].h1 == 0 && old[i].h2 == 0) continue;
        uint64_t s = mix64(old[i].h1) & SEEN_MASK;
        while (!(SEEN[s].h1 == 0 && SEEN[s].h2 == 0)) s = (s + 1) & SEEN_MASK;
        SEEN[s] = old[i];
    }
    free(old);
}

/* ------------------------------------------------------------------ closure */

static uint32_t *ELEMBUF;

static uint32_t closure(const uint32_t *gens, int ng, uint32_t *elems,
                        uint8_t *mark, uint64_t *h1, uint64_t *h2)
{
    uint32_t cnt = 0; uint64_t a = 0, b = 0;
    elems[cnt++] = 0; mark[0] = 1; a += R1[0]; b += R2[0];
    for (uint32_t head = 0; head < cnt; head++) {
        perm_t x = PERM[elems[head]];
        for (int s = 0; s < ng; s++) {
            uint32_t j = perm_index(pmul(x, PERM[gens[s]]));
            if (!mark[j]) { mark[j] = 1; elems[cnt++] = j; a += R1[j]; b += R2[j]; }
        }
    }
    for (uint32_t i = 0; i < cnt; i++) mark[elems[i]] = 0; /* leave mark[] clean */
    *h1 = a; *h2 = b;
    return cnt;
}

/* ------------------------------------------------------------------ classes */

typedef struct {
    uint32_t gens[MAXGENS];
    uint8_t  ngens;
    uint32_t order;
    uint64_t class_size;
    uint8_t  abelian;
    int8_t   elab_rank;   /* -1 if not elementary abelian of exponent <= 2 */
} class_t;

static class_t *CLS; static uint64_t NCLS, CLSCAP;

int main(int argc, char **argv)
{
    if (argc < 2) { fprintf(stderr, "usage: %s <n> [out.json]\n", argv[0]); return 2; }
    N = atoi(argv[1]);
    if (N < 1 || N > 9) { fprintf(stderr, "n must be in 1..9\n"); return 2; }
    clock_t t0 = clock();

    NF = 1; for (int i = 2; i <= N; i++) NF *= i;

    PERM = malloc((size_t)NF * sizeof(perm_t));
    { uint8_t a[16]; for (int i = 0; i < 16; i++) a[i] = (uint8_t)i;
      for (uint32_t c = 0; c < NF; c++) {
          memcpy(&PERM[c], a, 16);
          int i = N - 2; while (i >= 0 && a[i] >= a[i + 1]) i--;
          if (i < 0) break;
          int j = N - 1; while (a[j] <= a[i]) j--;
          uint8_t t = a[i]; a[i] = a[j]; a[j] = t;
          for (int l = i + 1, r = N - 1; l < r; l++, r--) { t = a[l]; a[l] = a[r]; a[r] = t; } } }

    IDXMASK = 1; while (IDXMASK < NF * 4u) IDXMASK <<= 1; IDXMASK -= 1;
    IDXTAB = calloc((size_t)IDXMASK + 1, 4);
    for (uint32_t i = 0; i < NF; i++) {
        uint64_t h = mix64(PERM[i].lo ^ mix64(PERM[i].hi));
        uint32_t s = (uint32_t)(h & IDXMASK);
        while (IDXTAB[s]) s = (s + 1) & IDXMASK;
        IDXTAB[s] = i + 1;
    }

    R1 = malloc((size_t)NF * 8); R2 = malloc((size_t)NF * 8);
    { uint64_t st = 0x9E3779B97F4A7C15ULL;
      for (uint32_t i = 0; i < NF; i++) {
          st += 0x9E3779B97F4A7C15ULL; R1[i] = mix64(st) | 1ULL;
          st += 0x9E3779B97F4A7C15ULL; R2[i] = mix64(st) | 1ULL; } }

    seen_init(16);
    CLSCAP = 256; CLS = malloc(CLSCAP * sizeof(class_t)); NCLS = 0;

    /* two generators of S_n: the transposition (0 1) and the n-cycle */
    uint32_t sgen[2]; int nsgen = 0;
    if (N >= 2) {
        uint8_t a[16]; perm_t p;
        for (int i = 0; i < 16; i++) a[i] = (uint8_t)i;
        a[0] = 1; a[1] = 0; memcpy(&p, a, 16); sgen[nsgen++] = perm_index(p);
        for (int i = 0; i < 16; i++) a[i] = (uint8_t)i;
        for (int i = 0; i < N; i++) a[i] = (uint8_t)((i + 1) % N);
        memcpy(&p, a, 16); sgen[nsgen++] = perm_index(p);
    }

    uint8_t *mark = calloc(NF, 1);
    uint8_t *cover = calloc(NF, 1);
    uint32_t *helem = malloc((size_t)NF * 4);
    uint32_t *kelem = malloc((size_t)NF * 4);
    ELEMBUF = malloc((size_t)NF * 4);

    /* orbit BFS scratch */
    uint32_t *oq = malloc((size_t)NF * MAXGENS * 4);

    /* register a class given a generating set; returns 1 if new */
    /* (written inline as a lambda-ish block via a helper macro is awkward in C,
       so it is a plain function using file-scope state) */

    /* --- seed: trivial subgroup --- */
    {
        uint64_t h1, h2;
        uint32_t ord = closure(NULL, 0, kelem, mark, &h1, &h2);
        seen_add(h1, h2);
        CLS[0].ngens = 0; CLS[0].order = ord; CLS[0].class_size = 1;
        CLS[0].abelian = 1; CLS[0].elab_rank = 0; NCLS = 1;
    }

    for (uint64_t ci = 0; ci < NCLS; ci++) {
        class_t H = CLS[ci];
        uint64_t hh1, hh2;
        uint32_t ho = closure(H.gens, H.ngens, helem, mark, &hh1, &hh2);
        if (ho != H.order) { fprintf(stderr, "internal: order mismatch\n"); return 1; }

        memset(cover, 0, NF);
        for (uint32_t i = 0; i < ho; i++) cover[helem[i]] = 1;

        uint32_t newgens[MAXGENS];
        memcpy(newgens, H.gens, sizeof(uint32_t) * H.ngens);
        if (H.ngens + 1 > MAXGENS) { fprintf(stderr, "MAXGENS too small\n"); return 1; }

        for (uint32_t g = 0; g < NF; g++) {
            if (cover[g]) continue;
            /* mark the whole right coset Hg: <H,g> = <H,hg> */
            perm_t pg = PERM[g];
            for (uint32_t i = 0; i < ho; i++)
                cover[perm_index(pmul(PERM[helem[i]], pg))] = 1;

            newgens[H.ngens] = g;
            uint64_t k1, k2;
            uint32_t ko = closure(newgens, H.ngens + 1, kelem, mark, &k1, &k2);
            if (!seen_add(k1, k2)) continue;   /* class already known */

            /* --- new class: record invariants, then walk the orbit --- */
            if (NCLS == CLSCAP) { CLSCAP *= 2; CLS = realloc(CLS, CLSCAP * sizeof(class_t)); }
            class_t *S = &CLS[NCLS];
            memcpy(S->gens, newgens, sizeof(uint32_t) * (H.ngens + 1));
            S->ngens = H.ngens + 1; S->order = ko;

            int ab = 1;
            for (int x = 0; x < S->ngens && ab; x++)
                for (int y = x + 1; y < S->ngens && ab; y++)
                    if (!peq(pmul(PERM[S->gens[x]], PERM[S->gens[y]]),
                             pmul(PERM[S->gens[y]], PERM[S->gens[x]]))) ab = 0;
            S->abelian = (uint8_t)ab;
            S->elab_rank = -1;
            if ((ko & (ko - 1)) == 0) {
                int e2 = 1;
                for (uint32_t t = 0; t < ko && e2; t++)
                    if (!peq(pmul(PERM[kelem[t]], PERM[kelem[t]]), PERM[0])) e2 = 0;
                if (e2) { int r = 0; while ((1u << r) < ko) r++; S->elab_rank = (int8_t)r; }
            }

            uint64_t size = 1, qh = 0, qt = 0;
            memcpy(oq + qt * MAXGENS, S->gens, sizeof(uint32_t) * S->ngens); qt++;
            while (qh < qt) {
                uint32_t cur[MAXGENS];
                memcpy(cur, oq + qh * MAXGENS, sizeof(uint32_t) * S->ngens); qh++;
                for (int s = 0; s < nsgen; s++) {
                    perm_t x = PERM[sgen[s]], xi = pinv(x);
                    uint32_t cg[MAXGENS];
                    for (int t = 0; t < S->ngens; t++)
                        cg[t] = perm_index(pmul(pmul(xi, PERM[cur[t]]), x));
                    uint64_t c1, c2;
                    closure(cg, S->ngens, ELEMBUF, mark, &c1, &c2);
                    if (seen_add(c1, c2)) {
                        size++;
                        memcpy(oq + qt * MAXGENS, cg, sizeof(uint32_t) * S->ngens); qt++;
                    }
                }
            }
            S->class_size = size;
            NCLS++;
        }
        if ((ci & 63) == 0)
            fprintf(stderr, "\r  class %llu/%llu  subgroups %llu    ",
                    (unsigned long long)ci, (unsigned long long)NCLS,
                    (unsigned long long)SEEN_COUNT);
    }
    fprintf(stderr, "\r                                                        \r");

    uint64_t total = 0;
    for (uint64_t i = 0; i < NCLS; i++) total += CLS[i].class_size;
    printf("n=%d  classes=%llu  subgroups=%llu  (%.1fs)\n", N,
           (unsigned long long)NCLS, (unsigned long long)total,
           (double)(clock() - t0) / CLOCKS_PER_SEC);
    if (total != SEEN_COUNT) {
        fprintf(stderr, "internal: class sizes sum to %llu but table holds %llu\n",
                (unsigned long long)total, (unsigned long long)SEEN_COUNT);
        return 1;
    }

    if (argc >= 3) {
        uint64_t *oh = calloc((size_t)NF + 1, 8), *ah = calloc((size_t)NF + 1, 8);
        uint64_t *ch = calloc((size_t)NF + 1, 8);       /* classes by order */
        uint64_t er[16]; memset(er, 0, sizeof er);
        for (uint64_t i = 0; i < NCLS; i++) {
            oh[CLS[i].order] += CLS[i].class_size;
            ch[CLS[i].order] += 1;
            if (CLS[i].abelian) ah[CLS[i].order] += CLS[i].class_size;
            if (CLS[i].elab_rank >= 0) er[CLS[i].elab_rank] += CLS[i].class_size;
        }
        FILE *o = fopen(argv[2], "w");
        fprintf(o, "{\n \"n\": %d,\n \"total_subgroups\": %llu,\n \"total_classes\": %llu,\n",
                N, (unsigned long long)total, (unsigned long long)NCLS);
        const char *names[3] = {"order_histogram", "abelian_order_histogram", "class_order_histogram"};
        uint64_t *arrs[3] = {oh, ah, ch};
        for (int k = 0; k < 3; k++) {
            fprintf(o, " \"%s\": {", names[k]); int first = 1;
            for (uint32_t x = 1; x <= NF; x++) if (arrs[k][x]) {
                fprintf(o, "%s\"%u\": %llu", first ? "" : ", ", x,
                        (unsigned long long)arrs[k][x]); first = 0; }
            fprintf(o, "},\n");
        }
        fprintf(o, " \"elementary_abelian_rank\": [");
        for (int r = 0; r <= N / 2; r++)
            fprintf(o, "%s%llu", r ? ", " : "", (unsigned long long)er[r]);
        fprintf(o, "]\n}\n");
        fclose(o);
    }
    return 0;
}
