/*
 * census16.c -- the subgroups of R = prod_i R_{b_i} that are transitive on every block, for one
 * orbit type (b_1, ..., b_t), b_i in {2, 4, 8, 16}, sum b_i <= 16, each with the numbers
 * s_{b_i}(H|O_i) of Sylow 2-subgroups of Sym(O_i) that contain its projections.
 *
 * R_b is the Sylow 2-subgroup C2 wr ... wr C2 of S_b preserving the nest of consecutive pairs,
 * quadruples, ...; R_16 = (R_8 x R_8) : <sigma>, sigma = (0 8)(1 9)...(7 15).  It is the
 * C counterpart of count_type() in two_subgroups.py, written for the orbit types of k = 16,
 * whose R have up to 2^15 elements.  Method, as there:
 *
 *   - descent from R through maximal subgroups: the maximal subgroups of M are the preimages
 *     of the hyperplanes of M/Phi(M), Phi(M) = <squares>; a candidate that is intransitive on
 *     some block is dropped together with everything below it (a subgroup of an intransitive
 *     group is intransitive);
 *   - every subgroup is stored once, keyed exactly: a generating set and the order; a candidate
 *     equals a stored subgroup iff the orders agree and every stored generator lies in it
 *     (a 64-bit hash of the element bitset only picks the bucket);
 *   - s_b(P), for P transitive on b points, is the number of flags B_2 < B_4 < ... < B_{b/2}
 *     of blocks of P containing the first point (equivalently, of maximal chains from a point
 *     stabiliser up to P), found with minimal-block closures; a Sylow 2-subgroup of S_b is the
 *     stabiliser of such a nest, so this is the number of Sylow 2-subgroups containing P;
 *   - Frattini-closed (FC): (h^2)|O_i lies in H for every h in H and every block of size >= 4.
 *
 * Output (stdout): one line per class of (log2|H|, FC, elementary abelian, s_1, ..., s_t),
 *   "key j fc ea d s_1 ... s_t count",
 *   with j = log2|H| and d = d(H) = log2|H/Phi(H)|, then "listed <number of subgroups listed>" and
 *   "layer <j> <count>" lines.  Progress goes to stderr.  The Python driver
 *   src/two_subgroups16.py turns these into exact weighted counts.
 *
 * Usage: census16 b_1 b_2 ... b_t      e.g. census16 8 8, census16 16
 *
 * Exit status 0 on success, 2 on a usage error or if an invariant fails: |R_b| = 2^(b-1), the
 * root generators give all of R, a stored order equals the closure of the stored generators,
 * M/Phi(M) is elementary abelian, d(M) <= 8, and the fixed table sizes are not exceeded.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef uint16_t el;
#define MAXF 8
#define MAXR 8 /* d(H) <= n/2 = 8 for a 2-subgroup of S_16 */

/* ------------------------------------------------------------------ R_2, R_4, R_8 */
typedef struct {
  int b, n;
  uint8_t perm[128][8];
  uint8_t mul[128][128];
} Small;
static Small SM[4]; /* SM[1] = R_2, SM[2] = R_4, SM[3] = R_8 */

static void build_small(int lg) {
  Small *S = &SM[lg];
  int b = 1 << lg;
  S->b = b;
  uint8_t gens[3][8];
  int ng = 0;
  for (int s = 1; s < b; s *= 2, ng++) {
    for (int i = 0; i < b; i++) gens[ng][i] = (uint8_t)i;
    for (int i = 0; i < s; i++) { gens[ng][i] = (uint8_t)(i + s); gens[ng][i + s] = (uint8_t)i; }
  }
  for (int i = 0; i < b; i++) S->perm[0][i] = (uint8_t)i;
  S->n = 1;
  for (int f = 0; f < S->n; f++) {
    for (int g = 0; g < ng; g++) {
      uint8_t q[8];
      for (int i = 0; i < b; i++) q[i] = S->perm[f][gens[g][i]]; /* perm[f] after gen */
      int found = 0;
      for (int e = 0; e < S->n && !found; e++) found = memcmp(S->perm[e], q, (size_t)b) == 0;
      if (!found) memcpy(S->perm[S->n++], q, (size_t)b);
    }
  }
  if (S->n != (1 << (b - 1))) { fprintf(stderr, "R_%d has %d elements\n", b, S->n); exit(2); }
  for (int x = 0; x < S->n; x++)
    for (int y = 0; y < S->n; y++) {
      uint8_t q[8];
      for (int i = 0; i < b; i++) q[i] = S->perm[x][S->perm[y][i]]; /* x after y */
      int z = -1;
      for (int e = 0; e < S->n; e++)
        if (memcmp(S->perm[e], q, (size_t)b) == 0) { z = e; break; }
      S->mul[x][y] = (uint8_t)z;
    }
}

/* ------------------------------------------------------------------ the ambient R */
static int T, K, NBITS;           /* number of blocks, points, bits of an element index */
static int B[MAXF], OFF[MAXF], SH[MAXF], LG[MAXF];
static uint32_t MK[MAXF];
static uint32_t N;                 /* |R| = 2^NBITS */
static el *TAB;                    /* full multiplication table when N <= 16384 */
static el *SQ;                     /* squares */
static uint8_t *PERM;              /* PERM[x*16 + i] = image of point i under x */

static inline uint32_t mul16(uint32_t x, uint32_t y) {
  uint32_t a = x & 127, b = (x >> 7) & 127, s = x >> 14;
  uint32_t a2 = y & 127, b2 = (y >> 7) & 127, s2 = y >> 14;
  const Small *S = &SM[3];
  if (!s) return S->mul[a][a2] | (uint32_t)S->mul[b][b2] << 7 | s2 << 14;
  return S->mul[a][b2] | (uint32_t)S->mul[b][a2] << 7 | (s2 ^ 1u) << 14;
}

static inline el mul_slow(el x, el y) {
  uint32_t z = 0;
  for (int i = 0; i < T; i++) {
    uint32_t cx = (x >> SH[i]) & MK[i], cy = (y >> SH[i]) & MK[i], cz;
    cz = LG[i] == 4 ? mul16(cx, cy) : SM[LG[i]].mul[cx][cy];
    z |= cz << SH[i];
  }
  return (el)z;
}

static inline el mul(el x, el y) { return TAB ? TAB[(uint32_t)x * N + y] : mul_slow(x, y); }

static void comp_perm(int i, uint32_t c, uint8_t *out) { /* local images of block i under component c */
  int b = B[i];
  if (LG[i] < 4) { memcpy(out, SM[LG[i]].perm[c], (size_t)b); return; }
  uint32_t a = c & 127, bb = (c >> 7) & 127, s = c >> 14;
  for (int x = 0; x < 16; x++) {
    if (!s) out[x] = x < 8 ? SM[3].perm[a][x] : (uint8_t)(SM[3].perm[bb][x - 8] + 8);
    else out[x] = x < 8 ? (uint8_t)(SM[3].perm[bb][x] + 8) : SM[3].perm[a][x - 8];
  }
}

static void build_ambient(void) {
  NBITS = 0;
  K = 0;
  for (int i = 0; i < T; i++) {
    LG[i] = __builtin_ctz((unsigned)B[i]);
    OFF[i] = K;
    SH[i] = NBITS;
    MK[i] = (1u << (B[i] - 1)) - 1;
    K += B[i];
    NBITS += B[i] - 1;
  }
  N = 1u << NBITS;
  if (N <= 16384) {
    TAB = malloc((size_t)N * N * sizeof(el));
    for (uint32_t x = 0; x < N; x++)
      for (uint32_t y = 0; y < N; y++) TAB[(size_t)x * N + y] = mul_slow((el)x, (el)y);
  }
  SQ = malloc(N * sizeof(el));
  PERM = malloc((size_t)N * 16);
  for (uint32_t x = 0; x < N; x++) {
    SQ[x] = mul((el)x, (el)x);
    for (int i = 0; i < T; i++) {
      uint8_t loc[16];
      comp_perm(i, (x >> SH[i]) & MK[i], loc);
      for (int p = 0; p < B[i]; p++) PERM[x * 16 + OFF[i] + p] = (uint8_t)(OFF[i] + loc[p]);
    }
  }
}

/* ------------------------------------------------------------------ s_b by block flags */
static int uf[16];
static int ufind(int x) { while (uf[x] != x) x = uf[x] = uf[uf[x]]; return x; }

/* the smallest block of <gens> (on b points) containing the points of seed; returned as a mask */
static uint32_t minblock(uint8_t (*gens)[16], int ng, int b, uint32_t seed) {
  int qa[256], qb[256], qn = 0;
  for (int i = 0; i < b; i++) uf[i] = i;
  int s0 = __builtin_ctz(seed);
  for (int i = 0; i < b; i++)
    if ((seed >> i & 1) && i != s0) { uf[ufind(i)] = ufind(s0); qa[qn] = s0; qb[qn] = i; qn++; }
  for (int q = 0; q < qn; q++)
    for (int g = 0; g < ng; g++) {
      int x = ufind(gens[g][qa[q]]), y = ufind(gens[g][qb[q]]);
      if (x != y) { uf[x] = y; qa[qn] = x; qb[qn] = y; qn++; }
    }
  uint32_t m = 0;
  int r0 = ufind(0);
  for (int i = 0; i < b; i++)
    if (ufind(i) == r0) m |= 1u << i;
  return m;
}

static long sylow_count(uint8_t (*gens)[16], int ng, int b) {
  if (b == 2) return 1;
  uint32_t lev[5][1024];
  long ways[5][1024];
  int cnt[5];
  lev[0][0] = 1;
  ways[0][0] = 1;
  cnt[0] = 1;
  int L = 0;
  for (int m = 1; 2 * m < b; m *= 2, L++) {
    cnt[L + 1] = 0;
    for (int u = 0; u < cnt[L]; u++)
      for (int x = 0; x < b; x++) {
        if (lev[L][u] >> x & 1) continue;
        uint32_t c = minblock(gens, ng, b, lev[L][u] | 1u << x);
        if (__builtin_popcount(c) != 2 * m) continue;
        int seen = 0;
        for (int v = 0; v < cnt[L + 1]; v++)
          if (lev[L + 1][v] == c) { seen = 1; break; }
        if (!seen) {
          if (cnt[L + 1] == 1024) { fprintf(stderr, "too many blocks\n"); exit(2); }
          lev[L + 1][cnt[L + 1]++] = c;
        }
      }
    for (int v = 0; v < cnt[L + 1]; v++) {
      ways[L + 1][v] = 0;
      for (int u = 0; u < cnt[L]; u++)
        if ((lev[L][u] & ~lev[L + 1][v]) == 0) ways[L + 1][v] += ways[L][u];
    }
  }
  long s = 0;
  for (int v = 0; v < cnt[L]; v++) s += ways[L][v];
  return s;
}

/* ------------------------------------------------------------------ subgroup store */
typedef struct {
  uint8_t ng, lg;
  el g[15];
} Sub;
static Sub *SUBS;
static uint64_t *HASH;
static size_t NSUB, CAPSUB;
static uint32_t *HIDX; /* open addressing, entry index + 1 */
static size_t HCAP;

static uint64_t hash_bits(const uint64_t *w, size_t nw) {
  uint64_t h = 0x9E3779B97F4A7C15ull;
  for (size_t i = 0; i < nw; i++) {
    h ^= w[i] + 0x9E3779B97F4A7C15ull + (h << 6) + (h >> 2);
    h *= 0xff51afd7ed558ccdull;
  }
  h ^= h >> 33;
  h *= 0xc4ceb9fe1a85ec53ull;
  h ^= h >> 33;
  return h;
}

static void hgrow(void) {
  size_t nc = HCAP ? HCAP * 2 : (1u << 16);
  uint32_t *nh = calloc(nc, sizeof(uint32_t));
  for (size_t i = 0; i < HCAP; i++)
    if (HIDX[i]) {
      size_t p = HASH[HIDX[i] - 1] & (nc - 1);
      while (nh[p]) p = (p + 1) & (nc - 1);
      nh[p] = HIDX[i];
    }
  free(HIDX);
  HIDX = nh;
  HCAP = nc;
}

static int insert(const uint64_t *bits, size_t nw, int lg, const el *gens, int ng) {
  if ((NSUB + 1) * 2 > HCAP) hgrow();
  uint64_t h = hash_bits(bits, nw);
  size_t p = h & (HCAP - 1);
  while (HIDX[p]) {
    size_t e = HIDX[p] - 1;
    if (HASH[e] == h && SUBS[e].lg == lg) {
      int all = 1;
      for (int i = 0; i < SUBS[e].ng && all; i++) all = bits[SUBS[e].g[i] >> 6] >> (SUBS[e].g[i] & 63) & 1;
      if (all) return 0;
    }
    p = (p + 1) & (HCAP - 1);
  }
  if (NSUB == CAPSUB) {
    CAPSUB = CAPSUB ? CAPSUB * 2 : (1u << 16);
    SUBS = realloc(SUBS, CAPSUB * sizeof(Sub));
    HASH = realloc(HASH, CAPSUB * sizeof(uint64_t));
  }
  if (ng > 15) { fprintf(stderr, "too many generators\n"); exit(2); }
  SUBS[NSUB].ng = (uint8_t)ng;
  SUBS[NSUB].lg = (uint8_t)lg;
  memcpy(SUBS[NSUB].g, gens, (size_t)ng * sizeof(el));
  HASH[NSUB] = h;
  NSUB++;
  HIDX[p] = (uint32_t)NSUB;
  return 1;
}

/* ------------------------------------------------------------------ classes of the result */
#define MAXKEYS 4096
static int NKEY;
static int KEYJ[MAXKEYS], KEYFC[MAXKEYS], KEYEA[MAXKEYS], KEYD[MAXKEYS];
static long KEYS[MAXKEYS][MAXF];
static long KEYCNT[MAXKEYS];
static long LAYER[16];

static void tally(int j, int fc, int ea, int d, const long *s) {
  for (int k = 0; k < NKEY; k++)
    if (KEYJ[k] == j && KEYFC[k] == fc && KEYEA[k] == ea && KEYD[k] == d && memcmp(KEYS[k], s, (size_t)T * sizeof(long)) == 0) { KEYCNT[k]++; return; }
  if (NKEY == MAXKEYS) { fprintf(stderr, "too many keys\n"); exit(2); }
  KEYJ[NKEY] = j; KEYFC[NKEY] = fc; KEYEA[NKEY] = ea; KEYD[NKEY] = d;
  memcpy(KEYS[NKEY], s, (size_t)T * sizeof(long));
  KEYCNT[NKEY++] = 1;
}

/* cache of s_b for projections onto blocks of size <= 8, keyed by the projection as a subset of R_b */
typedef struct { uint64_t lo, hi; long s; } SCache;
static SCache SC[4][1024];
static int NSC[4];

static long s_of_block(int i, const el *gens, int ng) {
  int b = B[i];
  uint8_t pg[15][16];
  for (int g = 0; g < ng; g++)
    for (int p = 0; p < b; p++) pg[g][p] = (uint8_t)(PERM[(uint32_t)gens[g] * 16 + OFF[i] + p] - OFF[i]);
  if (LG[i] == 4) return sylow_count(pg, ng, b);
  /* projection as a subset of R_b: closure of the components */
  const Small *S = &SM[LG[i]];
  uint8_t in[128] = {0}, list[128];
  int n = 1;
  list[0] = 0;
  in[0] = 1;
  for (int q = 0; q < n; q++)
    for (int g = 0; g < ng; g++) {
      uint8_t y = S->mul[list[q]][(gens[g] >> SH[i]) & MK[i]];
      if (!in[y]) { in[y] = 1; list[n++] = y; }
    }
  uint64_t lo = 0, hi = 0;
  for (int x = 0; x < S->n; x++)
    if (in[x]) { if (x < 64) lo |= 1ull << x; else hi |= 1ull << (x - 64); }
  for (int c = 0; c < NSC[LG[i]]; c++)
    if (SC[LG[i]][c].lo == lo && SC[LG[i]][c].hi == hi) return SC[LG[i]][c].s;
  long s = sylow_count(pg, ng, b);
  SCache *e = &SC[LG[i]][NSC[LG[i]]++];
  e->lo = lo; e->hi = hi; e->s = s;
  return s;
}

/* ------------------------------------------------------------------ the descent */
static uint64_t *MEM, *PHI, *CB;
static el *ELS, *PHIEL;
static int16_t *COORD;
static size_t NW;

static uint32_t closure(const el *gens, int ng, uint64_t *mem, el *els) {
  memset(mem, 0, NW * 8);
  uint32_t n = 1;
  els[0] = 0;
  mem[0] |= 1;
  for (uint32_t q = 0; q < n; q++)
    for (int g = 0; g < ng; g++) {
      el y = mul(els[q], gens[g]);
      if (!(mem[y >> 6] >> (y & 63) & 1)) { mem[y >> 6] |= 1ull << (y & 63); els[n++] = y; }
    }
  return n;
}

static void process(size_t idx) {
  Sub M = SUBS[idx];
  uint32_t n = closure(M.g, M.ng, MEM, ELS);
  if ((1u << M.lg) != n) { fprintf(stderr, "order mismatch\n"); exit(2); }

  /* Phi(M) = <squares>, with a generating set chosen greedily */
  memset(PHI, 0, NW * 8);
  PHI[0] |= 1;
  PHIEL[0] = 0;
  uint32_t np = 1;
  el pg[15];
  int npg = 0;
  for (uint32_t u = 0; u < n; u++) {
    el s = SQ[ELS[u]];
    if (PHI[s >> 6] >> (s & 63) & 1) continue;
    pg[npg++] = s;
    uint32_t old = np;
    for (uint32_t v = 0; v < old; v++) {
      el y = mul(PHIEL[v], s);
      if (!(PHI[y >> 6] >> (y & 63) & 1)) { PHI[y >> 6] |= 1ull << (y & 63); PHIEL[np++] = y; }
    }
    for (uint32_t v = old; v < np; v++)
      for (int g = 0; g < npg; g++) {
        el y = mul(PHIEL[v], pg[g]);
        if (!(PHI[y >> 6] >> (y & 63) & 1)) { PHI[y >> 6] |= 1ull << (y & 63); PHIEL[np++] = y; }
      }
  }

  /* classification */
  int ea = np == 1;
  int fc = 1;
  int big = 0;
  for (int i = 0; i < T; i++) big += B[i] >= 4;
  if (big >= 2)
    for (uint32_t u = 0; u < n && fc; u++)
      for (int i = 0; i < T && fc; i++)
        if (B[i] >= 4) {
          el r = (el)(SQ[ELS[u]] & (MK[i] << SH[i]));
          fc = MEM[r >> 6] >> (r & 63) & 1;
        }

  /* coordinates on M/Phi(M) */
  for (uint32_t u = 0; u < n; u++) COORD[ELS[u]] = -1;
  for (uint32_t v = 0; v < np; v++) COORD[PHIEL[v]] = 0;
  /* reuse PHIEL as the list of assigned elements */
  uint32_t na = np;
  el basis[MAXR];
  int r = 0;
  for (uint32_t u = 0; u < n; u++) {
    el x = ELS[u];
    if (COORD[x] >= 0) continue;
    if (r == MAXR) { fprintf(stderr, "d(M) > %d\n", MAXR); exit(2); }
    basis[r] = x;
    uint32_t L = na;
    for (uint32_t v = 0; v < L; v++) {
      el z = mul(PHIEL[v], x);
      COORD[z] = (int16_t)(COORD[PHIEL[v]] | 1 << r);
      PHIEL[na++] = z;
    }
    r++;
  }
  if (na != n || (np << r) != n) { fprintf(stderr, "M/Phi(M) is not elementary abelian\n"); exit(2); }
  long s[MAXF];
  for (int i = 0; i < T; i++) s[i] = s_of_block(i, M.g, M.ng);
  tally(M.lg, fc, ea, r, s);
  LAYER[M.lg]++;
  if (r == 0) return;

  /* transitivity: hit[i][p] = coordinates of the elements taking the first point of block i to p */
  uint64_t hit[16][4];
  memset(hit, 0, sizeof hit);
  for (uint32_t u = 0; u < n; u++) {
    el x = ELS[u];
    int c = COORD[x];
    for (int i = 0; i < T; i++) {
      int p = PERM[(uint32_t)x * 16 + OFF[i]];
      hit[p][c >> 6] |= 1ull << (c & 63);
    }
  }
  int R2 = 1 << r;
  for (int l = 1; l < R2; l++) {
    uint64_t z[4] = {0, 0, 0, 0};
    for (int c = 0; c < R2; c++)
      if (!(__builtin_popcount((unsigned)(l & c)) & 1)) z[c >> 6] |= 1ull << (c & 63);
    int ok = 1;
    for (int p = 0; p < K && ok; p++) {
      int any = 0;
      for (int w = 0; w < 4; w++) any |= (hit[p][w] & z[w]) != 0;
      ok = any;
    }
    if (!ok) continue;
    memset(CB, 0, NW * 8);
    for (uint32_t u = 0; u < n; u++) {
      el x = ELS[u];
      if (!(__builtin_popcount((unsigned)(l & COORD[x])) & 1)) CB[x >> 6] |= 1ull << (x & 63);
    }
    /* generators: those of Phi(M) and lifts of a basis of the hyperplane l^perp */
    el gens[15];
    int ng = 0;
    for (int g = 0; g < npg; g++) gens[ng++] = pg[g];
    int i0 = __builtin_ctz((unsigned)l);
    for (int j = 0; j < r; j++) {
      if (j == i0) continue;
      gens[ng++] = (l >> j & 1) ? mul(basis[j], basis[i0]) : basis[j];
    }
    insert(CB, NW, M.lg - 1, gens, ng);
  }
}

int main(int argc, char **argv) {
  if (argc < 2) { fprintf(stderr, "usage: census16 b_1 ... b_t\n"); return 2; }
  T = argc - 1;
  if (T > MAXF) { fprintf(stderr, "at most %d blocks\n", MAXF); return 2; }
  int tot = 0;
  for (int i = 0; i < T; i++) {
    B[i] = atoi(argv[i + 1]);
    if (B[i] != 2 && B[i] != 4 && B[i] != 8 && B[i] != 16) { fprintf(stderr, "block sizes are 2, 4, 8 or 16\n"); return 2; }
    tot += B[i];
  }
  if (tot > 16) { fprintf(stderr, "at most 16 points\n"); return 2; }
  for (int lg = 1; lg <= 3; lg++) build_small(lg);
  build_ambient();
  NW = (N + 63) / 64;
  MEM = calloc(NW, 8); PHI = calloc(NW, 8); CB = calloc(NW, 8);
  ELS = malloc(N * sizeof(el)); PHIEL = malloc(N * sizeof(el));
  COORD = malloc(N * sizeof(int16_t));

  /* the root: R itself, generated by the wreath generators of each factor (and sigma for R_16) */
  el gens[15];
  int ng = 0;
  for (int i = 0; i < T; i++) {
    int lg = LG[i] < 4 ? LG[i] : 3, b = 1 << lg;
    for (int s = 1; s < b; s *= 2) {
      uint8_t want[8];
      for (int p = 0; p < b; p++) want[p] = (uint8_t)p;
      for (int p = 0; p < s; p++) { want[p] = (uint8_t)(p + s); want[p + s] = (uint8_t)p; }
      int c = -1;
      for (int e = 0; e < SM[lg].n; e++)
        if (memcmp(SM[lg].perm[e], want, (size_t)b) == 0) c = e;
      gens[ng++] = (el)((uint32_t)c << SH[i]); /* for R_16: (g, 1) on the first half */
    }
    if (LG[i] == 4) gens[ng++] = (el)(1u << 14 << SH[i]); /* sigma */
  }
  uint32_t n = closure(gens, ng, MEM, ELS);
  if (n != N) { fprintf(stderr, "the generators give %u elements, not %u\n", n, N); return 2; }
  insert(MEM, NW, NBITS, gens, ng);

  clock_t t0 = clock();
  for (size_t idx = 0; idx < NSUB; idx++) {
    process(idx);
    if ((idx & 0xFFFFF) == 0xFFFFF) fprintf(stderr, "  %zu processed, %zu stored, %.0fs\n", idx + 1, NSUB, (double)(clock() - t0) / CLOCKS_PER_SEC);
  }
  for (int k = 0; k < NKEY; k++) {
    printf("key %d %d %d %d", KEYJ[k], KEYFC[k], KEYEA[k], KEYD[k]);
    for (int i = 0; i < T; i++) printf(" %ld", KEYS[k][i]);
    printf(" %ld\n", KEYCNT[k]);
  }
  printf("listed %zu\n", NSUB);
  for (int j = 0; j < 16; j++)
    if (LAYER[j]) printf("layer %d %ld\n", j, LAYER[j]);
  printf("cpu_seconds %.1f\n", (double)(clock() - t0) / CLOCKS_PER_SEC);
  return 0;
}
