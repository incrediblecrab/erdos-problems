/* Cyclic k-AP-free 2-colourings of Z_p and Z_2p from quadratic residues.
 *
 * Rabung's method [Rabung 1979; Rabung-Lotts, EJC 19(2) (2012) #P35]: colour
 * x in Z_p* by its quadratic character, c(x) = 1 iff x is a non-residue.  A
 * k-AP a + jd with d != 0 equals d(y + j) for y = a/d, and c(dz) = c(z) XOR
 * c(d) for z != 0, so every progression is a rescaled window of k consecutive
 * residues.  Hence the colouring has no monochromatic k-AP in Z_p, degenerate
 * or not, iff
 *
 *   (i)  no k consecutive elements of 1..p-1 share a colour, and
 *   (ii) no window of k consecutive residues through 0 is constant on its
 *        non-zero elements.
 *
 * (ii) does not depend on c(0): both multiplier classes occur, so a window whose
 * non-zero part is constant can always be rescaled to match c(0).
 *
 * The "zipped" colouring of Z_2p [Herwig-Heule-van Lambalgen-van Maaren, EJC 14
 * (2007) R6] is, through Z_2p = Z_p x Z_2, simply
 *
 *   z(x) = c(x mod p) XOR [x even]        for x != 0, p,   z(0) != z(p).
 *
 * It is cyclic k-AP-free iff the Z_p colouring is (even d stays in one parity
 * class, where z is c up to a constant), no window avoiding {0,p} is constant
 * (odd d rescales to d = 1), no window through 0 or p is constant off it, and
 * z(0) != z(p) (d = p).
 *
 * These criteria only nominate candidates.  Nothing here is a certificate: the
 * words are re-checked by src/verify_zm, a plain double loop over all (a, d).
 *
 *   rabung search <kmin> <kmax> <pmax>    largest good p and 2p for each k
 *   rabung word <p> [zip]                 print the colouring as a 0/1 string
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static unsigned char *qr_colour(long p) {
    unsigned char *c = malloc(p);
    memset(c, 1, p);
    c[0] = 0;
    for (long x = 1; x <= p / 2; x++) c[(x * x) % p] = 0;
    return c;
}

static unsigned char *zip_colour(const unsigned char *c, long p) {
    long m = 2 * p;
    unsigned char *z = malloc(m);
    for (long x = 0; x < m; x++) z[x] = c[x % p] ^ (unsigned char)(x % 2 == 0);
    z[0] = 0;
    z[p] = 1;
    return z;
}

/* Longest run of equal colours among consecutive residues of Z_m, taken only
 * over windows that avoid the excluded positions (0, and p when zipped). */
static long longest_run(const unsigned char *c, long m, long skip) {
    long best = 0, run = 0;
    unsigned char prev = 2;
    for (long x = 1; x < m; x++) {
        if (x == skip) { run = 0; prev = 2; continue; }
        run = (c[x] == prev) ? run + 1 : 1;
        prev = c[x];
        if (run > best) best = run;
    }
    return best;
}

/* Does some window of k consecutive residues through position q have a
 * constant colour on its other k-1 elements? */
static int through_constant(const unsigned char *c, long m, long q, int k) {
    for (int a = 0; a < k; a++) {                 /* q sits at offset a */
        int first = -1, same = 1;
        for (int j = 0; j < k; j++) {
            if (j == a) continue;
            long x = ((q + j - a) % m + m) % m;
            if (first < 0) first = c[x];
            else if (c[x] != first) { same = 0; break; }
        }
        if (same) return 1;
    }
    return 0;
}

static int good_prime(const unsigned char *c, long p, int k) {
    return longest_run(c, p, -1) < k && !through_constant(c, p, 0, k);
}

static int good_zip(const unsigned char *z, long p, int k) {
    long m = 2 * p;
    return longest_run(z, m, p) < k && !through_constant(z, m, 0, k)
           && !through_constant(z, m, p, k);
}

int main(int argc, char **argv) {
    if (argc >= 3 && strcmp(argv[1], "word") == 0) {
        long p = atol(argv[2]);
        unsigned char *c = qr_colour(p);
        unsigned char *w = c;
        long m = p;
        if (argc >= 4 && strcmp(argv[3], "zip") == 0) { w = zip_colour(c, p); m = 2 * p; }
        for (long x = 0; x < m; x++) putchar('0' + w[x]);
        putchar('\n');
        return 0;
    }
    if (argc != 5 || strcmp(argv[1], "search") != 0) {
        fprintf(stderr, "usage: %s search <kmin> <kmax> <pmax> | word <p> [zip]\n", argv[0]);
        return 2;
    }
    int kmin = atoi(argv[2]), kmax = atoi(argv[3]);
    long pmax = atol(argv[4]);
    if (kmin < 3 || kmax < kmin || kmax > 64) { fprintf(stderr, "bad k range\n"); return 2; }

    unsigned char *composite = calloc(pmax + 1, 1);
    for (long i = 2; i * i <= pmax; i++)
        if (!composite[i]) for (long j = i * i; j <= pmax; j += i) composite[j] = 1;

    long bestp[65] = {0}, bestz[65] = {0}, nprime[65] = {0}, nzip[65] = {0};
    for (long p = 3; p <= pmax; p += 2) {
        if (composite[p] || p <= kmin) continue;
        unsigned char *c = qr_colour(p);
        unsigned char *z = NULL;
        for (int k = kmin; k <= kmax; k++) {
            if (p <= k) break;                    /* a window would meet 0 twice */
            if (!good_prime(c, p, k)) continue;
            bestp[k] = p; nprime[k]++;
            if (!z) z = zip_colour(c, p);
            if (good_zip(z, p, k)) { bestz[k] = p; nzip[k]++; }
        }
        free(c);
        free(z);
    }
    printf("primes 3..%ld; QR colouring of Z_p and zipped colouring of Z_2p\n", pmax);
    printf("  k  largest p  (#good p)   largest 2p  (#good zips)   best m   delta_k <= 1/(2(k-1)m)\n");
    for (int k = kmin; k <= kmax; k++) {
        long m = bestp[k] > 2 * bestz[k] ? bestp[k] : 2 * bestz[k];
        printf(" %2d  %9ld  (%6ld)   %10ld  (%6ld)   %7ld   1/%ld\n", k, bestp[k], nprime[k],
               2 * bestz[k], nzip[k], m, 2L * (k - 1) * m);
    }
    return 0;
}
