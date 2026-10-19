/* Independent brute-force count of Phi_k over Z_m.
 *
 * Shares no code with exact.py: it is a direct double loop over all (a,b) in
 * Z_m^2 with no weight table, no set deduplication and no rational arithmetic.
 *
 *     Phi_k(c) = (1/m^2) #{(a,b) : c(a) = c(a+b) = ... = c(a+(k-1)b)}
 *
 * Degenerate progressions are counted, including b = 0, which is Lu and Peng's
 * m_k(Z_n) convention: "Here we allow k-APs to be degenerated".
 *
 * usage:  verify_zm <k> <bits>          bits has length m
 * prints the monochromatic count, m^2, and a breakdown by b for small m.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static long gcd_l(long a, long b) { while (b) { long t = a % b; a = b; b = t; } return a; }

int main(int argc, char **argv) {
    if (argc != 3) { fprintf(stderr, "usage: %s <k> <bits>\n", argv[0]); return 2; }
    int k = atoi(argv[1]);
    const char *bits = argv[2];
    long m = (long)strlen(bits);
    if (k < 3 || m < 1) { fprintf(stderr, "bad k or empty bit string\n"); return 2; }
    unsigned char *c = malloc(m);
    for (long i = 0; i < m; i++) {
        if (bits[i] != '0' && bits[i] != '1') {
            fprintf(stderr, "non-binary character at position %ld\n", i); return 2;
        }
        c[i] = (unsigned char)(bits[i] - '0');
    }

    unsigned long long mono = 0, mono_nonzero = 0;
    for (long b = 0; b < m; b++) {
        unsigned long long here = 0;
        for (long a = 0; a < m; a++) {
            unsigned char c0 = c[a];
            long x = a;
            int j;
            for (j = 1; j < k; j++) {
                x += b; if (x >= m) x -= m;
                if (c[x] != c0) break;
            }
            if (j == k) here++;
        }
        mono += here;
        if (b != 0) mono_nonzero += here;
    }
    /* reduce count/m^2 to lowest terms */
    unsigned long long den = (unsigned long long)m * (unsigned long long)m;
    long g = gcd_l((long)mono, (long)den);
    printf("k=%d m=%ld  mono=%llu  m^2=%llu  Phi = %llu/%llu = %.15f\n",
           k, m, mono, den, mono / (unsigned long long)g, den / (unsigned long long)g,
           (double)mono / (double)den);
    printf("  b=0 contributes %ld, b!=0 contributes %llu\n", m, mono_nonzero);
    free(c);
    return 0;
}
