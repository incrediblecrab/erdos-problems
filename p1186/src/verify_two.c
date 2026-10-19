/* Independent brute-force check of the two-scale colouring bound.
 *
 * Shares no code with the search or with exact.py.  It builds the colouring of
 * {0,...,n-1} straight from the definition
 *
 *     colour(i) = C[ (i mod m) * B + floor(B*i/n) ]
 *
 * and counts, by direct enumeration, the monochromatic k-term arithmetic
 * progressions a, a+d, ..., a+(k-1)d with d >= 0 and a+(k-1)d <= n-1.
 * Psi_k is the limit of that count divided by n^2.
 *
 * usage:  verify_two <k> <m> <B> <colour bits, length m*B> <n> [n2 ...]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv) {
    if (argc < 6) {
        fprintf(stderr, "usage: %s <k> <m> <B> <bits> <n> [n2 ...]\n", argv[0]);
        return 2;
    }
    int k = atoi(argv[1]);
    int m = atoi(argv[2]);
    int B = atoi(argv[3]);
    const char *bits = argv[4];
    if ((int)strlen(bits) != m * B) {
        fprintf(stderr, "bit string has length %d, expected m*B = %d\n",
                (int)strlen(bits), m * B);
        return 2;
    }
    unsigned char *C = malloc(m * B);
    for (int i = 0; i < m * B; i++) {
        if (bits[i] != '0' && bits[i] != '1') {
            fprintf(stderr, "non-binary character at position %d\n", i);
            return 2;
        }
        C[i] = (unsigned char)(bits[i] - '0');
    }

    for (int arg = 5; arg < argc; arg++) {
        long n = atol(argv[arg]);
        unsigned char *col = malloc(n);
        for (long i = 0; i < n; i++) {
            long q = (long)((( __int128 )B * i) / n);      /* floor(B*i/n) */
            if (q >= B) q = B - 1;                          /* cannot trigger for i<n */
            col[i] = C[(i % m) * B + q];
        }
        unsigned long long mono = 0;
        long dmax = (n - 1) / (k - 1);
        for (long d = 0; d <= dmax; d++) {
            long amax = n - 1 - (k - 1) * d;
            for (long a = 0; a <= amax; a++) {
                unsigned char c0 = col[a];
                int j;
                for (j = 1; j < k; j++)
                    if (col[a + (long)j * d] != c0) break;
                if (j == k) mono++;
            }
        }
        free(col);
        double psi = (double)mono / ((double)n * (double)n);
        printf("k=%d m=%d B=%d n=%ld  mono=%llu  count/n^2 = %.12f\n",
               k, m, B, n, mono, psi);
        fflush(stdout);
    }
    free(C);
    return 0;
}
