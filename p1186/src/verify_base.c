/* Independent brute-force check of a wildcard-recursion base (NOTES.md s.3).
 *
 * Shares no code with cosetprobe.py, safesat.py or final_check.py.  The base
 * is a 0/1 word on Z_b with wildcard coset F = { x : x = r mod g }, g = b/f.
 * For every (a, d) with 0 <= a < b and 0 < d < b it walks a + jd (j < k) and
 * counts the progressions whose terms outside F are all one colour:
 *
 *   h1  no term in F                   (condition H1 forbids these)
 *   h2  terms both inside and outside  (condition H2 forbids these)
 *
 * Progressions lying entirely in F are ignored, as are the colours on F.
 *
 * usage:  verify_base <k> <f> <r> <bits>      exit 0 iff h1 = h2 = 0
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv) {
    if (argc != 5) {
        fprintf(stderr, "usage: %s <k> <f> <r> <bits>\n", argv[0]);
        return 2;
    }
    long k = atol(argv[1]), f = atol(argv[2]), r = atol(argv[3]);
    const char *bits = argv[4];
    long b = (long)strlen(bits);
    if (k < 2 || f < 1 || f >= b || b % f) {
        fprintf(stderr, "need k >= 2 and f a proper divisor of b = %ld\n", b);
        return 2;
    }
    long g = b / f;
    r = ((r % g) + g) % g;
    for (long i = 0; i < b; i++)
        if (bits[i] != '0' && bits[i] != '1') {
            fprintf(stderr, "bits must be 0/1\n");
            return 2;
        }
    long h1 = 0, h2 = 0;
    for (long a = 0; a < b; a++) {
        for (long d = 1; d < b; d++) {
            long x = a, in = 0, out = 0;
            int first = -1, constant = 1;
            for (long j = 0; j < k; j++) {
                if (x % g == r) {
                    in++;
                } else {
                    int c = bits[x] - '0';
                    if (first < 0) first = c;
                    else if (c != first) constant = 0;
                    out++;
                }
                x += d;
                if (x >= b) x -= b;
            }
            if (out && constant) {
                if (in) h2++;
                else h1++;
            }
        }
    }
    printf("k=%ld b=%ld f=%ld r=%ld  h1=%ld h2=%ld\n", k, b, f, r, h1, h2);
    return (h1 || h2) ? 1 : 0;
}
