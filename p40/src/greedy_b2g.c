/*
Greedy B2[g] set: a_1 = 1, and a_{k+1} is the least m > a_k such that every n has at most g representations n = a + b with a <= b and a, b in {a_1, ..., a_k, m}. With g = 1 this is the Mian-Chowla sequence (OEIS A005282).
Usage: greedy_b2g G K VMAX OUT. Writes a_1..a_K to OUT, one per line. Exit 0 when K terms were found, 4 when every m <= VMAX was tried first, 2 on bad arguments, 3 when memory is refused. Memory is 2 VMAX + 8 K bytes, fixed at start; time is bounded by VMAX candidates.
src/greedy_check.py checks the output without sharing any code with this file.
*/
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
  if (argc != 5) {
    fprintf(stderr, "usage: %s G K VMAX OUT\n", argv[0]);
    return 2;
  }
  long g = atol(argv[1]), K = atol(argv[2]), V = atol(argv[3]);
  if (g < 1 || g > 200 || K < 1 || V < 1 || V > 2000000000L) {
    fprintf(stderr, "bad arguments\n");
    return 2;
  }
  /* r[n] counts the pairs a <= b in the current set with a + b = n. */
  uint8_t *r = calloc(2 * (size_t)V + 1, 1);
  long *A = malloc(sizeof(long) * (size_t)K);
  if (!r || !A) {
    fprintf(stderr, "out of memory\n");
    return 3;
  }
  long k = 0;
  for (long m = 1; m <= V && k < K; m++) {
    /* Adding m creates the pairs (a, m) for a in the set, with distinct sums a + m, and (m, m) with sum 2m, which differs from every a + m. */
    if (r[2 * m] >= g) continue;
    long i = 0;
    while (i < k && r[A[i] + m] < g) i++;
    if (i < k) continue;
    for (i = 0; i < k; i++) r[A[i] + m]++;
    r[2 * m]++;
    A[k++] = m;
  }
  FILE *f = fopen(argv[4], "w");
  if (!f) {
    fprintf(stderr, "cannot open %s\n", argv[4]);
    return 2;
  }
  for (long i = 0; i < k; i++) fprintf(f, "%ld\n", A[i]);
  fclose(f);
  fprintf(stderr, "g=%ld: %ld terms, last %ld\n", g, k, k ? A[k - 1] : 0);
  return k == K ? 0 : 4;
}
