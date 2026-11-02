/*
 * twof.c -- floating-point Psi_k of a two-scale colouring of {1,...,n} with
 * free coarse breakpoints, for SEARCH ONLY (src/twofree.py).
 *
 * Colour i in {1..n} by C[i mod m][q], q the coarse block of i/n, where the
 * coarse blocks are [a_q, a_{q+1}), 0 = a_0 < ... < a_nb = 1, freely placed.
 * As in exact.two_scale_weights_int, the residues (alpha, beta) of a
 * progression's start and difference are asymptotically uniform on Z_m^2 and
 * independent of its scaled position (x, t) in T = {x, t >= 0,
 * x + (k-1) t <= 1}, so
 *
 *   Psi_k = (1/m^2) sum_{alpha,beta} area{(x,t) in T : the colours
 *            C[alpha + j beta][block(x + j t)], j < k, are all equal}.
 *
 * m = 1 is the block family of src/blocks.py.  The recursion clips T by the
 * strips a_q <= x + j t < a_{q+1}, as blocks.psi_blocks does.  Nothing here is
 * recorded as a result; a word that mattered would have to be re-derived
 * exactly (exact.psi_two_scale on an equal-block grid) and by src/verify_two.
 *
 * Build:  cc -O3 -shared -fPIC -o src/libtwof.so src/twof.c
 */

#include <stddef.h>

#define MAXV 64

typedef struct { int n; double x[MAXV], y[MAXV]; } poly;

static void clip(const poly *in, double n0, double n1, double c, poly *out) {
    out->n = 0;
    if (in->n == 0) return;
    int np = in->n;
    double px = in->x[np - 1], py = in->y[np - 1];
    double dp = n0 * px + n1 * py - c;
    for (int i = 0; i < np; i++) {
        double cx = in->x[i], cy = in->y[i];
        double dc = n0 * cx + n1 * cy - c;
        if (dc <= 0) {
            if (dp > 0) {
                double s = dp / (dp - dc);
                out->x[out->n] = px + s * (cx - px);
                out->y[out->n] = py + s * (cy - py);
                out->n++;
            }
            out->x[out->n] = cx; out->y[out->n] = cy; out->n++;
        } else if (dp <= 0) {
            double s = dp / (dp - dc);
            out->x[out->n] = px + s * (cx - px);
            out->y[out->n] = py + s * (cy - py);
            out->n++;
        }
        px = cx; py = cy; dp = dc;
    }
}

static double area(const poly *p) {
    if (p->n < 3) return 0.0;
    double s = 0.0;
    double px = p->x[p->n - 1], py = p->y[p->n - 1];
    for (int i = 0; i < p->n; i++) {
        s += px * p->y[i] - p->x[i] * py;
        px = p->x[i]; py = p->y[i];
    }
    return s >= 0 ? s / 2 : -s / 2;
}

static int K, NB;
static const double *A;
static const unsigned char *ROW[16];     /* colour row of residue alpha + j beta */
static double TOTAL;

static void rec(int j, const poly *p, int colour) {
    if (j == K) { TOTAL += area(p); return; }
    double lo = p->x[0] + j * p->y[0], hi = lo;
    for (int i = 1; i < p->n; i++) {
        double v = p->x[i] + j * p->y[i];
        if (v < lo) lo = v;
        if (v > hi) hi = v;
    }
    poly a, b;
    for (int q = 0; q < NB; q++) {
        if (ROW[j][q] != colour) continue;
        double left = A[q], right = A[q + 1];
        if (right <= lo || left >= hi) continue;
        const poly *r = p;
        if (left > lo) { clip(r, -1.0, -(double)j, -left, &a); if (a.n < 3) continue; r = &a; }
        if (right < hi) { clip(r, 1.0, (double)j, right, &b); if (b.n < 3) continue; r = &b; }
        rec(j + 1, r, colour);
    }
}

double psi_two(const unsigned char *C, int m, int nb, const double *bounds, int k) {
    K = k; NB = nb; A = bounds;
    double tot = 0.0;
    poly tri = {3, {0, 1, 0}, {0, 0, 1.0 / (k - 1)}}, p1, p2;
    for (int al = 0; al < m; al++) {
        for (int be = 0; be < m; be++) {
            for (int j = 0; j < k; j++) ROW[j] = C + (size_t)((al + (long)j * be) % m) * nb;
            TOTAL = 0.0;
            for (int q = 0; q < nb; q++) {
                const poly *p = &tri;
                if (bounds[q] > 0) { clip(p, -1.0, 0.0, -bounds[q], &p1); p = &p1; }
                if (bounds[q + 1] < 1) { clip(p, 1.0, 0.0, bounds[q + 1], &p2); p = &p2; }
                if (p->n >= 3) rec(1, p, ROW[0][q]);
            }
            tot += TOTAL;
        }
    }
    return tot / ((double)m * m);
}
