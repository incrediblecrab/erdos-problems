/*
 * circlef.c -- floating-point Phi_k of a circle colouring, for SEARCH ONLY.
 *
 * A line-by-line port of circle.phi_circle (src/circle.py) to C doubles, so
 * that src/circleopt.py can afford many starts.  Nothing it returns is ever
 * recorded as a result: every word circleopt.py stores is re-derived by
 * exact.phi in rational arithmetic, and circleopt.py checks this port against
 * the Python evaluator on random inputs before using it.
 *
 * bounds: 0 = a_0 < a_1 < ... < a_nb = 1, nb even; arc i = [a_i, a_{i+1})
 * has colour i mod 2.  Phi_k = area of {(x,d) in [0,1)^2 : every x+jd mod 1,
 * j < k, lies in an arc of one colour}, summed over both colours.
 *
 * Build:  cc -O3 -shared -fPIC -o src/libcirclef.so src/circlef.c
 */

#include <stdlib.h>

#define MAXV 64

typedef struct { int n; double x[MAXV], y[MAXV]; } poly;

/* Intersect with the half-plane n0*x + n1*y <= c (blocks.clip). */
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

static int K, NEXT;
static const double *EXT;
static double TOTAL;

/* first index g with ext[g] > v (bisect_right) */
static int bisect_right(double v) {
    int lo = 0, hi = NEXT;
    while (lo < hi) {
        int mid = (lo + hi) / 2;
        if (v < EXT[mid]) hi = mid; else lo = mid + 1;
    }
    return lo;
}

static void rec(int j, const poly *p, int parity) {
    if (j == K) { TOTAL += area(p); return; }
    double lo = p->x[0] + j * p->y[0], hi = lo;
    for (int i = 1; i < p->n; i++) {
        double v = p->x[i] + j * p->y[i];
        if (v < lo) lo = v;
        if (v > hi) hi = v;
    }
    int first = bisect_right(lo) - 1;
    if (first < 0) first = 0;
    int last = bisect_right(hi);
    if (last > NEXT - 1) last = NEXT - 1;
    poly a, b;
    for (int g = first; g < last; g++) {
        if ((g & 1) != parity) continue;
        double left = EXT[g], right = EXT[g + 1];
        if (right <= lo || left >= hi) continue;
        const poly *q = p;
        if (left > lo) {
            clip(q, -1.0, -(double)j, -left, &a);
            if (a.n < 3) continue;
            q = &a;
        }
        if (right < hi) {
            clip(q, 1.0, (double)j, right, &b);
            if (b.n < 3) continue;
            q = &b;
        }
        rec(j + 1, q, parity);
    }
}

double phi_circle(const double *bounds, int nb, int k) {
    if (nb < 2 || (nb & 1)) return -1.0;
    double *ext = malloc(sizeof(double) * (nb * k + 1));
    int n = 0;
    for (int w = 0; w < k; w++)
        for (int i = 0; i < nb; i++) ext[n++] = bounds[i] + w;
    ext[n++] = (double)k;
    K = k; NEXT = n; EXT = ext; TOTAL = 0.0;
    poly sq = {4, {0, 1, 1, 0}, {0, 0, 1, 1}}, p1, p2;
    for (int i0 = 0; i0 < nb; i0++) {
        const poly *p = &sq;
        if (bounds[i0] > 0) { clip(p, -1.0, 0.0, -bounds[i0], &p1); p = &p1; }
        if (bounds[i0 + 1] < 1) { clip(p, 1.0, 0.0, bounds[i0 + 1], &p2); p = &p2; }
        if (p->n >= 3) rec(1, p, i0 & 1);
    }
    free(ext);
    return TOTAL;
}

/* Phi of positive arc lengths L[0..nb-1] (normalised to sum 1). */
double phi_lengths(const double *L, int nb, int k) {
    double *b = malloc(sizeof(double) * (nb + 1));
    double s = 0;
    for (int i = 0; i < nb; i++) s += (L[i] < 0 ? -L[i] : L[i]);
    b[0] = 0;
    double acc = 0;
    for (int i = 0; i < nb; i++) { acc += (L[i] < 0 ? -L[i] : L[i]) / s; b[i + 1] = acc; }
    b[nb] = 1.0;
    double v = phi_circle(b, nb, k);
    free(b);
    return v;
}
