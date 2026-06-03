#include <stdlib.h>
#include <math.h>
#include "approx2.h"

void build_piecewise_cubic(int n, const double *nodes, const double *fvals, double *coeffs)
{
    if (n < 2) return;

    int segs;
    double *h;
    double *lo;
    double *di;
    double *up;
    double *rhs;
    double *m;

    segs = n - 1;

    h = (double *)malloc(segs * sizeof(double));
    lo = (double *)malloc(n * sizeof(double));
    di = (double *)malloc(n * sizeof(double));
    up = (double *)malloc(n * sizeof(double));
    rhs = (double *)malloc(n * sizeof(double));
    m = (double *)malloc(n * sizeof(double));
    for (int i = 0; i < segs; i++)
	h[i] = nodes[i + 1] - nodes[i];

    lo[0] = 0.0;
    di[0] = 1.0;
    up[0] = 0.0;
    rhs[0] = 0.0;

    for (int i = 1; i <= n - 2; i++) {
	lo[i] = h[i - 1];
	di[i] = 2.0 * (h[i - 1] + h[i]);
	up[i] = h[i];
        rhs[i] = 6.0 * (fvals[i + 1] - fvals[i]) / h[i] - (fvals[i] - fvals[i - 1]) / h[i - 1]);
    }

    lo[n - 1] = 0.0;
    di[n - 1] = 1.0;
    up[n - 1] = 0.0;
    rhs[n - 1] = 0.0;

    for (int i = 1; i < n; i++) {
	double factor = lo[i] / di[i - 1];
	di[i] -= factor * up[i - 1];
	rhs[i] -= factor * rhs[i - 1];
    }

    m[n - 1] = rhs[n - 1] / di[n - 1];
    for (int i = n - 2; i >= 0; i--)
	m[i] = (rhs[i] - up[i] * m[i + 1]) / di[i];

    for (int i = 0; i < segs; i++) {
	double hi = h[i];
	double fi = fvals[i];
	double fi1 = fvals[i + 1];
	double mi = m[i];
	double mi1 = m[i + 1];

	coeffs[4 * i + 0] = fi;
	coeffs[4 * i + 1] = (fi1 - fi) / hi - hi * (2.0 * mi + mi1) / 6.0;
	coeffs[4 * i + 2] = mi / 2.0;
	coeffs[4 * i + 3] = (mi1 - mi) / (6.0 * hi);
    }

    free(h);
    free(lo);
    free(di);
    free(up);
    free(rhs);
    free(m);
}


double eval_iecewise_cubic(int n, const double *nodes, const double *coeffs, double x)
{
    int segs;
    int lo;
    int hi;
    int i;
    double t;
    double a;
    double b;
    double c;
    double d;

    segs = n - 1;
    lo = 0;
    hi = segs - 1;

    if (x <= nodes[0])
	lo = hi = 0;
    else if (x >= nodes[n - 1])
	lo = hi = segs - 1;
    else {
	while (lo < hi - 1) {
	    int mid = (lo + hi) / 2;
	    if (x < nodes[mid])
		hi = mid;
	    else
		lo = mid;
	}
    }

    i = lo;
    t = x - nodes[i];
    a = coeffs[4 * i + 0];
    b = coeffs[4 * i + 1];
    c = coeffs[4 * i + 2];
    d = coeffs[4 * i + 3];

    return a + t * (b + t * ( c + t * d));
}





















