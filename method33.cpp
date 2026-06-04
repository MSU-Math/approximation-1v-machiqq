#include "method33.h"
#include <math.h>
#include <stddef.h>

static double divdiff(double xi, double fi, double xi1, double fi1)
{
    return (fi1 - fi) / (xi1 - xi);
}

void method33_build(int n, const double *x, const double *f,
                    double *a, double *work)
{
    int i;
    double *xe = work;
    double *fe = work + (n + 2);

    for (i = 0; i < n; i++) {
        xe[i + 1] = x[i];
        fe[i + 1] = f[i];
    }

    xe[0] = x[0] - (x[1] - x[0]);
    fe[0] = f[0] - (f[1] - f[0]) / (x[1] - x[0]) * (x[1] - x[0]);
    fe[0] = 2.0 * f[0] - f[1];

    xe[n + 1] = x[n - 1] + (x[n - 1] - x[n - 2]);
    fe[n + 1] = 2.0 * f[n - 1] - f[n - 2];

    double *d = work + 2 * (n + 2);

    double dd_left, dd_right;

    dd_left = divdiff(xe[0], fe[0], xe[1], fe[1]);

    for (i = 0; i < n; i++) {
        dd_right = divdiff(xe[i + 1], fe[i + 1], xe[i + 2], fe[i + 2]);

        if ((dd_left >= 0.0) == (dd_right >= 0.0) &&
            !(dd_left == 0.0 && dd_right == 0.0)) {
            double sign = (dd_left >= 0.0) ? 1.0 : -1.0;
            double abs_left  = fabs(dd_left);
            double abs_right = fabs(dd_right);
            d[i] = sign * (abs_left < abs_right ? abs_left : abs_right);
        } else {
            d[i] = 0.0;
        }

        dd_left = dd_right;
    }

    for (i = 0; i < n - 1; i++) {
        double h = x[i + 1] - x[i];
        double dd_seg = (f[i + 1] - f[i]) / h;
        double c0 = f[i];
        double c1 = d[i];
        double c2 = (3.0 * dd_seg - 2.0 * d[i] - d[i + 1]) / h;
        double c3 = (d[i] + d[i + 1] - 2.0 * dd_seg) / (h * h);
        a[4 * i + 0] = c0;
        a[4 * i + 1] = c1;
        a[4 * i + 2] = c2;
        a[4 * i + 3] = c3;
    }
}

double method33_eval(double t, double a_left, double b_right,
                     int n, const double *x, const double *a)
{
    int lo, hi, mid, seg;
    (void)a_left;
    (void)b_right;

    if (n <= 1) return a[0];

    if (t <= x[0]) {
        seg = 0;
    } else if (t >= x[n - 1]) {
        seg = n - 2;
    } else {
        lo = 0;
        hi = n - 2;
        while (lo < hi) {
            mid = (lo + hi) / 2;
            if (x[mid + 1] <= t)
                lo = mid + 1;
            else
                hi = mid;
        }
        seg = lo;
    }

    double dt = t - x[seg];
    return a[4*seg+0] + dt*(a[4*seg+1] + dt*(a[4*seg+2] + dt*a[4*seg+3]));
}
