#include "method43.h"
#include <math.h>
#include <stddef.h>


void method43_build(int n, const double *x, const double *f,
                    double *a, double *work)
{
    int i;
    double *b = work;
    double *c = work + n;
    double *r = work + 2 * n;
    double *d = work + 3 * n;

    if (n == 1) {
        a[0] = f[0];
        a[1] = 0.0;
        a[2] = 0.0;
        a[3] = 0.0;
        return;
    }

    if (n == 2) {
        double h = x[1] - x[0];
        double dd = (f[1] - f[0]) / h;
        a[0] = f[0];
        a[1] = dd;
        a[2] = 0.0;
        a[3] = 0.0;
        return;
    }

    {
        double h0 = x[1] - x[0];
        double dd0 = (f[1] - f[0]) / h0;
        b[0] = 2.0;
        c[0] = 1.0;
        r[0] = 3.0 * dd0;
    }

    for (i = 1; i <= n - 2; i++) {
        double h_left  = x[i]     - x[i - 1];
        double h_right = x[i + 1] - x[i];
        double dd_left  = (f[i]     - f[i - 1]) / h_left;
        double dd_right = (f[i + 1] - f[i])     / h_right;
        d[i] = h_right;
        b[i] = 2.0 * (h_left + h_right);
        c[i] = h_left;
        r[i] = 3.0 * (dd_left * h_right + dd_right * h_left);
    }

    {
        double h_last = x[n - 1] - x[n - 2];
        double dd_last = (f[n - 1] - f[n - 2]) / h_last;
        d[n - 1] = 1.0;
        b[n - 1] = 2.0;
        c[n - 1] = 0.0;
        r[n - 1] = 3.0 * dd_last;
    }

    for (i = 1; i < n; i++) {
        double l_i = d[i];
        double factor = l_i / b[i - 1];
        b[i] -= factor * c[i - 1];
        r[i] -= factor * r[i - 1];
        d[i] = 0.0;
    }

    d[n - 1] = r[n - 1] / b[n - 1];
    for (i = n - 2; i >= 0; i--) {
        d[i] = (r[i] - c[i] * d[i + 1]) / b[i];
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

double method43_eval(double t, double a_left, double b_right,
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
