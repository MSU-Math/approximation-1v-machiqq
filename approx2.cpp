#include "approx2.h"
#include <stdlib.h>
#include <string.h>

static void divided_diff4(const double *z, const double *v, double *dd)
{
    double d[4];
    int i, j;
    for (i = 0; i < 4; i++) d[i] = v[i];
    dd[0] = d[0];
    for (j = 1; j < 4; j++) {
        for (i = 3; i >= j; i--) {
            d[i] = (d[i] - d[i-1]) / (z[i] - z[i-j]);
        }
        dd[j] = d[j];
    }
}

static void newton_to_power(const double *z, const double *dd,
                            double shift, double *out)
{
    double poly[4] = {dd[3], 0.0, 0.0, 0.0};
    double alpha;
    int j, k;
    double a1;
    double a2;

    alpha = shift - z[2];
    for (k = 3; k >= 1; k--)
        poly[k] = poly[k-1] + alpha * poly[k];
    poly[0] = alpha * poly[0];
    {
        double tmp[4] = {dd[3], 0.0, 0.0, 0.0};
        double a0 = shift - z[2];
        for (k = 3; k >= 1; k--)
            tmp[k] = tmp[k-1] + a0 * tmp[k];
        tmp[0] = a0 * tmp[0];
        tmp[0] += dd[2];
        a1 = shift - z[1];
        for (k = 3; k >= 1; k--)
            tmp[k] = tmp[k-1] + a1 * tmp[k];
        tmp[0] = a1 * tmp[0];
        tmp[0] += dd[1];
        a2 = shift - z[0];
        for (k = 3; k >= 1; k--)
            tmp[k] = tmp[k-1] + a2 * tmp[k];
        tmp[0] = a2 * tmp[0];
        tmp[0] += dd[0];

        for (j = 0; j < 4; j++) out[j] = tmp[j];
    }
    (void)poly;
}

void build_piecewise_cubic_m1(int n,
                           const double *x,
                           const double *f,
                           double *a)
{
    int i;
    double z[4], v[4], dd[4];

    if (n < 4) {
        memset(a, 0, sizeof(double) * 4 * (n > 1 ? n-1 : 1));
        return;
    }

    for (i = 0; i < n - 1; i++) {
        int i0;

        if (i == 0) {
            i0 = 0;
        } else if (i == n - 2) {
            i0 = n - 4;
        } else {
            i0 = i - 1;
        }

        z[0] = x[i0];     v[0] = f[i0];
        z[1] = x[i0+1];   v[1] = f[i0+1];
        z[2] = x[i0+2];   v[2] = f[i0+2];
        z[3] = x[i0+3];   v[3] = f[i0+3];

        divided_diff4(z, v, dd);

        newton_to_power(z, dd, x[i], &a[4*i]);
    }
}

double eval_piecewise_cubic_m1(int n,
                            const double *x,
                            const double *a,
                            double xval)
{
    int i;
    int lo;
    int hi;
    int mid;
    double t, val;

    if (xval <= x[0]) {
        i = 0;
    } else if (xval >= x[n-1]) {
        i = n - 2;
    } else {
        lo = 0;
        hi = n - 2;
        while (lo < hi - 1) {
            mid = (lo + hi) / 2;
            if (xval < x[mid]) hi = mid;
            else               lo = mid;
        }
        i = lo;
    }

    t = xval - x[i];

    val = a[4*i+3];
    val = val * t + a[4*i+2];
    val = val * t + a[4*i+1];
    val = val * t + a[4*i+0];

    return val;
}

