#include "method43.h"
#include <math.h>
#include <stddef.h>

void method43_build(int n, const double *x, const double *f,
                    double *a, double *work)
{
    int i;
    double *c = work;

    if (n == 1) {
        a[0] = f[0];
        a[1] = 0.0;
        a[2] = 0.0;
        return;
    }
    if (n == 2) {
        double h   = x[1] - x[0];
        double xi0 = 0.5 * (x[0] + x[1]);
        c[0]     = 0.0;
        a[0]     = 0.5 * (f[0] + f[1]);
        a[1]     = (f[1] - f[0]) / h;     
        a[2]     = 0.0;                    
        (void)xi0;
        return;
    }

    c[0] = 0.0;

    for (i = 1; i <= n - 2; i++) {
        double h_prev  = x[i]     - x[i - 1];
        double h_curr  = x[i + 1] - x[i];
        double b_prev  = (f[i]     - f[i - 1]) / h_prev;
        double b_curr  = (f[i + 1] - f[i])     / h_curr;
        double xi_prev = 0.5 * (x[i - 1] + x[i]);
        double xi_curr = 0.5 * (x[i]     + x[i + 1]);
        double gap     = xi_curr - xi_prev;

        if (i < n - 2) {
            double xi_next = 0.5 * (x[i + 1] + x[i + 2]);
            double span    = xi_next - xi_curr;
            c[i] = (b_curr - b_prev - 2.0 * c[i - 1] * gap)
                   / (2.0 * span);
        } else {
            c[i] = 0.0;
            (void)b_curr;
            (void)gap;
        }
    }

    for (i = 0; i < n - 1; i++) {
        double h    = x[i + 1] - x[i];
        double b_i  = (f[i + 1] - f[i]) / h;
        double a_i  = 0.5 * (f[i] + f[i + 1]) - c[i] * h * h * 0.25;

        a[3 * i + 0] = a_i;
        a[3 * i + 1] = b_i;
        a[3 * i + 2] = c[i];
    }
}

double method43_eval(double t, double a_left, double b_right,
                     int n, const double *x, const double *a)
{
    int lo, hi, mid, seg;
    double xi_seg, dt;
    (void)a_left;
    (void)b_right;

    if (n <= 1) return a[0];
    
    if (t <= x[0]) {
        seg = 0;
    } else if (t >= x[n - 1]) {
        seg = n - 2;
    } else {
        lo = 0; hi = n - 2;
        while (lo < hi) {
            mid = (lo + hi) / 2;
            if (x[mid + 1] <= t) lo = mid + 1;
            else                  hi = mid;
        }
        seg = lo;
    }

    xi_seg = 0.5 * (x[seg] + x[seg + 1]);
    dt     = t - xi_seg;

    return a[3*seg+0] + dt * (a[3*seg+1] + dt * a[3*seg+2]);
}
