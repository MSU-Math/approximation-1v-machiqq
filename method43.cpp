#include "method43.h"

#include <cmath>
#include <cstdlib>


static inline int A_IDX(int i, int comp)
{
    return 4 * i + comp;
}

void method43_build(int n, const double *x, const double *f,
                     double *a, double *workspace)
{
    if (n < 3) {
        return;
    }

    double *xi = workspace;
    double *low = workspace + (n + 1);
    double *diag = workspace + 2 * (n + 1);
    double *up = workspace + 3 * (n + 1);
    double *rhs = workspace + 4 * (n + 1);

    double x0 = x[0] - (x[1] - x[0]);
    double xn1 = x[n - 1] + (x[n - 1] - x[n - 2]);

    xi[0] = (x0 + x[0]) / 2.0;
    for (int j = 1; j < n; ++j) {
        xi[j] = (x[j - 1] + x[j]) / 2.0;
    }
    xi[n] = (x[n - 1] + xn1) / 2.0;

    double *u = (double *)malloc(sizeof(double) * (size_t)n);
    double *w = (double *)malloc(sizeof(double) * (size_t)n);
    double *D = (double *)malloc(sizeof(double) * (size_t)n);

    for (int i = 0; i < n; ++i) {
        u[i] = xi[i] - x[i];
        w[i] = xi[i + 1] - x[i];
        D[i] = u[i] * w[i] * (w[i] - u[i]);
    }

    low[0] = 0.0;
    diag[0] = -w[0];
    up[0] = u[0];
    rhs[0] = f[0] * (u[0] - w[0]);

    for (int j = 1; j < n; ++j) {
        int iL = j - 1;
        int iR = j;

        double uL = u[iL], wL = w[iL], DL = D[iL];
        double uR = u[iR], wR = w[iR], DR = D[iR];

        double coef_vjm1 = -(wL * wL) / DL;
        double coefA_vj = uL * (2.0 * wL - uL) / DL;
        double coefB_vj = wR * (wR - 2.0 * uR) / DR;
        double coef_vjp1 = uR * uR / DR;

        double const_left = f[iL] * (wL - uL) / (uL * wL);
        double const_right = -f[iR] * (wR - uR) / (uR * wR);

        low[j] = coef_vjm1;
        diag[j] = coefA_vj - coefB_vj;
        up[j] = -coef_vjp1;
        rhs[j] = const_right - const_left;
    }

    low[n] = -w[n - 1];
    diag[n] = u[n - 1];
    up[n] = 0.0;
    rhs[n] = f[n - 1] * (u[n - 1] - w[n - 1]);

    for (int j = 1; j <= n; ++j) {
        double m = low[j] / diag[j - 1];
        diag[j] -= m * up[j - 1];
        rhs[j] -= m * rhs[j - 1];
    }

    rhs[n] = rhs[n] / diag[n];
    for (int j = n - 1; j >= 0; --j) {
        rhs[j] = (rhs[j] - up[j] * rhs[j + 1]) / diag[j];
    }

    double *v = rhs;

    for (int i = 0; i < n; ++i) {
        double vi = v[i];
        double vip1 = v[i + 1];
        double fi = f[i];

        double c2 = (vi - fi) * w[i] * w[i] - (vip1 - fi) * u[i] * u[i];
        c2 /= D[i];

        double c3 = (vip1 - fi) * u[i] - (vi - fi) * w[i];
        c3 /= D[i];

        a[A_IDX(i, 0)] = x[i];
        a[A_IDX(i, 1)] = fi;
        a[A_IDX(i, 2)] = c2;
        a[A_IDX(i, 3)] = c3;
    }

    free(u);
    free(w);
    free(D);
}

static int find_piece(double t, int n, const double *x)
{
    if (t <= x[0]) {
        return 0;
    }
    if (t >= x[n - 1]) {
        return n - 1;
    }

    int lo = 0;
    int hi = n - 1;
    while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        if (x[mid] <= t) {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    double mid_point = (x[lo] + x[lo + 1]) / 2.0;
    if (t < mid_point) {
        return lo;
    } else {
        return lo + 1;
    }
}

double method43_value(double t, double a_left, double b_right,
                       int n, const double *x, const double *a)
{
    (void)a_left;
    (void)b_right;

    if (n < 3) {
        return 0.0;
    }

    int i = find_piece(t, n, x);

    double xi_center = a[A_IDX(i, 0)];
    double fi = a[A_IDX(i, 1)];
    double c2 = a[A_IDX(i, 2)];
    double c3 = a[A_IDX(i, 3)];

    double dx = t - xi_center;

    return fi + c2 * dx + c3 * dx * dx;
}
