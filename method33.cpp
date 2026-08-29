#include "method33.h"

#include <cmath>
#include <cstdlib>

static inline int A_IDX(int i, int comp)
{
    return 4 * i + comp;
}

void method33_build(int n, const double *x, const double *f,
                     double *a, double *d_workspace)
{
    if (n < 3) {
        return;
    }

    double *X = (double *)malloc(sizeof(double) * (size_t)(n + 2));
    double *F = (double *)malloc(sizeof(double) * (size_t)(n + 2));
    double *fd = (double *)malloc(sizeof(double) * (size_t)(n + 1));

    for (int j = 0; j < n; ++j) {
        X[j + 1] = x[j];
        F[j + 1] = f[j];
    }

    X[0] = x[0] - (x[1] - x[0]);
    F[0] = 2.0 * f[0] - f[1];

    X[n + 1] = x[n - 1] + (x[n - 1] - x[n - 2]);
    F[n + 1] = 2.0 * f[n - 1] - f[n - 2];

    for (int j = 0; j <= n; ++j) {
        fd[j] = (F[j + 1] - F[j]) / (X[j + 1] - X[j]);
    }

    double *d = d_workspace;
    for (int i = 1; i <= n; ++i) {
        double left = fd[i - 1];
        double right = fd[i];

        double sign_left = (left > 0.0) - (left < 0.0);
        double sign_right = (right > 0.0) - (right < 0.0);

        if (sign_left == sign_right && sign_left != 0.0) {
            double abs_left = std::fabs(left);
            double abs_right = std::fabs(right);
            double m = (abs_left < abs_right) ? abs_left : abs_right;
            d[i - 1] = sign_right * m;
        } else {
            d[i - 1] = 0.0;
        }
    }

    for (int i = 0; i < n - 1; ++i) {
        double h = x[i + 1] - x[i];
        double fdiv = (f[i + 1] - f[i]) / h; /* f[x_i;x_{i+1}] */

        double a1 = f[i];
        double a2 = d[i];
        double a3 = (fdiv - d[i]) / h;
        double a4 = (d[i] + d[i + 1] - 2.0 * fdiv) / (h * h);

        a[A_IDX(i, 0)] = a1;
        a[A_IDX(i, 1)] = a2;
        a[A_IDX(i, 2)] = a3;
        a[A_IDX(i, 3)] = a4;
    }

    free(X);
    free(F);
    free(fd);
}

static int find_piece(double t, int n, const double *x)
{
    if (t <= x[0]) {
        return 0;
    }
    if (t >= x[n - 1]) {
        return n - 2;
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
    return lo;
}

double method33_value(double t, double a_left, double b_right,
                       int n, const double *x, const double *a)
{
    (void)a_left;
    (void)b_right;

    if (n < 3) {
        return 0.0;
    }

    int i = find_piece(t, n, x);

    double a1 = a[A_IDX(i, 0)];
    double a2 = a[A_IDX(i, 1)];
    double a3 = a[A_IDX(i, 2)];
    double a4 = a[A_IDX(i, 3)];

    double dx = t - x[i];
    double dxn = t - x[i + 1];

    return a1 + a2 * dx + a3 * dx * dx + a4 * dx * dx * dxn;
}
