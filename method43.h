#ifndef METHOD43_H
#define METHOD43_H

#ifdef __cplusplus
extern "C" {
#endif

void method43_build(int n, const double *x, const double *f,
                    double *a, double *work);

double method43_eval(double t, double a_left, double b_right,
                     int n, const double *x, const double *a);

#ifdef __cplusplus
}
#endif

#endif
