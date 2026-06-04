#ifndef METHOD33_H
#define METHOD33_H

#ifdef __cplusplus
extern "C" {
#endif

void method33_build(int n, const double *x, const double *f,
                    double *a, double *work);


double method33_eval(double t, double a_left, double b_right,
                     int n, const double *x, const double *a);

#ifdef __cplusplus
}
#endif

#endif
