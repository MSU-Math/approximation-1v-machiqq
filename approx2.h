#ifndef APPROX2_H
#define APPROX2_H

void build_piecewise_cubic(int n,
                           const double *x,
                           const double *f,
                           double *a);

double eval_piecewise_cubic(int n,
                            const double *x,
                            const double *a,
                            double xval);

#endif
