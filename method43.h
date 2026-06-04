#ifndef METHOD43_H
#define METHOD43_H


void build_piecewise_cubic(int n,
                           const double *nodes,
                           const double *fvals,
                           double *coeffs);

double eval_piecewise_cubic(int n,
                            const double *nodes,
                            const double *coeffs,
                            double x);
#endif
