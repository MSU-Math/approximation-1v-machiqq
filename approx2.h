#ifndef APPROX2_H
#define APPROX2_H

void build_piecewise_cubic_m1(int n,
                           const double *nodes,
                           const double *fvals,
                           double *coeffs);

void build_piecewise_cubic(int n,
                           const double *nodes,
                           const double *fvals,
                           double *coeffs);


double eval_piecewise_cubic_m1(int n,
                            const double *nodes,
                            const double *coeffs,
                            double x);


double eval_piecewise_cubic(int n,
                            const double *nodes,
                            const double *coeffs,
                            double x);


#endif
