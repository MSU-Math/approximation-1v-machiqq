#ifndef FUNCTIONS_H
#define FUNCTIONS_H


#ifdef __cplusplus
extern "C" {
#endif

int f_count(void);

double f_value(int k, double x);

double f_deriv1(int k, double x);

double f_deriv2(int k, double x);

const char *f_name(int k);

#ifdef __cplusplus
}
#endif

#endif
