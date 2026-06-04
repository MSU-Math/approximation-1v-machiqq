#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#ifdef __cplusplus
extern "C" {
#endif

double func(int k, double x);

double func_deriv(int k, double x);

const char* func_name(int k);

#define FUNC_COUNT 7

#ifdef __cplusplus
}
#endif

#endif
