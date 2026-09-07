#include "functions.h"

#include <cmath>

static const int FUNCTIONS_COUNT = 7;

int f_count(void)
{
    return FUNCTIONS_COUNT;
}

double f_value(int k, double x)
{
    switch (k) {
    case 0:
        return 1.0;
    case 1:
        return x;
    case 2:
        return x * x;
    case 3:
        return x * x * x;
    case 4:
        return x * x * x * x;
    case 5:
        return std::exp(x);
    case 6:
        return 1.0 / (25.0 * x * x + 1.0);
    default:
        return 0.0;
    }
}

double f_deriv1(int k, double x)
{
    switch (k) {
    case 0:
        return 0.0;
    case 1:
        return 1.0;
    case 2:
        return 2.0 * x;
    case 3:
        return 3.0 * x * x;
    case 4:
        return 4.0 * x * x * x;
    case 5:
        return std::exp(x);
    case 6: {
        double g = 25.0 * x * x + 1.0;
        return (-50.0 * x) / (g * g);
    }
    default:
        return 0.0;
    }
}

double f_deriv2(int k, double x)
{
    switch (k) {
    case 0:
        return 0.0;
    case 1:
        return 0.0;
    case 2:
        return 2.0;
    case 3:
        return 6.0 * x;
    case 4:
        return 12.0 * x * x;
    case 5:
        return std::exp(x);
    case 6: {
        double g = 25.0 * x * x + 1.0;
        return (3750.0 * x * x - 50.0) / (g * g * g);
    }
    default:
        return 0.0;
    }
}

const char *f_name(int k)
{
    switch (k) {
    case 0:
        return "f(x)=1";
    case 1:
        return "f(x)=x";
    case 2:
        return "f(x)=x^2";
    case 3:
        return "f(x)=x^3";
    case 4:
        return "f(x)=x^4";
    case 5:
        return "f(x)=e^x";
    case 6:
        return "f(x)=1/(25x^2+1)";
    default:
        return "f(x)=?";
    }
}
