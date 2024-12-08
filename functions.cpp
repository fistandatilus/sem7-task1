#include "functions.h"

double rho_test(double t, double x) {
    return exp(t)*(cos(3*M_PI*x) + 1.5);
}

double u_test(double t, double x) {
    return cos(2*M_PI*t)*sin(4*M_PI*x);
}

double f_test_0(double t, double x, double mu) {
    return exp(t)*(cos(3*M_PI*x) + 1.5)*(-2*M_PI*sin(2*M_PI*t)*sin(4*M_PI*x) + 4*M_PI*cos(2*M_PI*t)*sin(4*M_PI*x)*cos(2*M_PI*t)*cos(4*M_PI*x)) + (-rho_test(t, x))*3*M_PI*exp(t)*sin(3*M_PI*x) - mu*16*M_PI*M_PI*cos(2*M_PI*t)*sin(4*M_PI*x);
}
double f_test_1(double t, double x, double mu) {
    return exp(t)*(cos(3*M_PI*x) + 1.5)*(-2*M_PI*sin(2*M_PI*t)*sin(4*M_PI*x) + 4*M_PI*cos(2*M_PI*t)*sin(4*M_PI*x)*cos(2*M_PI*t)*cos(4*M_PI*x)) + (-10*(rho_test(t, x)))*3*M_PI*exp(t)*sin(3*M_PI*x) - mu*16*M_PI*M_PI*cos(2*M_PI*t)*sin(4*M_PI*x);
}
double f_test_2(double t, double x, double mu) {
    return exp(t)*(cos(3*M_PI*x) + 1.5)*(-2*M_PI*sin(2*M_PI*t)*sin(4*M_PI*x) + 4*M_PI*cos(2*M_PI*t)*sin(4*M_PI*x)*cos(2*M_PI*t)*cos(4*M_PI*x)) + (-100*(rho_test(t, x)))*3*M_PI*exp(t)*sin(3*M_PI*x) - mu*16*M_PI*M_PI*cos(2*M_PI*t)*sin(4*M_PI*x);
}
double f_test_poc(double t, double x, double mu) {
    return exp(t)*(cos(3*M_PI*x) + 1.5)*(-2*M_PI*sin(2*M_PI*t)*sin(4*M_PI*x) + 4*M_PI*cos(2*M_PI*t)*sin(4*M_PI*x)*cos(2*M_PI*t)*cos(4*M_PI*x)) + (-log(GAMMA)*pow(rho_test(t, x), GAMMA))*3*M_PI*exp(t)*sin(3*M_PI*x) - mu*16*M_PI*M_PI*cos(2*M_PI*t)*sin(4*M_PI*x);
}

double f_0_test(double t, double x) {
    return exp(t)*cos(2*M_PI*t)*M_PI*(4*cos(4*M_PI*x)*(cos(3*M_PI*x) + 1.5)- 3*sin(3*M_PI*x)*sin(4*M_PI*x));
}
