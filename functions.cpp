#include "functions.h"
/*
double rho_test(double t, double x) {
    return exp(t)*(cos(3*M_PI*x) + 1.5);
}

double u_test(double t, double x) {
    return cos(2*M_PI*t)*sin(4*M_PI*x);
}

double f_test_0(double t, double x, double mu) {
    return -2*M_PI*sin(2*M_PI*t)*sin(4*M_PI*x) + u_test(t, x)*cos(2*M_PI*t)*4*M_PI*cos(4*M_PI*x) - 1 * 3*M_PI*sin(3*M_PI*x)/(cos(3*M_PI*x) + 1.5) - mu/rho_test(t, x)*(-4*4*M_PI*M_PI*u_test(t, x));
}
double f_test_1(double t, double x, double mu) {
    return -2*M_PI*sin(2*M_PI*t)*sin(4*M_PI*x) + u_test(t, x)*cos(2*M_PI*t)*4*M_PI*cos(4*M_PI*x) - 10 * 3*M_PI*sin(3*M_PI*x)/(cos(3*M_PI*x) + 1.5) - mu/rho_test(t, x)*(-4*4*M_PI*M_PI*u_test(t, x));
}
double f_test_2(double t, double x, double mu) {
    return -2*M_PI*sin(2*M_PI*t)*sin(4*M_PI*x) + u_test(t, x)*cos(2*M_PI*t)*4*M_PI*cos(4*M_PI*x) - 100 * 3*M_PI*sin(3*M_PI*x)/(cos(3*M_PI*x) + 1.5) - mu/rho_test(t, x)*(-4*4*M_PI*M_PI*u_test(t, x));
}
double f_test_poc(double t, double x, double mu) {
    return -2*M_PI*sin(2*M_PI*t)*sin(4*M_PI*x) + u_test(t, x)*cos(2*M_PI*t)*4*M_PI*cos(4*M_PI*x) - GAMMA * exp(t)*3*M_PI*sin(3*M_PI*x) * pow(rho_test(t, x), GAMMA - 2) - mu/rho_test(t, x)*(-4*4*M_PI*M_PI*u_test(t, x));
}

double f_0_test(double t, double x) {
    return rho_test(t,x) + (-u_test(t,x)*exp(t)*3*M_PI*sin(3*M_PI*x) + rho_test(t, x)*cos(2*M_PI*t)*4*M_PI*cos(4*M_PI*x));
}
*/

double rho_1(double /*t*/, double x) {
    return (x < 4.5 || x > 5.5) ? 1 : 2;
}

double u_1(double /*t*/, double /*x*/) {
    return 0;
}

double f(double /*t*/, double /*x*/, double /*mu*/) {
    return 0;
}
/*
double f_test_1(double t, double x, double mu) {
    return 10/(x+1);//10/(x+t+1);//((x+1)*(x*x-x) + (t+1)*(t+1)*(x*x*x-x)*(2*x-1) + 10 - mu*2)/x+1;
}
double f_test_2(double t, double x, double mu) {
    return 100/(x+1);//((x+1)*(x*x-x) + (t+1)*(t+1)*(x*x*x-x)*(2*x-1) + 100 - mu*2)/(x+1);
}
double f_test_poc(double t, double x, double mu) {
    return GAMMA*pow(x+1, GAMMA-2);
}

double f_0_test(double t, double x) {
    return 0;
}
*/
