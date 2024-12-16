#include "functions.h"
/*
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
*/

double rho_test(double t, double x) {
    return x+t+1;
}

double u_test(double t, double x) {
    return 0;
}

double f_test_0(double t, double x, double mu) {
    return 1/(x+t+1);//M_PI*sin(M_PI*x)*cos(M_PI*x) + M_PI*M_PI*sin(M_PI*x)*mu;
}
double f_test_1(double t, double x, double mu) {
    return 10/(x+t+1);//((x+1)*(x*x-x) + (t+1)*(t+1)*(x*x*x-x)*(2*x-1) + 10 - mu*2)/x+1;
}
double f_test_2(double t, double x, double mu) {
    return 100/(x+t+1);//((x+1)*(x*x-x) + (t+1)*(t+1)*(x*x*x-x)*(2*x-1) + 100 - mu*2)/(x+1);
}
double f_test_poc(double t, double x, double mu) {
    return GAMMA*pow(x+t+1, GAMMA-2);
}

double f_0_test(double t, double x) {
    return 1;
}
