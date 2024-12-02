#include "gas_one.h"


int progonka(int n, double *a, double *b, double *c, const double *f) {
    double pa = a[0], pb = b[0];
        
    for (int i = 1; i < n; i++){
        double ca, cb, denominator;
        denominator = c[i-1] + pa*a[i-1];
        if (fabs(denominator) < SMALL_NUMBER)
            return -1;
        ca = -pb/denominator;
        cb = (f[i-1] - pa*b[i-1])/denominator;
        pa = a[i];
        pb = b[i];
        a[i] = ca;
        b[i] = cb;
    }

    double denominator = c[n-1] + pa*a[n-1];
    if (fabs(denominator) < SMALL_NUMBER)
        return -1;
    c[n-1] = (f[n-1] - pa*b[n-1])/denominator;

    for (int i = n-2; i >= 0; i--)
        c[i] = a[i+1]*c[i+1] + b[i+1];

    return 0;
}


void solve(const P_gas &p_gas, const P_she &p_she, double *v, double *h, double *buf) {
    double *a, *b, *ch, *cv, *f;
    double tau = p_she.tau;
    double h_x = p_she.h_x;
    double mu = p_gas.mu;
    int M = p_she.M_x;
    //один шаг
    //заполнение для V
    a = buf;
    b = buf + M;
    cv = buf + 2*M;
    ch = buf + 3*M;
    f = buf + 4*M;

    double mu_loc = p_gas.mu;
    double frac = tau*mu_loc/h_x/h_x;

    a[0] = 0; 
    b[0] = tau/(6*h_x)*v[1] - frac;
    cv[0] = 1 + frac;
    f[0] = tau*(p(h[1])/(2*h_x*h[0]) - (mu_loc - mu/h[0])*(v[1] - 2*v[0])/h_x/h_x + func(0))
    b[M] = tau/(6*h_x)*v[1] - frac;
    a[M] = 0; 
    cv[M] = 1 + frac;
    f[M] = tau*(p(h[1])/(2*h_x*h[0]) - (mu_loc - mu/h[0])*(v[1] - 2*v[0])/h_x/h_x + func(0))
    

}
