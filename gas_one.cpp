#include "string.h"

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
    auto func = p_gas.f;
    auto otlad_func = p_gas.f_0;
    auto p = p_gas.p;
    double mu_loc;

    a = buf;
    b = buf + M;
    cv = buf + 2*M;
    ch = buf + 3*M;
    f = buf + 4*M;
    
    {
        //один шаг
        //нахождение мю с волной
        mu_loc = mu/h[0];
        for (int i = 1; i <= M; i++) {
            double sr = mu/h[i];
            if (mu_loc < sr) mu_loc = sr;
        }

            //заполнение для V
            double frac = tau*mu_loc/h_x/h_x;
        a[0] = 0;
        double cososim = tau/(6*h_x)*v[1];
        b[0] = cososim - frac;
        a[1] = -cososim - frac;
        cv[0] = 1 + 2*frac;
        f[0] = tau*(-p(h[1])/(2*h_x*h[0]) - (mu_loc - mu/h[0])*(v[1] - 2*v[0])/h_x/h_x + func(0));
        cososim = tau/(6*h_x)*v[M-1];
        a[M] = -cososim - frac;
        b[M-1] = cososim - frac;
        b[M] = 0; 
        cv[M] = 1 + 2*frac;
        f[M] = tau*(p(h[M-1])/(2*h_x*h[M]) - (mu_loc - mu/h[M])*(v[M-1] - 2*v[M])/h_x/h_x + func(p_gas.Segm_X));
        for (int i = 1; i <= M-2; i++) {
            cv[i] = 1 + 2*frac;
            cososim = tau/(6*h_x)*(v[i] + v[i+1]);
            b[i] = cososim - frac;
            a[i+1] = -cososim - frac;
            f[i] = v[i] + tau*(-(p(h[i+1]) - p(h[i-1]))/2/h_x/h[i] - (mu_loc - mu/h[i])*(v[i+1] - 2*v[i] + v[i-1])/h_x/h_x + func(double(i)/M*p_gas.Segm_X));
        }
        f[M-1] = v[M-1] + tau*(-(p(h[M]) - p(h[M-2]))/2/h_x/h[M-1] - (mu_loc - mu/h[M-1])*(v[M] - 2*v[M-1] + v[M-2])/h_x/h_x + func(double(M-1)/M*p_gas.Segm_X)); //не заполнится в циклe
        progonka(M+1, a, b, cv, f);

        //заполнение для H
        ch[0] = 1;
        b[0] = tau/2/h_x*cv[1];
        a[0] = 0;
        a[1] = -tau/4/h_x*(cv[0] + cv[1]);
        f[0] = h[0] - tau/2/h_x*(h[0]*cv[1] - h[2]*v[2] + 2*h[1]*v[1] + 0.5*h[3]*v[3] - h[2]*v[2] + 0.5*h[1]*v[1]
            - h[0]*v[2] + 2*h[0]*v[1] + 0.5*h[0]*v[3] - h[0]*v[2] + 0.5*h[0]*v[1]) + tau*otlad_func(0);
        ch[M] = 1;
        a[M] = tau/2/h_x*v[M-1];
        b[M] = 0;
        b[M-1] = tau/4/h_x*(cv[M-1] + cv[M]);
        f[M] = h[M] - tau/2/h_x*(h[M]*cv[M-1] - h[M-2]*v[M-2] + 2*h[M-1]*v[M-1] + 0.5*h[M-3]*v[M-3] - h[M-2]*v[M-2] + 0.5*h[M-1]*v[M-1] 
            - h[M]*v[M-2] + 2*h[M]*v[M-1] + 0.5*h[M]*v[M-3] - h[M]*v[M-2] + 0.5*h[M]*v[M-1]) + tau*otlad_func(p_gas.Segm_X);
        for (int i = 1; i <= M-2; i++) {
            ch[i] = 1;
            cososim = tau/4/h_x*(cv[i] + cv[i+1]);
            b[i] = cososim;
            a[i+1] = -cososim;
            f[i] = h[i] - tau/8/h_x/h_x*h[i]*(cv[i+1] - cv[i-1]) + tau*otlad_func(double(i)/M*p_gas.Segm_X);
        }
        f[M-1] = h[M-1] - tau/8/h_x/h_x*h[M-1]*(cv[M] - cv[M-2]) + tau*otlad_func(double(M-1)/M*p_gas.Segm_X);
        progonka(M+1, a, b, ch, f);

        memcpy(v, cv, (M+1)*sizeof(double));
        memcpy(h, ch, (M+1)*sizeof(double));
    }

}
