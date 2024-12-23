#include <string.h>
#include <stdio.h>

#include "functions.h"
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


void solve(const P_gas &p_gas, const P_she &p_she, double *res, double *buf) {
    double *a, *b, *ch, *cv, *f, *v, *h, *ph;
    double tau = p_she.tau;
    double h_x = p_she.h_x;
    double mu = p_gas.mu;
    int M = p_she.M_x;
    int N = p_she.N;
    int dim = p_she.Dim;
    auto func = p_gas.f;
    double mu_loc;
    double *res2 = new double[2*(M+1)];

    v = res;
    h = v + dim;

    a = buf;
    b = buf + M + 1;
    f = buf + 2*(M + 1);
    ph = buf + 3*(M + 1);
    cv = buf + 4*(M + 1);
    ch = buf + 5*(M + 1);

    //начальные данные
    //пока захардкожу
    for (int i = 0; i <= M; i++) {
        v[i] = u_1(0, i*h_x);
        h[i] = rho_1(0, i*h_x);
    }
    
    for(int n = 1; n <= N; n++) {
        //один шаг
        //нахождение мю с волной
        mu_loc = mu/h[0];
        for (int i = 1; i <= M; i++) {
            double sr = mu/h[i];
            if (mu_loc < sr) mu_loc = sr;
        }
        //вычисление p(h)
        if (p_gas.p_mode < 3) {
            for (int i = 0; i <= M; i++) {
                ph[i] = p_gas.p_ro*h[i];
            }
        }
        else {
            for (int i = 0; i <= M; i++) {
                ph[i] = pow(h[i], p_gas.p_gamma);
            }
        }

        v[0] = 0;
        v[M] = 0;

        //заполнение для V
        double frac = tau*mu_loc/h_x/h_x;
        double cososim;
        for (int i = 1; i <= M-1; i++) {
            cv[i] = 1 + 2*frac;
            a[i] = -tau/6/h_x*(v[i-1] + v[i]) - frac;
            b[i] =  tau/6/h_x*(v[i] + v[i+1]) - frac;
            f[i] = v[i] + tau*(-(ph[i+1] - ph[i-1])/2/h_x/h[i] - (mu_loc - mu/h[i])*(v[i+1] - 2*v[i] + v[i-1])/h_x/h_x + func(tau*(n-1), i*h_x, mu));
        }
        cv[0] = 1;
        a[0] = 0;
        b[0] = 0;
        f[0] = 0;
        cv[M] = 1;
        a[M] = 0;
        b[M] = 0;
        f[M] = 0;
        progonka(M+1, a, b, cv, f);
        cv[0] = 0;
        cv[M] = 0;
        //for (int i = 0; i <= M; i++) printf("cv = %le, f = %le %d\n", cv[i], f[i], i);

        //заполнение для H
        ch[0] = 1 - (tau/2/h_x)*cv[0];
        b[0] = tau/2/h_x*cv[1];
        a[0] = 0;
        a[1] = -tau/4/h_x*(cv[0] + cv[1]);
        f[0] = h[0] - tau/2/h_x*(h[0]*cv[1] - h[0]*cv[0] - (2*h[0]*v[0] - 2.5*h[1]*v[1] + 2*h[2]*v[2] - 0.5*h[3]*v[3] - 2.5*h[0]*v[1] + 2*h[0]*v[2] - 0.5*h[0]*v[3]));
        for (int i = 1; i <= M-1; i++) {
            ch[i] = 1;
            cososim = tau/4/h_x*(cv[i] + cv[i+1]);
            b[i] = cososim;
            a[i+1] = -cososim;
            f[i] = h[i] - tau/4/h_x*h[i]*(cv[i+1] - cv[i-1]);
        }
        ch[M] = 1 + (tau/2/h_x)*cv[M];
        a[M] = -tau/2/h_x*v[M-1];
        b[M] = 0;
        f[M] = h[M] - tau/2/h_x*(-h[M]*cv[M-1] + h[M]*cv[M] + 2*h[M]*v[M] - 2.5*h[M-1]*v[M-1] + 2*h[M-2]*v[M-2] - 0.5*h[M-3]*v[M-3] - 2.5*h[M]*v[M-1] + 2*h[M]*v[M-2] - 0.5*h[M]*v[M-3]);
        progonka(M+1, a, b, ch, f);
        //printf("h[0] = %le, h[M] = %le\n", ch[0], ch[M]);
        /*
        double *cv2 = res2;
        double *ch2 = res2 + (M+1);
        for (int j = 0; j <= M; j++) {
            cv2[j] = u_1(n*tau, double(j)/M);
            ch2[j] = rho_1(n*tau, double(j)/M);
        }

        double c_norm = C_norm(p_she, cv, res2);
        double l_norm = L_norm(p_she, cv, res2);
        double w_norm = W_norm(p_she, cv, res2);
        printf("step = %d, c_norm = %le, l_norm = %le, w_norm = %le\n", n, c_norm, l_norm, w_norm);
        */
        memcpy(v, cv, (M+1)*sizeof(double));
        memcpy(h, ch, (M+1)*sizeof(double));
    }
    delete[] res2;
}
