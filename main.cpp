#include <stdio.h>
#include "gas_one.h"

int main(int argc, char *argv[]) {
    int n, m;
    double mu;
    if (!(argc == 4 && sscanf(argv[1], "%d", &n) == 1 && sscanf(argv[2], "%d", &m) == 1 && sscanf(argv[3], "%lf", &mu) == 1)) {
        printf("program usage: %s N M mu\n", argv[0]);
        return 1;
    }

    double *res1, *res2, *buf;
    res1 = new double[(n+1)*2*(m+1)];
    res2 = new double[(1u << (4+4))*(n+1)*2*(m+1)];
    buf = new double[5*(1u << (4+4))*(m+1)];

    if (!(res1 && res2 && buf)) {
        printf("Memory error");
        if (res1) delete[] res1;
        if (res2) delete[] res2;
        if (buf) delete[] buf;
        return -2;
    }

    P_gas p_gas;
    P_she p_she;
    p_gas.Segm_T = 1;
    p_gas.Segm_X = 1;
    p_gas.mu = mu;
    p_gas.f_0 = f_0_test;

    for(int i = 0; i <= n; i++) {
        double *cv = res2 + 2*i*(m+1);
        double *ch = res2 + (2*i + 1)*(m+1);
        for (int j = 0; j <= m; j++) {
            cv[i] = u_test(double(i)/n, double(j)/m);
            ch[i] = rho_test(double(i)/n, double(j)/m);
        }

    }

    for (int i = 0; i < 4; i++) {
        p_gas.p_mode = (i < 3);
        switch (i) {
            case 0:
                p_gas.p_ro = 1;
                p_gas.f = f_test_0;
                break;
            case 1:
                p_gas.p_ro = 10;
                p_gas.f = f_test_1;
                break;
            case 2:
                p_gas.p_ro = 100;
                p_gas.f = f_test_2;
                break;
            case 3:
                p_gas.p_gamma = 1.4;
                p_gas.f = f_test_poc;
                break;
        };

        p_she.M_x = m;
        p_she.N = n;
        p_she.Dim = m + 1;
        p_she.h_x = 1./m;
        p_she.tau = 1./n;
        solve(p_gas, p_she, res1, buf);
        double norm_diff = C_norm(p_she, res1, res2);
        printf("mode = %d, norm = %le\n", i, norm_diff);

/*
        for (int j = 1; j <= 4; j++) {
            p_she.M_x = m*(1u << i);
            p_she.N = n*(1u << i);
            p_she.Dim = m*(1u << i) + 1;
            p_she.h_x = 1./(m*(1u << i));
            p_she.tau = 1./(n*(1u << i));
            solve(p_gas, p_she, res2, buf);
            double norm_diff = C_norm(p_she, res1, res2);
            printf("mode = %d, k = %d, norm = %le\n", i, j, norm_diff);
        }
        */
    }
    delete[] res1;
    delete[] res2;
    delete[] buf;
    return 0;
}
