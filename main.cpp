#include <stdio.h>
#include <time.h>
#include <fenv.h>
#include "gas_one.h"

int main(int argc, char *argv[]) {
//    feenableexcept(FE_ALL_EXCEPT ^ FE_INEXACT);
    int n, m;
    double mu;
    if (!(argc == 4 && sscanf(argv[1], "%d", &n) == 1 && sscanf(argv[2], "%d", &m) == 1 && sscanf(argv[3], "%lf", &mu) == 1)) {
        printf("program usage: %s N M mu\n", argv[0]);
        return 1;
    }

    double *res1, *res2, *buf;
    res1 = new double[2*(m+1)];
    res2 = new double[(1u << (4+4))*2*(m+1)];
    buf = new double[6*(1u << (4+4))*(m+1)];

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

   
    for (int i = 0; i < 4; i++) {
        double *cv = res2;
        double *ch = res2 + (m+1);
        for (int j = 0; j <= m; j++) {
            cv[j] = u_test(1, double(j)/m);
            ch[j] = rho_test(1, double(j)/m);
        }

        p_gas.p_mode = i;
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
                p_gas.p_gamma = GAMMA;
                p_gas.f = f_test_poc;
                break;
        };

        p_she.M_x = m;
        p_she.N = n;
        p_she.Dim = m + 1;
        p_she.h_x = 1./m;
        p_she.tau = 1./n;
        double time = clock();
        solve(p_gas, p_she, res1, buf);
        double c_norm_diff = C_norm(p_she, res1, res2);
        double l_norm_diff = L_norm(p_she, res1, res2);
        double w_norm_diff = W_norm(p_she, res1, res2);
        time = (clock() - time)/CLOCKS_PER_SEC;
        printf("mode = %d, c_norm = %le, l_norm = %le, w_norm = %le, time = %.2f\n", i, c_norm_diff, l_norm_diff, w_norm_diff, time);

        for (int j = 1; j <= 4; j++) {
            int scale = (1u << j);
            p_she.M_x = m*scale;
            p_she.N = n*scale;
            p_she.Dim = m*scale + 1;
            p_she.h_x = 1./(m*scale);
            p_she.tau = 1./(n*scale);
            time = clock();
            solve(p_gas, p_she, res2, buf);
            p_she.M_x = m;
            p_she.N = n;
            p_she.Dim = m + 1;
            p_she.h_x = 1./m;
            p_she.tau = 1./n;
            double c_norm = C_norm(p_she, res1, res2, scale);
            double l_norm = L_norm(p_she, res1, res2, scale);
            double w_norm = W_norm(p_she, res1, res2, scale);
            time = (clock() - time)/CLOCKS_PER_SEC;
            printf("mode = %d, k = %d, c_norm = %le, l_norm = %le, w_norm = %le, time = %.2f\n", i, j, c_norm, l_norm, w_norm, time);
        }
    }
    delete[] res1;
    delete[] res2;
    delete[] buf;
    return 0;
}
