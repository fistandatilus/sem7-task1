#include <cstdio>
#include <time.h>
#include <fenv.h>
#include "gas_one.h"

int main(int argc, char *argv[]) {
    //feenableexcept(FE_ALL_EXCEPT ^ FE_INEXACT);
    int n, m, mode;
    double mu, tau, eps;
    if (!(argc >= 7 && sscanf(argv[1], "%d", &n) == 1 && sscanf(argv[2], "%d", &m) == 1 
        && sscanf(argv[3], "%lf", &tau) == 1 && sscanf(argv[4], "%lf", &mu) == 1 && sscanf(argv[5], "%lf", &eps) == 1 
        && sscanf(argv[6], "%d", &mode) == 1 && mode > 0 && mode <= 4)) {
        printf("program usage: %s N M tau mu eps mode [u_file rho_file]\n", argv[0]);
        return 1;
    }

    mode--;

    FILE *fp_u = nullptr, *fp_rho = nullptr;

    if (argc == 9) {
        fp_u = fopen(argv[7], "w");
        fp_rho = fopen(argv[8], "w");
        if (!fp_u || !fp_rho) {
            printf("cannot open file!\n");
            if (fp_u) fclose(fp_u);
            if (fp_rho) fclose(fp_rho);
        }
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
        if (fp_u) fclose(fp_u);
        if (fp_rho) fclose(fp_rho);
        return -2;
    }

    P_gas p_gas;
    P_she p_she;
    p_gas.Segm_T = 100;
    p_gas.Segm_X = 10;
    p_gas.mu = mu;
    p_gas.f = f;
    p_gas.p_mode = mode;
    switch (mode) {
        case 0:
            p_gas.p_ro = 1;
            break;
        case 1:
            p_gas.p_ro = 10;
            break;
        case 2:
            p_gas.p_ro = 100;
            break;
        case 3:
            p_gas.p_gamma = GAMMA;
            break;
    };

    p_she.M_x = m;
    p_she.N = n;
    p_she.Dim = m + 1;
    p_she.h_x = p_gas.Segm_X/m;
    p_she.tau = tau;
    p_she.eps = eps;
    double time = clock();
    int n_st = 0;
    solve(p_gas, p_she, n_st, res1, buf, 0, fp_u, fp_rho);
    p_she.N = n_st;
    solve(p_gas, p_she, n_st, res1, buf, 1);
    double stab_norm = stabilization_norm(res1 + m + 1, res1, m);
    time = (clock() - time)/CLOCKS_PER_SEC;
    printf("mode = %d, n_st = %d, T_st = %f stab_norm = %le, time = %.2f\n", mode, n_st, n_st*tau, stab_norm, time);
   /*
    int n_st_plug = n_st;
    p_she.N = n_st/10;
    solve(p_gas, p_she, n_st_plug, res1, buf, 0);
    P_she p_she_other(p_she);

    for (int j = 1; j <= 4; j++) {
        int scale = (1u << j);
        p_she_other.M_x = m*scale;
        p_she_other.N = n_st*scale/10;
        p_she_other.Dim = m*scale + 1;
        p_she_other.h_x = p_gas.Segm_X/(m*scale);
        p_she_other.tau = tau/scale;
        time = clock();
        solve(p_gas, p_she_other, n_st_plug, res2, buf, 0);
        double c_norm = C_norm(p_she, res1, res2, scale);
        double l_norm = L_norm(p_she, res1, res2, scale);
        double w_norm = W_norm(p_she, res1, res2, scale);
        time = (clock() - time)/CLOCKS_PER_SEC;
        printf("mode = %d, k = %d, c_norm = %le, l_norm = %le, w_norm = %le, time = %.2f\n", mode, j, c_norm, l_norm, w_norm, time);
        printf("$\\|v-v^{%d}\\|_{C_h} = %e$\n", j, c_norm);
    }
*/
    delete[] res1;
    delete[] res2;
    delete[] buf;
    if (fp_u) fclose(fp_u);
    if (fp_rho) fclose(fp_rho);
    return 0;
}
