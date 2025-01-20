#include <pthread.h>
#include <cstdio>
#include <sys/sysinfo.h>
#include <sys/time.h>
#include <fenv.h>
#include "gas_one.h"

struct arguments {
    P_gas *p_gas;
    P_she *p_she;
    double t_st = -1; //return 
    double stab_norm = -1; //return
    int k;
    int thread;

    void set(P_gas *p_gas, P_she *p_she, int k, int thread) {
    this->p_gas = p_gas;
    this->p_she = p_she;
    this->k = k;
    this->thread = thread;
    }
};

void *thread_main(void *args_ptr);
double get_full_time();

double get_full_time() {
    timeval buf;
    gettimeofday(&buf, 0);
    return buf.tv_sec + buf.tv_usec / 1e6;
}

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

    P_gas *p_gas_ptr = new P_gas;
    P_she *p_she_ptr = new P_she;
    P_gas &p_gas = *p_gas_ptr;
    P_she &p_she = *p_she_ptr;
    p_gas.Segm_X = 1;
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

    pthread_t tid[19];
    arguments arg[19];
    int thread_count = 0;
    for (int i = 1; i <= 10; i++) { 
        arg[thread_count].set(p_gas_ptr, p_she_ptr, i, thread_count);
        if (i > 1) pthread_create(tid+thread_count, 0, thread_main, arg+thread_count);
        thread_count++;
    }
    for (int i = 1; i < 10; i++) { 
        arg[thread_count].set(p_gas_ptr, p_she_ptr, 10 + i*m/10, thread_count);
        pthread_create(tid+thread_count, 0, thread_main, arg+thread_count);
        thread_count++;
    }
    thread_main(arg+0);

    for (int i = 1; i < 19; i++) {
        pthread_join(tid[i], 0);
    }

    FILE *fp = fopen("task3_table_2.tex", "a");
    if (fp != nullptr) 
    { 
        fprintf(fp, "$\\begin{array}{c}\\mu = %.3f\\\\p(\\rho) = %s\\\\\\tau = %f\\\\h = %f\\end{array}$", mu, (mode == 0 ? "\\rho" : mode == 1 ? "10\\rho" : mode == 2 ? "100\\rho" : "\\rho^{1,4}"), tau, p_she.h_x);
        for (int i = 0; i < 19; i++) {
            fprintf(fp, "&$%f$ ", arg[i].t_st);
        }
        fprintf(fp, "\\\\\n\\hline\n");

        fclose(fp);
    }
    
/*   
    p_she.N = n_st;
    solve(p_gas, p_she, n_st, res1, buf, 1);
    time = (clock() - time)/CLOCKS_PER_SEC;
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
    delete p_gas_ptr;
    delete p_she_ptr;
    if (fp_u) fclose(fp_u);
    if (fp_rho) fclose(fp_rho);
    return 0;
}



void *thread_main(void *args_ptr){
    arguments &args = *((arguments *)args_ptr);

    P_gas &p_gas = *args.p_gas;
    P_she &p_she = *args.p_she;

    int m = p_she.M_x;

    double *res1, *buf;
    res1 = new double[2*(m+1)];
    buf = new double[6*(m+1)];
    cpu_set_t cpu;
    CPU_ZERO(&cpu);
    int nproc = get_nprocs();
    int cpu_id = nproc - 1 - args.thread%nproc;
    CPU_SET(cpu_id, &cpu);
    pthread_t tid = pthread_self();
    pthread_setaffinity_np(tid, sizeof(cpu_set_t), &cpu);

    double time = get_full_time();
    int n_st = 0;
    solve(p_gas, p_she, n_st, res1, buf, 0, args.k);
    double stab_norm = stabilization_norm(res1 + m + 1, res1, m);
    time = get_full_time() - time;
    printf("mode = %d, n_st = %d, T_st = %f stab_norm = %le, time = %.2f\n", p_gas.p_mode, n_st, n_st*p_she.tau, stab_norm, time);

    args.t_st = p_she.tau*n_st;
    args.stab_norm = stab_norm;
    return nullptr;
}
