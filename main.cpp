#include <pthread.h>
#include <sys/sysinfo.h>
#include <stdio.h>
#include <time.h>
#include <fenv.h>
#include "gas_one.h"
#include "table_print.h"

#define MAX_THREADS 8

struct arguments {
    double *data;
    double mu;
    int thread;
    int mode;

    void set(double *data, double mu, int thread, int mode) {
        this->data = data;
        this->mu = mu;
        this->thread = thread;
        this->mode = mode;
    }
};

void *thread_main(void *args_ptr);

int main(void) {
    int table_size = 3*4*5;
    double *data = new double[12*table_size];
    double mu = 0.1;
    int thread_count = 0;
    pthread_t tid[12];
    arguments arg[12];
    for (int mu_count = 0; mu_count < 3; mu_count++, mu *= 0.1) {
        for(int mode = 0; mode < 4; mode++) {
            arg[thread_count].set(data + table_size*thread_count, mu, thread_count, mode);
            pthread_create(tid+thread_count, 0, thread_main, arg+thread_count);
            thread_count++;
        }
    }
    FILE *fp = fopen("nested_tabels.tex", "w");
    for (int i = 0; i < thread_count; i++) {
        pthread_join(tid[i], 0);
    }
    mu = 0.1;
    for (int mu_count = 0; mu_count < 3; mu_count++, mu *= 0.1) {
        for(int mode = 0; mode < 4; mode++) {
            int min_n, min_m, max_n, max_m, line_length;
            if ((mu > 0.09 && mode <= 1 ) || (mu > 0.009 && mu < 0.02 && mode == 0)) {
                min_n = min_m = 10;
                max_n = max_m = 10000;
                line_length = 4;
            }
            else if (mu > 0.09 && mode == 3) {
                min_n = min_m = 100;
                max_n = max_m = 10000;
                line_length = 3;
            }
            else if (mu > 0.009 && mu < 0.02 && mode == 3) {
                min_n = min_m = 1000;
                max_n = max_m = 10000;
                line_length = 2;
            }
            else if ((mu > 0.009 && mu < 0.02 && mode == 2) || (mu < 0.002 && mode == 2)) {
                min_n = 100000;  min_m = 100;
                max_n = 1000000; max_m = 1000;
                line_length = 2;
            }
            else {
                min_n = 10000;  min_m = 100;
                max_n = 100000; max_m = 1000;
                line_length = 2;
            }
            char title[1234];
            if (mode < 3){
                int c = 1;
                for(int i = 0; i < mode; i++, c*=10);
                sprintf(title, "$\\mu = %.*f, p(\\rho) = %d\\rho$", mu_count+1, mu, c);
            }
            else {
                sprintf(title, "$\\mu = %.*f, p(\\rho) = \\rho^{1.4}$", mu_count+1, mu);
            }
            char (headers_or[4])[1234];
            char *(headers[4]);
            for (int i = 0; i < 4; i++) headers[i] = &*(headers_or[i]);

            headers[0][0] = 0;
            for (int n = min_n, m = min_m, i = 1; n <= max_n; n*=10, m*=10, i++) {
                sprintf(headers[i], "$\\tau = %f, h = %f$", 1./n, 1./m);
            }
            print_header(title, headers, line_length+1, fp);
            for (int i = 0; i < 4; i++) {
                fprintf(fp, "$v-v^{%d}$", i+1);
                print_row(data + table_size*(mu_count*4 + mode) + i*3*line_length, line_length+1, 3, fp);
            }
                fprintf(fp, "$v-u$"); print_row(data + table_size*(mu_count*4 + mode) + 4*3*line_length, line_length+1, 3, fp);
            fprintf(fp, "\\end{tabular}\n");
        }
    }
    
    delete[] data;
    return 0;
}

void *thread_main(void *args_ptr) {
    //feenableexcept(FE_ALL_EXCEPT ^ FE_INEXACT);
    arguments &args = *((arguments *)args_ptr);

    int thread = args.thread;
    cpu_set_t cpu;
    CPU_ZERO(&cpu);
    int nproc = get_nprocs();
    int cpu_id = nproc - 1 - thread%nproc;
    CPU_SET(cpu_id, &cpu);
    pthread_t tid = pthread_self();
    pthread_setaffinity_np(tid, sizeof(cpu_set_t), &cpu);

    double *data = args.data;

    int min_n, min_m, max_n, max_m, i = args.mode, mode = args.mode;
    int line_length;
    double mu = args.mu;
    if ((mu > 0.09 && mode <= 1 ) || (mu > 0.009 && mu < 0.02 && mode == 0)) {
        min_n = min_m = 10;
        max_n = max_m = 10000;
        line_length = 4;
    }
    else if (mu > 0.09 && mode == 3) {
        min_n = min_m = 100;
        max_n = max_m = 10000;
        line_length = 3;
    }
    else if (mu > 0.009 && mu < 0.02 && mode == 3) {
        min_n = min_m = 1000;
        max_n = max_m = 10000;
        line_length = 2;
    }
    else if ((mu > 0.009 && mu < 0.02 && mode == 2) || (mu < 0.002 && mode == 2)) {
        min_n = 100000;  min_m = 100;
        max_n = 1000000; max_m = 1000;
        line_length = 2;
    }
    else {
        min_n = 10000;  min_m = 100;
        max_n = 100000; max_m = 1000;
        line_length = 2;
    }
    (void)max_n;

    double *res1, *res2, *buf;
    res1 = new double[2*(max_m+1)];
    res2 = new double[(1u << (4+4))*2*(max_m+1)];
    buf = new double[6*(1u << (4+4))*(max_m+1)];

    if (!(res1 && res2 && buf)) {
        printf("Memory error");
        if (res1) delete[] res1;
        if (res2) delete[] res2;
        if (buf) delete[] buf;
        return nullptr;
    }

    P_gas p_gas;
    P_she p_she;
    p_gas.Segm_T = 1;
    p_gas.Segm_X = 1;
    p_gas.mu = mu;
    p_gas.f_0 = f_0_test;
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

    for (int m = min_m, n = min_n, mc = 0; mc < line_length; m *= 10, mc++, n *= 10) {
        double *cv = res2;
        double *ch = res2 + (m+1);
        for (int j = 0; j <= m; j++) {
            cv[j] = u_test(1, double(j)/m);
            ch[j] = rho_test(1, double(j)/m);
        }

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
        data[line_length*4*3 + mc*3] = c_norm_diff;
        data[line_length*4*3 + mc*3 + 1] = l_norm_diff;
        data[line_length*4*3 + mc*3 + 2] = w_norm_diff;
        printf("mu = %e, mode = %d, c_norm = %le, l_norm = %le, w_norm = %le, time = %.2f\n", mu, i, c_norm_diff, l_norm_diff, w_norm_diff, time);

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
            data[line_length*(j-1)*3 + mc*3] = c_norm;
            data[line_length*(j-1)*3 + mc*3 + 1] = l_norm;
            data[line_length*(j-1)*3 + mc*3 + 2] = w_norm;
            time = (clock() - time)/CLOCKS_PER_SEC;
            printf("mu = %e, mode = %d, k = %d, c_norm = %le, l_norm = %le, w_norm = %le, time = %.2f\n", mu, i, j, c_norm, l_norm, w_norm, time);
        }
    }

    delete[] res1;
    delete[] res2;
    delete[] buf;
    return nullptr;
}
