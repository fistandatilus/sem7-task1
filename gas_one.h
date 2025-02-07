#include <cstdio>
#include <math.h>

#include "functions.h"

//Структуры из пособия

//Параметры дифференциальной задачи
struct P_gas {
    double Segm_T;
    double Segm_X;
    double p_ro;
    double p_gamma;
    int p_mode;
    double mu;
    //граничные условия в 4 задаче
    double u_left;
    double rho_left;

    double (*f)(double, double, double);   //правая часть
    double (*f_0)(double, double); //дополнительная правая часть для отладочного теста
    int k;

    P_gas() = default;
    P_gas& operator=(P_gas &a) {
    Segm_T = a.Segm_T;
    Segm_X = a.Segm_X;
    p_ro = a.p_ro;
    p_gamma = a.p_gamma;
    p_mode = a.p_mode;
    mu = a.mu;
    u_left = a.u_left;
    rho_left = a.rho_left;
    f = a.f;
    f_0 = a.f_0;
    k = a.k;
    return *this;
    }
};

//Параметры схемы
struct P_she {
    int M_x = 0;
    //Для задачи 1 N - число разбиений, для отсальных - максимум итераций
    int N = 0;
    int Dim = 0;
    int stab_const = 1;
    double h_x = 0;
    double tau = 0;
    double eps = 0;

    P_she() = default;
    P_she(P_she &a) {
        M_x = a.M_x;
        N = a.N;
        Dim = a.Dim;
        stab_const = a.stab_const;
        h_x = a.h_x;
        tau = a.tau;
        eps = a.eps;
    }
};


#define SMALL_NUMBER (1e-200)

//метод прогонки
int progonka(int n, double *a, double *b, double *c, const double *f);
void check_matrix(int n, const double *a, const double *b, const double *c);

//схема
void solve(const P_gas &p_gas, const P_she &p_she, int &n, double *res, double *buf, double &stab_norm, int print, FILE *fp_u=nullptr, FILE *fp_rho=nullptr);

//нормы
double C_norm(const P_she &p_she, const double *res1, const double *res2, const int scale=1);
double L_norm(const P_she &p_she, const double *res1, const double *res2, const int scale=1);
double W_norm(const P_she &p_she, const double *res1, const double *res2, const int scale=1);

//норма для стабилизации
double stabilization_norm(const double *h, const double *v, int M);
