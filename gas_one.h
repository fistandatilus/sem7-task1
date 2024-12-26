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
    double (*f)(double, double, double);   //правая часть
    double (*f_0)(double, double); //дополнительная правая часть для отладочного теста
};

//Параметры схемы
struct P_she {
    int M_x = 0;
    //Для задачи 1 N - число разбиений, для отсальных - максимум итерация и слой стабилизации
    int N = 0;
    int Dim = 0;
    double h_x = 0;
    double tau = 0;
    double eps = 0;

    P_she() = default;
    P_she(P_she &a) {
        M_x = a.M_x;
        N = a.N;
        Dim = a.Dim;
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
void solve(const P_gas &p_gas, const P_she &p_she, int &n, double *res, double *buf, int print, FILE *fp_u=nullptr, FILE *fp_rho=nullptr);

//нормы
double C_norm(const P_she &p_she, const double *res1, const double *res2, const int scale=1);
double L_norm(const P_she &p_she, const double *res1, const double *res2, const int scale=1);
double W_norm(const P_she &p_she, const double *res1, const double *res2, const int scale=1);

//норма для стабилизации
double stabilization_norm(const double *h, const double *v, int M);
