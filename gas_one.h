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
    int M_x;
    int N;
    int Dim;
    double h_x;
    double tau;
    double eta;
};


#define SMALL_NUMBER (1e-200)

//метод прогонки
int progonka(int n, double *a, double *b, double *c, const double *f);

//схема 1
void solve(const P_gas &p_gas, const P_she &p_she, double *res, double *buf);

//нормы
double C_norm(const P_she &p_she, const double *res1, const double *res2, const int scale=1);
double L_norm(const P_she &p_she, const double *res1, const double *res2, const int scale=1);
double W_norm(const P_she &p_she, const double *res1, const double *res2, const int scale=1);
