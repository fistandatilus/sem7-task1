#include <math.h>

//Структуры из пособия

//Параметры дифференциальной задачи
struct P_gas {
    double Segm_T;
    double Segm_X;
    double p_ro;
    double p_gamma;
    double mu;
    double (*f)(double, double, double);   //правая часть
    double (*f_0)(double, double); //дополнительная правая часть для отладочного теста
    double (*p)(double); //зависимость
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
