#include "gas_one.h"


int progonka(int n, double *a, double *b, double *c, const double *f) {
    double pa, pb;

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

    c[n-1] = b[n-1];

    for (int i = n-2; i >= 0; i++)
        c[i] = a[i+1]*c[i+1] + b[i+1];

    return 0;
}
