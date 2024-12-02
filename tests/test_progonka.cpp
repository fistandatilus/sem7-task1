#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "../gas_one.h"

int main(int argc, char **argv) {
    srand(clock());
    size_t n = 0;
    if (!(argc == 2 && sscanf(argv[1], "%lu", &n) == 1 && n >= 2)) {
        printf("Usage: %s n\nn -size of matrix\n", argv[0]);
        return -1;
    }

    double *a = new double[n];
    double *b = new double[n];
    double *c = new double[n];
    double *f = new double[n];
    c[0] = 1;
    a[0] = 0;
    b[0] = double(rand())/RAND_MAX;
    f[0] = b[0];
    c[n-1] = 1;
    a[n-1] = double(rand())/RAND_MAX;
    b[n-1] = 0;
    f[n-1] = ((n - 1) & 1u) + (a[n-1] + b[n-1]) * (1 - ((n - 1) & 1u));
    for (size_t i = 1; i < n-1; i++) {
        c[i] = 1;
        a[i] = .5*double(rand())/RAND_MAX;
        b[i] = .5*double(rand())/RAND_MAX;
        f[i] = (i & 1u) + (a[i] + b[i]) * (1 - (i & 1u));
    }

    double t = clock();
    progonka(n, a, b, c, f);
    t = (clock() - t)/CLOCKS_PER_SEC;
    double res = 0;
    for (size_t i = 0; i < n; i++) {
        printf("%e\n", c[i]);
        res += fabs(c[i] - double(i&1u));
    }

    printf("time = %.2f, residual = %e\n", t, res);

    delete[] a;
    delete[] b;
    delete[] c;
    delete[] f;

    return 0;
}
