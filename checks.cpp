#include <stdio.h>
#include "gas_one.h"

void check_matrix(int n, const double *a, const double *b, const double *c) {
    int s = 0;
    for (int i = 0; i < n; i++) {
        if (a[i] + b[i] > c[i])
            s++;
    }
    //if (s) printf("bad matrix! s = %d\n", s);

}
