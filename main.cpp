#include <stdio.h>

int main(int argc, char *argv[]) {
    int p_mode;
    double tau, h, mu, p_coeff;
    if (!(argc == 6 && sscanf(argv[1], "%lf", &tau) == 1 && sscanf(argv[2], "%lf", &h) == 1 && sscanf(argv[3], "%lf", &mu) == 1 && sscanf(argv[1], "%d", &p_mode) == 1 && sscanf(argv[1], "%lf", &p_coeff) == 1 )) {
        printf("program usage: %s tau h mu p_mode p_coeff\np_mode 0 -- linear correlation\n       1 -- nonlinear\n", argv[0]);
        return 1;
    }
    return 0;
}
