#include "gas_one.h"

double C_norm(P_she &p_she, double *res1, double *res2) {
    double norm = 0;
    for (int i = 0; i <= p_she.N; i++) {
        double *cv1 = res1 + 2*i*p_she.Dim;
        double *ch1 = res1 + (2*i + 1)*p_she.Dim;
        double *cv2 = res2 + 2*i*p_she.Dim;
        double *ch2 = res2 + (2*i + 1)*p_she.Dim;
        for (int j = 0; j <= p_she.M_x; j++) {
            double s = sqrt((cv1[j] - cv2[j])*(cv1[j] - cv2[j]) + (ch1[j] - ch2[j])*(ch1[j] - ch2[j]));
            if (s > norm) norm = s;
        }
    }
    return norm;
}
