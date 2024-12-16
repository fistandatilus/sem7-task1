#include "gas_one.h"

double C_norm(const P_she &p_she, const double *res1, const double *res2, const int scale) {
    double norm = 0;
    const double *cv1 = res1;
    const double *ch1 = res1 + p_she.Dim;
    const double *cv2 = res2;
    const double *ch2 = res2 + p_she.M_x*scale + 1;
    for (int j = 0; j <= p_she.M_x; j++) {
        double s = sqrt((cv1[j] - cv2[j*scale])*(cv1[j] - cv2[j*scale]) + (ch1[j] - ch2[j*scale])*(ch1[j] - ch2[j*scale]));
        if (s > norm) norm = s;
    }
    return norm;
}

double L_norm(const P_she &p_she, const double *res1, const double *res2, const int scale) {
    double norm1 = 0;
    double norm2 = 0;
    const double *cv1 = res1;
    const double *ch1 = res1 + p_she.Dim;
    const double *cv2 = res2;
    const double *ch2 = res2 + p_she.M_x*scale + 1;
    int M = p_she.M_x;
    double h = p_she.h_x;
    for (int i = 1; i < M; i++) {
        norm1 += (cv1[i] - cv2[scale*i])*(cv1[i] - cv2[scale*i]);
        norm2 += (ch1[i] - ch2[scale*i])*(ch1[i] - ch2[scale*i]);
    }
    norm1 += 0.5*(cv1[0] - cv2[scale*0])*(cv1[0] - cv2[scale*0]);
    norm2 += 0.5*(ch1[0] - ch2[scale*0])*(ch1[0] - ch2[scale*0]);
    norm1 += 0.5*(cv1[M] - cv2[scale*M])*(cv1[M] - cv2[scale*M]);
    norm2 += 0.5*(ch1[M] - ch2[scale*M])*(ch1[M] - ch2[scale*M]);
    norm1 = sqrt(h*norm1);
    norm2 = sqrt(h*norm2);
    return norm1 + norm2;
}

double W_norm(const P_she &p_she, const double *res1, const double *res2, const int scale) {
    double norm1 = 0;
    double norm2 = 0;
    double norm1w = 0;
    double norm2w = 0;
    const double *cv1 = res1;
    const double *ch1 = res1 + p_she.Dim;
    const double *cv2 = res2;
    const double *ch2 = res2 + p_she.M_x*scale + 1;
    int M = p_she.M_x;
    double h = p_she.h_x;
    for (int i = 0; i < M; i++) {
        norm1 += (cv1[i] - cv2[scale*i])*(cv1[i] - cv2[scale*i]);
        norm2 += (ch1[i] - ch2[scale*i])*(ch1[i] - ch2[scale*i]);
        norm1w += (cv1[i+1] - cv1[i] - cv2[scale*(i+1)] + cv2[scale*i])*(cv1[i+1] - cv1[i] - cv2[scale*(i+1)] + cv2[scale*i]);
        norm2w += (ch1[i+1] - ch1[i] - ch2[scale*(i+1)] + ch2[scale*i])*(ch1[i+1] - ch1[i] - ch2[scale*(i+1)] + ch2[scale*i]);
    }
    norm1 -= 0.5*(cv1[0] - cv2[scale*0])*(cv1[0] - cv2[scale*0]);
    norm2 -= 0.5*(ch1[0] - ch2[scale*0])*(ch1[0] - ch2[scale*0]);
    norm1 += 0.5*(cv1[M] - cv2[scale*M])*(cv1[M] - cv2[scale*M]);
    norm2 += 0.5*(ch1[M] - ch2[scale*M])*(ch1[M] - ch2[scale*M]);
    norm1 = sqrt(h*norm1 + norm1w/h);
    norm2 = sqrt(h*norm2 + norm2w/h);
    return norm1 + norm2;
}
