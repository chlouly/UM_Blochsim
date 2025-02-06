#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "blochsim_core.h"

#define RKADD_REP(M, k, M_out, scale) \
    (M_out)[X] = (M)[X] + ((k)[X] * (scale)); \
    (M_out)[Y] = (M)[Y] + ((k)[Y] * (scale)); \
    (M_out)[Z] = (M)[Z] + ((k)[Z] * (scale)); \

#define RKADD(k, M_out, scale) \
    (M_out)[X] += (k)[X] * (scale); \
    (M_out)[Y] += (k)[Y] * (scale); \
    (M_out)[Z] += (k)[Z] * (scale); \

void c_blochsim_rk4(//double *Mi,
		double *B,
    double *M,
		double T1,
		double T2,
    int nstep,
		double dt
)   { 
    /* dt in ms (B timestep, not M. M timestep is h) */
    double h = 2 * dt;
    printf("h:%lf\n", h);
    int niter = nstep - (nstep % 2);

    printf("T1: %lf   T2: %lf\n\n", T1, T2);

    // T1 = dt / T1;
    // T2 = (1. - dt / T2);

    M[X]=0.0;    //Mi[0];
    M[Y]=0.0;    //Mi[1];
    M[Z]=1.0;    //Mi[2];

    // temp vectors for the RK method
    double k[3];
    double M_tmp[3];

    // Loop initialization
    //FILL(M_tmp, M)

    PRINT_VEC(M_tmp);

    int n;
    for (n = 0; n < niter; n += 2) { // += 2 because we simulate t+2 and interpolate t+1
        FILL(M_tmp, M + off(n))
        //printf("\n%d",n);
        //PRINT_VEC(B + off(n));
        //printf("m_tmp");
        //PRINT_VEC(M + off(n));
        c_bloch_step(B + off(n), M + off(n), M_tmp, k, T1, T2, h, (double)(1./6.), (double)(0.5));
        c_bloch_step(B + off(n + 1), M + off(n), M_tmp, k, T1, T2, h, (double)(1./3.), (double)(0.5));
        c_bloch_step(B + off(n + 1), M + off(n), M_tmp, k, T1, T2, h, (double)(1./3.), (double)(1.0));
        c_bloch_step(B + off(n + 2), M + off(n), M_tmp, k, T1, T2, h, (double)(1./6.), (double)(0.0));

        // Add the previous Magnetization values
        ADD_VEC(M + off(n + 2), M + off(n));

        // Interpollate to find the M value between the current and next
        MIDPT(M + off(n), M + off(n + 2), M + off(n + 1));
        //if (n > 20) break;
    }
}

void c_bloch_step(
    double * B,
    double * M,
    double * M_tmp,
    double * k,
    double T1,
    double T2,
    double h,
    double scale,
    double scale_tmp
)   {
    // Bloch Calculation
    CROSS(M_tmp, B, k, GAMMA);

    // Decay
    k[X] -= M_tmp[X] / T2;
    k[Y] -= M_tmp[Y] / T2;
    k[Z] += (1. - M_tmp[Z]) / T1;

    //PRINT_VEC(k);

    // Add to M_tmp
    RKADD_REP(M, k, M_tmp, h * scale_tmp);

    // Add to the output M
    RKADD(k, M + off(2), h * scale);
}
