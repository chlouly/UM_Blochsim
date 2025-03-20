#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "blochsim_core.h"
#include "blochsim_3d.h"


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
    for (n = 0; n < niter - 2; n += 2) { // += 2 because we simulate t+2 and interpolate t+1
        //printf("==== %d / %d ====\n", n, niter);
        FILL(M_tmp, M + off(n))

        c_bloch_step(B + off(n), M + off(n), M_tmp, k, T1, T2, h, (double)(1./6.), (double)(0.5));
        c_bloch_step(B + off(n + 1), M + off(n), M_tmp, k, T1, T2, h, (double)(1./3.), (double)(0.5));
        c_bloch_step(B + off(n + 1), M + off(n), M_tmp, k, T1, T2, h, (double)(1./3.), (double)(1.0));
        c_bloch_step(B + off(n + 2), M + off(n), M_tmp, k, T1, T2, h, (double)(1./6.), (double)(0.0));

        // Add the previous Magnetization values
        ADD_VEC(M + off(n + 2), M + off(n));

        // Interpollate to find the M value between the current and next
        MIDPT(M + off(n), M + off(n + 2), M + off(n + 1));
    }

    // Handling final case
    // This means we missed the very last entry of M
    if (n == (niter - 2)) {
        EXTRAP(M + off(niter - 3), M + off(niter - 2), M + off(niter - 1));
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
