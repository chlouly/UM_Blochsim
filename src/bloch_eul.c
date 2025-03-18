#include <math.h>
#include <stdio.h>
#include "blochsim_3d.h"
#include "blochsim_core.h"


void c_blochsim_eul(//double *Mi,
		long double *beff,
    long double *M,
		long double T1,
		long double T2,
		int nstep,
		long double dt  
) {
    M[X]=0.0;    //Mi[0];
    M[Y]=0.0;    //Mi[1];
    M[Z]=1.0;    //Mi[2];

    for(int lp=1;lp<nstep; lp++) {
        CROSS(M + off(lp - 1), beff + off(lp - 1), M + off(lp), dt * GAMMA);
        M[x(lp)] -= M[x(lp - 1)] * dt / T2;
        M[y(lp)] -= M[y(lp - 1)] * dt / T2;
        M[z(lp)] += (M[z(0)] - M[z(lp - 1)]) * dt / T1;
        ADD_VEC(M + off(lp), M + off(lp - 1));
    }
}
