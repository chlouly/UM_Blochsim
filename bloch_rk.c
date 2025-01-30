#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "blochsim_core.h"

void c_blochsim_rk4(//double *Mi,
		     double *B,
             double *M,
		     double T1,
		     double T2,
		     int nstep,
		     double dt)
{ 

        printf("HELLO\n\n");
    /* dt in ms */
    T1 = dt / T1;
    T2 = (1. - dt / T2);

    M[X]=0.0;    //Mi[0];
    M[Y]=0.0;    //Mi[1];
    M[Z]=1.0;    //Mi[2];

    double h = 2 * dt;

    double k_1[3] = {0.0, 0.0, 0.0};
    double k_2[3] = {0.0, 0.0, 0.0};
    double k_3[3] = {0.0, 0.0, 0.0};
    double k_4[3] = {0.0, 0.0, 0.0};

    for (int n=1; n < nstep; n++) {
        //c_bloch_step(B + IND(n - 1), M + IND(n - 1), NULL, k_1, T1, T2);      // Find k_1
        ADD_SCALED_VEC(M + IND(n), k_1, h / 6);                         // Add k_1

        //c_bloch_step();


        ADD_SCALED_VEC(M + IND(n), k_2, h / 3);
        ADD_SCALED_VEC(M + IND(n), k_3, h / 3);
        ADD_SCALED_VEC(M + IND(n), k_4, h / 6);
        ADD_VEC(M + IND(n), M + IND(n - 1))
    }

}


void c_bloch_step(
    double * B,
    double * M,
    double * k_in,
    double * k_out,
    double T1,
    double T2,
    double scale
) {
      
      double Bx = B[X] * GAMMA; 
      double By = B[Y] * GAMMA; 
      double Bz = B[Z] * GAMMA;

      /*
	Compute sines and cosines of field angles:
	Theta = angle w.r.t. positive z axis
	Phi = angle w.r.t. positive x axis
	Psi = angle w.r.t. transformed positive x axis
      */

      double Bmag = sqrt(Bx*Bx+By*By+Bz*Bz);  /* Magnitude of applied field */
      double Btrans = sqrt(Bx*Bx+By*By);      /* Magnitude of transverse applied field */

      double ct = 1;
      if(Bmag > 0)
	    ct = Bz/Bmag;  /* cos(theta) */

      double st = sqrt(1 - ct*ct);  /* sin(theta) > 0 */

      double cphi = 1;
      if(Btrans > 0)
	    cphi = Bx/Btrans;  /* cos(phi) */
      
      double sphi; 
      if(By < 0.)
	    sphi = sqrt(1 - cphi*cphi)*(-1.);
      else
	    sphi = sqrt(1 - cphi*cphi);

      double cpsi = cos(Bmag);  /* cos(psi) */
      double spsi = sin(Bmag);  /* sin(psi) */
      
      double Mx1, My1, Mz1;
      double Mx0 = 0;
      double My0 = 0;
      double Mz0 = 0;

    if (k_in) {
        Mx0 += k_in[X] * scale;
        My0 += k_in[Y] * scale;
        Mz0 += k_in[Z] * scale;
    }
      
    if(Bmag > 0)
	{
	  Mx0 += M[X];
	  My0 += M[Y];
	  Mz0 += M[Z];
	  
	  Mx1 = cphi*(ct*(cpsi*(ct*(sphi*My0+cphi*Mx0)-st*Mz0) + spsi*(cphi*My0-sphi*Mx0))+st*(ct*Mz0+st*(sphi*My0+cphi*Mx0))) - sphi*(-1.*spsi*(ct*(sphi*My0+cphi*Mx0)-st*Mz0) + cpsi*(cphi*My0-sphi*Mx0));
	  My1 = sphi*(ct*(cpsi*(ct*(sphi*My0+cphi*Mx0)-st*Mz0) + spsi*(cphi*My0-sphi*Mx0))+st*(ct*Mz0+st*(sphi*My0+cphi*Mx0))) + cphi*(-1.*spsi*(ct*(sphi*My0+cphi*Mx0)-st*Mz0) + cpsi*(cphi*My0-sphi*Mx0));
	  Mz1 = ct*(ct*Mz0+st*(sphi*My0+cphi*Mx0)) - st*(cpsi*(ct*(sphi*My0+cphi*Mx0)-st*Mz0) + spsi*(cphi*My0-sphi*Mx0));

	} else {
	  
	  Mx1 = M[X];
	  My1 = M[Y];
	  Mz1 = M[Z];
	  
	}
      
      /* relaxation effects: "1" in Mz since Mo=1 by assumption */
      
      k_out[X] = Mx1*T2;
      k_out[Y] = My1*T2;
      k_out[Z] = Mz1 + (1 - Mz1)*T1;
}