#include <math.h>
#include "blochsim_core.h"


void c_blochsim_eul(//double *Mi,
		     double *beff,
             double *M,
		     double T1,
		     double T2,
		     int nstep,
		     double dt)
{
  /* dt in ms */
  T1 = dt / T1;
  T2 = (1. - dt / T2);

  M[X]=0.0;    //Mi[0];
  M[Y]=0.0;    //Mi[1];
  M[Z]=1.0;    //Mi[2];

  for(int lp=1;lp<nstep; lp++)
    {   
      double Bx = beff[x(lp-1)] * dt * GAMMA; 
      double By = beff[y(lp-1)] * dt * GAMMA; 
      double Bz = beff[z(lp-1)] * dt * GAMMA;

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
      
      double Mx1,My1,Mz1;
      
    if(Bmag > 0)
	{
	  double Mx0 = M[x(lp-1)];
	  double My0 = M[y(lp-1)];
	  double Mz0 = M[z(lp-1)];
	  
	  Mx1 = cphi*(ct*(cpsi*(ct*(sphi*My0+cphi*Mx0)-st*Mz0) + spsi*(cphi*My0-sphi*Mx0))+st*(ct*Mz0+st*(sphi*My0+cphi*Mx0))) - sphi*(-1.*spsi*(ct*(sphi*My0+cphi*Mx0)-st*Mz0) + cpsi*(cphi*My0-sphi*Mx0));
	  My1 = sphi*(ct*(cpsi*(ct*(sphi*My0+cphi*Mx0)-st*Mz0) + spsi*(cphi*My0-sphi*Mx0))+st*(ct*Mz0+st*(sphi*My0+cphi*Mx0))) + cphi*(-1.*spsi*(ct*(sphi*My0+cphi*Mx0)-st*Mz0) + cpsi*(cphi*My0-sphi*Mx0));
	  Mz1 = ct*(ct*Mz0+st*(sphi*My0+cphi*Mx0)) - st*(cpsi*(ct*(sphi*My0+cphi*Mx0)-st*Mz0) + spsi*(cphi*My0-sphi*Mx0));

	} else {
	  
	  Mx1 = M[x(lp-1)];
	  My1 = M[y(lp-1)];
	  Mz1 = M[z(lp-1)];
	  
	}
      
      /* relaxation effects: "1" in Mz since Mo=1 by assumption */
      
      M[x(lp)] = Mx1*T2;
      M[y(lp)] = My1*T2;
      M[z(lp)] = Mz1 + (1 - Mz1)*T1;
      
    }
}
