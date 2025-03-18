#include <math.h>
#include <stdio.h>
#include "blochsim_core.h"
#include "blochsim_ljn.h"

void bloch_LJN(
    //INPUT ARGS
    double * M,
    double * B,
    int n_time,
    double dt,
    double T1s,
    double T2s,
    double T1f,
    double T2f
) {
    
    // M_a is the same as M_d from the previous itteration
    double M_b[4];
    // M_c is our output data, since that's what we keep
    double M_d[4];

    for (int i = 1; i < n_time; i++) {
        // M_b[i] = R_i M_d[i - 1]
        // Equivalent to RF pulses and gradients.
        LJN_RF_excite(M + off(i - 1), M_b, B + off(i), dt);

        // out = A(z_i) C(z_i) E(z_i) M_b[i] + (I - A(z_i) C(z_i) E(z_i)) D()
        // TODO

        // M_d[i] = A(t_i - z_i) C(t_i - z_i) E(t_i - z_i) M_b[i] + (I - A(t_i - z_i) C(t_i - _i) E(t_i - z_i)) D()
        // TODO
    }
}

void LJN_RF_excite(double * M_in, double * M_out, double * B_eff, double dt) {
    // This function essentially acts as left multiplication by a rotation matrix for the
    // first three dimensions, and a scaling of the fourth dimension. Be careful, for speed
    // there is no input validation, so make sure the values passed to the whole bloch sim
    // are within normal ranges.
    //
    // rot_ax should be a unit-norm 3-vector that describes the axis to be rotated about

    double B_mag = L2_NORM(B_eff);

    // Rotations (theta !~= 0)
    if (isnaprox(B_mag, 0.0)) {
        double theta = GAMMA * dt * B_mag;
        double u[3];
        NORMALIZE(B_eff, B_mag, u);

        double c_t = cos(theta);
        double s_t = sin(theta);
    
        double u_x = u[X];
        double u_y = u[Y];
        double u_z = u[Z];

        M_out[X] = 
            (c_t + pow(u_x, 2.0) * (1 - c_t)) * (M_in[X]) +
            (u_x * u_y * (1 - c_t) - u_z * s_t) * (M_in[Y]) + 
            (u_x * u_z * (1 - c_t) + u_y * s_t) * (M_in[Z]);
        M_out[Y] = 
            (u_x * u_y * (1 - c_t) + u_z * s_t) * (M_in[X]) + 
            (c_t + pow(u_y, 2.0) * (1 - c_t)) * (M_in[Y]) +
            (u_y * u_z * (1 - c_t) + u_x * s_t) * (M_in[Z]);
        M_out[Z] = 
            (u_x * u_z * (1 - c_t) - u_y * s_t) * (M_in[X]) +
            (u_z * u_y * (1 - c_t) + u_x * s_t) * (M_in[Y]) + 
            (c_t + pow(u_z, 2.0) * (1 - c_t)) * (M_in[Z]);
    }

    // Scaling the Semisoid component
    // TODO - Absorption linshape of the semisoid pool
    M_out[S] = M_in[S] * exp(-PI * pow(GAMMA * B_mag, 2.0) * dt);
}
