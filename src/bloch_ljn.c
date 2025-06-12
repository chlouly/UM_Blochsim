#include <math.h>
#include <stdio.h>
#include "blochsim_core.h"
#include "blochsim_ljn.h"

void c_blochsim_ljn(
    //INPUT ARGS
    double * M,             // Acts as function output
    double * B,
    double * s,             // s(t) values
    double * M_start,
    int * crush_inds,    // Indices of any crushers
    int num_crushers,       // Number of Crushers in the sequence
    int n_time,             // Number of time Indices
    double dt,              // tau
    double f,
    double ks,
    double kf,
    double R1f_app,
    double R2f_app,
    double R1s_app,
    double M0_s,
    double M0_f,
    double absorp,
    double s_sat
) {
    /* -- Precalculations -- */
    // double T1_app = (1 / T1f) + (F / lam);
    // double T2_app = (1 / T2f) + (F / lam);
    double T1f_app = 1.0 / R1f_app;
    double T1s = 1.0 / R1s_app;

    // Now we calculate the necessary 6 values of the ACE Matrix
    // (this matrix is constant given a constanat timestep)
    double exp_val = exp(-(f + 1.0) * ks * dt);

    double ace_11_22 = exp(-R2f_app * dt) * exp_val;
    double ace_33 = (1.0 + f * exp_val) * exp(-R1f_app * dt) / (1.0 + f);
    double ace_34 = (1.0 - exp_val) * exp(-R1s_app * dt) / (1.0 + f);
    double ace_43 = (f - f * exp_val) * exp(-R1f_app * dt) / (1.0 + f);
    double ace_44 = (f + exp_val) * exp(-R1s_app * dt) / (1.0 + f);





    // // M_a is the same as M_d from the previous itteration
    // double M_b[4];
    // // M_c is our output data, since that's what we keep
    // double M_d[4];

    // Initializing M
    ASSIGN_VEC(M_start, M);
    //ASSIGN_VEC(M_start, M_d);

    int cur_crush = 0;

    for (int i = 1; i < n_time; i++) {
        // M_b[i] = R_i M_d[i - 1]
        // Equivalent to RF pulses and gradients.
        LJN_RF_excite(M + off(i - 1), M + off(i), B + off3(i), dt, absorp, s_sat);

        if ((crush_inds) && (crush_inds[cur_crush] == i)) {
            CRUSH(M + off(i));      // Crush the transverse magnetization
            cur_crush += 1;
            cur_crush %= num_crushers;
        }

        // out = A(z_i) C(z_i) E(z_i) M_b[i] + (I - A(z_i) C(z_i) E(z_i)) D()
        LJN_decay_and_transfer(
            M + off(i),
            M + off(i),
            dt,
            s[i],
            T1f_app,
            T1s,
            kf,
            ks,
            M0_f,
            M0_s,
            ace_11_22,
            ace_33,
            ace_34,
            ace_43,
            ace_44
        );

        // M_d[i] = A(t_i - z_i) C(t_i - z_i) E(t_i - z_i) M_b[i] + (I - A(t_i - z_i) C(t_i - _i) E(t_i - z_i)) D()
        // LJN_decay_and_transfer(
        //     M + off(i),
        //     M_d,
        //     dt,
        //     s[i],
        //     T1_app,
        //     T2_app,
        //     T1s,
        //     kf,
        //     ks,
        //     1,   // M0z
        //     1,   // M0s
        //     f
        // );
    }
}

void LJN_RF_excite(
    double * M_in, 
    double * M_out, 
    double * B_eff, 
    double dt,
    double absorp,
    double s_sat
) {
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
            (c_t + pow(u_x, 2.0) * (1.0 - c_t)) * (M_in[X]) +
            (u_x * u_y * (1.0 - c_t) - u_z * s_t) * (M_in[Y]) + 
            (u_x * u_z * (1.0 - c_t) + u_y * s_t) * (M_in[Z]);
        M_out[Y] = 
            (u_x * u_y * (1.0 - c_t) + u_z * s_t) * (M_in[X]) + 
            (c_t + pow(u_y, 2.0) * (1.0 - c_t)) * (M_in[Y]) +
            (u_y * u_z * (1.0 - c_t) - u_x * s_t) * (M_in[Z]);
        M_out[Z] = 
            (u_x * u_z * (1.0 - c_t) - u_y * s_t) * (M_in[X]) +
            (u_z * u_y * (1.0 - c_t) + u_x * s_t) * (M_in[Y]) + 
            (c_t + pow(u_z, 2.0) * (1.0 - c_t)) * (M_in[Z]);

        // Scaling the Semisoid component
        // TODO - Absorption linshape of the semisoid pool
        // We're supposed to assume the semisoid pool gets saturated
        // NOTE: should there be a dt here?
        M_out[S] = M_in[S] * exp(-PI * pow(GAMMA * B_mag, 2.0) * absorp);
    } else {
        // Otherwise we just copy the vector
        ASSIGN_VEC(M_in, M_out);
    }

    if (isnaprox(s_sat, 0.0)) {
        M_out[S] *= exp(-PI * s_sat);
    }
}

void LJN_decay_and_transfer(
    double * M_in, 
    double * M_out, 
    double dt, 
    double s_t,
    double T1f_app,
    //double T2_app,
    double T1_s,
    double kf,
    double ks,
    double M_0f,
    double M_0s,
    double ace_11_22,
    double ace_33,
    double ace_34,
    double ace_43,
    double ace_44
    //double f
) {
    // if (isaprox(dt, 0.0)) {
    //     // For short time, this operation approaches I
    //     // so M_out = M_in. THis means we can skip all
    //     // the math below.
    //     ASSIGN_VEC(M_out, M_in);
    //     return;
    // }

    //D vector Definition (Only the last two components 
    // are nonzero)
    double denom = 1 + T1f_app * kf + T1_s * ks;
    double D_z = ((1 + T1_s * ks) / denom) * (M_0f + s_t * T1f_app) + ((T1f_app * ks) / denom) * M_0s;
    double D_s = ((T1_s * kf) / denom) * (M_0f + s_t * T1f_app) + ((1 + T1f_app * kf) / denom) * M_0s;

    // Finding relevant A matrix values
    // double exp_val = exp(-(f + 1) * ks * dt);
    // double a_33 = (1 + f * exp_val) / (1 + f);
    // double a_34 = (1 - exp_val) / (1 + f);
    // double a_43 = (f - f * exp_val) / (1 + f);
    // double a_44 = (f + exp_val) / (1 + f);

    // This comes from simplifying the operation 
    // shown in the paper
    // M_out = ACE M_in + (I - ACE)D
    //       = ACE(M_in - D) + D
    M_out[X] = ace_11_22 * M_in[X];
    M_out[Y] = ace_11_22 * M_in[Y];

    double mz_temp = ace_33 * (M_in[Z] - D_z) + ace_43 * (M_in[S] - D_s) + D_z;
    double ms_temp = ace_34 * (M_in[Z] - D_z) + ace_44 * (M_in[S] - D_s) + D_s;
    M_out[Z] = mz_temp;
    M_out[S] = ms_temp;

    return;
}
