#include <math.h>
#include <stdio.h>
#include "blochsim_core.h"

void bloch_LJN(
    //INPUT ARGS
    int n_time,
    double ** out
) {
    
    // M_a is the same as M_d from the previous itteration
    double M_b[4];
    // M_c is our output data, since that's what we keep
    double M_d[4];

    for (int i = 0; i < n_time; i++) {
        // M_b[i] = R_i M_d[i - 1]
        // Equivalent to RF pulses and gradients.
        // TODO

        // out = A(z_i) C(z_i) E(z_i) M_b[i] + (I - A(z_i) C(z_i) E(z_i)) D()
        // TODO

        // M_d[i] = A(t_i - z_i) C(t_i - z_i) E(t_i - z_i) M_b[i] + (I - A(t_i - z_i) C(t_i - _i) E(t_i - z_i)) D()
        // TODO
    }
}

void LJN_RF_excite(double ** M_in, double ** M_out) {

}