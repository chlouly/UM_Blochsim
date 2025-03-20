#ifndef __BLOCHSIM_CORE__
#define __BLOCHSIM_CORE__

/*
    Definitions shared across the entire program.
    Macros, and function declarations specific to 
    certain programs can be found in other header 
    files of the form "incl/blochsim_*.h".
*/

#define PI 3.141592653589
#define GAMBAR 42577.478461     // gamma/2pi [kHz/T]
#define GAMMA GAMBAR * 2 * PI    // gamma [rad/us/T]
#define TOLERANCE 1e-12

#define isnaprox(a, b) ((((a) - (b)) > TOLERANCE) || (((a) - (b)) < -TOLERANCE))
#define isaprox(a, b) (!isnaprox(a, b))

#endif /* __BLOCHSIM_CORE__ */
