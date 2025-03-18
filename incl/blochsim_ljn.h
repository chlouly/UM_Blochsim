#ifndef __BLOCHSIM_LJN__
#define __BLOCHSIM_LJN__

/*
    This file contains macros and function
    declarations relevant to the 4-dimensional
    semisoid and flow formalism as outlined in
    the 2020 paper titled "Numerical approximation 
    to the general kinetic model for ASL 
    quantification", by Lee, Javed, Jao, and Nayak.
    Source code that uses this file includes:
        - src/bloch_ljn.c
*/

void bloch_LJN(/* ADD INPUT TYPES */);

#define X 0
#define Y 1
#define Z 2
#define S 3

#define x(n) (4 * (n))          // Gets the x index of a given timept
#define y(n) (4 * (n) + 1)      // Gets the y index of a given timept
#define z(n) (4 * (n) + 2)      // Gets the z index of a given timept
#define s(n) (4 * (n) + 3)      // Gets the s index of a given timept
#define off(n) (4 * (n))        // Gets the index offset of a 4-vector at a given timept
#define off3(n) (3 * (n))       // Gets the index offset of a 3-vector at a given timept

#define L2_SQRD(A) pow((A)[X], 2) + pow((A)[Y], 2) + pow((A)[Z], 2)
#define L2_NORM(A) sqrt(L2_SQRD(A))
#define NORMALIZE(A, norm, out) \
    (out)[X] = (A)[X] / (norm); \
    (out)[Y] = (A)[Y] / (norm); \
    (out)[Z] = (A)[Z] / (norm); \

#endif /* __BLOCHSIM_LJN__ */
