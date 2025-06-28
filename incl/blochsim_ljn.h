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

void c_blochsim_ljn(
    double *,
    double *,
    double *,
    double *,
    int *,
    int,
    int,
    double,
    double,
    double,
    double,
    double,
    double,
    double,
    double,
    double,
    double,
    double
);

void c_blochsim_ljn_dyntime(
    double *,
    double *,
    double *,
    double *,
    double *,    
    int *,
    int,
    int,
    double,
    double,
    double,
    double,
    double,
    double,
    double,
    double,    
    double,
    double
);

void LJN_RF_excite(
    double *,
    double *,
    double *,
    double,
    double,
    double
);

void LJN_decay_and_transfer(
    double *,
    double *,
    double,
    double,
    double,
    double,
    double,
    double,
    double,
    double,
    double,
    double,
    double,
    double,
    double
);

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

#define ASSIGN_VEC(A, B) \
    (B)[X] = (A)[X]; \
    (B)[Y] = (A)[Y]; \
    (B)[Z] = (A)[Z]; \
    (B)[S] = (A)[S]; \

#define CRUSH(M) \
    (M)[X] = 0.0; \
    (M)[Y] = 0.0; \

#define PRINTVEC(A) printf("[ %f   %f   %f   %f ]\n", (A)[X], (A)[Y], (A)[Z], (A)[S]);
#define PRINTVEC3(A) printf("[ %f   %f   %f ]\n", (A)[X], (A)[Y], (A)[Z]);



#endif /* __BLOCHSIM_LJN__ */
