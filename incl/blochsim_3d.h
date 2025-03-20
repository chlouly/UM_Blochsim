#ifndef __BLOCHSIM_3D__
#define __BLOCHSIM_3D__

/*
    This header file contains macros and 
    function declarations specific to the
    3-dimensional Bloch simulation formalisms
    This includes:
        - src/bloch_eul.c
        - src/bloch_rk4.c
*/

void c_blochsim(double *, double *, double, double, int, double);
void c_blochsim_rk4(double *, double *, double, double, int, double);
void c_bloch_step(double *, double *, double *, double *, double, double, double, double, double);

#define X 0
#define Y 1
#define Z 2

#define x(n) (3 * (n))
#define y(n) (3 * (n) + 1)
#define z(n) (3 * (n) + 2)
#define off(n) (3 * (n))

#define PRINT_VEC(A) printf("[ %f   %f   %f ]\n", (A)[X], (A)[Y], (A)[Z]);

#define FILL(A, B) \
    (A)[X] = (B)[X]; \
    (A)[Y] = (B)[Y]; \
    (A)[Z] = (B)[Z]; \

#define CROSS(A, B, C, scale) \
    (C)[X] = ((A)[Y] * (B)[Z] - (A)[Z] * (B)[Y]) * (scale); \
    (C)[Y] = ((A)[Z] * (B)[X] - (A)[X] * (B)[Z]) * (scale); \
    (C)[Z] = ((A)[X] * (B)[Y] - (A)[Y] * (B)[X]) * (scale); \

#define SCALE_VEC(A, s) \
    (A)[X] *= (s); \
    (A)[Y] *= (s); \
    (A)[Z] *= (s); \

#define ADD_SCALED_REPLACE(A, B, C, s) \
    (C)[X] = ((A)[X] + (B)[X]) * (s); \
    (C)[Y] = ((A)[Y] + (B)[Y]) * (s); \
    (C)[Z] = ((A)[Z] + (B)[Z]) * (s); \

#define ADD_SCALED_VEC(A, B, s) \
    (A)[X] += ((B)[X] * (s)); \
    (A)[Y] += ((B)[Y] * (s)); \
    (A)[Z] += ((B)[Z] * (s)); \

#define ADD_VEC(A, B) ADD_SCALED_VEC((A), (B), 1)

#define ADD_VEC_REPLACE(A, B, C) ADD_SCALED_REPLACE((A), (B), (C), 1)

#define MIDPT(A, B, out) \
    (out)[X] = ((A)[X] + (B)[X]) / 2; \
    (out)[Y] = ((A)[Y] + (B)[Y]) / 2; \
    (out)[Z] = ((A)[Z] + (B)[Z]) / 2; \

#define EXTRAP(A, B, out) \
    (out)[X] = (2 * (B)[X] - (A)[X]); \
    (out)[Y] = (2 * (B)[Y] - (A)[Y]); \
    (out)[Z] = (2 * (B)[Z] - (A)[Z]); \

#define ZERO(M) \
    (M)[X] = 0.0; \
    (M)[Y] = 0.0; \
    (M)[Z] = 0.0; \

#define RKADD_REP(M, k, M_out, scale) \
    (M_out)[X] = (M)[X] + ((k)[X] * (scale)); \
    (M_out)[Y] = (M)[Y] + ((k)[Y] * (scale)); \
    (M_out)[Z] = (M)[Z] + ((k)[Z] * (scale)); \

#define RKADD(k, M_out, scale) \
    (M_out)[X] += (k)[X] * (scale); \
    (M_out)[Y] += (k)[Y] * (scale); \
    (M_out)[Z] += (k)[Z] * (scale); \

#endif /* __BLOCHSIM_3D__ */
