#ifndef __BLOCHSIM__
#define __BLOCHSIM__

void c_blochsim(double *, double *, double, double, int, double);
void c_blochsim_rk4(double *, double *, double, double, int, double);
void c_bloch_step(double *, double *, double *, double *, double, double, double, double, double);

#define PI 3.141592653589
#define GAMBAR 42577.478461     // gamma/2pi [kHz/T]
#define GAMMA GAMBAR * 2 * PI    // gamma [rad/us/T]

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

#define ZERO(M) \
    (M)[X] = 0.0; \
    (M)[Y] = 0.0; \
    (M)[Z] = 0.0; \

#endif /* __BLOCHSIM__ */
