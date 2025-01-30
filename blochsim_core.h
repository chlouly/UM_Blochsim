#ifndef __BLOCHSIM__
#define __BLOCHSIM__

void c_blochsim(double *, double *, double, double, int, double);
void c_blochsim_rk4(double *, double *, double, double, int, double);
void c_bloch_step(double *, double *, double *, double *, double, double, double);

#define PI 3.141592653589
#define GAMBAR 42577.478461     // gamma/2pi [kHz/T]
#define GAMMA GAMBAR * 2 * PI    // gamma [rad/us/T]

#define IND(n) (3 * (n))

#define X 0
#define Y 1
#define Z 2

#define x(n) (3 * (n))
#define y(n) (3 * (n) + 1)
#define z(n) (3 * (n) + 2)

#define SCALE_VEC(A, s) \
    (A)[X] *= (s); \
    (A)[Y] *= (s); \
    (A)[Z] *= (s); \

#define ADD_SCALED_VEC(A, B, s) \
    (A)[X] += ((B)[X] * (s)); \
    (A)[Y] += ((B)[Y] * (s)); \
    (A)[Z] += ((B)[Z] * (s)); \

#define ADD_VEC(A, B) ADD_SCALED_VEC((A), (B), 1)

#endif /* __BLOCHSIM__ */
