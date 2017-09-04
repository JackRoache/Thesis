#ifndef TYPES_H
#define TYPES_H
#include <complex.h>
#include <af/array.h>

#define MIN(x, y) (x < y ? x : y)
#define MAX(x, y) (x < y ? y : x)

typedef std::complex<float> comp;
typedef float real;

#define acomp c32
#define areal f32

typedef struct Pos {
    float x;
    float y;
} Position;

#endif // TYPES_H
