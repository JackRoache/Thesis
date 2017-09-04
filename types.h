#ifndef TYPES_H
#define TYPES_H
#include <complex.h>
#include <arrayfire.h>

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

struct carray
{
    carray();
    carray(int dim);

    carray(int dim1, int dim2);
    void refresh();
    void resize(int dim);
    void resize(int dim1, int dim2);

    af::array a;
    af::array r;
    af::array i;
};
#endif // TYPES_H
