#ifndef CLBESSEL_H
#define CLBESSEL_H
#include "arrayfire.h"

void initKernels(af::array dummy);
void bessj0(af::array &in, af::array &out);
void bessj1(af::array &in, af::array &out);
void bessy0(af::array &in, af::array &out);

#endif // CLBESSEL_H
