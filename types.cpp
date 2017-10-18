#include "types.h"

carray::carray(){}
carray::carray(int dim):
    r(dim),
    i(dim)
{
    refresh();
}

carray::carray(int dim1, int dim2):
    r(dim1, dim2),
    i(dim1, dim2)
{
    refresh();
}
void carray::refresh()
{
    a = af::complex(r, i);
}
void carray::resize(int dim)
{
    r = af::array(dim, TYPE_R);
    i = af::array(dim, TYPE_R);
    refresh();
}
void carray::resize(int dim1, int dim2){
    r = af::array(dim1, dim2, TYPE_R);
    i = af::array(dim1, dim2, TYPE_R);
    refresh();
}
