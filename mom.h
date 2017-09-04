#ifndef MOM_H
#define MOM_H

#include <af/array.h>

#include "types.h"

struct ImagingSpace {

    af::array x;
    af::array y;
    af::array Er;

    real dx {0};
    real dy {0};

    std::vector<Position> probes;
};


class MOM
{
public:
    MOM();





};

#endif // MOM_H
