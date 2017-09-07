#ifndef MOM_H
#define MOM_H

#include <af/array.h>

#include "types.h"

struct ImagingSpace {

    af::array x;
    af::array y; //corresponding positions of Er
    carray Er; //space that will be simulated

    real dx {0}; //voxel size
    real dy {0};

    real lx {0}; //dimensions of space
    real ly {0};

    std::vector<Position> probes; //probes positions to simulate
    std::vector<float> freqs; //frequencies to simulate
};

struct RunInfo {
    std::string name; //used for image names
    int iterations; //Iterations to run for

};

class MOM
{
public:
    MOM();

    void setImagingSpace(ImagingSpace *space);
    void setIterations(RunInfo *info);

    void run();

private:

    void mom(int probenum, float freq, bool simulate);
    void inverseBuilder(carray &Efunc, carray &C, real k);

    void simulateSpace();
    void initialGuess();
    void iterateMom();

    void spaceToImage();

    real wavenumber(real freq);
    void pinv(carray &A, carray &Ai);

    RunInfo *info {0};
    ImagingSpace *space {0};
    carray Es, Et;
    carray Ez; //Simulated data
};

#endif // MOM_H
