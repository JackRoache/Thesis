#ifndef MOM_H
#define MOM_H

#include <af/array.h>

#include "types.h"

struct ImagingSpace {

    af::array x;
    af::array y; //corresponding positions of Er
    carray Er; //space that will be simulated

    carray initalGuess;

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
    float lambda;

};

class MOM
{
public:
    MOM();

    void setImagingSpace(ImagingSpace *space);
    void setIterations(RunInfo *info);

    typedef std::function<void(RunInfo *, ImagingSpace *, carray &, int)> imageCB;
    void setCallBack(imageCB cb);

    void run();

private:

    void mom(int probenum, float freq, bool simulate);
    void inverseBuilder(carray &Efunc, carray &C, real k);

    void simulateSpace();
    void iterateMom();

    void spaceToImage();

    real wavenumber(real freq);
    void pinv(carray &A, carray &Ai);

    RunInfo *info {0};
    ImagingSpace *space {0};
    imageCB cb {0};
    carray Es, Et;
    carray Ez; //Simulated data

};

#endif // MOM_H
