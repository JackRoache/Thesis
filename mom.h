#ifndef MOM_H
#define MOM_H

#include <af/array.h>

#include "types.h"

struct ImagingSpace {

    af::array x;
    af::array y; //corresponding positions of Er
    carray Er; //space that will be simulated
    af::cfloat medium;

    double minReal, maxReal;
    double minImag, maxImag; //for imaging purposes

    carray initalGuess;

    double dx {0}; //voxel size
    double dy {0};

    double lx {0}; //dimensions of space
    double ly {0};

    std::vector<Position> probes; //probes positions to simulate
    std::vector<float> freqs; //frequencies to simulate
};

struct RunInfo {
    std::string name; //used for image names
    int iterations; //Iterations to run for
    float lambda;
    af::Window *window;

};

class MOM
{
public:
    MOM();

    void setImagingSpace(ImagingSpace *space);
    void setIterations(RunInfo *info);

    typedef void (*imageCB)(RunInfo *, ImagingSpace *, carray &, int);
    void setCallBack(imageCB cb);

    void run();

private:

    void mom(int probenum, af::cfloat freq, bool simulate, carray &Er, carray &Et, carray &Es);
    void inverseBuilder(carray &Efunc, carray &C, af::cfloat k);

    void simulateSpace();
    void iterateMom();

    void spaceToImage();

    af::cfloat wavenumber(float freq, af::cfloat em);
    void pinv(carray &A, carray &Ai);
    void least_squares(af::array &A, af::array &b, af::array &x, double alpha);

    RunInfo *info {0};
    ImagingSpace *space {0};
    imageCB cb {0};
//    carray Es, Et;
    carray Ez; //Simulated data

};

#endif // MOM_H
