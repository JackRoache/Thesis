#ifndef BIM_H
#define BIM_H

#include <af/array.h>

#include "types.h"

struct ImagingSpace {

    af::array x;
    af::array y; //corresponding positions of Er
    carray Er; //space that will be simulated
//    af::cfloat medium;

    double medium_es;
    double medium_cond;

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
    bool slow;
};

class BIM
{
public:
    BIM();

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

    af::cfloat wavenumber(double freq, double es, double cond);
    void pinv(af::array &A, carray &Ai);
    void slow_hankel(af::array &in, carray &out);

    RunInfo *info {0};
    ImagingSpace *space {0};
    imageCB cb {0};
//    carray Es, Et;
    carray Ez; //Simulated data

    void fast_hankel(af::array &in, carray &out);

    void tikhonov_reg(af::array &A, af::array &b, carray &out, double lambda);
};

#endif // BIM_H
