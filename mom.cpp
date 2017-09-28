#include "mom.h"

#include <assert.h>
#include <arrayfire.h>
#include <algorithm>

#include "bessel.h" //CPU implementation
#include "clbessel.h" //OpenCL implementation

#define CE0     (float)(8.854e-12)
#define CU0     (float)(12.56e-7)
#define PI      (float)(3.14159265)
#define CC      3e8
#define CEM     1


MOM::MOM()
{

}

void MOM::setImagingSpace(ImagingSpace *space)
{
    assert(space);
    this->space = space;
}

void MOM::setIterations(RunInfo *info)
{
    assert(info);
    this->info = info;
}

void MOM::setCallBack(MOM::imageCB fn)
{
    cb = fn;
}

void MOM::run()
{
    simulateSpace();

    for (int i = 0; i < info->iterations; i++){
        std::cout << "Iteration " << i + 1  << std::endl;
        iterateMom();

        if(cb)
            cb(info, space, space->Er, i);
    }
}

void MOM::mom(int probenum, float k, bool simulate)
{
    assert(space->x.elements() == space->y.elements());
    assert(space->x.elements() == space->Er.a.elements());
    int N = space->x.elements();

    if (Et.a.elements() != N)
        Et.resize(N);
    if (Es.a.elements() != space->probes.size())
        Es.resize(space->probes.size());


    af::array pho(N), d(N), Er_n, b(N); //intermidiearies
    carray bh(N);
    af::array  p(N, N), c(N, N); //main matrix

    //offset matrix by 1 relative permability
    Er_n = space->Er.a - 1.0;

    //calculate constants
    float a = sqrt(space->dx*space->dy / PI);
    af::cfloat scale(0, PI * k * a * (float)0.5);
    float bj = bessj(1, k*a);

    //Diaganol of matrix to solve
    af::cfloat D = af::cfloat(0,0.5)*(PI*k*a* af::cfloat(bessj(1, a*k), -1 * bessy(1, a*k)) - af::cfloat(0,2));
    d = D * (space->Er.a - 1) + 1;

    //Bulk of matrix to solve
    p = af::pow(af::tile(af::transpose(space->x), N, 1) - af::tile(space->x, 1, N), 2);
    p = p + af::pow(af::tile(af::transpose(space->y), N, 1) - af::tile(space->y, 1, N), 2);
    p = af::sqrt(p) * k;

    bessj0(p, bh.r);
    bessy0(p, bh.i);
    bh.i = bh.i * -1;
    bh.refresh();

    c = bj * af::tile(af::transpose(Er_n, false), N, 1);
    c = c * bh.a;
    c = c * scale;

    //combine c and d, d is the diaganol of c
    for (int i = 0; i < N; i++){
        c(i,i) = d(i);
    }

    float probeX = space->probes[probenum].x;
    float probeY = space->probes[probenum].y;

    pho = af::pow(space->x - probeX, 2) + af::pow(space->y - probeY, 2);
    pho = af::sqrt(pho) * k;

    af::array phobj, phoby;

    bessj0(pho, phobj);
    bessy0(pho, phoby);
    b = af::cfloat(0, 1) * k * 0.25 * (phobj - af::cfloat(0,1) * phoby);
    Et.a = af::solve(c,b);

    //Simulate recieved EM for the probes.
    if (simulate){
        for (size_t i = 0; i < space->probes.size(); i++){
            float x0 = space->probes[i].x;
            float y0 = space->probes[i].y;

            af::array dis;
            carray esb;

            dis = af::sqrt(af::pow(x0 - space->x, 2) + af::pow(y0 - space->y, 2));
            dis = dis * k;
            bessj0(dis, esb.r);
            bessy0(dis, esb.i);
            esb.i = esb.i * -1;
            esb.refresh();
            af::cfloat cons = af::cfloat(0,-1) * PI * k/ 2.0 * a * bessj(1, k*a);
            Es.a(i) = af::sum(cons * (space->Er.a - 1) * Et.a * esb.a);
        }
    }
}

void MOM::inverseBuilder(carray &Efunc, carray &C, float k)
{
    assert(space->x.elements() == space->y.elements());

    int M = space->probes.size();
    int N = space->x.elements();
    C.resize(M, N);
    af::array p, x2, y2;

    float a = sqrt(space->dx*space->dy / PI);

    af::array probeX(M), probeY(M);

    for(int i = 0; i < M; i++){
        probeX(i) = space->probes[i].x;
        probeY(i) = space->probes[i].y;
    }

    x2 = af::pow(af::tile(af::transpose(space->x), M, 1) - af::tile(probeX, 1, N), 2);
    y2 = af::pow(af::tile(af::transpose(space->y), M, 1) - af::tile(probeY, 1, N), 2);
    p = af::sqrt(x2 + y2);
    p = p * k;
    af::array bj0, by0;
    bessj0(p, bj0);
    bessy0(p, by0);

    af::cfloat cons = af::cfloat(0, -PI * k * 0.5) * a * bessj(1, k*a);
    C.a = cons * af::tile(af::transpose(Efunc.a),M,1);
    C.a = C.a * (bj0 - af::cfloat(0, 1) * by0);
}

void MOM::simulateSpace()
{
    Ez.resize(space->probes.size() * space->probes.size() * space->freqs.size());
    int probesSize = space->probes.size();
    for (size_t l = 0; l < space->freqs.size(); l++){
        float k = wavenumber(space->freqs[l]);
        for (size_t i = 0; i < space->probes.size(); i++){
            mom(i, k, true);

            af::seq s(i*probesSize, (i+1)*probesSize - 1);
            Ez.a(s) = Es.a;
        }
    }
}

void MOM::iterateMom()
{
    std::vector<carray> computations;
    carray Ereg, Treg;

    for (size_t f = 0; f < space->freqs.size(); f++){
        float k = wavenumber(space->freqs[f]);
        for (size_t i = 0; i < space->probes.size(); i++){
            mom(i, k, false);
            computations.push_back(Et);
        }
    }

    std::vector<carray> matrices;
    carray m, B2;
    af::array B;
    for (size_t f = 0; f < space->freqs.size(); f++){
        float k = wavenumber(space->freqs[f]);
        for (size_t i = 0; i < space->probes.size(); i++){
            inverseBuilder(computations[i+space->probes.size() * f], m, k);
            matrices.push_back(m);
        }
    }

    B =  matrices.at(0).a;

    for (int i = 1; i < matrices.size(); i++){
        B = af::join(0, B, matrices[i].a);
    }

    float lambda = info->lambda;

    int L = B.col(0).elements(); //rows
    int N = B.row(0).elements(); //columns
    B2.resize(N,N);

    B2.i = af::constant(0, B2.i.dims());
    B2.r = af::constant(0, B2.r.dims());

    //set diaganol to lambda
    for (int i = 0; i < N; i++)
        B2.r(i,i) = lambda;

    B2.refresh();
    B = af::join(0, B, B2.a);

    L = Ez.a.elements();
    Ereg.resize(N + L);
    Ereg.a(af::seq(0, L-1)) = Ez.a;
    Ereg.a(af::seq(L, L+N-1)) = 0;

    //solve the system
    pinv(Ereg, Treg);

    Treg.a = af::matmul(Treg.a, B);

    af::cfloat n(1);
    space->Er.a = af::transpose(Treg.a) + n; //correct for offset
}

void MOM::pinv(carray &A, carray &Ai)
{
    int minDim = min(A.a.row(0).elements(), A.a.col(0).elements());
    af::array E(minDim, minDim);
    af::array u, vt;
    af::array s;

    af::svd(u, s, vt, A.a);

    u = u(af::span, af::seq(minDim));
    s = af::diag(s, 0, false).as(c32);
    vt = vt(af::seq(minDim), af::span);

    for (int i = 0; i < minDim; i++)
        s(i,i) = 1/s(i,i);

    Ai.a = af::matmul(vt.H(), s, u.H());
}

void MOM::spaceToImage()
{

}

float MOM::wavenumber(float freq)
{
    return sqrt(pow(2 * PI * freq, 2) * CU0 * CEM * CE0);
}
