#include "mom.h"

#include <complex>
#include <assert.h>
#include <arrayfire.h>
#include <algorithm>

#include "bessel.h" //CPU implementation
#include "clbessel.h" //OpenCL implementation

#define CE0     (float)(8.854e-12)
#define CU0     (float)(12.56e-7)
#define PI      (float)(3.14159265)
#define CC      3e8

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

    space->Er = space->initalGuess;

    for (int i = 0; i < info->iterations; i++){
        std::cout << "Iteration " << i + 1  << std::endl;
        iterateMom();

        if(cb)
            cb(info, space, space->Er, i);
    }
}
/* Solve the forward problem using Method of Moments
 * E_m^{inc} = E_m^t + (jk^2/4)\sum_{n=1}^N(\epsilon_r-1)E_n\int\int_{area\ n}H_0^{(2)}(k\rho)dx'dy'
 */
void MOM::mom(int probenum, af::cfloat k, bool simulate, carray &Er, carray &Et, carray &Es)
{
    assert(space->x.elements() == space->y.elements());
    assert(space->x.elements() == Er.a.elements());
    int N = space->x.elements();

    if (Et.a.elements() != N)
        Et.resize(N);
    if (Es.a.elements() != space->probes.size())
        Es.resize(space->probes.size());


    af::array pho(N), d(N), Er_n, b(N); //intermidiearies
    carray bh(N);
    af::array  p(N, N), c(N, N); //main matrix

    //offset matrix by 1 relative permability
    Er_n = Er.a - 1.0;

    //calculate constants
    float a = sqrt(space->dx*space->dy / PI);
    af::cfloat scale(0, PI * k.real * a * (float)0.5);
    float bj = bessj(1, k.real*a);

    //Diaganol of matrix to solve
    af::cfloat D = af::cfloat(0,0.5)*(PI*k.real*a* af::cfloat(bessj(1, a*k.real), -1 * bessy(1, a*k.real)) - af::cfloat(0,2));
    d = D * (Er_n) + 1;

    //Bulk of matrix to solve
    p = af::pow(af::tile(af::transpose(space->x), N, 1) - af::tile(space->x, 1, N), 2);
    p = p + af::pow(af::tile(af::transpose(space->y), N, 1) - af::tile(space->y, 1, N), 2);
    p = af::sqrt(p) * k.real;

    bessj0(p, bh.r);
    bessy0(p, bh.i);
    bh.i = bh.i * -1;
    bh.refresh();

    c = bj * af::tile(af::transpose(Er_n, false), N, 1);
    c = c * bh.a;
    c = c * scale;

    //combine c and d, d is the diaganol of c
    assert(c.dims()[0] == N);
    for (int i = 0; i < N; i++){
        c(i,i) = d(i);
    }

    float probeX = space->probes[probenum].x;
    float probeY = space->probes[probenum].y;

    pho = af::pow(space->x - probeX, 2) + af::pow(space->y - probeY, 2);
    pho = af::sqrt(pho) * k.real;

    af::array phobj, phoby;

    bessj0(pho, phobj);
    bessy0(pho, phoby);
    b = af::cfloat(0, 1) * k.real * 0.25 * (phobj - af::cfloat(0,1) * phoby);

    Et.a = af::solve(c,b);

    //Simulate recieved EM for the probes.
    if (simulate){
        for (size_t i = 0; i < space->probes.size(); i++){
            float x0 = space->probes[i].x;
            float y0 = space->probes[i].y;

            af::array dis;
            carray esb;

            dis = af::sqrt(af::pow(x0 - space->x, 2) + af::pow(y0 - space->y, 2));
            dis = dis * k.real;
            bessj0(dis, esb.r);
            bessy0(dis, esb.i);
            esb.i = esb.i * -1;
            esb.refresh();
            af::cfloat cons = af::cfloat(0,-1) * PI * k.real/ 2.0 * a * bessj(1, k.real*a);
            af::array sum = af::sum(cons * (Er_n) * Et.a * esb.a);
            Es.a(i) = sum;
        }
    }
}
/*
 * This function determines the inverting matrix for the born method
 *
 * b_{i,j} = - \frac{j}{4}\int_{s_j} \bm{E}_{inc}(x^{'},y^{'})  H_0^{(2)}(k_m\rho(x,y,x^{'},y^{'}))dx^{'}dy^{'}
 */
void MOM::inverseBuilder(carray &Efunc, carray &B, af::cfloat k)
{
    assert(space->x.elements() == space->y.elements());

    int M = space->probes.size();
    int N = space->x.elements();
    B.resize(M, N);
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
    p = p * k.real;
    af::array bj0, by0;
    bessj0(p, bj0);
    bessy0(p, by0);

    af::cfloat cons = af::cfloat(0, -PI * k.real * 0.5) * a * bessj(1, k.real*a);
    B.a = cons * af::tile(af::transpose(Efunc.a),M,1);
    B.a = B.a * (bj0 - af::cfloat(0, 1) * by0);
}

void MOM::simulateSpace()
{
    Ez.resize(space->probes.size() * space->probes.size() * space->freqs.size());
    assert(space->initalGuess.a.elements() == space->Er.a.elements());
    int probesSize = space->probes.size();
    for (size_t l = 0; l < space->freqs.size(); l++){
        af::cfloat k = wavenumber(space->freqs[l], space->medium);
        for (size_t i = 0; i < space->probes.size(); i++){
            carray Et, Es;
            mom(i, k, true, space->Er, Et, Es);
            assert(Es.a.elements() == probesSize);

            af::seq s(i*probesSize, (i+1)*probesSize - 1);
            assert(s.size == probesSize);
            Ez.a(s) = Es.a;
        }
    }
}

void MOM::iterateMom()
{
    std::vector<carray> computations;
    carray Ereg, Treg, Inverse;

    for (size_t f = 0; f < space->freqs.size(); f++){
        af::cfloat k = wavenumber(space->freqs[f], space->medium);
        for (size_t i = 0; i < space->probes.size(); i++){
            carray Et, Es;
            computations.push_back(Et);
            mom(i, k, false, space->Er, computations.at(computations.size()-1), Es);
        }
    }

    std::vector<carray> matrices;
    carray m, B2;
    af::array B;
    for (size_t f = 0; f < space->freqs.size(); f++){
        af::cfloat k = wavenumber(space->freqs[f], space->medium);
        for (size_t i = 0; i < space->probes.size(); i++){
            inverseBuilder(computations[i+space->probes.size() * f], m, k);
            matrices.push_back(m);
        }
    }

    B =  matrices.at(0).a;

    for (int i = 1; i < matrices.size(); i++){
        B = af::join(0, B, matrices[i].a);
    }

    Ereg.a = Ez.a;

    /* Solve the system
     * BO^{'}=E_{scat}
     */
    least_squares(Ereg.a, B, Treg.a, info->lambda);

    af::cfloat n(1,0);
    space->Er.a = af::transpose(Treg.a + n, true) ; //correct for offset and apply complex conjugate / transpose
//    space->Er.a = af::transpose(Treg.a) + n; //correct for offset
}

void MOM::pinv(carray &A, carray &Ai)
{
    int minDim = min(A.a.row(0).elements(), A.a.col(0).elements());
    af::array u, vt;
    af::array s;

    af::svd(u, s, vt, A.a);

    //prevents divide by zero problems
    af::array index = (af::iszero(s) - 1) * -1 * af::seq(s.elements());
    s(index) = 1/s(index);

    u = u(af::span, af::seq(minDim));
    s = af::diag(s, 0, false).as(c32);
    vt = vt(af::seq(minDim), af::span);

//    for (int i = 0; i < minDim; i++){
//        s(i,i) = 1/s(i,i);
//    }

    Ai.a = af::matmul(vt.H(), s, u.H());
}

void MOM::least_squares(af::array &A, af::array &b, af::array &x, double alpha)
{
    int minDim = min(A.row(0).elements(), A.col(0).elements());

    af::array u, vt, s, d;
    af::svd(u, s, vt, A);

    //remove trailing zeroes
    u = u(af::span, af::seq(minDim));
    vt = vt(af::seq(minDim), af::span);

    //D_{ii} = sigma_i / (sigma_i^2 + alpha^2)
    d = s / (af::pow(s,2) + alpha * alpha);
    d = af::diag(d, 0, false).as(c32);

    u = af::transpose(u, false);
    x = af::matmul(vt.H(), d, u, b);
}

void MOM::spaceToImage()
{

}

af::cfloat MOM::wavenumber(float freq, af::cfloat em)
{
    std::complex<float> complex(0.0,1.0);

    float real = pow(2 * PI * freq, 2)* CU0 * em.real * CE0;
    float imag = 2 * PI * freq * CU0 * em.imag / CE0;

    std::complex<float> num = sqrt(real - complex * imag);
    af::cfloat out(num.real(), num.imag());
    return out;
}
