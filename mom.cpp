#include "mom.h"

#include <complex>
#include <assert.h>
#include <arrayfire.h>
#include <algorithm>

#include "bessel.h" //CPU implementation
#include "clbessel.h" //OpenCL implementation

#include "Complex_Bessel/complex_bessel.h"

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
        std::cout << "Iteration " << i + 1 << " of " << info->iterations << std::endl;
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
    double a = sqrt(space->dx*space->dy / PI);
    af::cfloat scale;
    af::cfloat bj;
    if (k.imag == 0 && !info->slow){
        scale = af::cfloat(0, PI * a * (float)0.5) * k.real;
        bj = bessj(1, k.real*a);
    } else
    {
        std::complex<double> comp = sp_bessel::besselJ(1, std::complex<double>(k.real, k.imag));
        bj = af::cfloat(comp.real(), comp.imag());
        scale = af::cfloat(0, PI * a * (float)0.5) * k;
    }

    //Diaganol of matrix to solve
    af::cfloat D;
    if (k.imag == 0 && !info->slow){
        D = af::cfloat(0,0.5)*(PI*k.real*a* af::cfloat(bessj(1, a*k.real), -1 * bessy(1, a*k.real)) - af::cfloat(0,2));
    } else
    {
        std::complex<double> comp(k.real, k.imag);
        comp = sp_bessel::besselJ(1, a*comp) - std::complex<double>(0,1) * sp_bessel::besselY(1, a*comp);
        af::cfloat h(comp.real(), comp.imag());
        D = af::cfloat(0,0.5)*(PI* k * a * h - af::cfloat(0,2));
    }

    d = D * (Er_n) + 1;

    //Bulk of matrix to solve
    p = af::pow(af::tile(af::transpose(space->x), N, 1) - af::tile(space->x, 1, N), 2);
    p = p + af::pow(af::tile(af::transpose(space->y), N, 1) - af::tile(space->y, 1, N), 2);
    p = af::sqrt(p);
    if (k.imag == 0 && !info->slow){
        p *= k.real;
    } else {
        p *= k;
    }

    if (k.imag == 0 && !info->slow){
        fast_hankel(p, bh);
    } else {
        slow_hankel(p, bh);
    }

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

    carray phoHank;
    if (k.imag == 0 && !info->slow){
        pho = af::sqrt(pho) * k.real;
        fast_hankel(pho, phoHank);
        b = af::cfloat(0, 1) * k.real * 0.25 * phoHank.a;
    } else {
        pho = af::sqrt(pho) * k;
        slow_hankel(pho, phoHank);
        b = af::cfloat(0, 1) * k * 0.25 * phoHank.a;
    }

    Et.a = af::solve(c,b);

    //Simulate recieved EM for the probes.
    if (simulate){
        for (size_t i = 0; i < space->probes.size(); i++){
            float x0 = space->probes[i].x;
            float y0 = space->probes[i].y;

            af::array dis;
            carray esb;

            dis = af::sqrt(af::pow(x0 - space->x, 2) + af::pow(y0 - space->y, 2));
            if (k.imag == 0 && !info->slow){
                dis = dis * k.real;
                fast_hankel(dis, esb);
            } else {
                dis = dis * k;
                slow_hankel(dis, esb);
            }

            af::cfloat cons;
            if (k.imag == 0 && !info->slow){
                cons = af::cfloat(0,-1) * PI * k.real/ 2.0 * a * bessj(1, k.real*a);
            }else{
                std::complex<double> num(k.real, k.imag);
                num = sp_bessel::besselJ(1, num * a);
                af::cfloat f(num.real(), num.imag());
                cons = af::cfloat(0,-1) * PI * k/ 2.0 * a * f;
            }

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

    double a = sqrt(space->dx*space->dy / PI);

    af::array probeX(M), probeY(M);

    for(int i = 0; i < M; i++){
        probeX(i) = space->probes[i].x;
        probeY(i) = space->probes[i].y;
    }

    x2 = af::pow(af::tile(af::transpose(space->x), M, 1) - af::tile(probeX, 1, N), 2);
    y2 = af::pow(af::tile(af::transpose(space->y), M, 1) - af::tile(probeY, 1, N), 2);
    p = af::sqrt(x2 + y2);
    carray bj;
    af::cfloat cons;
    if (k.imag == 0 && !info->slow){
        p = p * k.real;
        fast_hankel(p, bj);
        cons = af::cfloat(0, -PI * k.real * 0.5) * a * bessj(1, k.real*a);
    } else {
        p = p * k;
        slow_hankel(p, bj);
        std::complex<double> num(k.real, k.imag);
        num = sp_bessel::besselJ(1, num * a);
        af::cfloat besselj(num.real(), num.imag());
        cons = af::cfloat(0, -PI * 0.5) * k * a * besselj;
    }

    B.a = cons * af::tile(af::transpose(Efunc.a),M,1);
    B.a = B.a * bj.a;
}

void MOM::simulateSpace()
{
    Ez.resize(space->probes.size() * space->probes.size() * space->freqs.size());
    assert(space->initalGuess.a.elements() == space->Er.a.elements());
    int probesSize = space->probes.size();
    for (size_t l = 0; l < space->freqs.size(); l++){
        af::cfloat k = wavenumber(space->freqs[l], space->medium_es, space->medium_cond);
        for (size_t i = 0; i < space->probes.size(); i++){
            std::cout << "MoM " << i + 1 << " of " << space->probes.size() << std::endl;
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
        af::cfloat k = wavenumber(space->freqs[f], space->medium_es, space->medium_cond);
        for (size_t i = 0; i < space->probes.size(); i++){
            carray Et, Es;
            std::cout << "MoM " << i + 1 << " of " << space->probes.size() << std::endl;
            mom(i, k, false, space->Er, Et, Es);
            computations.push_back(Et);
        }
    }

    std::vector<carray> matrices;
    carray m, B2;
    af::array B;
    for (size_t f = 0; f < space->freqs.size(); f++){
        af::cfloat k = wavenumber(space->freqs[f], space->medium_es, space->medium_cond);
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
    for (int i = 0; i < N; i++){
        B2.r(i,i) = lambda;
    }
//    af::array a = af::constant(lambda, N);
//    B2.r = af::diag(a, 0, false);

    B = af::join(0, B, B2.a);

    L = Ez.a.elements();
    Ereg.resize(N + L);
    Ereg.a(af::seq(0, L-1)) = Ez.a;
    Ereg.a(af::seq(L, L+N-1)) = 0;

    /* Solve the system
     * BO'=E_{scat}
     */
    pinv(Ereg, Inverse);

    Treg.a = af::matmul(Inverse.a, B);

    af::cfloat n(1,0);
//    space->Er.a = af::transpose(Treg.a + n, true) ; //correct for offset and apply complex conjugate / transpose
    space->Er.i = af::imag(af::transpose(Treg.a + n, true));
    space->Er.r = af::real(af::transpose(Treg.a + n, true));
    space->Er.refresh();
//    space->Er.a = af::transpose(Treg.a) + n; //correct for offset
}

void MOM::pinv(carray &A, carray &Ai)
{
    int minDim = std::min(A.a.row(0).elements(), A.a.col(0).elements());
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
void MOM::fast_hankel(af::array &in, carray &out)
{
    bessj0(in, out.r);
    bessy0(in, out.i);
    out.i = out.i * -1;
    out.refresh();
}

#if 1
#include <QtConcurrent/QtConcurrentMap>
void map_function(std::complex<double> &num)
{
    num = sp_bessel::hankelH2(0, num);
}
#endif

void MOM::slow_hankel(af::array &in, carray &out)
{
#if 0
    float *_real = af::real(in).host<float>();
    float *_imag = af::imag(in).host<float>();
    out.a = af::array(in.dims());

    for (int i = 0; i < in.elements(); i++){
        std::complex<double> num(_real[i], _imag[i]);
        num = sp_bessel::hankelH2(0, num);
        _real[i] = num.real();
        _imag[i] = num.imag();
    }
    out.r = af::array(in.dims(), _real);
    out.i = af::array(in.dims(), _imag);
    out.refresh();

    af::freeHost(_real);
    af::freeHost(_imag);
#else
    float *_real = af::real(in).host<float>();
    float *_imag = af::imag(in).host<float>();
    out.a = af::array(in.dims());

    QVector<std::complex<double>> nums;
    for (int i = 0; i < in.elements(); i++){
        nums.append(std::complex<double>(_real[i], _imag[i]));
    }

    QFuture<void> future = QtConcurrent::map(nums, map_function);
    future.waitForFinished();

    for (int i = 0; i < in.elements(); i++){
        _real[i] = nums[i].real();
        _imag[i] = nums[i].imag();
    }

    out.r = af::array(in.dims(), _real);
    out.i = af::array(in.dims(), _imag);
    out.refresh();

    af::freeHost(_real);
    af::freeHost(_imag);

#endif
}

void MOM::spaceToImage()
{

}

af::cfloat MOM::wavenumber(double freq, double es, double cond)
{
    std::complex<float> complex(0.0,1.0);

    float real = std::pow(2 * PI * freq, 2)* CU0 * es * CE0;
    float imag = 2 * PI * freq * CU0 * cond / CE0;

    std::complex<float> num = sqrt(real - complex * imag);
    af::cfloat out(num.real(), num.imag());
    return out;
}
