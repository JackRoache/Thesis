#include "mom.h"

#include <assert.h>
#include <arrayfire.h>

#include "bessel.h" //CPU implementation
#include "clbessel.h" //OpenCL implementation

#define CE0     (real)(8.854e-12)
#define CU0     (real)(12.56e-7)
#define PI      (real)(3.14159265)
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

void MOM::run()
{
    simulateSpace();

    initialGuess();


    for (int i = 0; i < info->iterations; i++){
        std::cout << "Iteration " << i + 1  << std::endl;
        iterateMom();
        af::print("Er",  space->Er.a);
    }
}

void MOM::mom(int probenum, real k)
{
    //void mom(af::array &space->x, af::array &space->y, real dx, real dy, carray &Er, real k, std::vector<Position> &probes, int prob_num, carray &Es, carray &Et)

    assert(space->x.elements() == space->y.elements());
    assert(space->x.elements() == space->Er.a.elements());
    comp com(0,1);
    int N = space->x.elements();

    Et.resize(N);
    Es.resize(space->probes.size());

    af::array p(N), x2(N), y2(N), pho(N);
    carray bh(N), Er_n, d(N), bj(N), b(N), distance(N,N);
    carray c(N, N);

    float *x = space->x.host<float>();
    float *y = space->y.host<float>();
    real *Err = space->Er.r.host<float>();
    real *Eri = space->Er.i.host<float>();

    //offset matrix by 1 relative permability
    Er_n.r = space->Er.r - 1;
    Er_n.i = space->Er.i;
    Er_n.refresh();

    real a = sqrt(space->dx*space->dy / PI);
    comp scale = com * PI * k * a * (real)0.5;

    for (int i = 0; i < N; i++){
        comp Er0 = comp(Err[i], Eri[i]);

        //generate distances along x and y axis
        x2 = af::tile(space->x(i), space->x.elements()) - space->x;
        y2 = af::tile(space->y(i), space->y.elements()) - space->y;

        //generate actual distances
        p = af::pow(x2, 2) + af::pow(y2, 2);
        p = af::sqrt(p);

        //scale
        p = p * k;

        float b_j = bessj(1, k*a);
        bj.r = b_j;
        bj.i = 0;
        bj.refresh();

        //generate bessel fields
        bessj0(p, bh.r);
        bessy0(p, bh.i);
        bh.i = bh.i * -1;
        bh.refresh();

        carray temp(bj.a.elements());
        af::cfloat s(scale.real(), scale.imag());
        temp.a = bj.a * Er_n.a;
        temp.a = temp.a * bh.a;
        temp.a = temp.a * s;

        c.a.row(i) = temp.a;

        //this can be done outside the loop if wanted
        comp n = comp(1,0) + (Er0 - comp(1)) * com / (real)2.0 * PI * k * a * comp(bessj(1, a*k) - com * bessy(1, a*k)) - comp(2.0) * com;
        d.r(i) = n.real();
        d.i(i) = n.imag();

    }

    d.refresh();

    //combine c and d, d is the diaganol of c
    for (int i = 0; i < N; i++){
        c.a(i,i) = d.a(i);
    }

    for (int i = 0; i <N; i++){
        pho(i) = pow(x[i] - space->probes[probenum].x, 2) + pow(y[i] - space->probes[probenum].y, 2);
    }
    pho = af::sqrt(pho);

    real *_p = pho.host<real>();
    for (int i = 0; i < N; i++){
        comp n = com * k * (real)0.25 * comp(bessj(0,k*_p[i]) - com * bessy(0, k * _p[i]));
        b.r(i) = n.real();
        b.i(i) = n.imag();
    }

    af::freeHost(_p);
    b.refresh();
    Et.a = af::solve(c.a,b.a);

    //DELETE
    //    real *_Err = space->Er.r.host<real>();
    //    real *_Eri = Er.i.host<real>();
    real *_Etr = Et.r.host<real>();
    real *_Eti = Et.i.host<real>();

    for (size_t i = 0; i < space->probes.size(); i++){
        float x0 = space->probes[i].x;
        float y0 = space->probes[i].y;
        comp total = 0;

        for (int j = 0; j < N; j++){
            float p = sqrt(pow(x0 - x[j], 2) + pow(y0 - y[j], 2));
            comp er = comp(Err[j], Eri[j]);
            comp et = comp(_Etr[j], _Eti[j]);
            total += -com * PI * k / (real)2.0 * (er - (real)1.0) * et * a * comp(bessj(1, k*a)) * comp(bessj(0, k*p) - com * bessy(0, k*p));
        }
        Es.r(i) = total.real();
        Es.i(i) = total.imag();


    }

    //DELETE
    //    af::freeHost(_Err);
    //    af::freeHost(_Eri);
    af::freeHost(_Etr);
    af::freeHost(_Eti);
    af::freeHost(x);
    af::freeHost(y);
    af::freeHost(Eri);
    af::freeHost(Err);

}

void MOM::inverseBuilder(carray &Efunc, carray &C, real k)
{
    //void inverseMatrixBuilder(af::array &xc, af::array &yc, real dx, real dy, std::vector<Position> &probes, real k ,real EM, carray &Efunc, carray &C)

    assert(space->x.elements() == space->y.elements());

    int M = space->probes.size();
    int N = space->x.elements();
    C.resize(M, N);
    af::array p(N);

    real a = sqrt(space->dx*space->dy / PI);

    for (int i = 0; i < M; i++){
        real x0 = space->probes[i].x;
        real y0 = space->probes[i].y;

        p = af::pow(space->x - x0, 2) + af::pow(space->y - y0, 2);
        p = af::sqrt(p) * k;
        float bj1 = bessj(1, k*a);
        af::array bj0, by0;
        bessj0(p, bj0);
        bessy0(p, by0);

        C.a.row(i) = Efunc.a * af::cfloat(0, -PI * k * 0.5) * a * bj1 * (bj0 + af::cfloat(-1) * by0);
    }


}

void MOM::simulateSpace()
{
    Ez.resize(space->probes.size() * space->probes.size() * space->freqs.size());

    for (size_t l = 0; l < space->freqs.size(); l++){
        float k = wavenumber(space->freqs[l]);
        for (size_t i = 0; i < space->probes.size(); i++){
            mom(i, k);
            real *_Esr = Es.r.host<float>();
            real *_Esi = Es.i.host<float>();
            for (size_t j = 0; j < space->probes.size(); j++){
                comp e(_Esr[i], _Esi[i]);
                Ez.r(i*space->probes.size() + j) = e.real();
                Ez.i(i*space->probes.size() + j) = e.imag();
            }
            af::freeHost(_Esr);
            af::freeHost(_Esi);

        }
    }
    Ez.refresh();
}

void MOM::initialGuess()
{
    //simple initial guess, modify space    
    space->Er.r = af::constant(1.0, space->Er.a.dims());
    space->Er.i = af::constant(0.0, space->Er.a.dims());
    space->Er.refresh();
}

void MOM::iterateMom()
{
    std::vector<carray> computations;
    carray Ereg, Treg;
    carray inv;

    for (size_t f = 0; f < space->freqs.size(); f++){
        real k = wavenumber(space->freqs[f]);
        for (size_t i = 0; i < space->probes.size(); i++){
            mom(i, k);
            computations.push_back(Et);
        }
    }

    std::vector<carray> matrices;
    carray m, B, B2;
    for (size_t f = 0; f < space->freqs.size(); f++){
        float k = wavenumber(space->freqs[f]);
        for (size_t i = 0; i < space->probes.size(); i++){
            inverseBuilder(computations[i+space->probes.size() * f], m, k);
            matrices.push_back(m);
        }
    }

    B.resize(pow(space->probes.size(), 2) * space->freqs.size(), space->x.elements());

    for (size_t i = 0; i < space->probes.size(); i++){
        for (size_t j = 0; j < space->probes.size(); j++){
            for (dim_t k = 0; k < space->y.elements(); k++){
                B.a(i*space->probes.size() + j, k) = matrices[i].a(j,k);
            }
        }
    }

    real lambda = 100;

    int L = B.a.col(0).elements(); //rows
    int N = B.a.row(0).elements(); //columns
    B2.resize(N,N);

    B2.i = af::constant(0, B2.i.dims());
    B2.r = af::constant(0, B2.r.dims());

    for (int i = 0; i < N; i++)
        B2.r(i,i) = lambda;

    B2.refresh();


    B.a = af::join(0, B.a, B2.a);

    L = Ez.a.elements();
    Ereg.resize(N + L);
    for (int i = 0; i < L; i++){
        Ereg.a(i) = Ez.a(i);
    }
    for (int i = L; i < L + N; i++){
        Ereg.a(i) = 0;
    }

    //solve the system
    pinv(Ereg, Treg);

    Treg.a = af::matmul(Treg.a, B.a);

    af::cfloat n(1);
    space->Er.a = Treg.a + n; //correct for offset
}

void MOM::pinv(carray &A, carray &Ai)
{
    int minDim = std::min(A.a.row(0).elements(), A.a.col(0).elements());
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

real MOM::wavenumber(real freq)
{
    return sqrt(pow(2 * PI * freq, 2) * CU0 * CEM * CE0);
}
