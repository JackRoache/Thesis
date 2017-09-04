#include <QCoreApplication>
#include <QImage>
#include <QTime>
#include <QFile>

#include <stdio.h>
#include <algorithm>
#include <complex>
#include <cmath>
#include <string>
#include <sstream>
#include <iostream>
#include <assert.h>

#include "arrayfire.h"

#include "bessel.h"
#include "clbessel.h"

#include "types.h"



#define CE0     (real)(8.854e-12)
#define CU0     (real)(12.56e-7)
#define PI      (real)(3.14159265)

#define CLX     0.015
#define CLY     CLX
#define CDX     2.5e-3
#define CDY     CDX
#define CDX2    CDX/2
#define CDY2    CDY/2
#define CNX     (int)(CLX/CDX)
#define CNY     (int)(CLY/CDY)
#define CEM     1
#define CC      3e8
#define CFREQ   5e9

#define RX_XRADIUS  0.1
#define RX_YRADIUS  0.1
#define RX_ANGLE    30

#define NITER   5



struct carray
{
    carray(){}
    carray(int dim):
        r(dim),
        i(dim)
    {
        refresh();
    }

    carray(int dim1, int dim2):
        r(dim1, dim2),
        i(dim1, dim2)
    {
        refresh();
    }
    void refresh()
    {
        a = af::complex(r, i);
    }
    void resize(int dim)
    {
        r = af::array(dim);
        i = af::array(dim);
        refresh();
    }
    void resize(int dim1, int dim2){
        r = af::array(dim1, dim2);
        i = af::array(dim1, dim2);
        refresh();
    }


    af::array a;
    af::array r;
    af::array i;
};

void mom(af::array &_x, af::array &_y, real dx, real dy, carray &Er, real k, std::vector<Position> &probes, int prob_num, carray &Es, carray &Et)
{
    assert(_x.elements() == _y.elements());
    assert(_x.elements() == Er.a.elements());
    comp com(0,1);
    int N = _x.elements();

    Et.resize(N);
    Es.resize(probes.size());

    af::array p(N), x2(N), y2(N), pho(N);
    carray bh(N), Er_n, d(N), bj(N), b(N), distance(N,N);
    carray c(N, N);

    float *x = _x.host<float>();
    float *y = _y.host<float>();
    real *Err = Er.r.host<float>();
    real *Eri = Er.i.host<float>();

    //offset matrix by 1 relative permability
    Er_n.r = Er.r - 1;
    Er_n.i = Er.i;
    Er_n.refresh();

    real a = sqrt(dx*dy / PI);
    comp scale = com * PI * k * a * (real)0.5;

    for (int i = 0; i < N; i++){
        comp Er0 = comp(Err[i], Eri[i]);

        //generate distances along x and y axis
        x2 = af::tile(_x(i), _x.elements()) - _x;
        y2 = af::tile(_y(i), _y.elements()) - _y;

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
        pho(i) = pow(x[i] - probes[prob_num].x, 2) + pow(y[i] - probes[prob_num].y, 2);
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

    real *_Err = Er.r.host<real>();
    real *_Eri = Er.i.host<real>();
    real *_Etr = Et.r.host<real>();
    real *_Eti = Et.i.host<real>();

    for (size_t i = 0; i < probes.size(); i++){
        float x0 = probes[i].x;
        float y0 = probes[i].y;
        comp total = 0;

        for (int j = 0; j < N; j++){
            float p = sqrt(pow(x0 - x[j], 2) + pow(y0 - y[j], 2));
            comp er = comp(_Err[j], _Eri[j]);
            comp et = comp(_Etr[j], _Eti[j]);
            total += -com * PI * k / (real)2.0 * (er - (real)1.0) * et * a * comp(bessj(1, k*a)) * comp(bessj(0, k*p) - com * bessy(0, k*p));
        }
        Es.r(i) = total.real();
        Es.i(i) = total.imag();


    }

    af::freeHost(_Err);
    af::freeHost(_Eri);
    af::freeHost(_Etr);
    af::freeHost(_Eti);
    af::freeHost(x);
    af::freeHost(y);
    af::freeHost(Eri);
    af::freeHost(Err);
}

void surfToImage(const af::array &x, const af::array &y, const carray &Er, const QString &name)
{
    QImage imReal(CNX, CNY, QImage::Format_RGB32);
    QImage imImag(CNX, CNY, QImage::Format_RGB32);
    float minR = 1000, maxR = -1000, minI = 1000, maxI = -1000;

//    real * dataR = Er.r.host<real>();
//    real * dataI = Er.i.host<real>();

    real * dataR = af::real(Er.a).host<real>();
    real * dataI = af::imag(Er.a).host<real>();

    for (dim_t i = 0; i < Er.a.elements(); i++){
        comp num = comp(dataR[i], dataI[i]);
        minR = MIN(minR, num.real());
        maxR = MAX(maxR, num.real());
        minI = MIN(minI, num.imag());
        maxI = MAX(maxI, num.imag());
    }

    uint32_t * dReal = (uint32_t *)imReal.bits();
    uint32_t * dImag  = (uint32_t*)imImag.bits();

    for (dim_t i = 0; i < Er.a.elements(); i++){
        uint32_t im,re;
        comp num = comp(dataR[i], dataI[i]);
        im = 255 * (num.imag() - minI) / (maxI - minI);
        re = 255 * (num.real() - minR) / (maxR - minR);
        im |= im << 8;
        im |= im << 8;
        re |= re << 8;
        re |= re << 8;

        *dReal++ = re;
        *dImag++ = im;
    }

    af::freeHost(dataR);
    af::freeHost(dataI);

    QFile f1(name + QString("_imag.bmp"));
    f1.open(QIODevice::ReadWrite);
    imImag.save(&f1);

    QFile f2(name + QString("_real.bmp"));
    f2.open(QIODevice::ReadWrite);
    imReal.save(&f2);
}

void inverseMatrixBuilder(af::array &xc, af::array &yc, real dx, real dy, std::vector<Position> &probes, real k ,real EM, carray &Efunc, carray &C)
{
    assert(xc.elements() == yc.elements());

    int M = probes.size();
    int N = xc.elements();
    C.resize(M, N);
    af::array p(N);

    real a = sqrt(dx*dy / PI);

    for (int i = 0; i < M; i++){
        real x0 = probes[i].x;
        real y0 = probes[i].y;

        p = af::pow(xc - x0, 2) + af::pow(yc - y0, 2);
        p = af::sqrt(p) * k;
        float bj1 = bessj(1, k*a);
        af::array bj0, by0;
        bessj0(p, bj0);
        bessy0(p, by0);

        C.a.row(i) = Efunc.a * af::cfloat(0, -PI * k * 0.5) * a * bj1 * (bj0 + af::cfloat(-1) * by0);
    }

}

void pinv(carray &A, carray &Ai)
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

float waveNumber(float freq)
{
    return sqrt(pow(2 * PI * freq, 2) * CU0 * CEM * CE0);
}

void testBessel()
{

    int size = 100;
    af::array x(size);
    af::array y1(size), y2(size);
    af::array j1(size), j2(size);
    af::array k1(size), k2(size);

    for (int i = 1; i < size + 1; i++)
        x(i - 1) = i;

    float *_x = x.host<float>();

    for (int i = 0; i < size; i++){
        y2(i) = bessj(0, _x[i]);
    }

    bessj0(x, y1);

    float sum = af::sum<float>(af::abs(y1 - y2));
    std::cout << "bessel j0 sum " << sum << std::endl;

    for (int i = 0; i < size; i++){
        j2(i) = bessj(1, _x[i]);
    }

    bessj1(x, j1);

    af::sync();

    sum = af::sum<float>(af::abs(j1 - j2));
    std::cout << "bessel j1 sum " << sum << std::endl;

    for (int i = 0; i < size; i++){
        k2(i) = bessy(0, _x[i]);
    }

    bessy0(x, k1);

    sum = af::sum<float>(af::abs(k1-k2));
    std::cout << "bessel y0 sum " << sum << std::endl;

    af::freeHost(_x);
}

void testpinv()
{
    carray a(2,3);
    carray ai;
    for (int i = 0; i < 3; i++){
        a.r(0,i) = i + 1;
        a.i(0,i) = i + 2;
    }

    for (int i = 0; i < 3; i++){
        a.r(1,i) = i + 4;
        a.i(1,i) = i + 6;
    }

    a.refresh();
    assert(a.a.iscomplex());

    af::print("a", a.a);
    pinv(a, ai);
    af::print("ai", ai.a);
    af::print("AAA", af::matmul(a.a, ai.a, a.a));
}

int main()
{
    af::info();

    carray Er(CNX * CNY);

    af::array cx(CNX * CNY);
    af::array cy(CNX * CNY);
    carray Et, Es;

    std::vector<Position> probes;
    std::vector<real> freqs = {5e9, 50e9, 100e9};

    initKernels(cx);

    //    Build up an imaging space
    int index = 0;
    for (int i = 0; i <CNX; ++i){
        for (int j = 0; j < CNY; ++j){
            real CX = CDX2 + i* CDX-CLX/2;
            real CY = CDY2 + j * CDY - CLY/2;

            cx(index) = CX;
            cy(index) = CY;

            if ((CX > -0.02 + CDX2)&&(CX < 0.02 - CDX2) && (CY > -0.02 + CDY2) && (CY < 0.02 - CDY2)){
                Er.r(index) = 1.5 * CEM;
                Er.i(index) = -0.2;
            } else {
                Er.r(index) = 1;
                Er.i(index) = 0;
            }
            index++;
        }
    }

    //Generate Probes
    int numProbes = 360 / RX_ANGLE;

    for (int i = 0; i < numProbes; i++){
        Position p;
        p.x = RX_XRADIUS * sin(i * RX_ANGLE * PI / 180);
        p.y = RX_YRADIUS * cos(i * RX_ANGLE * PI / 180);
        probes.push_back(p);
    }
    Er.refresh();
    surfToImage(cx, cy, Er, "source");


    //Generate initial space / What would be recorded
    carray Ez(pow(probes.size(), 2) * freqs.size());

    af::timer tim1= af::timer::start();
    std::cout << "start mom 1" << std::endl;
    for (size_t l = 0; l < freqs.size(); l++){
        float k = waveNumber(freqs[l]);
        for (size_t i = 0; i < probes.size(); i++){
            mom(cx, cy, CDX, CDY, Er, k, probes, i, Es, Et);
            real *_Esr = Es.r.host<float>();
            real *_Esi = Es.i.host<float>();
            for (size_t j = 0; j < probes.size(); j++){
                comp e(_Esr[i], _Esi[i]);
                Ez.r(i*probes.size() + j) = e.real();
                Ez.i(i*probes.size() + j) = e.imag();
            }
            af::freeHost(_Esr);
            af::freeHost(_Esi);
            Ez.refresh();
        }
    }
    std::cout << "end mom 1" << af::timer::stop(tim1) << " s" << std::endl;

    //Initial Guess of space. Assume all 1
    af::cfloat c(1);
    Er.a = af::constant(c, Er.a.dims());

    surfToImage(cx, cy, Er, "Empty space");

    for (int iter = 0; iter < NITER; iter++){
        std::cout << "Iteration " << iter + 1  << std::endl;
        af::timer::start();
        std::vector<carray> computations;
        carray Ereg, Treg;
        carray inv;

        for (size_t f = 0; f < freqs.size(); f++){
            real k = waveNumber(freqs[f]);
            for (size_t i = 0; i < probes.size(); i++){
                mom(cx, cy, CDX, CDY, Er, k, probes, i, Es, Et);
                computations.push_back(Et);
            }
        }

        std::vector<carray> matrices;
        carray m, B, B2;
        for (size_t f = 0; f < freqs.size(); f++){
            float k = waveNumber(freqs[f]);
            for (size_t i = 0; i < probes.size(); i++){
                inverseMatrixBuilder(cx, cy, CDX, CDY, probes, k, CEM, computations[i+probes.size() * f], m);
                matrices.push_back(m);
            }
        }

        B.resize(pow(probes.size(), 2) * freqs.size(), cx.elements());

        for (size_t i = 0; i < probes.size(); i++){
            for (size_t j = 0; j < probes.size(); j++){
                for (dim_t k = 0; k < cy.elements(); k++){
                    B.a(i*probes.size() + j, k) = matrices[i].a(j,k);
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
        Er.a = Treg.a + n; //correct for offset

        std::cout << "Time " << af::timer::stop() << std::endl;
        char str[40];
        sprintf(str, "it %d", iter);
        surfToImage(cx, cy, Er, str);
    }

    surfToImage(cx, cy, Er, "final");

    return 1;
}




/*
 * S parameter data
 * Data from  a full wave simulator
 *
 * Pre seeding data
 *
 * Can it make it better
 * Does it screw it up
 * Iteration speed up??
 * realistic seeding
 * average dialetric or approx internal structure??
 * How accurate does the interanl seeding need to be??
 *
 * How fine of a resolution do you need?? (Voxels)
 */


























