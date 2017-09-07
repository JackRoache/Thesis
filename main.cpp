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
#include "mom.h"


#define CLX     0.075
#define CLY     CLX
#define CDX     2.5e-3
#define CDY     CDX
#define CDX2    CDX/2
#define CDY2    CDY/2
#define CNX     (int)(CLX/CDX)
#define CNY     (int)(CLY/CDY)
#define CFREQ   5e9

#define RX_XRADIUS  0.1
#define RX_YRADIUS  0.1
#define RX_ANGLE    30

#define PI      (real)(3.14159265)

#define NITER   5


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

//void testBessel()
//{

//    int size = 100;
//    af::array x(size);
//    af::array y1(size), y2(size);
//    af::array j1(size), j2(size);
//    af::array k1(size), k2(size);

//    for (int i = 1; i < size + 1; i++)
//        x(i - 1) = i;

//    float *_x = x.host<float>();

//    for (int i = 0; i < size; i++){
//        y2(i) = bessj(0, _x[i]);
//    }

//    bessj0(x, y1);

//    float sum = af::sum<float>(af::abs(y1 - y2));
//    std::cout << "bessel j0 sum " << sum << std::endl;

//    for (int i = 0; i < size; i++){
//        j2(i) = bessj(1, _x[i]);
//    }

//    bessj1(x, j1);

//    af::sync();

//    sum = af::sum<float>(af::abs(j1 - j2));
//    std::cout << "bessel j1 sum " << sum << std::endl;

//    for (int i = 0; i < size; i++){
//        k2(i) = bessy(0, _x[i]);
//    }

//    bessy0(x, k1);

//    sum = af::sum<float>(af::abs(k1-k2));
//    std::cout << "bessel y0 sum " << sum << std::endl;

//    af::freeHost(_x);
//}

//void testpinv()
//{
//    carray a(2,3);
//    carray ai;
//    for (int i = 0; i < 3; i++){
//        a.r(0,i) = i + 1;
//        a.i(0,i) = i + 2;
//    }

//    for (int i = 0; i < 3; i++){
//        a.r(1,i) = i + 4;
//        a.i(1,i) = i + 6;
//    }

//    a.refresh();
//    assert(a.a.iscomplex());

//    af::print("a", a.a);
//    pinv(a, ai);
//    af::print("ai", ai.a);
//    af::print("AAA", af::matmul(a.a, ai.a, a.a));
//}

int main()
{
    af::info();

    carray Er(CNX * CNY);

    af::array cx(CNX * CNY);
    af::array cy(CNX * CNY);
    carray Et, Es;

    std::vector<Position> probes;
    std::vector<real> freqs = {5e9};

    initKernels(cx);

    //    Build up an imaging space
    int index = 0;
    for (int i = 0; i <CNX; ++i){
        for (int j = 0; j < CNY; ++j){
            real CX = CDX2 + i* CDX-CLX/2;
            real CY = CDY2 + j * CDY - CLY/2;

            cx(index) = CX;
            cy(index) = CY;

            if ((CX > -0.02)&& (CX < 0.02) && (CY > -0.02)&&(CY < 0.02)){
                Er.r(index) = 1.5 ;
                Er.i(index) = -0.1;
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

    float _t[] = {1,2,3};
    af::array t(3, 1, _t);
    af::print("t", t);
    af::print("tiled", af::tile(t, 1, 3));
    af::print("transposed", af::tile(af::transpose(t), 3, 1));
    af::print("operation", af::tile(af::transpose(t), 3, 1) - af::tile(t, 1, 3));

    //init mom structs

    ImagingSpace space;
    space.x = cx;
    space.y = cy;
    space.Er = Er;
    space.dx = CDX;
    space.dy = CDY;
    space.lx = CLX;
    space.ly - CLY;
    space.probes = probes;
    space.freqs = freqs;

    RunInfo info;
    info.name = "Basic";
    info.iterations = 3;

    surfToImage(cx, cy, space.Er, "Start");
    MOM mom;

    mom.setImagingSpace(&space);
    mom.setIterations(&info);

    mom.run();
    surfToImage(cx, cy, space.Er, "Done");

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
