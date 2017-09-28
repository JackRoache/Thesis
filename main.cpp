#include <QCoreApplication>
#include <QImage>
#include <QTime>
#include <QFile>
#include <QDir>

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

#include "phantomfile.h"
#include "tissueproperties.h"


#define CLX     0.20
#define CLY     CLX
#define CDX     2.5e-3
#define CDY     CDX
#define CDX2    CDX/2
#define CDY2    CDY/2
#define CNX     (int)(CLX/CDX)
#define CNY     (int)(CLY/CDY)
#define CFREQ   5e9

#define RX_XRADIUS  0.25
#define RX_YRADIUS  0.25
#define RX_ANGLE    30

#define PI      (float)(3.14159265)

#define NITER   5

int linear(int a, int b, float t)
{
    return a * (1 - t) + b * t;
}
QColor interpolate(QColor a, QColor b, float t)
{
    Q_ASSERT(t >= 0 && t <= 1);
    QColor final;

    int h = linear(a.hue(), b.hue(), t);
    int s = linear(a.saturation(), b.saturation(), t);
    int v = linear(a.value(), b.value(), t);

    final.setHsv(h,s,v);

    return final;
}
void surfToImage(const af::array &x, const af::array &y, const carray &Er, const QString &name, int cnx, int cny)
{
    QImage imReal(cnx, cny, QImage::Format_RGB32);
    QImage imImag(cnx, cny, QImage::Format_RGB32);
    float minR = 1000, maxR = -1000, minI = 1000, maxI = -1000;

    Q_ASSERT(cnx * cny == Er.a.elements());

    float * dataR = af::real(Er.a).host<float>();
    float * dataI = af::imag(Er.a).host<float>();

    for (dim_t i = 0; i < Er.a.elements(); i++){
        comp num = comp(dataR[i], dataI[i]);
        minR = MIN(minR, num.real());
        maxR = MAX(maxR, num.real());
        minI = MIN(minI, num.imag());
        maxI = MAX(maxI, num.imag());
    }

    uint32_t * dReal = (uint32_t *)imReal.bits();
    uint32_t * dImag  = (uint32_t*)imImag.bits();

    QColor first("yellow");
    QColor second("blue");

    for (dim_t i = 0; i < Er.a.elements(); i++){
        comp num = comp(dataR[i], dataI[i]);
        QColor im,re;

        if (minR != maxR){
            float sr = (num.real() - minR) / (maxR - minR);
            re = interpolate(first, second, sr);
            *dReal++ = re.rgb();
        }

        if (minI != maxI){
            float si =  (num.imag() - minI) / (maxI - minI);
            im = interpolate(first, second, si);
            *dImag++ = im.rgb();
        }
    }

//    af::freeHost(dataR);
//    af::freeHost(dataI);

    if (minI != maxI){
        QFile f1(name + QString("_imag.bmp"));
        f1.open(QIODevice::ReadWrite);
        imImag.save(&f1);
    }

    if (minR != maxR){
        QFile f2(name + QString("_real.bmp"));
        f2.open(QIODevice::ReadWrite);
        imReal.save(&f2);
    }
}

void iterationImage(RunInfo *info, ImagingSpace *space, carray &field, int iteration)
{
    std::ostringstream name;
    name << info->name << "/" << iteration;
    int nx = space->lx / space->dx;
    int ny = space->ly / space->dy;
    surfToImage(space->x, space->y, field, name.str().c_str(), nx, ny);
}

void generatePhantom(PhantomFile &file, int slice, int downsample, ImagingSpace &space, QMap<int, float> &permitivity, QMap<int, float> &conductivity, af::cfloat medium, float freq)
{
    std::cout << "Generating Phantom" << std::endl;
    QByteArray data = file.getSlice(slice);
    float x = space.lx / -2;
    float y = space.ly / -2;
    int size = file.getFrameSize() / downsample /downsample;
    space.x = af::array(size);
    space.y = af::array(size);
    space.Er = carray(size);

    int index = 0;
    for (int i = 0; i < file.getWidth(); i++){
        for (int j = 0; j < file.getHeight(); j++){\
            //relying on int rounding to get this right
            int h1 = i / downsample;
            int h2 = file.getWidth() / downsample;
            int l1 = j / downsample;
            int next = h1 * h2 + l1;

            space.x(next) = x;
            space.y(next) = y;

            int id = data[index];
            float val = permitivity.value(id, -1);
            Q_ASSERT(val != -1);
            if (val == 1)
                space.Er.r(next) = space.Er.r(next) + medium.real;
            else
                space.Er.r(next) = space.Er.r(next) + val;

            val = conductivity.value(id, -1);
            Q_ASSERT(val != -1);
            if (val == 0)
                space.Er.i(next) = space.Er.i(next) + medium.imag;
            else
                space.Er.i(next) = space.Er.i(next) + val;

            y+=space.dy / downsample;
            index++;
        }
        y = space.ly / -2;
        x += space.dx / downsample;
    }

//    af::print("x", space.x);
//    af::print("y", space.y);

    space.Er.refresh();
    float omega = 2 * PI * freq;
    space.Er.i = space.Er.i * -1 / freq;
    space.Er.a = space.Er.a / (downsample * downsample); //correct for downsampling
    space.Er.a = space.Er.a / medium;
    std::cout << "Finished Generating Phantom" << std::endl;
}

void generateProbes(ImagingSpace &space, float rx, float ry, int num)
{
    float angleDelta = 360 / num;

    for (int i = 0; i < num; i++){
        Position p;
        p.x = rx * sin(i * angleDelta * PI / 180);
        p.y = ry * cos(i * angleDelta * PI / 180);
        space.probes.push_back(p);
    }
}

void runBim(ImagingSpace &space, RunInfo &info)
{
    QDir dir;
    dir.mkdir(QString(info.name.c_str()));

    int nx = space.lx / space.dx;
    int ny = space.ly / space.dy;
    QString startName = QString(info.name.c_str()) + QString("/Start");
    QString endName = QString(info.name.c_str()) + QString("/Done");
    QString initialGuess = QString(info.name.c_str()) + QString("/Guess");
    surfToImage(space.x, space.y, space.Er, startName.toUtf8().data(), nx, ny);
    surfToImage(space.x, space.y, space.initalGuess, initialGuess.toUtf8().data(), nx, ny);
    std::cout << "Start BIM" << std::endl;
    MOM mom;
    mom.setCallBack(iterationImage);
    mom.setImagingSpace(&space);
    mom.setIterations(&info);

    mom.run();
    surfToImage(space.x, space.y, space.Er, endName.toUtf8().data(), nx, ny);
}

void simple_phantom()
{
    af::cfloat medium(250, 0);

    ImagingSpace space;
    space.dx = 1.1e-3 * 4;
    space.dy = 1.1e-3 * 4;
    space.lx = space.dx * 256 / 4;
    space.ly = space.dy * 256 / 4;
    space.freqs.push_back(1e9); //1GHz

    RunInfo info;
    info.iterations = 50;
    info.lambda = 0;

    PhantomFile phant("det_head_u2.med", 256, 256, 128);
    TissueProperties props;

    generatePhantom(phant, 64, 4, space, props.permitivity1ghz(), props.conductivity1ghz(), medium, 1e9);

    generateProbes(space, 0.2, 0.2, 12);
    carray phantom = space.Er;

    //simple initial guess
    space.initalGuess.r = af::constant(1, space.Er.a.dims());
    space.initalGuess.i = af::constant(0, space.Er.a.dims());
    space.initalGuess.refresh();
    info.name = "Phantom Constant";
    runBim(space, info);


    //Circular initial guess

    int nx = space.lx / space.dx;
    int ny = space.ly / space.ly;

    int index = 0;
    for (int i = 0; i < nx; i++){
        for (int j = 0; j < ny; j++){
            float x = i * space.lx / nx - space.lx / 2;
            float y = j * space.ly / ny - space.lx / 2;
            float r = std::sqrt(x*x + y*y);

            if (r < 7){
                space.initalGuess.r(index) = 860;
                space.initalGuess.i(index) = 0;
            } else {
                space.initalGuess.r(index) = 250;
                space.initalGuess.i(index) = 0;
            }
            index++;
        }
    }

    space.Er = phantom;
    space.initalGuess.refresh();
    info.name = "Phantom Round";
    runBim(space, info);


    //Perfect guess
    space.Er = phantom;
    space.initalGuess = phantom;
    info.name = "Phantom Perfect";
    runBim(space, info);
}

int main()
{
    af::setBackend(AF_BACKEND_OPENCL);
    af::info();
    std::cout << "Current Device " << af::getDevice() << std::endl;

    carray Er(CNX * CNY);

    af::array cx(CNX * CNY);
    af::array cy(CNX * CNY);
    carray Et, Es;

    std::vector<Position> probes;
    std::vector<float> freqs = {5e9};

    initKernels(cx);

    simple_phantom();
    return 0;

    //    Build up an imaging space
    int index = 0;
    for (int i = 0; i <CNX; ++i){
        for (int j = 0; j < CNY; ++j){
            float CX = CDX2 + i* CDX-CLX/2;
            float CY = CDY2 + j * CDY - CLY/2;

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

    surfToImage(cx, cy, space.Er, "Start", space.dx, space.dy);
    MOM mom;

    mom.setImagingSpace(&space);
    mom.setIterations(&info);

    mom.run();
    surfToImage(cx, cy, space.Er, "Done", space.dx, space.dy);

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
