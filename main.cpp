#include <QCoreApplication>
#include <QImage>
#include <QTime>
#include <QFile>
#include <QDir>
#include <QDateTime>
#include <QDebug>

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

#include "Complex_Bessel/complex_bessel.h"

#define CE0     (8.854e-12)

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
void surfToImage(const af::array &x, const af::array &y, const carray &Er, const QString &name, int cnx, int cny, ImagingSpace *space, bool firstRun)
{
    QImage imReal(cnx, cny, QImage::Format_RGB32);
    QImage imImag(cnx, cny, QImage::Format_RGB32);
    float minR = 1000, maxR = -1000, minI = 1000, maxI = -1000;

    Q_ASSERT(cnx * cny == Er.a.elements());
    uint32_t * dReal = (uint32_t *)imReal.bits();
    uint32_t * dImag  = (uint32_t*)imImag.bits();
    float * dataR = af::real(Er.a).host<float>();
    float * dataI = af::imag(Er.a).host<float>();

    //    if (firstRun){
    for (dim_t i = 0; i < Er.a.elements(); i++){
        comp num = comp(dataR[i], dataI[i]);
        minR = MIN(minR, num.real());
        maxR = MAX(maxR, num.real());
        minI = MIN(minI, num.imag());
        maxI = MAX(maxI, num.imag());
    }
    //        space->maxImag = maxI;
    //        space->minImag = minI;
    //        space->maxReal = maxR;
    //        space->minReal = minR;
    //    } else {
    //        minR = space->minReal;
    //        maxR = space->maxReal;
    //        minI = space->minImag;
    //        maxI = space->maxImag;
    //    }

    QList<QColor> map;


    QColor first("black");
    QColor second("red");
    QColor third("white");
    map << first << second << third;

    for (dim_t i = 0; i < Er.a.elements(); i++){
        comp num = comp(dataR[i], dataI[i]);
        QColor im,re;

        if (minR != maxR){
            double sr = (num.real() - minR) / (maxR - minR);
            sr = std::min(std::max(0.0, sr), 1.0);
            int lower = sr * map.size();
            lower = std::min(std::max(lower, 0), map.size()-2);
            re = interpolate(map[lower], map[lower+1], sr);
            *dReal++ = re.rgb();
        }

        if (minI != maxI){
            double si =  (num.imag() - minI) / (maxI - minI);
            si = std::min(std::max(0.0, si), 1.0);
            int lower = si * map.size();
            lower = std::min(std::max(lower, 0), map.size()-2);
            im = interpolate(map[lower], map[lower+1], si);
            *dImag++ = im.rgb();
        }
    }

    //    af::freeHost(dataR);
    //    af::freeHost(dataI);

    std::cout << "Real " << minR << " " << maxR << std::endl;
    std::cout << "Imag " << minI << " " << maxI << std::endl;

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
    surfToImage(space->x, space->y, field, name.str().c_str(), nx, ny, space, false);

    //    info->window->surface(space->x, space->y, af::imag(field.a), name.str().c_str());
}

void generatePhantom(PhantomFile &file, int slice, int downsample, int border, ImagingSpace &space, TissueProperties &props, float freq)
{
    std::cout << "Generating Phantom" << std::endl;
    QByteArray data = file.getSlice(slice);
    float omega = 2 * PI * freq;
    float x = space.lx / -2;
    float y = space.ly / -2;
    int size = (file.getWidth() - 2 * border) * (file.getHeight() - 2 * border) / downsample / downsample;
    space.x = af::array(size);
    space.y = af::array(size);
    space.Er = carray(size);
    space.Er.r = af::constant(0, space.Er.a.dims());
    space.Er.i = af::constant(0, space.Er.a.dims());

    af::cfloat medium;
    {
        std::complex<double> cond = 0;
        std::complex<double> imag(0,1);
        cond = space.medium_cond / (omega * CE0 * imag);
        medium = af::cfloat(space.medium_es + cond.real(), -space.medium_es + cond.imag());
        medium.imag = 0;//40.0;
    }

    for (int i = border; i < file.getWidth() - border; i++){
        for (int j = border; j < file.getHeight() - border; j++){\
            //relying on int rounding to get this right
            int h1 = (i - border) / downsample;
            int h2 = (file.getWidth() - 2 * border) / downsample;
            int l1 = (j - border) / downsample;
            int next = h1 * h2 + l1;
            int index = i * file.getWidth() + j;

            space.x(next) = x;
            space.y(next) = y;

            double Einf =0;
            double Es_cole[4] = {0};
            double Tau[4] = {0};
            double alpha[4] = {0};
            double Cond_cole = 0;
            int K = 0;

            int id = data[index];
            switch (id) {
            //Background, handled below
            case 0:
                break;

                //Skin
            case 1:
                Einf = 4;
                Es_cole[0] = 32;
                Es_cole[1] = 1100;
                Tau[0] = 7.23e-12;
                Tau[1] = 32.48e-9;
                alpha[0] = 0;
                alpha[1] = 0.2;
                Cond_cole = 0.0002;
                K  = 2;
                break;

                //Skull
            case 4:
            case 5:
            case 9:
            case 70:
            case 71:
            case 76:
            case 81:
            case 100:
            case 102:
            case 125:
                Einf = 4;
                Es_cole[0] = 8;
                Es_cole[1] = 4;
                Tau[0] = 17e-12;
                Tau[1] = 0.46e-9;
                alpha[0] = 0;
                alpha[1] = 0;
                Cond_cole = 0.065;
                K  = 2;
                break;

                // Fat
            case 22:
            case 98:
                Einf  = 3;
                Es_cole[0] = 1.42;
                Es_cole[1] = 1.87;
                Tau[0] = 13e-12;
                Tau[1] = 0.651e-9;
                alpha[0] =0;
                alpha[1] = 0;
                Cond_cole = 0.026;
                K = 2;
                break;

                //Dura
            case 75:
            case 113:
                Einf  = 4;
                Es_cole[0] = 40;
                Es_cole[1] = 200;
                Es_cole[2] = 1e4;
                Es_cole[3] = 1e6;
                Tau[0] =7.958e-12;
                Tau[1]=7.958e-9;
                Tau[2]=159.15e-6;
                Tau[3]=15.915e-3;
                alpha[0] =0.15;
                alpha[1] =0.1;
                alpha[2] = 0.2;
                alpha[3] =0;
                Cond_cole = 0.5;
                K = 4;
                break;

                //Brain (Gray Matter)
            case 118:
            case 101:
            case 89:
            case 120:
            case 96:
            case 103:
            case 95:
            case 117:
            case 124:
            case 105:
            case 112:
            case 114:
            case 109:
            case 99:
                Einf = 1;
                Es_cole[0] =49;
                Es_cole[1] = 55;
                Tau[0] = 10e-12;
                Tau[1] = 1.3e-9;
                Cond_cole = 0.5;
                alpha[0] =0;
                alpha[1] = 0;
                K  = 2;
                break;
                //White Matter
            case 111:
            case 107:
            case 108:
            case 83:
                Einf = 8;
                Es_cole[0] = 29;
                Es_cole[1] = 26;
                Tau[0] = 13e-12;
                Tau[1] = 0.94e-9;
                alpha[0] = 0;
                alpha[1] = 0;
                Cond_cole = 0.3;
                K = 2;
                break;
                //CSF
            case 122:
            case 2:
            case 115:
            case 123:
            case 92:
                Einf = 1;
                Es_cole[0] = 66;
                Es_cole[1] = 17;
                Tau[0] =6.9e-12;
                Tau[1] = 1.7e-9;
                alpha[0] =0;
                alpha[1] = 0;
                Cond_cole = 2.2;
                K = 2;
                break;
                //Blood
            case 23:
            case 84:
                Einf = 8;
                Es_cole[0] = 47;
                Es_cole[1] = 23;
                Tau[0] =10.8e-12;
                Tau[1] = 1.2e-9;
                alpha[0] =0;
                alpha[1] = 0;
                Cond_cole = 1.5;
                K = 2;
                break;
                //            case 9:
                //                Einf = 4;
                //                Es_cole[0] =50;
                //                Es_cole[1] = 7000;
                //                Es_cole[2] = 1.2e6;
                //                Es_cole[3] = 2.5e7;
                //                Tau[0] = 7.234e-12;
                //                Tau[1] = 353.678e-9;
                //                Tau[2] = 313.310e-6;
                //                Tau[3] = 2.264e-3;
                //                alpha[0] = 0.1;
                //                alpha[1] = 0.1;
                //                alpha[2] = 0.1;
                //                alpha[3] = 0;
                //                Cond_cole = 0.2;
                //                K =4 ;
                //                break

                //Bone Marrow
            case 26:
                Einf = 2.5;
                Es_cole[0] = 9;
                Es_cole[1] = 80;
                Es_cole[2] = 1e4;
                Es_cole[3] = 2e6;
                Tau[0] = 14.469e-12;
                Tau[1] = 15.91e-9;
                Tau[2] = 1591.549e-6;
                Tau[3] = 15.915e-3;
                alpha[0] = 0.2;
                alpha[1] = 0.1;
                alpha[2] = 0.1;
                alpha[3] = 0.1;
                Cond_cole = 0.1;
                K =4 ;
                break;

            default:
                std::cout << "ID " << id << std::endl;
                assert(false);
            }

            std::complex<double> temp = 0;
            std::complex<double> cond = 0;
            std::complex<double> imag(0,1);
            for (int j = 0; j < K; j++){
                temp += Es_cole[j] / (1.0 + omega * std::pow(Tau[j], 1-alpha[j]) * imag);
            }

            cond = Cond_cole / (omega * CE0 * imag);

            if (id != 0){
                space.Er.r(next) += Einf + temp.real() + cond.real();
                space.Er.i(next) += temp.imag() + cond.imag();
            } else {
                space.Er.r(next) += medium.real;
                space.Er.i(next) += medium.imag;
            }


            y+=space.dy / downsample;
        }
        y = space.ly / -2;
        x += space.dx / downsample;
    }


    space.Er.refresh();
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
    surfToImage(space.x, space.y, space.Er, startName.toUtf8().data(), nx, ny, &space, true);
    surfToImage(space.x, space.y, space.initalGuess, initialGuess.toUtf8().data(), nx, ny, &space, false);
    std::cout << "Start BIM " << info.name << std::endl;
    MOM mom;
    mom.setCallBack(iterationImage);
    mom.setImagingSpace(&space);
    mom.setIterations(&info);

    mom.run();
    surfToImage(space.x, space.y, space.Er, endName.toUtf8().data(), nx, ny, &space, false);
}

void simple_phantom(int freq, double lambda, af::Window &window)
{
    int border = 40;
    int downsample = 4;

    ImagingSpace space;
    space.dx = 1.1e-3 * downsample;
    space.dy = 1.1e-3 * downsample;
    space.lx = space.dx * (256 - border * 2) / downsample;
    space.ly = space.dy * (256 - border * 2) / downsample;
    space.freqs.push_back(freq);
//    space.medium = af::cfloat(30, 0);
    space.medium_es = 40;
    space.medium_cond = 0;

    RunInfo info;
    info.iterations = 5;
    info.lambda = lambda;
    info.window = &window;
    info.slow = false;

    PhantomFile phant("det_head_u2.med", 256, 256, 128);
    TissueProperties props;

    generatePhantom(phant, 30, downsample, border, space, props, freq);

    generateProbes(space, 0.11, 0.09, 36);
    carray phantom = space.Er;

    //simple initial guess
    space.initalGuess.r = af::constant(1, space.Er.a.dims());
    space.initalGuess.i = af::constant(0, space.Er.a.dims());
    space.initalGuess.refresh();
    {
        std::stringstream ss;
        ss << "Phantom Constant " << freq << "hz " << lambda << " lambda";
        info.name = ss.str();
    }
    runBim(space, info);

    //Circular initial guess
    int nx = space.lx / space.dx;
    int ny = space.ly / space.dy;

    int index = 0;
    for (int i = 0; i < nx; i++){
        for (int j = 0; j < ny; j++){
            float x = i * space.lx / nx - space.lx / 2;
            float y = j * space.ly / ny - space.lx / 2;
            float r = std::sqrt(x*x + y*y);

            if (r < 0.07){
                space.initalGuess.r(index) = 1;
                space.initalGuess.i(index) = 0;
            } else {
                space.initalGuess.r(index) = 0;
                space.initalGuess.i(index) = 0;
            }
            index++;
        }
    }

    space.Er = phantom;
    space.initalGuess.refresh();
    {
        std::stringstream ss;
        ss << "Phantom Round " << freq << "hz " << lambda <<" lambda";
        info.name = ss.str();
    }
//    runBim(space, info);


    //Perfect guess
    space.Er = phantom;
    space.initalGuess = phantom;
    {
        std::stringstream ss;
        ss << "Phantom Perfect "<< freq << "hz " << lambda <<" lambda";
        info.name = ss.str();
    }
//    runBim(space, info);
}

void very_simple()
{
    ImagingSpace space;
    space.dx = 5e-3;
    space.dy = 5e-3;
    space.lx = 1e-3 * 75;
    space.ly = 1e-3 * 75;
    space.freqs.push_back(5e9);
//    space.medium = af::cfloat(1,0);
    space.medium_es = 1;
    space.medium_cond = 0;

    RunInfo info;
    info.iterations = 5;
    info.lambda = 0.1;
    info.name  = "very simple";
    info.slow = false;

    int nx = space.lx / space.dx;
    int ny = space.ly / space.dy;

    space.x = af::array(nx * ny);
    space.y = af::array(nx * ny);
    space.Er.resize(nx * ny);

    std::cout << nx << " " << ny << std::endl;

    int index = 0;
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++){
            double Xc = space.dx/2 + i * space.dx - space.lx / 2;
            double Yc = space.dy/2 + j * space.dy - space.ly / 2;
            space.x(index) = Xc;
            space.y(index) = Yc;

            if ((Xc > -0.02) && (Xc < 0.02) && (Yc > -0.02) && (Yc < 0.02)){
                space.Er.r(index) = 1.5;
                space.Er.i(index) = -0.1;
            } else {
                space.Er.r(index) = 1;
                space.Er.i(index) = 0;
            }
            index++;
        }
    }
    assert(index == nx * ny);
    space.Er.refresh();
    af::print("Er", space.Er.a);
    generateProbes(space, 0.05, 0.05, 12);
    space.initalGuess.r = af::constant(1.0, space.Er.a.dims());
    space.initalGuess.i = af::constant(0.0, space.Er.a.dims());
    space.initalGuess.refresh();
    runBim(space, info);
}

void bessel()
{
    std::complex<double> comp(0,1);
    std::cout << sp_bessel::besselJ(0, comp) << std::endl;
    std::cout << sp_bessel::besselY(0, comp) << std::endl;
    std::cout << sp_bessel::hankelH1(0, comp) << std::endl;
}

int main()
{
    af::setBackend(AF_BACKEND_OPENCL);
    af::info();
    std::cout << "Current Device " << af::getDevice() << std::endl;
    af::array dummy(1);
    initKernels(dummy);

    af::Window window(512, 512, "Window!");

    bessel();

    //    very_simple();

//    simple_phantom(500000000, 0.07, window);
//    simple_phantom(500000000, 0.08, window);
//    simple_phantom(500000000, 0.09, window);
//    simple_phantom(500000000, 0.11, window);
//    simple_phantom(500000000, 0.12, window);
//    simple_phantom(500000000, 0.13, window);

    qDebug() << QDateTime::currentDateTime();
    uint64_t epoch = QDateTime::currentMSecsSinceEpoch();
//    simple_phantom(500000000, 0.1, window);
//    simple_phantom(500000000, 0.01, window);
//    simple_phantom(500000000, 0, window);
//    simple_phantom(500000000, 0.2, window);
//    simple_phantom(500000000, 0.5, window);
//    simple_phantom(500000000, 1, window);

//    simple_phantom(300000000, 0.1, window);
//    simple_phantom(400000000, 0.1, window);
//    simple_phantom(600000000, 1.1, window);
    simple_phantom(700000000, 1.1, window);
    simple_phantom(800000000, 0.1, window);
    simple_phantom(900000000, 0.1, window);
    simple_phantom(1000000000, 0.1, window);
    simple_phantom(1100000000, 0.1, window);
    simple_phantom(1200000000, 0.1, window);
    simple_phantom(1300000000, 0.1, window);
    simple_phantom(1400000000, 0.1, window);
    simple_phantom(1500000000, 0.1, window);
    simple_phantom(1600000000, 0.1, window);
    simple_phantom(1700000000, 0.1, window);
    simple_phantom(1800000000, 0.1, window);
    simple_phantom(1900000000, 0.1, window);
    simple_phantom(2000000000, 0.1, window);
    qDebug() << QDateTime::currentDateTime();
    qDebug() << epoch - QDateTime::currentMSecsSinceEpoch();


    //    simple_phantom(500000000, 0.2);
    //    simple_phantom(600000000, 0.2);

    //    simple_phantom(500000000, 0.5);
    //    simple_phantom(600000000, 0.5);

    //    simple_phantom(500000000, 0.01);
    //    simple_phantom(600000000, 0.01);

    return 0;
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
