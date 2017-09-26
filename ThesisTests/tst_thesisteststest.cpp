#include <QString>
#include <QtTest>
#include <QImage>
#include <QFile>

#include "../types.h"
#include "../mom.h"
#include "../clbessel.h"

#define PI      (real)(3.14159265)

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

void surfToImage(const af::array &x, const af::array &y, const carray &Er, const QString &name, int CNX, int CNY)
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

    af::freeHost(dataR);
    af::freeHost(dataI);

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

void iterationImage(RunInfo *info, ImagingSpace *space, carray field, int iteration)
{
    std::ostringstream name;
    name << info->name << "/" << iteration;
    int nx = space->lx / space->dx;
    int ny = space->ly / space->dy;
    surfToImage(space->x, space->y, field, name.str().c_str(), nx, ny);
}

void runBIM(RunInfo &info, ImagingSpace &space, int CNX, int CNY)
{
    QDir ().mkdir(info.name.c_str());

    surfToImage(space.x, space.y, space.Er, (info.name + "/start").c_str(), CNX, CNY);
    surfToImage(space.x, space.y, space.initalGuess, (info.name + "/intial guess").c_str(), CNX, CNY);
    MOM mom;

    mom.setImagingSpace(&space);
    mom.setIterations(&info);
    mom.setCallBack(iterationImage);

    mom.run();

    surfToImage(space.x, space.y, space.Er, (info.name + "/Done").c_str(), CNX, CNY);
}

class ThesisTestsTest : public QObject
{
    Q_OBJECT

public:
    ThesisTestsTest();

private Q_SLOTS:
    void testCase1();
    void testCase2();
};

ThesisTestsTest::ThesisTestsTest()
{
    af::array dummy(10);
    af::info();
    initKernels(dummy);
}

void ThesisTestsTest::testCase1()
{
    real CLX = 0.075;
    real CLY = CLX;
    real CDX = 0.0025;
    real CDY = CDX;
    real CDX2 = CDX/2;
    real CDY2 = CDY/2;
    int CNX = CLX/CDX;
    int CNY = CLY/CDY;

    real RX_XRADIUS = 0.1;
    real RX_YRADIUS = 0.1;
    real RX_ANGLE =30;

    carray Er(CNX * CNY);
    af::array cx(CNX * CNY);
    af::array cy(CNX * CNY);

    std::vector<Position> probes;
    std::vector<real> freqs = {5e9};

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

    //init mom structs
    ImagingSpace space;
    space.x = cx;
    space.y = cy;
    space.Er = Er;
    space.dx = CDX;
    space.dy = CDY;
    space.lx = CLX;
    space.ly = CLY;
    space.probes = probes;
    space.freqs = freqs;

    space.initalGuess.r = af::constant(1, Er.a.dims());
    space.initalGuess.i = af::constant(0, Er.a.dims());
    space.initalGuess.refresh();

    RunInfo info;
    info.name = "Basic";
    info.iterations = 10;
    info.lambda = 10;

    runBIM(info, space, CNX, CNY);
}

void ThesisTestsTest::testCase2()
{
    std::string name("arm");
    real CLX = 0.1;
    real CLY = CLX;
    real CDX = 0.0025;
    real CDY = CDX;
    real CDX2 = CDX/2;
    real CDY2 = CDY/2;
    int CNX = CLX/CDX;
    int CNY = CLY/CDY;

    real RX_XRADIUS = 0.1;
    real RX_YRADIUS = 0.1;
    real RX_ANGLE = 30;

    carray Er(CNX * CNY);
    af::array cx(CNX * CNY);
    af::array cy(CNX * CNY);

    std::vector<Position> probes;
    std::vector<real> freqs = {1e9};

    int index = 0;
    //generate imaging space
    for (int i = 0; i <CNX; ++i){
        for (int j = 0; j < CNY; ++j){
            real CX = CDX2 + i* CDX-CLX/2;
            real CY = CDY2 + j * CDY - CLY/2;

            real distance = std::sqrt(std::pow(CX, 2) + std::pow(CY, 2));

            Er.r(index) = 1;
            Er.i(index) = 0;
            //Fat
            if (distance <= 0.045){
                Er.r(index) = 27.222;
                Er.i(index) = 0.025079;
            }

            //Muscle
            if (distance <= 0.040){
                Er.r(index) = 1836.4;
                Er.i(index) = 0.50268;
            }

            //fat slice
            real angle = std::atan2(CX, CY) * 180 / PI;
            if (distance < 0.040 && angle > -120  && angle < -60){
                Er.r(index) = 27.222;
                Er.i(index) = 0.025079;
            }

            //Bone
            if (distance <= 0.01){
                Er.r(index) = 144.51;
                Er.i(index) = 0.024353;
            }

            //artery
            //            if (CX == 0 && CY == 0.020){
            //                Er.r(index) = 43;
            //                Er.i(index) = 1.85;
            //            }

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

    //init mom structs
    ImagingSpace space;
    space.x = cx;
    space.y = cy;
    space.Er = Er;
    space.dx = CDX;
    space.dy = CDY;
    space.lx = CLX;
    space.ly = CLY;
    space.probes = probes;
    space.freqs = freqs;

    space.initalGuess.r = af::constant(1, Er.a.dims());
    space.initalGuess.i = af::constant(0, Er.a.dims());

    index = 0;
    for (int i = 0; i <CNX; ++i){
        for (int j = 0; j < CNY; ++j){
            real CX = CDX2 + i* CDX-CLX/2;
            real CY = CDY2 + j * CDY - CLY/2;

            real distance = std::sqrt(std::pow(CX, 2) + std::pow(CY, 2));
            if (distance < 0.045){
                space.initalGuess.r(index) = 1800;
            }
            index++;
        }
    }

    space.initalGuess.refresh();

//    space.initalGuess = space.Er;

    RunInfo info;
    info.name = "arm";
    info.iterations = 10;
    info.lambda = 0.1;

    runBIM(info, space, CNX, CNY);
}

QTEST_APPLESS_MAIN(ThesisTestsTest)

#include "tst_thesisteststest.moc"
