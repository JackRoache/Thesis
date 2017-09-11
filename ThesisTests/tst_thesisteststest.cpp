#include <QString>
#include <QtTest>
#include <QImage>
#include <QFile>

#include "../types.h"
#include "../mom.h"
#include "../clbessel.h"

#define PI      (real)(3.14159265)

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
    std::string name("simple");
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

    RunInfo info;
    info.name = "Basic";
    info.iterations = 50;
    info.lambda = 10;

    surfToImage(cx, cy, space.Er, (name + "Start").c_str(), CNX, CNY);
    MOM mom;

    mom.setImagingSpace(&space);
    mom.setIterations(&info);

    mom.run();

    surfToImage(cx, cy, space.Er, (name + "Done").c_str(), CNX, CNY);
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

    real RX_XRADIUS = 0.2;
    real RX_YRADIUS = 0.2;
    real RX_ANGLE =30;

    carray Er(CNX * CNY);
    af::array cx(CNX * CNY);
    af::array cy(CNX * CNY);

    std::vector<Position> probes;
    std::vector<real> freqs = {2.45e9};

    int index = 0;
    for (int i = 0; i <CNX; ++i){
        for (int j = 0; j < CNY; ++j){
            real CX = CDX2 + i* CDX-CLX/2;
            real CY = CDY2 + j * CDY - CLY/2;

            cx(index) = CX;
            cy(index) = CY;
            Er.r(index) = 1;
            Er.i(index) = 0;
            index++;
        }
    }
    //generate Fat
    index = 0;
    for (int i = 0; i <CNX; ++i){
        for (int j = 0; j < CNY; ++j){
            real CX = CDX2 + i* CDX-CLX/2;
            real CY = CDY2 + j * CDY - CLY/2;

            real distance = std::sqrt(std::pow(CX, 2) + std::pow(CY, 2));
            //Fat
            if (distance <= 0.045){
                Er.r(index) = 21;
                Er.i(index) = 0.02;
            }

            //Muscle
            if (distance <= 0.040){
                Er.r(index) = 639.6;
                Er.i(index) = 0.55;
            }

            //fat slice
            real angle = std::atan2(CX, CY) * 180 / PI;
            if (distance < 0.040 && angle > -120  && angle < -60){
                Er.r(index) = 21;
                Er.i(index) = 0022;
            }

            //Bone
            if (distance <= 0.01){
                Er.r(index) = 164;
                Er.i(index) = 0.00021;
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

    RunInfo info;
    info.name = "Basic";
    info.iterations = 100;
    info.lambda = 100;

    surfToImage(cx, cy, space.Er, (name + "Start").c_str(), CNX, CNY);

    MOM mom;

    mom.setImagingSpace(&space);
    mom.setIterations(&info);

    mom.run();

    surfToImage(cx, cy, space.Er, (name + "Done").c_str(), CNX, CNY);
}

QTEST_APPLESS_MAIN(ThesisTestsTest)

#include "tst_thesisteststest.moc"
