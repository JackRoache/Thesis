#ifndef TISSUEPROPERTIES_H
#define TISSUEPROPERTIES_H
#include "QMap"

class TissueProperties
{
public:
    TissueProperties();

    QMap<int, float> permitivity100mhz() const{
        return p100mhz;
    }
    QMap<int, float> conductivity100mhz() const{
        return c100mhz;
    }

    QMap<int, float> permitivity500mhz() const{
        return p500mhz;
    }
    QMap<int, float> conductivity500mhz() const{
        return c500mhz;
    }

    QMap<int, float> permitivity1ghz() const{
        return p1ghz;
    }
    QMap<int, float> conductivity1ghz() const{
        return c1ghz;
    }

    QMap<int, float> permitivity2ghz() const {
        return p2ghz;
    }
    QMap<int, float> conductivity2ghz() const {
        return c2ghz;
    }

    QMap<int, float> permitivity5ghz() const {
        return p5ghz;
    }
    QMap<int, float> conductivity5ghz() const {
        return c5ghz;
    }

private:
    QMap<int, float> p100mhz;
    QMap<int, float> c100mhz;

    QMap<int, float> p500mhz;
    QMap<int, float> c500mhz;

    QMap<int, float> p1ghz;
    QMap<int, float> c1ghz;

    QMap<int, float> p2ghz;
    QMap<int, float> c2ghz;

    QMap<int, float> p5ghz;
    QMap<int, float> c5ghz;
};

#endif // TISSUEPROPERTIES_H
