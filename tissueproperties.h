#ifndef TISSUEPROPERTIES_H
#define TISSUEPROPERTIES_H
#include "QMap"

class TissueProperties
{
public:
    TissueProperties();

    QMap<int, float> permitivity(int freq) const;
    QMap<int, float> conductivity(int freq) const;

    QMap<int, float> lossTangent(int freq) const;


private:
    QMap<int, float> p300mhz;
    QMap<int, float> c300mhz;
    QMap<int, float> lt300mhz;

    QMap<int, float> p400mhz;
    QMap<int, float> c400mhz;
    QMap<int, float> lt400mhz;

    QMap<int, float> p500mhz;
    QMap<int, float> c500mhz;
    QMap<int, float> lt500mhz;

    QMap<int, float> p600mhz;
    QMap<int, float> c600mhz;
    QMap<int, float> lt600mhz;

    QMap<int, float> p700mhz;
    QMap<int, float> c700mhz;
    QMap<int, float> lt700mhz;

    QMap<int, float> p800mhz;
    QMap<int, float> c800mhz;
    QMap<int, float> lt800mhz;

    QMap<int, float> p900mhz;
    QMap<int, float> c900mhz;
    QMap<int, float> lt900mhz;

    QMap<int, float> p1000mhz;
    QMap<int, float> c1000mhz;
    QMap<int, float> lt1000mhz;
};

#endif // TISSUEPROPERTIES_H
