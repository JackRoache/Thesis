#include "tissueproperties.h"
#include <QFile>
#include <QTextStream>
/*
 * Data was collected from
 */


TissueProperties::TissueProperties()
{
    QFile file("body_tissue.csv");
    Q_ASSERT(file.open(QFile::ReadOnly));
    QTextStream stream(&file);

    //first 2 lines are the headers
    stream.readLine();
    stream.readLine();


    //read line by line
    while (!stream.atEnd()){
        QString line = stream.readLine();

        QStringList list = line.split(",");

        if (list.size() != 27)
            continue;

        int id = list.at(0).toInt();
        //1 tissue
        //2 organ

        //300mhz

        c300mhz.insert(id, list.at(3).toDouble());
        p300mhz.insert(id, list.at(4).toDouble());
        lt300mhz.insert(id, list.at(5).toDouble());

        //400mhz
        c400mhz.insert(id, list.at(6).toDouble());
        p400mhz.insert(id, list.at(7).toDouble());
        lt400mhz.insert(id, list.at(8).toDouble());

        //500mhz
        c500mhz.insert(id, list.at(9).toDouble());
        p500mhz.insert(id, list.at(10).toDouble());
        lt500mhz.insert(id, list.at(11).toDouble());

        //600mhz
        c600mhz.insert(id, list.at(12).toDouble());
        p600mhz.insert(id, list.at(13).toDouble());
        lt600mhz.insert(id, list.at(14).toDouble());

        //700mhz
        c700mhz.insert(id, list.at(15).toDouble());
        p700mhz.insert(id, list.at(16).toDouble());
        lt700mhz.insert(id, list.at(17).toDouble());

        //800mhz
        c800mhz.insert(id, list.at(18).toDouble());
        p800mhz.insert(id, list.at(19).toDouble());
        lt800mhz.insert(id, list.at(20).toDouble());

        //900mhz
        c900mhz.insert(id, list.at(21).toDouble());
        p900mhz.insert(id, list.at(22).toDouble());
        lt900mhz.insert(id, list.at(23).toDouble());

        //1000mhz
        c1000mhz.insert(id, list.at(24).toDouble());
        p1000mhz.insert(id, list.at(25).toDouble());
        lt1000mhz.insert(id, list.at(26).toDouble());
    }

}

QMap<int, float> TissueProperties::permitivity(int freq) const
{
    switch (freq){
    case 300000000:
        return p300mhz;
    case 400000000:
        return p400mhz;
    case 500000000:
        return p500mhz;
    case 600000000:
        return p600mhz;
    case 700000000:
        return p700mhz;
    case 800000000:
        return p800mhz;
    case 900000000:
        return p900mhz;
    case 1000000000:
        return p1000mhz;
    default:
        Q_ASSERT(false);
    }
    //compiler warning
    return QMap<int, float>();
}

QMap<int, float> TissueProperties::conductivity(int freq) const
{
    switch (freq){
    case 300000000:
        return c300mhz;
    case 400000000:
        return c400mhz;
    case 500000000:
        return c500mhz;
    case 600000000:
        return c600mhz;
    case 700000000:
        return c700mhz;
    case 800000000:
        return c800mhz;
    case 900000000:
        return c900mhz;
    case 1000000000:
        return c1000mhz;
    default:
        Q_ASSERT(false);
    }
    //compiler warning
    return QMap<int, float>();
}

QMap<int, float> TissueProperties::lossTangent(int freq) const
{
    switch (freq){
    case 300000000:
        return lt300mhz;
    case 400000000:
        return lt400mhz;
    case 500000000:
        return lt500mhz;
    case 600000000:
        return lt600mhz;
    case 700000000:
        return lt700mhz;
    case 800000000:
        return lt800mhz;
    case 900000000:
        return lt900mhz;
    case 1000000000:
        return lt1000mhz;
    default:
        Q_ASSERT(false);
    }
    //compiler warning
    return QMap<int, float>();
}
