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

        if (list.size() != 13)
            continue;

        int id = list.at(0).toInt();
        //1 tissue
        //2 organ

        //100mhz
        p100mhz.insert(id, list.at(3).toDouble());
        c100mhz.insert(id, list.at(4).toDouble());

        //500mhz
        p500mhz.insert(id, list.at(5).toDouble());
        c500mhz.insert(id, list.at(6).toDouble());

        //1ghz
        p1ghz.insert(id, list.at(7).toDouble());
        c1ghz.insert(id, list.at(8).toDouble());

        //2ghz
        p2ghz.insert(id, list.at(9).toDouble());
        c2ghz.insert(id, list.at(10).toDouble());

        //5ghz
        p5ghz.insert(id, list.at(11).toDouble());
        c5ghz.insert(id, list.at(12).toDouble());


    }

}
