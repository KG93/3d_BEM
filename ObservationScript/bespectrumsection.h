#ifndef BESPECTRUMSECTION_H
#define BESPECTRUMSECTION_H

#include <QString>

#include "global.h"
#include "bespectrumitem.h"
#include "globallogstrings.h"


class BESpectrumSection
{
public:
    BESpectrumSection();
    BESpectrumSection(QString name){
        this->name=name;
    }
    bool handleScriptLine(const QString&);
    bool containsNoElements();
    QString name;
    double Range;
    double Range_min;
    double Range_max;
    QString GraphHeader;
    double PolarRangeStartAngle;
    double PolarRangeEndAngle;
    quint64 PolarRangeNumberofPoints;
    double Distance=1;
    QString RefNodes;
    QList<BESpectrumItem> bESpectrumItems;
    QList<BESpectrumPoint> bESpectrumPoints;

private:
    bool handleSpectrumItemListLine(const QString& currentLine);
    bool handleParameterSectionLine(const QString& currentLine);
    bool readingSpectrumItemList = false;
    bool validLine;

};

#endif // BESPECTRUMSECTION_H
