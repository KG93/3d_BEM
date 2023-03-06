#ifndef CONTROLSPECTRUMSECTION_H
#define CONTROLSPECTRUMSECTION_H
#include <QString>
#include "global.h"
//#include "globallogstrings.h"


class ControlSpectrumSection
{
public:
    ControlSpectrumSection();
    bool handleScriptLine(const QString& line);

    QString Name;
    double f1;
    double f2;
    double MeshFrequency;
    double c;
    double Altitude;
    double StaticPressure;
    double Temperature;
    double rho;

    quint64 NumFrequencies;
    quint64 abscissa;
    enum abscissaEnum{log=0,lin};

private:
    bool validLine;
    double Tiny = 0.0000000000000001;


};

#endif // CONTROLSPECTRUMSECTION_H
