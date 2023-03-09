#ifndef CONTROLFIELDSECTION_H
#define CONTROLFIELDSECTION_H

#include <QString>

#include "global.h"
//#include "globallogstrings.h"

class ControlFieldSection
{
public:
    ControlFieldSection();

    /*!
    Handles a script line for the control field section.
    @param line The line to be handled.
    @return True if the line was valid, false otherwise.
    */
    bool handleScriptLine(const QString& line);
    double f1;
    double f2;
    double MeshFrequency;
    double EdgeLength;
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
};

#endif // CONTROLFIELDSECTION_H
