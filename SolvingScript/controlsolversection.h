#ifndef CONTROLSOLVERSECTION_H
#define CONTROLSOLVERSECTION_H_H
#include "global.h"
#include <iostream>
#include <QObject>
#include <QRegularExpressionMatch>
#include <QRegularExpression>
#include "globallogstrings.h"



class ControlSolverSection
{
//    Q_OBJECT
public:

    ControlSolverSection();

    bool handleScriptLine(const QString line);

    double f1=0; /*!< Lowest frequency for the BEM simulation. */
    double f2=0; /*!< Highest frequency for the BEM simulation. */
    double MeshFrequency=0; /*!< MeshFrequency controls the mesh refinement. (max triangle edge < wavelength/2)*/
    double c = 0; /*!< Speed of sound in m/s. */
    double Altitude = 0;  /*!< Altitute in meters above sea level. */
    double StaticPressure = 101325; /*!< The static pressure in the acoustic medium in Pa. */
    double Temperature = 25; /*!< Temperature of the acoustic medium in degrees celsius. */
    double rho = 1.2041; /*!< Density of the acoustic medium. */

    quint64 NumFrequencies=0; /*!< Number of frequencies to be solved for in the interval [f1, f2]. */
    quint64 abscissa;
    enum abscissaEnum{log=0,lin};  /*!< If log; simulation frequencies are chosen equidistantly in [log(f1), log(f2)]. If lin; simulation frequencies are chosen equidistantly in [f1, f2]. */
    bool allMandatoryParametersSet(); /*!< Check whether all mandatory parameters of the solver are set. */
    QVector<double> getFrequencies(); /*!< Return the simulation frequencies. */
    double getWavespeed(); /*!< Return the wavespeed. */
    double getEdgelength(); /*!< Return the maximum admissible mesh edgelength. */
    double getTemperatureInC(); /*!< Return the temperature of the acoustic medium in degrees celsius. */
    double getDensity(); /*!< Return the density of the acoustic medium. */
    void clear();

private:
    bool validLine = false;
};

#endif // CONTROLSOLVERSECTION_H_H
