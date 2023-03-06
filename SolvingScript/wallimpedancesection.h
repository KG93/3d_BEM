#ifndef WALLIMPEDANCESECTION_H
#define WALLIMPEDANCESECTION_H
#include "global.h"
#include "SolvingScript/wallimpedanceelement.h"
#include "globallogstrings.h"

enum ElTypeEnum {Reflection=0, Damping, Impedance, Admittance};

/**
 * \class WallImpedanceSection
 * \brief This class represents a solving script section to set impedance boundary conditions to boundary elements.
 *
 * A WallImpedanceSection refers to an ElementSection. It can be used to control the boundary conditions of the boundary elements in the elements section.
 * Parameters to control reflection, damping, impedance and admittance are available.
 */
class WallImpedanceSection
{
public:
    WallImpedanceSection(){}

    /**
    * \brief Constructor for WallImpedanceSection
    * \param name The name of the wall impedance section
    */
    WallImpedanceSection(QString name){this->name=name;}

    QString name;
    ElTypeEnum ImpType = Reflection;
    bool Normalized = true;
    QString RefElements;
    std::complex<double> Value = 1.0;
    RobinBoundaryCondition robinBoundaryCondition;
    double ft;
    QString DataFileAlias;
    QVector<WallImpecanceElement> wallImpecanceElements;

    /**
    * \brief Handles a solving script line that defines properties of the wall impedance section.
    * \param currentLine The script line to handle.
    * \return True if the script line was valid and was handled successfully, false otherwise.
    */
    bool handleScriptLine(const QString& line);
//    bool checkIdentifiersAllUnique(QVector<Node> nodes);
    bool checkIdentifiersAllUnique();

    /**
    * \brief Set the wavenumber for the wall impedance section.
    * \param wavenumber The wavenumber to use.
    */
    void setWavenumer(std::complex<double> wavenumber){this->wavenumber=wavenumber;}

    /**
    * \brief Set the medium impedance for the wall impedance section.
    * \param mediumImpedance The medium impedance to use.
    */
    void setMediumImpedance(std::complex<double> mediumImpedance){this->mediumImpedance=mediumImpedance;}

    /**
    * \brief Set up the impedance values for the wall impedance section.
    */
    /*bool*/ void setUpValues();

private:
    bool handleWallImpedancesListLine(const QString& line);
    bool handleParameterSectionLine(const QString& line);
    /*const*/ std::complex<double> imaginaryUnit = {0, 1};

    std::complex<double> wavenumber;
    std::complex<double> mediumImpedance;

    bool readingWallImpedancesList = false;
    bool validLine = false;
};

#endif // WALLIMPEDANCESECTION_H
