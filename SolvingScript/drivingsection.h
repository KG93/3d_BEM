#ifndef DRIVINGSECTION_H
#define DRIVINGSECTION_H

#include "globallogstrings.h"
#include "global.h"
#include "drivingelement.h"

/**
 * \class DrivingSection
 * \brief The DrivingSection represents a solving script section to set radiation boundary conditions to boundary elements.
 *
 * A DrivingSection refers to an ElementSection. It can be used to control the radiation conditions of the boundary elements in the elements section.
 */
class DrivingSection
{
public:
    DrivingSection();
    DrivingSection(QString name){this->name=name;}
    QString name;
    int ImpType;
    bool Normalized;
    QString RefElements;
    QString RefNodes;
    QString MeshFileAlias;
    int DrvGroup;
    Eigen::Vector3d Direction;
    bool DirIsAxial=false;
    std::complex<double> DrvWeight = 0;
    double DrvDelay = 0;
    QVector<DrivingElement> drivingElements;

    /**
    * \brief Handles a solving script line that defines properties of the driving section.
    * \param currentLine The script line to handle.
    * \return True if the script line was valid and was handled successfully, false otherwise.
    */
    bool handleScriptLine(const QString& line);
//    bool checkIdentifiersAllUnique(QVector<Node> nodes);

//    /**
//     * \brief Check that all identifiers in the driving section are unique.
//     * \return True if all identifiers are unique, false otherwise.
//     */
//    bool checkIdentifiersAllUnique();
//    void setWavenumer(std::complex<double> wavenumber){this->wavenumber=wavenumber;}
//    void setMediumImpedance(std::complex<double> mediumImpedance){this->mediumImpedance=mediumImpedance;}
//    bool setUpValues();

private:
    bool handleDrivingElementsListLine(const QString& line);
    bool handleParameterSectionLine(const QString& line);
//    const std::complex<double> imaginaryUnit = std::complex<double>(0.0,1.0);
//    std::complex<double> wavenumber;
//    std::complex<double> mediumImpedance;

    bool readingDrivingElementsList=false;
    bool validLine = false;
};

#endif // DRIVINGSECTION_H
