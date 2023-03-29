#ifndef SOLVINGSCRIPTREADER_H
#define SOLVINGSCRIPTREADER_H

#include "boundaryelements.h"
#include "LinearElements/linearboundaryelements.h"
#include "nodesSection.h"
#include "elementSection.h"
#include "meshfilepropertiessection.h"
#include "controlsolversection.h"
#include "subdomainpropertiessection.h"
#include "SolvingScript/wallimpedancesection.h"
#include "SolvingScript/pressurepointssection.h"
#include "SolvingScript/drivingsection.h"
#include "SolvingScript/infinitebafflesection.h"

#include "pointsource.h"
#include "meshfunctions.h"
#include "globallogstrings.h"

#include <iostream>
#include <QTextStream>
#include <iostream>
#include <fstream>
#include <QFile>
#include <QVector>
#include <QStringView>
#include <QRegularExpressionMatch>
#include <QRegularExpression>
#include <QtMath>

/**
* \class SolvingScriptReader
* \brief This class reads the solving script file and extracts information from it.
*
* The SolvingScriptReader class reads a solving script file, extracts the simulation solver parameters, the boundary elements and their attributes, and handles the mesh refinement.
* The extracted information is used to create input for the solver class BoundaryElementSolver.
*/
class SolvingScriptReader: public QObject
{
    Q_OBJECT
public:
    explicit SolvingScriptReader(QObject* parent = nullptr);

    /**
    * \brief Set the mesh files and their aliases for the current project.
    * \param[in] meshFiles The list of mesh files.
    * \param[in] meshFileAlias The list of mesh file aliases.
    * \return True if the mesh files and aliases were set successfully, false otherwise.
    */
    bool setMeshFilesAndAliases( const QStringList &meshFiles, const QStringList &meshFileAlias);

    /**
    * \brief Read the solving script with the given file name.
    * \param[in] fileName The file name of the solving script.
    */
    void readSolvingScript(const QString &fileName);

    /**
    * \brief Clears the data stored by this SolvingScriptReader object.
    */
    void clear();

    /**
    * \brief Removes invalid elements from the solving script.
    * \return True if any invalid elements were removed, false otherwise.
    */
    bool removeInvalidElements();

    /**
    * \brief Sets up the boundary elements in the solving script.
    * \return True if the boundary elements were set up successfully, false otherwise.
    */
    bool setupConstBoundaryElements();

    /**
    * \brief Sets up the boundary elements in the solving script.
    * \return True if the boundary elements were set up successfully, false otherwise.
    */
    bool setupLinearBoundaryElements();

    /**
    * \brief Sets up the boundary conditions in the solving script.
    */
    void setupBoundaryConditions();

    /**
    * \brief Sets up the driving/rediating elements in the solving script.
    */
    void setupDrivingElements();

    /**
    * \brief Sets up the pressure points in the solving script.
    */
    void setupPressurePoints();

    /**
    * \brief Gets the boundary elements from the solving script.
    * \return The boundary elements from the solving script.
    */
    BoundaryElements getBoundaryElements();

    /**
    * \brief Gets the point sources from the solving script.
    * \return The point sources from the solving script.
    */
    QVector<PointSource> getPointSources(){return pointSources;}

    /**
    * \brief Checks the elements identifiers and sorts them.
    */
    void checkNodesIdentifiersAndSort();
    void checkElementsIdentifiersAndSort();

    /**
    * \brief Refines the triangle and quadrilateral mesh so that the maximum element dimension is smaller than half the wavelength.
    */
    void refineElements();

    /**
    * \brief Gets the center of mass of the boundary elements in the solving script.
    * \return The center of mass of the boundary elements in the solving script.
    */
    Eigen::Vector3d getCenterOfMass(){return centerOfMass;}
    double getContainingRadius(){return containingSphereRadius;}
    double getGlobalEdgelength(){return globalEdgelength;}
    QVector<double> getFrequencies(){return frequencies;}
    double getWavespeed(){return c;}
    double getAirDensity(){return airDensity;}
    void calculateCenterOfMass();
    QList<ElementSection> getElementsSections(){return elementsSections;}
    ControlSolverSection getControlSolverSection(){return controlSolverSection;}
    QStringList getElementSectionsNames();

private:
//    void setUpBoundaryElementsParameters(BoundaryElements& boundaryElements,const ElementSection& elementSection);
    bool controlSolverSectionDeclared;
    bool newSectionDeclared(const QString &Line);
    QRegularExpression sectionIdentifierRegEx(const QString &identifier);
    void setSolverParameters();
    void handleLineAccordingToCurrentSection(const QString &Line);
    enum sections {none, Control_Solver, Driving, Subdomain_Properties, MeshFile_Properties, Elements, Nodes, WallImpedance, Pressure_Points, Infinite_Baffle};
    int currentSection = none;
    QStringList getNodesSectionsNames();
    QStringList getPressurePointsSectionsNames();

    ControlSolverSection controlSolverSection;
    QVector<NodesSection> nodesSections;
//    QStringList nodesNames;

    QList<ElementSection> elementsSections;
    mshReader meshFileReader;

//    QStringList elementsNames;
    QVector<SubdomainPropertiesSection> subdomainPropertiesSections;
    QVector<WallImpedanceSection> wallImpedanceSections;
    QVector<PressurePointsSection> pressurePointsSections;
    QVector<InfiniteBaffleSection> infiniteBaffleSections;
    QVector<MeshFilePropertiesSection> meshFilePropertiesSections;
    QVector<DrivingSection> drivingSections;
    QVector<PointSource> pointSources;
    QVector<BoundaryElements> boundaryElements;
    QVector<double> frequencies;
    double staticPressure = 101325; //Pa at 20°C and sea level
    double airDensity =1.2041;
//    const double PI = std::numbers::pi;
    std::complex<double> mediumImpedance;

    bool meshfilesAndAliasesSet;
    QStringList meshFiles;
    QStringList meshFileAlias;

    Eigen::Vector3d centerOfMass;
    double objectVolume;
    double containingSphereRadius;
    double globalEdgelength;
    double c;

signals:
    void errorMessage(QString err);
    void logMessage(QString log);

};

#endif // SOLVINGSCRIPTREADER_H
