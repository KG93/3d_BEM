#ifndef INFINITEBAFFLESECTION_H
#define INFINITEBAFFLESECTION_H
#include "globallogstrings.h"
#include "global.h"

#include <QString>
#include <QVector>
#include <QHash>
#include <QMatrix4x4>


/**
* \class InfiniteBaffleSection
* \brief The InfiniteBaffleSection class represents an infinite baffle section in a solving script for a sound simulation.
*
* An infinite baffle section represents a flat, infinite surface that absorbs/reflects/impedes sound waves. It is used to model the absorption/reflection/impedance of sound by walls, floors, and ceilings in a sound simulation.
*/
class InfiniteBaffleSection
{
public:
    /**
    *
    * \brief Constructs an InfiniteBaffleSection with default values.
    */
    InfiniteBaffleSection();

    /**
    *
    * \brief Constructs an InfiniteBaffleSection with the given name.
    * \param name The name of the infinite baffle section.
    */
    InfiniteBaffleSection(QString name/*, const QVector<QString>& MeshFiles, const QVector<QString>& MeshFileAliases*/){
        this->name=name;
//        this->meshFiles=MeshFiles,
//        this->meshFileAliases=MeshFileAliases;
//        Scale={1,1,1};
    }

    /**
    *
    * \brief Handles a solving script line that defines properties of the infinite baffle section.
    * \param currentLine The script line to handle.
    * \return True if the script line was valid and was handled successfully, false otherwise.
    */
    bool handleScriptLine(const QString& currentLine);

    QString name; /*!< The name of the infinite baffle section. */
    QString RefNodes;  /*!< The reference nodes of the infinite baffle section. */
    QString MeshFileAlias;
    unsigned int Subdomain = 1;
    Eigen::Vector3d planePoint = Eigen::Vector3d(0,0,0);
    Eigen::Vector3d planeNormal = Eigen::Vector3d(0,0,1);
    double distanceToParallelPlane = 0;

    Eigen::Vector3d Scale = Eigen::Vector3d(1,1,1);
    Eigen::Vector3d Shift = Eigen::Vector3d(0,0,0);
    Eigen::Vector3d Rotate = Eigen::Vector3d(0,0,0);
    QMatrix4x4 transformationMatrix;

private:
    bool validLine = false;

};

#endif // INFINITEBAFFLESECTION_H
