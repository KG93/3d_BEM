#ifndef OBSERVATIONFIELD_H
#define OBSERVATIONFIELD_H
#include "global.h"

/**
* \class ObservationField
* \brief The class represents an observation field as it is used in the BoundaryElementSolver.
*
* The field is in triangulated form and the boundary element solver calculates the solution for all the triangle midpoints.
*/
class ObservationField
{
public:
    ObservationField(){}
    ObservationField(QString name, QVector<VectorTriangle> triangles)
    {
        this->name=name;
        this->triangles=triangles;
    }

    /**
    * \brief The triangles that make up the observation field.
    */
    QVector<VectorTriangle> triangles;

    /**
    * \brief The colors for the triangles. The colors visualize the solution on the field.
    */
    QVector<Eigen::Vector4d> triangleColors;
    QVector<Eigen::Vector3d> triangleMidPoints;

    /**
    * \brief The potential at the triangles' midpoints.
    */
    Eigen::VectorXcd phiSolution;

    /**
    * \brief The sound pressure at the triangles' midpoints.
    */
    Eigen::VectorXcd soundPressure;

    /**
    * \brief The name of the sound field.
    */
    QString name;

    /**
    * \brief Calculate the geometric midpoint of the triangles.
    */
    void calculateTrianglesMidpoints();
};

#endif // OBSERVATIONFIELD_H
