#ifndef LINEARBOUNDARYELEMENTS_H
#define LINEARBOUNDARYELEMENTS_H

#include "lineartrianglenode.h"
#include "lineartriangle.h"
#include "impedanceplane.h"
#include "global.h"

#include <QVector>


class LinearBoundaryElements
{
public:
    LinearBoundaryElements();

    double frequency; /*!< \brief The frequency to which the solution/pressure corresponds. */

    QVector<LinearTriangleNode> nodes; /*!< \brief A list of triangle nodes. */
    QVector<LinearTriangle> triangles; /*!< \brief A list of triangles that make up the discretization. */

    /**
    * @brief The average triangle dimension.
    */
    double averageTriangleDim = 0;

    QVector<ImpedancePlane> impedancePlanes; /*!< \brief A list of impedance planes. */

    Eigen::VectorXd trianglesArea;

    QVector<Eigen::Vector3d> collocationPoints;


private:
    void calculateAssociatedTriangles(); /*!< \brief Add the indices of the triangles that are associated with a node to the node. */

    /**
    * @brief Calculate the collocation points for the BEM.
    *
    * This method calculates the geometric midpoint of each triangle in the discretization and stores
    * the resulting points in the #collocationPoints member. These points are used as the collocation
    * points in the BEM.
    *
    */
    void calculateCollocationPoints();

    /**
    * @brief Calculate the area of each triangle in the discretization.
    *
    * This method calculates the area of each triangle in the discretization and stores the results in
    * the #trianglesArea member.
    *
    */
    void calculateTrianglesArea();

    /**
    * @brief Calculate the area of each triangle and estimate the average triangle dimension.
    */
    void calculateTrianglesAreaAndAverageDimension();
};

#endif // LINEARBOUNDARYELEMENTS_H
