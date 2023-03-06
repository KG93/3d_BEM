#ifndef BOUNDARYELEMENTS_H
#define BOUNDARYELEMENTS_H

#include "triangle.h"
#include "vectortriangle.h"
#include "vectorquadrilateral.h"
#include "robinboundarycondition.h"
#include "global.h"
#include "meshfunctions.h"
#include "impedanceplane.h"

#include <QString>
#include <QVector>

/**
*\class BoundaryElements
*\brief This class represents a boundary element method (BEM) discretization of a 3D surface.
*
*This class is responsible for storing the BEM discretization of a 3D surface, which consists of a list of
*triangles, the areas of the triangles, the collocation points for the BEM, and the solution to the BEM
*equations (the potential and its normal derivative at the collocation points). This class also provides
*methods for calculating the collocation points, the triangle areas, and the average triangle dimension.
*/

class BoundaryElements
{
public:
    BoundaryElements();
    BoundaryElements(QString name){this->name=name;}
//    QString elementsGroupName;

    /**
    * @brief A list of triangles that make up the discretization.
    */
    QVector<VectorTriangle> triangles;

    /**
    * @brief The areas of the triangles.
    */
    Eigen::VectorXd trianglesArea;

    /**
    * @brief The collocation points for the BEM.
    */
    QVector<Eigen::Vector3d> collocationPoints;
//    QVector<VectorQuadrilateral> elementsVectorQuadrilaterals;
//    QVector<RobinBoundaryCondition> robinBoundaryConditions;
    /**
    * @brief The color of each triangle.
    */
    QVector<Eigen::Vector4d> triangleColors;

    /**
    * @brief The solution to the BEM equations (the potential and its normal derivative) at the collocation points.
    */
    Eigen::VectorXcd phiSolution;
    Eigen::VectorXcd dPhiSolution;
    Eigen::VectorXcd soundPressure;

    /**
     * @brief The name of the boundary element group.
     */
    QString name;

    /**
    * @brief The frequency to which the solution/pressure corresponds.
    */
    double frequency;

    /**
    * @brief A flag indicating whether this is an exterior or interior problem.
    */
    bool exteriorProblem;

    /**
    * @brief The subdomain to which this boundary element group belongs.
    */
    unsigned long subdomain = 1;

    /**
    * @brief The average triangle dimension.
    */
    double averageTriangleDim = 0;

    /**
    * @brief A list of impedance planes.
    */
    QVector<ImpedancePlane> impedancePlanes;

//    std::vector<triangle> getTriangles();
//    void addTriangles(QString elementsGroupName, std::vector<triangle> triangleVector);
//    void addTriangles(std::vector<Triangle> triangleVector);

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
    * @brief Calculate the geometric midpoint of each triangle.
    */
    void calculateTriangleMidPoints();

    /**
    * @brief Calculate the area of each triangle and estimate the average triangle dimension.
    */
    void calculateTrianglesAreaAndAverageDimension();

    /**
    * @brief Flip the normals of the triangles.
    */
    void flipNormals();

    /**
    * @brief Swap the order of the elements at the given indexes.
    * @param index1 - Index of the first element to be swapped.
    * @param index2 - Index of the second element to be swapped.
    */
    void swapElementIndexes(int index1, int index2);

    /**
    * Reflect the geometry across the given reflection plane.
    * @param reflectionPlane - The plane across which to reflect the geometry.
    */
    void reflectGeometry(ImpedancePlane reflectionPlane);

private:
    bool geometryIsReflected = false;
};

#endif // BOUNDARYELEMENTS_H
