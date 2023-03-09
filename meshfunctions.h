#ifndef MESHFUNCTIONS_H
#define MESHFUNCTIONS_H

#include "global.h"
#include <QVector4D>
#include <QMatrix4x4>

/**
* \class MeshFunctions
* \brief The class contains functions for manipulating and refining triangular and quadrilateral meshes.
*
* The class contains functions for swapping triangle and quadrilateral normals, for converting quadrilaterals
* to triangles, and for refining triangular and quadrilateral meshes to smaller edge lengths. The class
* also contains functions for transforming triangles and quadrilaterals with an Eigen::Matrix4d and
* for reflecting triangles over a plane.
*/
class MeshFunctions
{
public:
    MeshFunctions();

    /**
    \brief Reverse the direction of the triangle normals.
    \param elementsVectorTriangles The triangles for which to reverse the normals.
    */
    static void swapTriangleNormals(QVector<VectorTriangle>& elementsVectorTriangles);

    /**
    \brief Reverse the direction of the quadrilateral normals.
    \param elementsVectorQuadrilaterals The quadrilaterals for which to reverse the normals.
    */
    static void swapQuadrilateralNormals(QVector<VectorQuadrilateral>& elementsVectorQuadrilaterals);

    /**
    \brief Convert a vector of quadrilaterals to a vector of triangles.
    \param quadrilaterals The quadrilaterals to convert.
    \return A vector of triangles that make up the quadrilaterals.
    */
    static QVector<VectorTriangle> quadrilateralsToTriangles(const QVector<VectorQuadrilateral>& quadrilaterals); /*!< Convert a vector of quadrilaterals to a vector of triangles. */

    /**
    \brief Refine quadrilaterals until all edges are smaller than EdgeLength.
    \param elementsVectorQuadrilaterals The quadrilaterals to refine.
    \param EdgeLength The maximum allowed edge length.
    \return A vector of triangles that make up the refined quadrilaterals.
    */
    static QVector<VectorTriangle> refineQuadrilaterals(const QVector<VectorQuadrilateral> &elementsVectorQuadrilaterals, const double EdgeLength); /*!< Refine quadrilaterals untill all edges are smaller than EdgeLength. */

    /**
    \brief Recursively refine a single quadrilateral.
    \param parentQuadrilateral The quadrilateral to refine.
    \param longestEdgeIndex The index of the longest edge in the quadrilateral.
    \param vectorToStoreQuadrilateral The vector to store the refined quadrilaterals in.
    \param EdgeLength The maximum allowed edge length.
    */
    static void recursiveRefineQuadrilateral(VectorQuadrilateral parentQuadrilateral, int longestEdgeIndex, QVector<VectorQuadrilateral>& vectorToStoreQuadrilateral, const double EdgeLength);

    /**
    \brief Refine triangles until all edges are smaller than EdgeLength.
    \param elementsVectorTriangles The triangles to refine.
    \param EdgeLength The maximum allowed edge length.
    \return A vector of triangles that make up the refined triangles.
    */
    static QVector<VectorTriangle> refineTriangles(const QVector<VectorTriangle> &elementsVectorTriangles, const double EdgeLength); /*!< Refine triangles untill all edges are smaller than EdgeLength. */

    /**
    \brief Recursively refine a single triangle.
    \param vectorTriangle The triangle to refine.
    \param longestEdgeIndex The index of the longest edge in the triangle.
    \param vectorToStoreTriangles The vector to store the refined triangles in.
    \param EdgeLength The maximum allowed edge length.
    */
    static void recursiveRefineTriangle(VectorTriangle vectorTriangle, int longestEdgeIndex, QVector<VectorTriangle>& vectorToStoreTriangles, const double EdgeLength);

    /**
    @brief Calculates the middle point of a line.
    @param point1 The first point of the line.
    @param point2 The second point of the line.
    @return The middle point of the line.
    */
    static Eigen::Vector3d middleOfLine(Eigen::Vector3d point1,Eigen::Vector3d point2);

    /**
    @brief Finds the maximum of four values.
    @param a The first value.
    @param b The second value.
    @param c The third value.
    @param d The fourth value.
    @return A pair containing the maximum value and its index.
    */
    static QPair<double,quint64> maxOfFour(double a, double b, double c,double d);

    /**
    @brief Converts a 3d vector to a 4d homogeneous vector.
    @param vector The 3d vector to convert.
    @return The 4d homogeneous vector.
    */
    static Eigen::Vector4d homogVec(const Eigen::Vector3d vector){return Eigen::Vector4d(vector(0),vector(1),vector(2),1);}

    /**
    @brief Converts a 4d homogeneous vector to a 3d vector.
    @param vector The 4d homogeneous vector to convert.
    @return The 3d vector.
    */
    static Eigen::Vector3d eigenVec(const Eigen::Vector4d vector){return Eigen::Vector3d(vector(0)/vector(3),vector(1)/vector(3),vector(2)/vector(3));}

    /**
    @brief Transforms a vector of triangles using a 4d transformation matrix.
    @param elementsVectorTriangles The vector of triangles to transform.
    @param transformationMatrix The transformation matrix to use.
    */
    static void transformTriangles(QVector<VectorTriangle> &elementsVectorTriangles,const Eigen::Matrix4d transformationMatrix);

    /**
    @brief Transforms a vector of quadrilaterals using a 4d transformation matrix.
    @param elementsVectorQuadrilaterals The vector of quadrilaterals to transform.
    @param transformationMatrix The transformation matrix to use.
    */
    static void transformQuadrilaterals(QVector<VectorQuadrilateral> &elementsVectorQuadrilaterals,const Eigen::Matrix4d transformationMatrix);

    /**
    @brief Reflects a vector of triangles by a plane.
    @param elementsVectorTriangles The vector of triangles to reflect.
    @param planePoint A point on the reflection plane.
    @param planeNormal The normal of the reflection plane.
    */
    static void reflectTriangles(QVector<VectorTriangle> &elementsVectorTriangles, Eigen::Vector3d planePoint, Eigen::Vector3d planeNormal); // reflect triangles by a plane

    /**
    @brief Converts a QT 4x4 matrix to an Eigen 4x4 matrix.
    @param qtMatrix The QT matrix to convert.
    @return The equivalent Eigen matrix.
    */
    static Eigen::Matrix4d qtToEigenMatrix(const QMatrix4x4& qtMatrix);
};

#endif // MESHFUNCTIONS_H
