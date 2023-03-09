#ifndef LINEARTRIANGLENODE_H
#define LINEARTRIANGLENODE_H

#include <eigen3/Eigen/Geometry>
#include <QVector>
#include "robinboundarycondition.h"

/**
*\brief A node for a 3d linear triangle.
*/
class LinearTriangleNode
{
public:
    LinearTriangleNode();

    LinearTriangleNode(Eigen::Vector3d coordinates){
        this->coordinates = coordinates;
    }


    Eigen::Vector3d coordinates = {0,0,0}; /*!< \brief The 3D coordinates of the node. */

    std::complex<double> phi = 0; /*!< \brief The value of the potential (phi) at the triangle node. */

    std::complex<double> dPhi = 0; /*!< \brief The value of the normal derivative (dPhi) at the triangle node. */
    QVector<long> associatedTriangles;

    RobinBoundaryCondition robinBoundaryCondition;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#endif // LINEARTRIANGLENODE_H
