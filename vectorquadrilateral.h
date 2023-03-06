#ifndef VECTORQUADRILATERAL_H
#define VECTORQUADRILATERAL_H

#include "robinboundarycondition.h"

#include <QtGlobal>
#include <eigen3/Eigen/Geometry>

/*!
 * \brief The VectorQuadrilateral class represents a simple quadrilateral in 3D space.
 *
 * This class is used in the Boundary Element Method to represent the boundary of a domain.
 */
class VectorQuadrilateral
{
public:
    VectorQuadrilateral(){
        node1 = {0,0,0};
        node2 = {0,0,0};
        node3 = {0,0,0};
        node4 = {0,0,0};
    }

    /*!
     * \brief Constructs a quadrilateral with the given nodes.
     *
     * \param Node1 The first node of the quadrilateral.
     * \param Node2 The second node of the quadrilateral.
     * \param Node3 The third node of the quadrilateral.
     * \param Node4 The fourth node of the quadrilateral.
     */
    VectorQuadrilateral(Eigen::Vector3d Node1, Eigen::Vector3d Node2, Eigen::Vector3d Node3, Eigen::Vector3d Node4){
            this->node1 = Node1;
            this->node2 = Node2;
            this->node3 = Node3;
            this->node4 = Node4;
    }

    /*!
     * \brief Constructs a quadrilateral with the given index and nodes.
     *
     * \param elementIndex The index of the quadrilateral.
     * \param Node1 The first node of the quadrilateral.
     * \param Node2 The second node of the quadrilateral.
     * \param Node3 The third node of the quadrilateral.
     * \param Node4 The fourth node of the quadrilateral.
     */
    VectorQuadrilateral(qint64 elementIndex, Eigen::Vector3d Node1, Eigen::Vector3d Node2, Eigen::Vector3d Node3, Eigen::Vector3d Node4){
            this->elementIndex = elementIndex;
            this->node1 = Node1;
            this->node2 = Node2;
            this->node3 = Node3;
            this->node4 = Node4;
    }

    /*!
    * \brief The index of the quadrilateral.
    */
    qint64 elementIndex;

    /*!
     * \brief The first node of the quadrilateral.
     */
    Eigen::Vector3d node1;

    /*!
     * \brief The second node of the quadrilateral.
     */
    Eigen::Vector3d node2;

    /*!
     * \brief The third node of the quadrilateral.
     */
    Eigen::Vector3d node3;

    /*!
     * \brief The fourth node of the quadrilateral.
     */
    Eigen::Vector3d node4;

    /*!
     * \brief The color of the quadrilateral.
     *
     * This property is used for visualizing the quadrilateral in a graphics view.
     */
    Eigen::Vector3d color = {0.5, 0.5, 1};

    /*!
     * \brief The Robin boundary condition of the quadrilateral.
     */
    RobinBoundaryCondition robinBoundaryCondition;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};
#endif // VECTORQUADRILATERAL_H
