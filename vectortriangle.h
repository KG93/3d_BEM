#ifndef VECTORTRIANGLE_H
#define VECTORTRIANGLE_H

#include "robinboundarycondition.h"

#include <QtGlobal>
#include <eigen3/Eigen/Geometry>

/*!
 * \brief The VectorTriangle class represents a simple triangle in 3D space.
 */
class VectorTriangle
{
public:
    //! Constructs a default VectorTriangle.
    VectorTriangle();

    /*!
     * \brief Constructs a VectorTriangle with the given nodes. The node ordering is counter-clockwise for the triangle front side.
     * \param Node1 the first node of the triangle
     * \param Node2 the second node of the triangle
     * \param Node3 the third node of the triangle
     */
    VectorTriangle(Eigen::Vector3d Node1, Eigen::Vector3d Node2, Eigen::Vector3d Node3){
        this->node1 = Node1;
        this->node2 = Node2;
        this->node3 = Node3;
    }
    /*!
     * \brief Constructs a VectorTriangle with the given element index and nodes.
     * \param elementIndex the index of the element
     * \param Node1 the first node of the triangle
     * \param Node2 the second node of the triangle
     * \param Node3 the third node of the triangle
     */
    VectorTriangle(long elementIndex, Eigen::Vector3d Node1, Eigen::Vector3d Node2, Eigen::Vector3d Node3){
        this->elementIndex = elementIndex;
        this->node1 = Node1;
        this->node2 = Node2;
        this->node3 = Node3;
    }

    /*!
     * \brief Constructs a VectorTriangle with the given element index, nodes, and normal vector.
     * \param elementIndex the index of the element
     * \param Node1 the first node of the triangle
     * \param Node2 the second node of the triangle
     * \param Node3 the third node of the triangle
     * \param normal the normal vector of the triangle
     */
    VectorTriangle(long elementIndex, Eigen::Vector3d Node1, Eigen::Vector3d Node2, Eigen::Vector3d Node3, Eigen::Vector3d normal){
        this->elementIndex = elementIndex;
        this->node1 = Node1;
        this->node2 = Node2;
        this->node3 = Node3;
        this->normal = normal;
    }

    /*!
     * \brief Constructs a VectorTriangle with the given element index, nodes, normal vector, and Robin boundary condition.
     * \param elementIndex the index of the element
     * \param Node1 the first node of the triangle
     * \param Node2 the second node of the triangle
     * \param Node3 the third node of the triangle
     * \param normal the normal vector of the triangle
     * \param robinBoundaryCondition the Robin boundary condition of the triangle
     */
    VectorTriangle(long elementIndex, Eigen::Vector3d Node1, Eigen::Vector3d Node2, Eigen::Vector3d Node3, Eigen::Vector3d normal, RobinBoundaryCondition robinBoundaryCondition){
        this->elementIndex = elementIndex;
        this->node1 = Node1;
        this->node2 = Node2;
        this->node3 = Node3;
        this->normal = normal;
        this->robinBoundaryCondition = robinBoundaryCondition;
    }

    //! The index of the element.
    long elementIndex = 0;

    //node ordering is counter clockwise for triangle front side
    //! The first node of the triangle.
    Eigen::Vector3d node1;

    //! The second node of the triangle.
    Eigen::Vector3d node2;

    //! The third node of the triangle.
    Eigen::Vector3d node3;

    //! The normal vector of the triangle. The normal vector points from the boundary into the domain.
    Eigen::Vector3d normal;

    //! The midpoint of the triangle.
    Eigen::Vector3d triangleMidpoint;
    Eigen::Vector3d color = {0.5, 0.5, 1};

    //! The Robin boundary condition of the triangle.
    RobinBoundaryCondition robinBoundaryCondition;

    /*!
     * \brief phi
     * The value of the solution (phi) at the triangle.
     */
    std::complex<double> phi;

    /*!
     * \brief dPhi
     * The value of the normal derivative (dPhi)
     */
    std::complex<double> dPhi;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#endif // VECTORTRIANGLE_H
