#ifndef LINEARTRIANGLE_H
#define LINEARTRIANGLE_H
#include <eigen3/Eigen/Geometry>
#include "../robinboundarycondition.h"

/*!
 * \brief The LinearTriangle class represents a simple triangle in 3D space with linear shape functions.
 */
class LinearTriangle
{
public:
    LinearTriangle();

    /*!
     * \brief Constructs a VectorTriangle with the given nodes. The node ordering is counter-clockwise for the triangle front side.
     * \param node1Index the first node of the triangle
     * \param node2Index the second node of the triangle
     * \param node3Index the third node of the triangle
     */
    LinearTriangle(long node1Index,long node2Index, long node3Index){
        this->node1Index = node1Index;
        this->node2Index = node2Index;
        this->node3Index = node3Index;
    }

    /*!
     * \brief Constructs a VectorTriangle with the given element index, nodes, normal vector, and Robin boundary condition.
     * \param elementIndex the index of the element
     * \param node1Index the first node of the triangle
     * \param node2Index the second node of the triangle
     * \param node3Index the third node of the triangle
     * \param normal the normal vector of the triangle
     * \param robinBoundaryCondition the Robin boundary condition of the triangle
     */
    LinearTriangle(long elementIndex, long node1Index,long node2Index, long node3Index, Eigen::Vector3d normal, RobinBoundaryCondition robinBoundaryCondition){
        this->elementIndex = elementIndex;
        this->node1Index = node1Index;
        this->node2Index = node2Index;
        this->node3Index = node3Index;
        this->normal = normal;
        this->robinBoundaryCondition = robinBoundaryCondition;
    }

    //! The index of the element.
    long elementIndex = 0;

    //node ordering is counter clockwise for triangle front side
    //! The first node of the triangle.
    long node1Index;

    //! The second node of the triangle.
    long node2Index;

    //! The third node of the triangle.
    long node3Index;

//    //! The first node of the triangle.
//    Eigen::Vector3d node1;

//    //! The second node of the triangle.
//    Eigen::Vector3d node2;

//    //! The third node of the triangle.
//    Eigen::Vector3d node3;

    //! The normal vector of the triangle. The normal vector points from the boundary into the domain.
    Eigen::Vector3d normal;

//    //! The midpoint of the triangle.
//    Eigen::Vector3d triangleMidpoint;
//    Eigen::Vector3d color = {0.5, 0.5, 1};

    //! The Robin boundary condition of the triangle.
    RobinBoundaryCondition robinBoundaryCondition;

//    /*!
//     * \brief phi
//     * The value of the solution (phi) at the triangle.
//     */
//    std::complex<double> phi;

//    /*!
//     * \brief dPhi
//     * The value of the normal derivative (dPhi)
//     */
//    std::complex<double> dPhi;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#endif // LINEARTRIANGLE_H
