#ifndef IMPEDANCEPLANE_H
#define IMPEDANCEPLANE_H

#include <eigen3/Eigen/Geometry>
#include "robinboundarycondition.h"

/**
\brief A plane in 3d space that bounds the acoustic domain with an impedance condition.
*/
class ImpedancePlane
{
public:
    ImpedancePlane();

    /**
    \brief Point on the plane that bounds the acoustic domain.
    */
    Eigen::Vector3d halfSpacePlanePoint = Eigen::Vector3d(0,0,0);

    /**
    \brief Normal vector to the plane that bounds the acoustic domain.
    */
    Eigen::Vector3d halfSpacePlaneNormal = Eigen::Vector3d(0,0,0);

    /**
    \brief Robin boundary condition on the plane that bounds the acoustic domain.
    */
    RobinBoundaryCondition halfSpaceBoundaryCondition;
};

#endif // IMPEDANCEPLANE_H
