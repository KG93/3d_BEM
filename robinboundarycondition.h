#ifndef ROBINBOUNDARYCONDITION_H
#define ROBINBOUNDARYCONDITION_H
//#include <eigen3/Eigen/Geometry>
#include <complex>

/**
* \brief The class represents Robin (radiation) boundary conditions.
*/
class RobinBoundaryCondition
{
public:

    // a*U+b*dU = g
    /**
    * \brief Constructs a default RobinBoundaryCondition.
    */
    RobinBoundaryCondition(){}

    /**
    * \brief Constructs a RobinBoundaryCondition with the given parameters a phi + b dPhi = g.
    * \param a Factor for the solution (phi).
    * \param b Factor for the normal derivative (dPhi).
    * \param g The right-hand side.
    */
    RobinBoundaryCondition( std::complex<double> a,  std::complex<double>b,  std::complex<double>g)
    {
        this->a = a;
        this->b = b;
        this->g = g;
    }

    /**
    * \brief Factor for the solution (phi).
    */
    std::complex<double> a = 0.0;

    /**
    * \brief Factor for the normal derivative (dPhi).
    */
    std::complex<double> b = 1.0;

    /**
    * \brief The right-hand side.
    */
    std::complex<double> g = 0.0;
};

#endif // ROBINBOUNDARYCONDITION_H


