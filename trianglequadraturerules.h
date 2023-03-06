#ifndef TRIANGLEQUADRATURERULES_H
#define TRIANGLEQUADRATURERULES_H

#include <eigen3/Eigen/Geometry>

//The quadrature rules are implemented after the text
//"Quadrature Formulas in Two Dimensions" by Dr. Shaozhong Deng
//Math 5172 - Finite Element Method
//Section 001, Spring 2010

/**
* \class triangleQuadratureRules
* \brief The class contains weights and abscissas for gauss quadrature with different orders for a triangular domain.
*
* The class contains weights and abscissas for gauss quadrature rules of varying orders for a triangular domain.
* The weights and abscissas are returned as an Eigen::MatrixX3d, where the first two column contain the abscissas and the second column contains the corresponding weights.
*/
class triangleQuadratureRules
{
public:
    triangleQuadratureRules();
    static Eigen::MatrixX3d weightsandAbscissa(int order); /*!< \brief Get weights and abscissas that correspond to the order. Implemented up to order 15. */
};

#endif // TRIANGLEQUADRATURERULES_H
