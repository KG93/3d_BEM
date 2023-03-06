#ifndef LINEQUADRATURERULES_H
#define LINEQUADRATURERULES_H
#include <eigen3/Eigen/Geometry>

/**
* @class lineQuadratureRules
* @brief The class contains weights and abscissas for gauss quadrature with different orders for an [-1, 1]-interval domain.
*
* The class contains weights and abscissas for gauss quadrature rules of varying orders for the domain [-1, 1].
* The weights and abscissas are returned as an Eigen::MatrixX2d, where the first column contains the abscissas and the second column contains the corresponding weights.
*/
class lineQuadratureRules
{
public:
    lineQuadratureRules();

    /**
    * @brief Get weights and abscissas for a gauss quadrature rule.
    * This method returns the weights and abscissas for a gauss quadrature rule of the specified order for the domain [-1, 1].
    * The weights and abscissas are returned as an Eigen::MatrixX2d, where the first column contains the abscissas and the second column contains the corresponding weights.
    * @param order The order of the gauss quadrature rule.
    * @return Eigen::MatrixX2d with the abscissas and weights of the quadrature rule.
    */
    static Eigen::MatrixX2d weightsandAbscissa(int order);
};

#endif // LINEQUADRATURERULES_H
