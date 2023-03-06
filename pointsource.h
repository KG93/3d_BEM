#ifndef POINTSOURCE_H
#define POINTSOURCE_H

#include "global.h"
#include <eigen3/Eigen/Geometry>

/*!<
* \brief The class describes an acoustic point source in 3d space.
*/
class PointSource
{
public:
    PointSource(){}
    PointSource(int index,double weight,double delay,Eigen::Vector3d position)
    {
        this->index = index;
        this->weight = weight;
        this->delay = delay;
        this->position = position;
    }

    int index;
    double weight = 1;
    double delay = 0;
    Eigen::Vector3d position;
    std::complex<double> observedPhiPlusDPhi(const Eigen::Vector3d listeningPosition, const Eigen::Vector3d up , const std::complex<double> wavenumber, const std::complex<double> couplingParameter);
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
     double PI4 = 4.0 * global::PI;
     std::complex<double> imaginaryUnit = std::complex<double>(0.0,1.0);
};

#endif // POINTSOURCE_H
