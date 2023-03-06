#include "pointsource.h"

//PointSource::PointSource()
//{
//}

std::complex<double> PointSource::observedPhiPlusDPhi(const Eigen::Vector3d listeningPosition, const Eigen::Vector3d up , const std::complex<double> wavenumber,const std::complex<double> couplingParameter)
{
    Eigen::Vector3d rVector=listeningPosition-position; //p-q
//    double upNq=up.dot(boundaryElements.triangles.at(column).normal);
    double r=rVector.norm();
//    double rr=r*r;
//            double rrr=rr*r;
//    double rnq=-1.0*(rVector.dot(boundaryElements.triangles.at(column).normal))/r;
    double rup=-1.0*(rVector.dot(up))/r; //???
    std::complex<double> ikr=imaginaryUnit*wavenumber*r;
    std::complex<double> greensFunction=weight*(std::exp(ikr))/(PI4*r);
    std::complex<double> dGreensDR=(greensFunction/r)*(ikr-1.0);
    return greensFunction+couplingParameter*dGreensDR*rup;
}
