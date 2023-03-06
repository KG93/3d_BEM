#ifndef OBSERVATIONPOINT_H
#define OBSERVATIONPOINT_H
#include"global.h"

class ObservationPoint
{
public:
    ObservationPoint();
    ObservationPoint(QString name,int index, Eigen::Vector3d coordinates)
    {
        this->name=name;
        this->index=index;
        this->coordinates=coordinates;
    }
    QString name;
    int index;
    std::complex<double> solution;
    Eigen::Vector3d coordinates;

};

#endif // OBSERVATIONPOINT_H
