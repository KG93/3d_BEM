#ifndef DRIVINGELEMENT_H
#define DRIVINGELEMENT_H

#include <QString>
#include <complex>

class DrivingElement
{
public:
    DrivingElement(){}
    DrivingElement(quint64 index, qint64 elementIndex, QString RefElements, int DrvGroup, std::complex<double> Weight, double Delay){
        this->index=index;
        this->elementIndex=elementIndex;
        this-> RefElements=RefElements;
        this-> DrvGroup=DrvGroup;
        this-> Weight=Weight;
        this-> Delay=Delay;
    }

    quint64 index;
    qint64 elementIndex;
    QString RefElements;
    int DrvGroup;
    std::complex<double> Weight;
    double Delay;
};
#endif // DRIVINGELEMENT_H
