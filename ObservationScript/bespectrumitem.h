#ifndef BESPECTRUMITEM_H
#define BESPECTRUMITEM_H
#include <QString>
#include <QVector>
#include <eigen3/Eigen/Geometry>


class BESpectrumItem{
public:
    BESpectrumItem(){}
    BESpectrumItem(int index, int nodeIndex, QString RefNodes, /*QVector<int> DrvGroups,*/ QString ID)
    {
        this->index=index;
        this->nodeIndex=nodeIndex;
        this->RefNodes=RefNodes;
//        this->DrvGroups=DrvGroups;
        this->ID=ID;
    }
    int index;
    int nodeIndex;
    QString RefNodes;
//    QVector<int> DrvGroups;
    QString ID;
};

class BESpectrumPoint{
public:
    BESpectrumPoint(){}
    BESpectrumPoint(int index, Eigen::Vector3d coordinates)
    {
        this->index=index;
        this->coordinates=coordinates;
    }

    int index;
    Eigen::Vector3d coordinates;

};

#endif // BESPECTRUMITEM_H
