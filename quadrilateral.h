#ifndef QUADRILATERAL_H
#define QUADRILATERAL_H

#include <QtGlobal>
#include <QString>

/**
* \class Quadrilateral
* \brief  The class represents a quadrilateral in 3d space.
*
* The Quadrilateral class is used in an ElementSection to define quadrilaterals. The nodes are references to a NodesSection.
*/
class Quadrilateral
{
public:
    Quadrilateral();
    Quadrilateral(qint64 elementNumber/*, quint64 numberOfPhysicalGroup*/, qint64 node1, qint64 node2, qint64 node3, qint64 node4,const QString& nodesReference){
        this->node1 = node1;
        this->node2 = node2;
        this->node3 = node3;
        this->node4 = node4;
        this->elementIndex = elementNumber;
//        this->numberOfPhysicalGroup = numberOfPhysicalGroup;
        this-> nodesReference = nodesReference;

    }

    qint64 elementIndex;
//    quint64 numberOfPhysicalGroup;
    qint64 node1;
    qint64 node2;
    qint64 node3;
    qint64 node4;
    QString nodesReference;
};

#endif // QUADRILATERAL_H
