#ifndef TRIANGLE_H
#define TRIANGLE_H

#include <QtGlobal>
#include <QString>

/**
* \class Triangle
* \brief  The class represents a triangle in 3d space.
*
* The Triangle class is used in an ElementSection to define a triangle. The nodes are references to a NodesSection.
*/
class Triangle
{
public:
    Triangle();
//    ~triangle();
    Triangle(qint64 elementIndex/*, quint64 numberOfPhysicalGroup*/, qint64 node1, qint64 node2, qint64 node3, const QString &nodesReference){
        this->node1 = node1;
        this->node2 = node2;
        this->node3 = node3;
        this->elementIndex = elementIndex;
        this-> nodesReference = nodesReference;

//        this->numberOfPhysicalGroup = numberOfPhysicalGroup;
    }
    //node ordering is counter clockwise for triangle front side
    qint64 elementIndex;
//    quint64 numberOfPhysicalGroup;
    qint64 node1;
    qint64 node2;
    qint64 node3;
    QString nodesReference;
};

#endif // TRIANGLE_H
