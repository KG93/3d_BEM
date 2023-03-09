#ifndef FIELD_H
#define FIELD_H

#include "global.h"

/**
* \brief An observation field for a sound field calculation defined by four corner nodes.
*/
class NodesField
{
public:
    NodesField(){}
    NodesField(quint64 index, quint64 node1, quint64 node2, quint64 node3, quint64 node4, QString RefNodes)
    {
        this-> index=index;
        this-> node1=node1;
        this-> node2=node2;
        this-> node3=node3;
        this-> node4=node4;
        this->RefNodes= RefNodes;
    }

    quint64 index;
    quint64 node1;
    quint64 node2;
    quint64 node3;
    quint64 node4;
    QString RefNodes;
};

#endif // FIELD_H
