#ifndef BOUNDARYFIELD_H
#define BOUNDARYFIELD_H

#include "global.h"

class BoundaryField
{
public:
    BoundaryField(){}
    BoundaryField(quint64 index, QString Refelements)
    {
        this-> index=index;
        this->Refelements= Refelements;
    }

    quint64 index;
    QString Refelements;

};

#endif // BOUNDARYFIELD_H
