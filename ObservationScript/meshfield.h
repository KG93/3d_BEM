#ifndef MESHFIELD_H
#define MESHFIELD_H

#include "global.h"

/**
* \brief An observation field thas is constructed from a provided mesh file.
*/
class MeshField
{
public:
    MeshField(){}
    MeshField(int index, QString MeshFileAlias, QStringList include, QStringList exclude, bool includeAll)
    {
        this-> index=index;
        this->MeshFileAlias= MeshFileAlias;
        this->include=include;
        this->exclude=exclude;
        this->includeAll=includeAll;
    }

    quint64 index;

    /**
    * \brief Alias of the mesh file from which to extract the observation field.
    */
    QString MeshFileAlias;

    /**
    * \brief ELements to be included from the mesh file.
    */
    QStringList include;

    /**
    * \brief ELements to be excluded from the mesh file.
    */
    QStringList exclude;

    /**
    * \brief If true, include all elements from the mesh file
    */
    bool includeAll;
};

#endif // MESHFIELD_H
