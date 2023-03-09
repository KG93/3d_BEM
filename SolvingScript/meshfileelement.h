#ifndef MESHFILEELEMENT_H
#define MESHFILEELEMENT_H

#include <QString>
#include <QVector>
#include <QStringList>

/**
* \class MeshFileElement
* \brief This class represents a boundary element that will be constructed from a provided meshfile.
*
* The MeshFileElement class stores information about the corresponding mesh file alias, whether to swap normals, the list of elements to include and the list of elements to exclude from the mesh file.
*/
class MeshFileElement
{
public:
    MeshFileElement();
    MeshFileElement(int index, QString MeshFileAlias, QStringList include, QStringList exclude, bool swapNormals, bool includeAll)
    {
        this->elementIndex=index;
        this->MeshFileAlias=MeshFileAlias;
        this->swapNormals=swapNormals;
        this->include=include;
        this->exclude=exclude;
        this->includeAll=includeAll;
    }

    int elementIndex;
    QString MeshFileAlias;
    bool swapNormals;
    bool includeAll;
    QStringList include;
    QStringList exclude;
};

#endif // MESHFILEELEMENT_H
