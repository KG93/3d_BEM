#ifndef FIELDSECTION_H
#define FIELDSECTION_H

#include "global.h"
#include "nodesfield.h"
#include "boundaryfield.h"
#include "meshfield.h"
#include "globallogstrings.h"

#include <QSet>
#include <QString>

/**
* \class FieldSection
* \brief The FieldSection class represents a field section in an observation script. It is used to set up a sound field calculation.
*/
class FieldSection
{
public:
    FieldSection();
    FieldSection(QString name){
        this->name=name;
    }
    bool handleScriptLine(const QString &currentLine);
    bool containsNoElements();

    QString name;
    double MeshFrequency=0;
    double EdgeLength=0;
    QString RefNodes;
    QString MeshFileAlias;
    double Alpha;
    //    QVector<int> DrvGroups;
    double Range;
    double Range_min;
    double Range_max;
    double StepSize;
    QVector<NodesField> nodesFields;
    QVector<BoundaryField> boundaryFields;
    QVector<MeshField> meshFields;
    QVector<VectorTriangle> fieldTriangles;
    QVector<VectorQuadrilateral> fieldQuadrilaterals;
private:
    bool handleFieldListLine(const QString& currentLine);
    bool handleParameterSectionLine(const QString& currentLine);

    bool validLine;
    bool readingFieldList = false;

};

#endif // FIELDSECTION_H
