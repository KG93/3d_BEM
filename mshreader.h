#ifndef MSHREADER_H
#define MSHREADER_H

#include "triangle.h"
#include "node.h"
#include "vectortriangle.h"
#include "vectorquadrilateral.h"
#include"global.h"

#include <limits>
#include <QHash>
#include <QFile>
#include <QObject>
#include <QTextStream>
#include <iostream>
#include <QVector>
#include <QString>
#include <QRegularExpressionMatch>
#include <QRegularExpression>

/**
* \class mshReader
* \brief The mshReader class is responsible for reading the .msh version 2 asci format for mesh files.
*
* This class reads .msh files and parses the elements and nodes into internal data structures. The class provides methods
* to filter the elements that are read based on their physical and geometrical group numbers.
*/

class mshReader : public QObject
{
    Q_OBJECT
public:
    /**
    \brief Constructor for the mshReader class.
    \param[in] parent The parent QObject.
    */
    explicit mshReader(QObject *parent = 0);

    /**
    \brief Read the .msh file specified by the given filename.
    This method reads the .msh file and parses the elements and nodes into internal data structures. The method also
    allows for filtering of the elements based on their physical and geometrical group numbers.
    \param[in] filename The path to the .msh file to be read.
    \param[in] elementIndex The index of the element to be read.
    \param[in] includeAll Whether to include all elements or to filter them based on their physical and geometrical
    group numbers.
    \param[in] include The physical and geometrical group numbers to include in the mesh.
    \param[in] exclude The physical and geometrical group numbers to exclude from the mesh.
    */
    void readMsh(const QString& filename, qint64 elementIndex,  bool includeAll,const  QStringList& include,const  QStringList& exclude);

    /**
     * \brief Reads the sections of the msh file.
     */
    void readSections();

    /**
     * \brief Sorts out the includes and excludes of the msh file.
     */
    void sortOutIncludesAndExcludes();

    /**
     * \brief Returns the list of triangles read from the msh file.
     * \return List of triangles in the msh file.
     */
    QVector<VectorTriangle> getTriangles(){return triangles;}

    /**
     * \brief Returns the list of quadrilaterals read from the msh file.
     * \return List of quadrilaterals in the msh file.
     */
    QVector<VectorQuadrilateral> getQuadrilaterals(){return quadrilaterals;}

private:
    bool elementIsIncluded(quint64 physicalGroupNumber, quint64 geometricalGroupNumber);
    QStringList includeList;
    QStringList excludeList;
    QHash<quint64,Node> nodes;
    QVector<VectorTriangle> triangles;
    QVector<VectorQuadrilateral> quadrilaterals;
    QStringList physicalGroupsNames;
    QVector<int> physicalGroupsNumbers;
    QVector<quint64> physicalGroupNumbersFilterList;
    QVector<QPair<quint64,quint64>> geometricalGroupNumbersFilterList;

    quint64 numberOfNodes;
    quint64 numberOfElements;
    quint64 numberOfPhysicalNames;
    QString mshFilename;
    qint64 elementIndex;
    bool validLine;
    bool includeAll;
};

#endif // MSHREADER_H
