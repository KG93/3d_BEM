#ifndef ELEMENTPARAMETERS_H
#define ELEMENTPARAMETERS_H

#include "mshreader.h"
#include "triangle.h"
#include "quadrilateral.h"
#include "SolvingScript/meshfileelement.h"
#include "global.h"
#include "meshfunctions.h"
#include "globallogstrings.h"

#include <math.h>
#include <QRegularExpressionMatch>
#include <QRegularExpression>
#include <QObject>
#include <QString>
#include <QVector>
#include <QSet>
#include <QMatrix4x4>

/**
* \class ElementSection
* \brief The ElementSection represents a section in a solving script that defines boundary elements used in a BEM simulation.
*
* The ElementSection class holds the boundary elements that are defined in a corresponding section in a solving script. The class contains information about the geometry and location of the elements, as well boundary conditions and other attributes.
*/
class ElementSection
{
public:
    ElementSection();

//    ElementSection(QString name, QVector<QString>& meshFiles,QVector<QString>& meshFileAlias){
//        this->name=name;
//        this->meshFiles=meshFiles;
//        this->meshFileAlias=meshFileAlias;
//    }
    ElementSection(QString name/*, const QVector<QString>& MeshFiles, const QVector<QString>& MeshFileAliases*/){
        this->name=name;
//        this->meshFiles=MeshFiles,
//        this->meshFileAliases=MeshFileAliases;
        Scale={1,1,1};
    }

    bool handleScriptLine(const QString& line);
    void sortElements(){std::sort( triangles.begin(), triangles.end(), lessThanTriangles );
                        std::sort( quadrilaterals.begin(), quadrilaterals.end(), lessThanQuadrilaterals );
                        std::sort( meshFileElements.begin(), meshFileElements.end(), lessThanMeshfileElements );
                       sorted=true;}
    bool getSorted(){return sorted;}
    bool containsNoElements();
    bool allMandatoryParametersSet();
    void swapTriangleNormals(QVector<VectorTriangle>& elementsVectorTriangles);
    void swapQuadrilateralNormals(QVector<VectorQuadrilateral>& elementsVectorQuadrilaterals);
    QVector<VectorTriangle> refineQuadrilaterals(const QVector<VectorQuadrilateral> &elementsVectorQuadrilaterals);
    void recursiveRefineQuadrilateral(VectorQuadrilateral parentQuadrilateral, int longestEdgeIndex, QVector<VectorQuadrilateral>& vectorToStoreQuadrilateral);
    QVector<VectorTriangle> refineTriangles(const QVector<VectorTriangle> &elementsVectorTriangles);
    void recursiveRefineTriangle(VectorTriangle vectorTriangle, int longestEdgeIndex, QVector<VectorTriangle>& vectorToStoreTriangles);
    void scaleAndTransformElements();
    void setTrianglesWallimpedance(int index, RobinBoundaryCondition robinBoundaryCondition);
    void setQuadrilateralsWallimpedance(int index, RobinBoundaryCondition robinBoundaryCondition);
    void setLocalWallimpedance(RobinBoundaryCondition robinBoundaryCondition);
    void setTrianglesDrivingWeight(int index, std::complex<double> weight);
    void setQuadrilateralsDrivingWeight(int index, std::complex<double> weight);
    void setLocalDrivingWeight(std::complex<double> weight);

    QString name;
    QString RefNodes;
    QString MeshFileAlias;
    unsigned int Subdomain=1;
    Eigen::Vector3d Scale=Eigen::Vector3d(1,1,1);
    Eigen::Vector3d Shift=Eigen::Vector3d(0,0,0);
    Eigen::Vector3d Rotate=Eigen::Vector3d(0,0,0);
    QMatrix4x4 transformationMatrix;
//    std::vector<int> index;
    QList<Triangle> triangles;
    QList<Quadrilateral> quadrilaterals;
    QList<MeshFileElement> meshFileElements;

    QVector<VectorTriangle> elementsVectorTriangles;
    QVector<VectorQuadrilateral> elementsVectorQuadrilaterals;

    bool ElType=0; //default interior
    enum ElTypeEnum {interior =0, exterior};
    double MeshFrequency=0;
    double EdgeLength=0;
    bool swapNormals=false;
    bool readingElementList=false;

private:
    bool handleParameterSectionLine(const QString& line);
    bool handleElementListLine(const QString& line);
    bool handleDirectElementDeclaration(const QString& line);
    bool handleMeshFileDeclaration(const QString& line);
    QRegularExpression regExWithOneValue(const QString& identifier);
    QRegularExpression regExWithOptionalValue(const QString& identifier);

    double lineLength(Eigen::Vector3d point1,Eigen::Vector3d point2);
    quint64 indexOfMaxOfThree(double a, double b, double c);
    QPair<double,quint64> maxOfFour(double a, double b, double c,double d);
    Eigen::Vector3d middleOfLine(Eigen::Vector3d point1,Eigen::Vector3d point2);
    QStringList stringToStringList(const QString& line);
    QStringList stringToStringListQuotesIntact(const QString& line);
    void removeQuotes(QString& line);
    bool validLine = false;
    bool sorted=false;

    static bool lessThanTriangles( const Triangle & triangle1, const Triangle & triangle2 )
     {
         if(triangle1.elementIndex<triangle2.elementIndex)
         {return true;}
         return false;
     }
    static bool lessThanQuadrilaterals( const Quadrilateral & quadrilateral1, const Quadrilateral & quadrilateral2 )
     {
         if(quadrilateral1.elementIndex<quadrilateral2.elementIndex)
         {return true;}
         return false;
     }
    static bool lessThanMeshfileElements( const MeshFileElement & meshElement1, const MeshFileElement & meshElement2 )
     {
         if(meshElement1.elementIndex<meshElement2.elementIndex)
         {return true;}
         return false;
     }
//    QVector<QString> meshFiles;
//    QVector<QString> meshFileAliases;

};

#endif // ELEMENTPARAMETERS_H
