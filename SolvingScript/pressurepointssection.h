#ifndef PRESSUREPOINTSSECTION_H
#define PRESSUREPOINTSSECTION_H
#include "globallogstrings.h"
#include "global.h"

#include <QString>
#include <QVector>
#include <QHash>

/**
* \brief The class describes a solving script declaration of an acoustic point source in 3d space.
*/
struct PressurePoint
{
public:
  int index,nodeIndex,DrvGroup;
  double weight,delay;
  QString RefNodes;
  PressurePoint(){}
//  ~PressurePoint();
  PressurePoint(int index, int nodeIndex, int DrvGroup, double weight, double delay, QString RefNodes)
  {
      this->index=index;
      this->nodeIndex=nodeIndex;
      this->DrvGroup=DrvGroup;
      this->weight=weight;
      this->delay=delay;
      this->RefNodes=RefNodes;
  }
};

/**
* \brief The class describes a section of a solving script for declaring acoustic point sources.
*/
class PressurePointsSection
{
public:

    PressurePointsSection();
    PressurePointsSection(QString name/*, const QVector<QString>& MeshFiles, const QVector<QString>& MeshFileAliases*/){
        this->name=name;
//        this->meshFiles=MeshFiles,
//        this->meshFileAliases=MeshFileAliases;
//        Scale={1,1,1};

    }

    bool handleScriptLine(const QString& currentLine);

    QString name;
    QString RefNodes;
    QString MeshFileAlias;
    quint64 Subdomain=1;
    int  DrvGroup;
    double DrvWeight=1;
    double DrvDelay=0;
    bool pointsOnCenterOrOnElements=true;
    double MeshMinDistance=0;
//    QHash<quint64,PressurePoint> pressurePoints;
    QVector<PressurePoint> pressurePoints;

    bool readingElementList=false;

private:
    bool handleParameterSectionLine(const QString& line);
    bool handlePressurePointsListLine(const QString& line);
    bool handleDirectPressurePointDeclaration(const QString& line);
    bool handleMeshFileDeclaration(const QString& line);
    bool validLine = false;
};

#endif // PRESSUREPOINTSSECTION_H
