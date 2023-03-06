#ifndef OBSERVATIONSCRIPTREADER_H
#define OBSERVATIONSCRIPTREADER_H
#include "nodesSection.h"
#include "meshfilepropertiessection.h"
#include "controlfieldsection.h"
#include "controlspectrumsection.h"
#include "fieldsection.h"
#include "bespectrumsection.h"
#include "global.h"
#include "mshreader.h"
#include "meshfunctions.h"
#include "observationpoint.h"
#include "observationfield.h"

#include <QTextStream>
#include <QFile>


class ObservationScriptReader: public QObject
{
    Q_OBJECT
public:

    ObservationScriptReader(QObject* parent=0);
    bool setMeshFilesAndAliases( const QStringList& meshFiles, const QStringList& meshFileAlias);
    void readObservationScript(const QStringList& fileName);
    void clear();
    QStringList getNodesSectionsNames();
    void checkNodesIdentifiersAndSort();
    bool removeInvalidElements();
    void setElementSectionsNames(QStringList ElementSectionsNames){
        this->elementSectionsNames=ElementSectionsNames;
    }
    void setupObservationElements();
    void setEdgeLength(double edgeLength){this->EdgeLength=edgeLength;}
    void setWavespeed(double waveSpeed){this->c=waveSpeed;}

    QVector<ObservationPoint> getObservationPoints();
    QVector<ObservationField> getObservationFields();




    ControlFieldSection controlFieldSection;
    ControlSpectrumSection controlSpectrumSection;
    QVector<NodesSection> nodesSections;
    QVector<FieldSection> fieldSections;
    QVector<BESpectrumSection> bESpectrumSections;
    QVector<MeshFilePropertiesSection> meshFilePropertiesSections;


    QStringList meshFiles;
    QStringList meshFileAlias;
    QStringList elementSectionsNames;

private:
    bool newSectionDeclared(const QString& Line);
    QRegularExpression sectionIdentifierRegEx(const QString& identifier);
    void handleLineAccordingToCurrentSection(const QString& Line);
    enum sections {none, Control_Spectrum, Control_Field, Field, BE_Spectrum, Nodes, Driving_Values, MeshFile_Properties, Radiation_Impedance};
    int currentSection=none;
    bool controlSpectrumSectionDeclared;
    bool controlFieldSectionDeclared;
    bool meshfilesAndAliasesSet;
    double EdgeLength;
    double c;

    mshReader meshFileReader;

signals:
    void errorMessage(QString err);
    void logMessage(QString log);
};

#endif // OBSERVATIONSCRIPTREADER_H
