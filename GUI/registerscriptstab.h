#ifndef REGISTERSCRIPTSTAB_H
#define REGISTERSCRIPTSTAB_H
#include <QFileInfo>
#include <QMenuBar>
#include <QWidget>
#include <QLabel>
#include <QLineEdit>
#include <QVBoxLayout>
#include <QListWidget>
#include <QMouseEvent>
#include <QFileDialog>
#include <QInputDialog>
#include <QDesktopServices>
#include <iostream>

/**
* \class RegisterScriptsTab
* \brief This class provides a widget to register solving and observation scripts and mesh files.
*
* The RegisterScriptsTab class provides a widget to register solving and observation scripts and mesh files. It also allows the user to edit the aliases of the mesh files.
*/
class RegisterScriptsTab: public QWidget
{
    Q_OBJECT
public:
//    registerScriptsTab();
    explicit RegisterScriptsTab(/*const QFileInfo &solvingScript,*/ QWidget *parent = nullptr);
    QFileInfo solvingScript;
    QVector<QFileInfo> observationScripts;
    QVector<QFileInfo> meshFiles;
    QStringList meshFileAliases;
//    void update();
    void setSolvScript(const QFileInfo& solvingScript);
    void setObsScripts(const QVector<QFileInfo> &obsFileInfos);
    void setMeshFiles(const QVector<QFileInfo>& meshFiles, const QStringList &meshFileAliases);
    QFileInfo getSolvingScript(){return solvingScript;}
    QVector<QFileInfo> getObservationScripts(){return observationScripts;}
    QVector<QFileInfo> getMeshFiles(){return meshFiles;}
    QStringList getMeshFileAliases(){return meshFileAliases;}
    QStringList getObservationScriptPaths();
    QStringList getMeshFilePaths();
    void clean();
private:
    QLabel *solvScriptLabel ;
    QLabel *pathObsValueLabel ;
    QLabel *meshFilesLabel ;
//    QListWidget *solvingScriptWidget;
    QListWidget *meshListWidget;
    QListWidget *observationListWidget;

    void updateSolvLabel();
    void updateMeshFiles();
    void updateObservationScripts();

private slots:
    void mousePressEvent(QMouseEvent *event);
    void showSolvContextMenu(const QPoint& pos);
    void showMeshContextMenu(const QPoint& pos);
    void showObsContextMenu(const QPoint& pos);
    void addSolvScript();
    void replaceSolvScript();
    void addMeshItem();
    void editMeshItemAlias();
    void eraseMeshItem();
    void addObservationItem();
    void replaceObservationItem();
    void eraseObservationItem();
    void solvScriptInEditor();
    void obsScriptInEditor();
    void meshFileInEditor();
signals:
    void scriptChanged();
    void obsScriptChanged();
};


#endif // REGISTERSCRIPTSTAB_H
