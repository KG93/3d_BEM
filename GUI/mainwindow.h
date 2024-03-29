#ifndef MAINWINDOW_H
#define MAINWINDOW_H

//#include "../mshreader.h"
#include "../projectfilehandler.h"
#include "openglwidget.h"
#include "../SolvingScript/solvingscriptreader.h"
#include "../ObservationScript/observationscriptreader.h"
//#include "../boundaryelements.h"
#include "../boundaryelementsolver.h"
#include "../bemsolvercontroller.h"
#include "registerscriptstab.h"
#include "logwidget.h"
#include "parameterdialog.h"
#include "freqlistwidget.h"

#include <QSize>
#include <QMainWindow>
#include <QTabWidget>
#include <iostream>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QFileInfo>
#include<QMessageBox>
#include<QSlider>
#include <QThread>

/**
* \class MainWindow
* \brief This class creates the main window of the application.

* The MainWindow class creates the main window of the application, which contains a menu bar and a tab widget.
* The tab widget is used to switch between the different application tabs. The tabs include a log window, a register scripts tab, and an OpenGL widget for visualizing the simulation results.
* The class also provides methods for setting up and handling the observation scripts, mesh files, and simulation solver.
*/

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT
public:
    explicit MainWindow(QWidget* parent = nullptr);
    ~MainWindow();

protected:
//        void resizeEvent(QResizeEvent *event) override;

private:
    Ui::MainWindow* ui;
    QTabWidget* tabWidget;
    LogWidget* generalLog;
    LogWidget* errorLog;
    QWidget* widget;
    QMenuBar* menuBar;
    RegisterScriptsTab* registerScriptsTab/*= new RegisterScriptsTab(this)*/;
    QVBoxLayout* mainLayout /*= new QVBoxLayout*/;
    openglWidget* openGlWidget /*= new openglWidget(this)*/;
    QHBoxLayout* openGlTabLayout;
    QSlider* maxSlider;
    QSlider* minSlider;
    ParameterDialog* parameterDialog;
    ProjectFileHandler* projFileHandler = new ProjectFileHandler;
    SolvingScriptReader* solvScriptReader;
    ObservationScriptReader* obsScriptReader;
    BemSolverController* bemSolverController = new BemSolverController;

    FreqListWidget* freqWidget;
    QVector<double> freq;
    int currentSelectedFreq = -1;
    QVector<BoundaryElements> boundaryElementsFreq;

    QString projectFileName;

    void setUpMenu(); /** Set up the GUI menu */
//    void updateMenu(); /** Update the GUI menu */
    void setUpObservationScripts(const QStringList &observationFiles); /** Set the list of observation scripts. */
    void setUpMeshFiles(const QStringList &meshFiles, const QStringList &meshFileAlias); /** Set the list of registered meshfiles. */
    QStringList getMeshFiles(); /** Get the list of observation scripts. */
    QStringList getObservationscripts(); /** Get the list of registered meshfiles. */

    void setUpBemSolverSignalsAndSlots();
    bool solutionMenuAlreadySetup = false;
    bool showSolutionAndFieldMenu = false;
    bool showCalculateFieldButton = false;
    bool showSolutionIfAvailable = true; /** If true, the OpenGl viewer shows the solution on the boundary, if it is already calculated*/


    QFileInfo solvinScript;
    QVector<QFileInfo> observationScripts;
    QVector<QFileInfo> meshFiles;
    QStringList meshFileAlias;
    bool hasBeenSolved;

signals:
    void slidersChanged(int min, int max); /** \brief Signals the OpenGl widget the changed range for the sound pressure visualization. */

private slots:
    void newProject(); /** Create a new project. */
    void openProject(); /** Open an existing project. */
    void loadProject(); /** Handles the actual work of (re-)loading a project. */
    void loadObservationScript(); /** (Re-)load the observation script. */
    void handle3dViewerLogic(); /** Handle the logic for the 3d viewer. */
    void calculateSolution(); /** Use the boundary element solver to calculate the solution. */
    void showSolverParameterDialog(); /** Set the parameters for the solver. */
    void showFreqSelector(); /** Open the frequency selector widget. */
    void showSolutionOnBoundaryElements(); /** Signal the opengl widget to show the solution on the boundary. */
    void hideSolutionOnBoundaryElements(); /** Signal the opengl widget to hide the solution on the boundary. */
    void calculateSolutionOnField(); /** Use the boundary element solver to calculate the field solution. */
    void showField(); /** Signal the opengl widget to show the field solution. */
    void hideField(); /** Signal the opengl widget to hide the field solution. */
    void showLegend(); /** Signal the opengl widget to show the legend. */
    void hideLegend(); /** Signal the opengl widget to hide the legend. */
    void drawPressure(); /** Signal the opengl widget to show the solution pressure. */
    void drawPhase(); /** Signal the opengl widget to show the solution phase. */
    void animatePhase(); /** Signal the opengl widget to animate the solution pressure. */
    void stopAnimation(); /** Signal the opengl widget to stop the pressure animation. */
    void useLocalPressureMaximaForSolutionColoring(); /** Signal the opengl widget to use the local pressure extrema for the solution vizualizaition. */
    void useGlobalPressureMaximaForSolutionColoring(); /** Signal the opengl widget to use the global pressure extrema for the solution vizualizaition. */
    void saveImage(); /** Signal the opengl widget capture a screenshot. */
    void updateScripts(); /** Update the script list here and in the project file. */
    void readSolvingScript(); /** Read and process the solving script. */
    void reloadProjectQuery(); /** Query in case the already solved project is to be reloaded by user input. */
    void preventSliderCrossingMax();
    void preventSliderCrossingMin();
};
#endif // MAINWINDOW_H
