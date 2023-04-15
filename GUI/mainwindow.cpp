#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget* parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    widget = new QWidget(this);
    mainLayout = new QVBoxLayout;
    tabWidget = new QTabWidget();
    menuBar = new QMenuBar();
    openGlWidget = new openglWidget();
    openGlTabLayout = new QHBoxLayout();
    QWidget* openQlContainerWidget = new QWidget(tabWidget);
    openQlContainerWidget->setLayout(openGlTabLayout);
    tabWidget->addTab(openQlContainerWidget,"3d Viewer");
    openGlTabLayout->addWidget(openGlWidget);
    maxSlider = new QSlider();
    maxSlider->setFixedWidth(6);
    maxSlider->setRange(1,100);
    maxSlider->setValue(100);
    connect(maxSlider, SIGNAL(valueChanged(int)), this, SLOT(preventSliderCrossingMax()));
    minSlider = new QSlider();
    minSlider->setFixedWidth(6);
    connect(minSlider, SIGNAL(valueChanged(int)), this, SLOT(preventSliderCrossingMin()));
    connect(this,SIGNAL(slidersChanged(int,int)),openGlWidget,SLOT(getMinAndMaxCut(int,int)));
    openGlTabLayout->addWidget(minSlider);
    openGlTabLayout->addWidget(maxSlider);
    openGlTabLayout->setContentsMargins(QMargins(0,0,0,0));
    registerScriptsTab = new RegisterScriptsTab(/*solvinScript,*/tabWidget);
    connect(registerScriptsTab, SIGNAL(scriptChanged()), this, SLOT(updateScripts()));
    connect(registerScriptsTab, SIGNAL(obsScriptChanged()), this, SLOT(updateScripts()));
    solvScriptReader = new SolvingScriptReader(this);
    obsScriptReader = new ObservationScriptReader(this);
    generalLog = new LogWidget;
    connect(solvScriptReader, SIGNAL(logMessage(QString)), generalLog, SLOT(writeOutLogString()));
    connect(obsScriptReader, SIGNAL(logMessage(QString)), generalLog, SLOT(writeOutLogString()));
    connect(bemSolverController, SIGNAL(updateLog()), generalLog, SLOT(writeOutLogString()));
    errorLog = new LogWidget;
    connect(solvScriptReader, SIGNAL(errorMessage(QString)), errorLog, SLOT(writeOutErrorLogString()));
    connect(obsScriptReader, SIGNAL(errorMessage(QString)), errorLog, SLOT(writeOutErrorLogString()));
    connect(bemSolverController, SIGNAL(updateLog()), errorLog, SLOT(writeOutErrorLogString()));

    tabWidget->addTab(registerScriptsTab,"Scripts");
    tabWidget->addTab(generalLog,"Log");
    tabWidget->addTab(errorLog,"Error log");
    mainLayout->setMenuBar(menuBar);
    mainLayout->addWidget(tabWidget);
    mainLayout->setContentsMargins(QMargins(3,3,3,3));
    widget->setLayout(mainLayout);
    setCentralWidget(widget);
    setUpMenu();

    setWindowTitle(tr("BEM"));
    widget->show();

    parameterDialog = new ParameterDialog(widget);
    freqWidget = new FreqListWidget();

    setUpBemSolverSignalsAndSlots();

    if(global::activeProgramming)
    {
//        QString filename = "/home/lrak/workspace/3d_BEM/Example/Loudspeaker/Project.abec"; //loudspeaker
        QString filename = "/home/lrak/workspace/3d_BEM/Example/Sphere/sphere.abec"; //sphere

        projFileHandler->readProjectFile(filename);
        solvinScript = QFileInfo(projFileHandler->getSolvingScript());
        setUpObservationScripts(projFileHandler->getObservScriptList());
        setUpMeshFiles(projFileHandler->getMeshFileList(),projFileHandler->getMeshFileAliasList());
        loadProject();
        BoundaryElements elements = bemSolverController->getBoundaryElements();
        openGlWidget->setBoundaryElements(elements);
        std::cout<<"Setting up boundary elements in gl widget."<< std::endl;
    }
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::setUpMenu()
{
    menuBar->clear();
    QMenu* menu = menuBar->addMenu(tr("Project"));
    QAction* newProjectAction = new QAction(tr("&New project"), this);
    QAction* openProjectAction = new QAction(tr("&Open project"), this);

    menu->addAction(newProjectAction);
    connect(newProjectAction, SIGNAL(triggered()), SLOT(newProject()));

    menu->addAction(openProjectAction);
    connect(openProjectAction, SIGNAL(triggered()), SLOT(openProject()));

    if(showSolutionAndFieldMenu)
    {
        QMenu* menu1 = menuBar->addMenu(tr("Viewer"));
        QAction* showFrequencySelectorAction = new QAction(tr("&Show frequency selector"), this);
        QAction* showSolutionAction = new QAction(tr("&Show solution on boundary"), this);
        QAction* hideSolutionAction = new QAction(tr("&Hide solution on boundary"), this);
        QAction* showFieldAction = new QAction(tr("&Show field"), this);
        QAction* hideFieldAction = new QAction(tr("&Hide field"), this);
        QAction* showLegendAction = new QAction(tr("&Show legend"), this);
        QAction* hideLegendAction = new QAction(tr("&Hide legend"), this);
        QAction* drawPhaseAction = new QAction(tr("&Show phase"), this);
        QAction* drawPressureAction = new QAction(tr("&Show pressure"), this);
        QAction* animatePhaseAction = new QAction(tr("&Animate phase"), this);
        QAction* stopAnimation = new QAction(tr("&Stop animation"), this);
        QAction* showGlobalDb = new QAction(tr("&Global extrema for coloring"), this);
        QAction* showLocalDb = new QAction(tr("&Local extrema for coloring"), this);

        QAction* saveImageAction = new QAction(tr("&Save image"), this);

        menu1->addAction(showFrequencySelectorAction);

        if(openGlWidget->getShowSolutionValue() && bemSolverController->getSolutionState() == AtLeastOneSolved)
        {
            menu1->addAction(hideSolutionAction);
        }
        else if (bemSolverController->getSolutionState() == AtLeastOneSolved)
        {
            menu1->addAction(showSolutionAction);
        }
        if(openGlWidget->getDrawField())
        {
            menu1->addAction(hideFieldAction);
        }
        else
        {
            menu1->addAction(showFieldAction);
        }
        if(openGlWidget->getDrawLegend())
        {
            menu1->addAction(hideLegendAction);
        }
        else
        {
            menu1->addAction(showLegendAction);
        }
        if((openGlWidget->getShowSolutionValue() || openGlWidget->getDrawField()) && openGlWidget->getShowPhaseValue() && bemSolverController->getSolutionState() == AtLeastOneSolved)
        {
            menu1->addAction(drawPressureAction);
            if(openGlWidget->getAnimatePhase())
            {
                menu1->addAction(stopAnimation);
            }
            else
            {
                menu1->addAction(animatePhaseAction);
            }
        }
        else if((openGlWidget->getShowSolutionValue() || openGlWidget->getDrawField()) && bemSolverController->getSolutionState() == AtLeastOneSolved)
        {
            menu1->addAction(drawPhaseAction);
            if(openGlWidget->getUseGlobalMinMaxDb())
            {
                menu1->addAction(showLocalDb);
            }
            else
            {
                menu1->addAction(showGlobalDb);
            }
        }
        menu1->addAction(saveImageAction);
        connect(showFrequencySelectorAction, SIGNAL(triggered()), SLOT(showFreqSelector()));
        connect(showSolutionAction, SIGNAL(triggered()), SLOT(showSolutionOnBoundaryElements()));
        connect(hideSolutionAction, SIGNAL(triggered()), SLOT(hideSolutionOnBoundaryElements()));
        connect(showFieldAction, SIGNAL(triggered()), SLOT(showField()));
        connect(hideFieldAction, SIGNAL(triggered()), SLOT(hideField()));
        connect(showLegendAction, SIGNAL(triggered()), SLOT(showLegend()));
        connect(hideLegendAction, SIGNAL(triggered()), SLOT(hideLegend()));
        connect(drawPressureAction, SIGNAL(triggered()), SLOT(drawPressure()));
        connect(drawPhaseAction, SIGNAL(triggered()), SLOT(drawPhase()));
        connect(animatePhaseAction, SIGNAL(triggered()), SLOT(animatePhase()));
        connect(stopAnimation, SIGNAL(triggered()), SLOT(stopAnimation()));
        connect(showLocalDb, SIGNAL(triggered()), SLOT(useLocalPressureMaximaForSolutionColoring()));
        connect(showGlobalDb, SIGNAL(triggered()), SLOT(useGlobalPressureMaximaForSolutionColoring()));
        connect(saveImageAction, SIGNAL(triggered()), SLOT(saveImage()));

        QMenu* menu2 = menuBar->addMenu(tr("Solving"));

        QAction* calculateSolutionAction = new QAction(tr("&Calculate solution on boundary"), this);
        QAction* setParametersAction = new QAction(tr("&Set solver parameters"), this);

        if(bemSolverController->getControllerState() == NotWorking)
        {
            menu2->addAction(calculateSolutionAction);
        }
        menu2->addAction(setParametersAction);

        connect(calculateSolutionAction, SIGNAL(triggered()), SLOT(calculateSolution()));
        connect(setParametersAction, SIGNAL(triggered()), SLOT(showSolverParameterDialog()));

        QMenu* menu3 = menuBar->addMenu(tr("Field"));
        if(showCalculateFieldButton && bemSolverController->getControllerState() == NotWorking)
        {
            QAction* calculateFieldAction = new QAction(tr("&Calculate solution on field"), this);
            menu3->addAction(calculateFieldAction);
            connect(calculateFieldAction, SIGNAL(triggered()), SLOT(calculateSolutionOnField()));
        }       

        QAction* reloadProjectAction = new QAction(tr("Reload Project"), this);
        connect(reloadProjectAction, SIGNAL(triggered()), SLOT(reloadProjectQuery()));
        QAction* reloadObsScript = new QAction(tr("Reload ObservationScript"), this);
        connect(reloadObsScript, SIGNAL(triggered()), SLOT(loadObservationScript()));
        menuBar->addAction(reloadProjectAction);
        menuBar->addAction(reloadObsScript);

        BemControllerStates controllerState = bemSolverController->getControllerState();
        if(controllerState == Working)
        {
            QAction* stopSolverThread = new QAction(tr("Stop Solver"), this);
            connect(stopSolverThread, SIGNAL(triggered()), bemSolverController, SLOT(terminateThread()));
            menuBar->addAction(stopSolverThread);
        }
    }
}

void MainWindow::newProject()
{
    QString filename = QFileDialog::getSaveFileName(this, "Creates new .abec project file.","",tr("Text files (*.abec)"));
                                                                                                    // filter, für nur .abec dateien
      if ( filename.isEmpty() )
      {
          std::cout << "Empty filename." << std::endl;
          return;
      }
      if(!filename.endsWith(".abec",Qt::CaseInsensitive))
      {
          filename.append(".abec");
      }
      bemSolverController->clear();

      showSolutionAndFieldMenu = true;
      showCalculateFieldButton = false;
      projFileHandler->createProjectFile(filename);
      setUpMenu();
}

void MainWindow:: openProject()
{
    QString filename = QFileDialog::getOpenFileName(this, "Select .abec project file.","",tr("Text files (*.abec)"));
                                                                                                    // filter, für nur .abec dateien
      if ( filename.isEmpty() )
      {
          std::cout << "Empty filename: " << std::endl;
          return;
      }
      bemSolverController->clear();
      if(global::activeProgramming)
      {
          std::cout << "Project path: " << filename.toStdString() << std::endl;
      }
      projFileHandler->readProjectFile(filename);

      solvinScript = QFileInfo(projFileHandler->getSolvingScript());
      setUpObservationScripts(projFileHandler->getObservScriptList());
      setUpMeshFiles(projFileHandler->getMeshFileList(), projFileHandler->getMeshFileAliasList());

      loadProject();
}

void MainWindow::loadProject()
{
    registerScriptsTab->clean();
    registerScriptsTab->setSolvScript(solvinScript);
    registerScriptsTab->setObsScripts(observationScripts);
    registerScriptsTab->setMeshFiles(meshFiles,meshFileAlias);
    solvScriptReader->setMeshFilesAndAliases(projFileHandler->getMeshFileList(), projFileHandler->getMeshFileAliasList());
    solvScriptReader->readSolvingScript(solvinScript.absoluteFilePath());
    solvScriptReader->checkNodesIdentifiersAndSort();
    solvScriptReader->removeInvalidElements();
    solvScriptReader->setupConstBoundaryElements();
    solvScriptReader->setupPressurePoints();
    solvScriptReader->refineElements();

    freq = solvScriptReader->getFrequencies();
    freqWidget->setFrequencies(freq);

    boundaryElementsFreq.resize(freq.length());
    for(int i=0; i<freq.length(); i++)
    {
        boundaryElementsFreq[i] = solvScriptReader->getBoundaryElements();
    }

    obsScriptReader->setElementSectionsNames(solvScriptReader->getElementSectionsNames());
    obsScriptReader->setMeshFilesAndAliases(projFileHandler->getMeshFileList(), projFileHandler->getMeshFileAliasList());
    obsScriptReader->setEdgeLength(solvScriptReader->getGlobalEdgelength());
    obsScriptReader->setWavespeed(solvScriptReader->getWavespeed());
    obsScriptReader->readObservationScript(projFileHandler->getObservScriptList());
    obsScriptReader->checkNodesIdentifiersAndSort();
    obsScriptReader->removeInvalidElements();
    obsScriptReader->setupObservationElements();
    bemSolverController->setBoundaryElements(solvScriptReader->getBoundaryElements());
    bemSolverController->setPointSources(solvScriptReader->getPointSources());
    bemSolverController->setObservationFields(obsScriptReader->getObservationFields());
//    bemSolverController->setObservationPoints(obsScriptReader->getObservationPoints());
//    boundaryElementSolver->setFrequency(solvScriptReader->getFrequencies());
//    boundaryElementSolver->prepareBoundaryElements();
    openGlWidget->setBoundaryElements(solvScriptReader->getBoundaryElements());
    openGlWidget->setCenterOfMass(solvScriptReader->getCenterOfMass());
    openGlWidget->setCameraRadius(solvScriptReader->getContainingRadius());
    openGlWidget->setElementSections(solvScriptReader->getElementsSections());
    openGlWidget->setPointSources(solvScriptReader->getPointSources());
    openGlWidget->setObservationFields(obsScriptReader->getObservationFields());
    openGlWidget->setObservationPoints(obsScriptReader->getObservationPoints());
    openGlWidget->setShowSolutionValue(false);

    showSolutionAndFieldMenu = true;
    showCalculateFieldButton = false;
    bemSolverController->clear();
    setUpMenu();
    update();
    openGlWidget->update();
}

void MainWindow::loadObservationScript()
{
    bemSolverController->clearFields();

    obsScriptReader->setElementSectionsNames(solvScriptReader->getElementSectionsNames());
    obsScriptReader->setMeshFilesAndAliases(projFileHandler->getMeshFileList(), projFileHandler->getMeshFileAliasList());
    obsScriptReader->setEdgeLength(solvScriptReader->getGlobalEdgelength());
    obsScriptReader->setWavespeed(solvScriptReader->getWavespeed());
    obsScriptReader->readObservationScript(projFileHandler->getObservScriptList());
    obsScriptReader->checkNodesIdentifiersAndSort();
    obsScriptReader->removeInvalidElements();
    obsScriptReader->setupObservationElements();

    bemSolverController->setObservationFields(obsScriptReader->getObservationFields());
//    bemSolverController->setObservationPoints(obsScriptReader->getObservationPoints());
    openGlWidget->setObservationFields(obsScriptReader->getObservationFields());
    openGlWidget->setObservationPoints(obsScriptReader->getObservationPoints());
    openGlWidget->update();
}

void MainWindow::handle3dViewerLogic()
{
    freqWidget->show();

    unsigned int currentSelectedFreqIndex = freqWidget->getCurrentSelectedFreqIndex();
    BemSolutionStates controllerState = bemSolverController->getSolutionState();

    switch(controllerState)
    {
        case NothingSolved:
            openGlWidget->setShowSolutionValue(false);
            openGlWidget->setDrawField(false);
            break;
        case AtLeastOneSolved:
            if(bemSolverController->isBoundarySolvedAtFreq(currentSelectedFreqIndex) && showSolutionIfAvailable)
            {
                openGlWidget->setBoundaryElements(bemSolverController->getBoundaryElements());
                BoundarySolution solution = bemSolverController->getBoundarySolutionAtFrequencyIndex(currentSelectedFreqIndex);
                openGlWidget->setNewBoundarySolution(solution);
                showCalculateFieldButton = true;
                openGlWidget->setShowSolutionValue(true);
            }
            else
            {
                openGlWidget->setShowSolutionValue(false);
            }
            if(bemSolverController->isFieldSolvedAtFreq(currentSelectedFreqIndex))
            {
                openGlWidget->setBoundaryElements(bemSolverController->getBoundaryElements());
                FieldSolutions solutions = bemSolverController->getFieldSolutionAtFrequencyIndex(currentSelectedFreqIndex);
                openGlWidget->setObservationFields(bemSolverController->getObservationFields());
                openGlWidget->setNewFieldSolution(solutions);
            }
            break;
    }
    setUpMenu();
    openGlWidget->update();
    freqWidget->show();
}

void MainWindow::calculateSolution()
{
    bemSolverController->setSolverObjective(SolveBoundary);
    bemSolverController->setBoundaryElements(solvScriptReader->getBoundaryElements());
    bemSolverController->setPointSources(solvScriptReader->getPointSources());
    bemSolverController->setBemSolverParameters(parameterDialog->getBemSolverParameters());
    BemParameters bemParameters = {0,solvScriptReader->getWavespeed(),solvScriptReader->getAirDensity()};
    bemSolverController->setBemParameters(bemParameters);
    bemSolverController->setFrequencies(freqWidget->getFrequencies());

    bemSolverController->start(); // start threaded BEM calculation
    openGlWidget->setBoundaryElements(bemSolverController->getBoundaryElements());
    setUpMenu();
}

void MainWindow::showSolverParameterDialog() // Set the parameters for the solver.
{
    parameterDialog->show();
}

void MainWindow::showFreqSelector()
{
    freqWidget->show();
}

void MainWindow::showSolutionOnBoundaryElements()
{
    showSolutionIfAvailable = true;
    handle3dViewerLogic();
    setUpMenu();
    openGlWidget->update();
}

void MainWindow::hideSolutionOnBoundaryElements()
{
    showSolutionIfAvailable = false;
    handle3dViewerLogic();
    setUpMenu();
    openGlWidget->update();
}

void MainWindow::calculateSolutionOnField()
{
    bemSolverController->setSolverObjective(SolveField);
    bemSolverController->setBemSolverParameters(parameterDialog->getBemSolverParameters());
    bemSolverController->setObservationFields(obsScriptReader->getObservationFields());
    bemSolverController->start(); // start threaded BEM calculation

    openGlWidget->setObservationFields(bemSolverController->getObservationFields());
    openGlWidget->setBoundaryElements(bemSolverController->getBoundaryElements());
    openGlWidget->setDrawField(true);

    setUpMenu();
}

void MainWindow::showField()
{
    openGlWidget->setDrawField(true);
    setUpMenu();
    openGlWidget->update();
}

void MainWindow::hideField()
{
    openGlWidget->setDrawField(false);
    setUpMenu();
    openGlWidget->update();
}

void MainWindow::showLegend()
{
    openGlWidget->setDrawLegend(true);
    setUpMenu();
    openGlWidget->update();
}

void MainWindow::hideLegend()
{
    openGlWidget->setDrawLegend(false);
    setUpMenu();
    openGlWidget->update();
}

void MainWindow::drawPressure()
{
    openGlWidget->setShowPhaseValue(false);
    setUpMenu();
    openGlWidget->update();
}

void MainWindow::drawPhase()
{
    openGlWidget->setShowPhaseValue(true);
    setUpMenu();
    openGlWidget->update();
}

void MainWindow::animatePhase()
{
    openGlWidget->setShowPhaseValue(true);
    openGlWidget->setAnimatePhase(true);
    setUpMenu();
    openGlWidget->update();
}

void MainWindow::stopAnimation()
{
    openGlWidget->setAnimatePhase(false);
    setUpMenu();
    openGlWidget->update();
}

void MainWindow::useLocalPressureMaximaForSolutionColoring()
{
    openGlWidget->setUseGlobalMinMaxDb(false);
    setUpMenu();
}

void MainWindow::useGlobalPressureMaximaForSolutionColoring()
{
    openGlWidget->setUseGlobalMinMaxDb(true);
    setUpMenu();
}

void MainWindow::saveImage()
{
    QString filename = QFileDialog::getSaveFileName(this, "Save 3d viewer image.","",tr("Image files (*.png)"));
                                                                                                // filter, für nur .png dateien
    if ( filename.isEmpty() )
    {
      std::cout <<"Empty filename."<<std::endl;
      return;
    }
    if(!filename.endsWith(".png",Qt::CaseInsensitive))
    {
      filename.append(".png");
    }

//    QSize windowSize = openGlWidget->size();
//    openGlWidget->resize(2*windowSize);
    openGlWidget->grabFramebuffer().save(filename,"PNG",100);
//    openGlWidget->resize(windowSize);
}

void MainWindow::setUpObservationScripts(const QStringList &observationFiles)
{
    observationScripts.resize(observationFiles.size());
    for(int i=0; i<observationFiles.size(); i++)
    {
        observationScripts[i] = QFileInfo(observationFiles.at(i));
    }
}

void MainWindow::setUpMeshFiles(const QStringList &MeshFiles,const QStringList &MeshFileAlias)
{
    if(MeshFiles.size() != MeshFileAlias.size())
    {
        std::cout<<"MeshFiles.size() != meshFileAlias.size() in MainWindow class."<<std::endl;
    }
    else
    {
        meshFiles.resize(MeshFiles.size());
        for(int i=0; i<meshFiles.size(); i++)
        {
            meshFiles[i] = QFileInfo(MeshFiles.at(i));
        }
        this->meshFileAlias = MeshFileAlias;
    }
}

void MainWindow::updateScripts()
{
    solvinScript = registerScriptsTab->getSolvingScript();
    observationScripts = registerScriptsTab->getObservationScripts();
    meshFiles = registerScriptsTab->getMeshFiles();
    meshFileAlias = registerScriptsTab->getMeshFileAliases();
    projFileHandler->setSolvingScript(solvinScript.absoluteFilePath());
    projFileHandler->setObservScriptList(getObservationscripts());
    projFileHandler->setMeshFileList(getMeshFiles());
    projFileHandler->setMeshFileAliasList(meshFileAlias);
    projFileHandler->writeProjectFile();
}

void MainWindow::readSolvingScript()
{
    solvScriptReader->setMeshFilesAndAliases(projFileHandler->getMeshFileList(), projFileHandler->getMeshFileAliasList());
    solvScriptReader->readSolvingScript(projFileHandler->getSolvingScript());
    solvScriptReader->checkNodesIdentifiersAndSort();
    solvScriptReader->removeInvalidElements();
    solvScriptReader->setupConstBoundaryElements();
    solvScriptReader->setupPressurePoints();
    solvScriptReader->refineElements();
}

QStringList MainWindow::getMeshFiles()
{
    QStringList meshfileList;
    for(int i=0; i<meshFiles.size(); i++)
    {
        meshfileList.append(meshFiles.at(i).absoluteFilePath());
    }
    return meshfileList;
}

QStringList MainWindow::getObservationscripts()
{
    QStringList obsScriptList;
    for(int i=0; i<observationScripts.size(); i++)
    {
        obsScriptList.append(observationScripts.at(i).absoluteFilePath());
    }
    return obsScriptList;
}

void MainWindow::setUpBemSolverSignalsAndSlots()
{
    connect(bemSolverController, SIGNAL(boundarySolvedForFreqIndex(int)), freqWidget, SLOT(colorSolvedFreq(int)));
//    connect(bemSolverController, SIGNAL(boundarySolved(Eigen::VectorXcd, Eigen::VectorXcd, Eigen::VectorXcd)), freqWidget, SLOT(colorFreqSolved(int)));
    connect(freqWidget, SIGNAL(currentRowChanged(int)), this, SLOT(handle3dViewerLogic()));
    connect(bemSolverController, SIGNAL(boundarySolvedForFreqIndex(int)), this, SLOT(handle3dViewerLogic()));
    connect(bemSolverController, SIGNAL(finished()), this, SLOT(handle3dViewerLogic()));
}

void MainWindow::reloadProjectQuery()
{
    if(bemSolverController->getSolutionState() == AtLeastOneSolved)
    {
        QMessageBox::StandardButton reply;
        reply = QMessageBox::question(this, "", "Project has already been solved.\r\n Reload Project?", QMessageBox::Yes|QMessageBox::No);
        if (reply == QMessageBox::Yes)
        {
            loadProject();
        }
    }
    else
    {
        loadProject();
    }
}

void MainWindow::preventSliderCrossingMax()
{
    if(maxSlider->value() <= minSlider->value())
    {
        if(minSlider->value() == minSlider->minimum())
        {
            maxSlider->setValue(minSlider->value()+1);
        }
        else
        {
            minSlider->setValue(maxSlider->value()-1);
        }
    }
    emit slidersChanged(minSlider->value(), maxSlider->value());
}

void MainWindow::preventSliderCrossingMin()
{
    if(maxSlider->value() <= minSlider->value())
    {
        if(maxSlider->value() == maxSlider->maximum())
        {
            minSlider->setValue(maxSlider->value()-1);
        }
        else
        {
            maxSlider->setValue(minSlider->value()+1);
        }
    }
     emit slidersChanged(minSlider->value(), maxSlider->value());
}
