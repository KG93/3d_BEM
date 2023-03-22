#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget* parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    widget = new QWidget(this);
    mainLayout = new QVBoxLayout;
    tabWidget =new QTabWidget();
    menuBar2 = new QMenuBar();
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
    solvScriptReader= new SolvingScriptReader(this);
    obsScriptReader= new ObservationScriptReader(this);
    generalLog= new LogWidget;
    connect(solvScriptReader, SIGNAL(logMessage(QString)), generalLog, SLOT(writeOutLogString()));
    connect(obsScriptReader, SIGNAL(logMessage(QString)), generalLog, SLOT(writeOutLogString()));
    connect(boundaryElementSolver, SIGNAL(updateLog()), generalLog, SLOT(writeOutLogString()));
    errorLog = new LogWidget;
    connect(solvScriptReader, SIGNAL(errorMessage(QString)), errorLog, SLOT(writeOutErrorLogString()));
    connect(obsScriptReader, SIGNAL(errorMessage(QString)), errorLog, SLOT(writeOutErrorLogString()));
    connect(boundaryElementSolver, SIGNAL(updateLog()), errorLog, SLOT(writeOutErrorLogString()));

    tabWidget->addTab(registerScriptsTab,"Scripts");
    tabWidget->addTab(generalLog,"Log");
    tabWidget->addTab(errorLog,"Error log");
    mainLayout->setMenuBar(menuBar2);
    mainLayout->addWidget(tabWidget);
    mainLayout->setContentsMargins(QMargins(3,3,3,3));
    widget->setLayout(mainLayout);
    setCentralWidget(widget);
    setUpMenu();

    setWindowTitle(tr("BEM"));
    widget->show();

    parameterDialog = new ParameterDialog(widget);
    freqWidget = new FreqListWidget();
    connect(freqWidget, SIGNAL(signalNewSelection(Eigen::VectorXcd, Eigen::VectorXcd, Eigen::VectorXcd)), openGlWidget, SLOT(getNewSolution(Eigen::VectorXcd, Eigen::VectorXcd, Eigen::VectorXcd)));
    connect(freqWidget, SIGNAL(signalNewSelection(Eigen::VectorXcd, Eigen::VectorXcd, Eigen::VectorXcd, QVector<Eigen::VectorXcd>, QVector<Eigen::VectorXcd>)), openGlWidget, SLOT(getNewSolutionField(Eigen::VectorXcd, Eigen::VectorXcd, Eigen::VectorXcd, QVector<Eigen::VectorXcd>, QVector<Eigen::VectorXcd>)));
//    void signalNewSelection(Eigen::VectorXcd phiSolution, Eigen::VectorXcd dPhiSolution, Eigen::VectorXcd soundPressure, QVector<Eigen::VectorXcd> phiSolutionField, QVector<Eigen::VectorXcd> soundPressureField);
//    void openglWidget::getNewSolutionField(Eigen::VectorXcd phiSolution, Eigen::VectorXcd dPhiSolution, Eigen::VectorXcd soundPressure, QVector<Eigen::VectorXcd> phiSolutionField, QVector<Eigen::VectorXcd> soundPressureField)
//    freqWidget->show();

    if(global::activeProgramming)
    {
//        QString filename = "/home/lrak/workspace/3d_BEM/Example/Loudspeaker/Project.abec"; //loudspeaker
        QString filename = "/home/lrak/workspace/3d_BEM/Example/Sphere/sphere.abec"; //sphere

        projFileHandler->readProjectFile(filename);

        solvinScript = QFileInfo(projFileHandler->getSolvingScript());
        setUpObservationScripts(projFileHandler->getObservScriptList());
        setUpMeshFiles(projFileHandler->getMeshFileList(),projFileHandler->getMeshFileAliasList());

        loadProject();
        BoundaryElements elements = boundaryElementSolver->getBoundaryElements();
//        elements.reflectGeometry();
        openGlWidget->setBoundaryElements(elements);
//        openGlWidget->setBoundaryElements(boundaryElementSolver->getBoundaryElements());
        std::cout<<"Setting up boundary elements in gl widget."<< std::endl;
//        openGlWidget->setupClustertree();
//        std::cout<<"Setting up clusterTree in gl widget."<< std::endl;
//        calculateSolution();
    }
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::setUpMenu()
{
    menuBar2->clear();
    QMenu* menu = menuBar2->addMenu(tr("Project"));
    QAction* newProjectAction = new QAction(tr("&New project"), this);
    QAction* openProjectAction = new QAction(tr("&Open project"), this);

    menu->addAction(newProjectAction);
    connect(newProjectAction, SIGNAL(triggered()), SLOT(newProject()));

    menu->addAction(openProjectAction);
    connect(openProjectAction, SIGNAL(triggered()), SLOT(openProject()));

    if(showSolutionAndFieldMenu)
    {
        QMenu* menu1 = menuBar2->addMenu(tr("Viewer"));
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

        if(openGlWidget->getShowSolutionValue() && hasBeenSolved)
        {
            menu1->addAction(hideSolutionAction);
        }
        else if (hasBeenSolved)
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
        if((openGlWidget->getShowSolutionValue() || openGlWidget->getDrawField()) && openGlWidget->getShowPhaseValue() && hasBeenSolved)
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
        else if((openGlWidget->getShowSolutionValue() || openGlWidget->getDrawField()) && hasBeenSolved)
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

        QMenu* menu2 = menuBar2->addMenu(tr("Solving"));
        QAction* calculateSolutionAction = new QAction(tr("&Calculate solution on boundary"), this);
        QAction* setParametersAction = new QAction(tr("&Set solver parameters"), this);
//        QAction* hSolvingAction = new QAction(tr("&Switch to H-Solver"), this);
//        QAction* setACArelError = new QAction(tr("&Set relative ACA  error"), this);
//        QAction* setPrecondRank = new QAction(tr("&Set local rank of preconditioner"), this);

//        QAction* regSolvingAction = new QAction(tr("&Switch to regular Solver"), this);
//        QAction* defaultCouplinAction = new QAction(tr("&No Coupling"), this);
//        QAction* burtonMillerCouplingAction = new QAction(tr("&Burton and Miller Coupling"), this);
//        QAction* KirkupCouplingAction = new QAction(tr("&Kirkup Coupling"), this);
//        QAction* flipCouplingSign = new QAction(tr("&Set sign of coupling parameter negative"), this);
//        QAction* setRegularQuadratureOrderAction = new QAction(tr("&Set regular quadrature order"), this);
//        QAction* setHigQuadratureOrderAction = new QAction(tr("&Set precision quadrature order"), this);
//        if(!couplingSignPos)
//        {
//            flipCouplingSign= new QAction(tr("&Set sign of coupling parameter positive"), this);
//        }
        menu2->addAction(calculateSolutionAction);
        menu2->addAction(setParametersAction);
//        if(boundaryElementSolver->getHSolving())
//        {
//            menu2->addAction(regSolvingAction);
//            menu2->addAction(setACArelError);
//            menu2->addAction(setPrecondRank);
//        }
//        else
//        {
//            menu2->addAction(hSolvingAction);
//        }
//        menu2->addAction(defaultCouplinAction);
//        menu2->addAction(burtonMillerCouplingAction);
//        menu2->addAction(KirkupCouplingAction);
//        menu2->addAction(flipCouplingSign);
//        menu2->addAction(setRegularQuadratureOrderAction);
//        menu2->addAction(setHigQuadratureOrderAction);
        connect(calculateSolutionAction, SIGNAL(triggered()), SLOT(calculateSolution()));
        connect(setParametersAction, SIGNAL(triggered()), SLOT(showSolverParameterDialog()));
//        connect(hSolvingAction, SIGNAL(triggered()), this,  SLOT(setHSolving()));
//        connect(setACArelError, SIGNAL(triggered()),this, SLOT(setACARelativeError()));
//        connect(setPrecondRank, SIGNAL(triggered()),this, SLOT(setPreconditionerRank()));

//        connect(regSolvingAction, SIGNAL(triggered()),this, SLOT(setRegularSolving()));
//        connect(defaultCouplinAction, SIGNAL(triggered()),boundaryElementSolver, SLOT(setNoCoupling()));
//        connect(burtonMillerCouplingAction, SIGNAL(triggered()),boundaryElementSolver, SLOT(setBurtonMillerCoupling()));
//        connect(KirkupCouplingAction, SIGNAL(triggered()),boundaryElementSolver, SLOT(setKirkupCoupling()));
//        connect(flipCouplingSign, SIGNAL(triggered()),this, SLOT(flipCouplingSign()));
//        connect(setRegularQuadratureOrderAction, SIGNAL(triggered()),this, SLOT(setRegularQuadrature()));
//        connect(setHigQuadratureOrderAction, SIGNAL(triggered()),this, SLOT(setHighOrderQuadrature()));

        QMenu* menu3 = menuBar2->addMenu(tr("Field"));
        if(showCalculateFieldButton)
        {
            QAction* calculateFieldAction = new QAction(tr("&Calculate solution on field"), this);
            menu3->addAction(calculateFieldAction);
            connect(calculateFieldAction, SIGNAL(triggered()), SLOT(calculateSolutionOnField()));
        }       

        QAction* reloadProjectAction = new QAction(tr("&Reload Project"), this);
        connect(reloadProjectAction, SIGNAL(triggered()), SLOT(reloadProjectQuery()));
        QAction* reloadObsScript = new QAction(tr("&Reload ObservationScript"), this);
        connect(reloadObsScript, SIGNAL(triggered()), SLOT(loadObservationScript()));
        menuBar2->addAction(reloadProjectAction);
        menuBar2->addAction(reloadObsScript);
    }
}

void MainWindow:: newProject()
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
      showSolutionAndFieldMenu = true;
      showCalculateFieldButton = false;
      hasBeenSolved = false;
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

//void MainWindow::loadProject()
//{
//    registerScriptsTab->clean();
//    registerScriptsTab->setSolvScript(solvinScript);
//    registerScriptsTab->setObsScripts(observationScripts);
//    registerScriptsTab->setMeshFiles(meshFiles,meshFileAlias);
//    solvScriptReader->setMeshFilesAndAliases(projFileHandler->getMeshFileList(), projFileHandler->getMeshFileAliasList());
//    solvScriptReader->readSolvingScript(solvinScript.absoluteFilePath());
//    solvScriptReader->checkNodesIdentifiersAndSort();
//    solvScriptReader->removeInvalidElements();
//    solvScriptReader->setupBoundaryElements();
//    solvScriptReader->setupPressurePoints();
//    solvScriptReader->refineElements();

//    obsScriptReader->setElementSectionsNames(solvScriptReader->getElementSectionsNames());
//    obsScriptReader->setMeshFilesAndAliases(projFileHandler->getMeshFileList(), projFileHandler->getMeshFileAliasList());
//    obsScriptReader->setEdgeLength(solvScriptReader->getGlobalEdgelength());
//    obsScriptReader->setWavespeed(solvScriptReader->getWavespeed());
//    obsScriptReader->readObservationScript(projFileHandler->getObservScriptList());
//    obsScriptReader->checkNodesIdentifiersAndSort();
//    obsScriptReader->removeInvalidElements();
//    obsScriptReader->setupObservationElements();
//    boundaryElementSolver->setBoundaryElements(solvScriptReader->getBoundaryElements());
//    boundaryElementSolver->setPointSources(solvScriptReader->getPointSources());
//    boundaryElementSolver->setObservationFields(obsScriptReader->getObservationFields());
//    boundaryElementSolver->setObservationPoints(obsScriptReader->getObservationPoints());
////    boundaryElementSolver->setFrequency(solvScriptReader->getFrequencies());
//    boundaryElementSolver->prepareBoundaryElements();
//    openGlWidget->setBoundaryElements(solvScriptReader->getBoundaryElements());
//    openGlWidget->setCenterOfMass(solvScriptReader->getCenterOfMass());
//    openGlWidget->setCameraRadius(solvScriptReader->getContainingRadius());
//    openGlWidget->setElementSections(solvScriptReader->getElementsSections());
//    openGlWidget->setPointSources(solvScriptReader->getPointSources());
//    openGlWidget->setObservationFields(obsScriptReader->getObservationFields());
//    openGlWidget->setObservationPoints(obsScriptReader->getObservationPoints());
//    openGlWidget->setShowSolutionValue(false);

//    showSolutionAndFieldMenu=true;
//    showCalculateFieldButton=false;
//    hasBeenSolved=false;
//    setUpMenu();
//    update();
//}

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
    boundaryElementSolver->setBoundaryElements(solvScriptReader->getBoundaryElements());
    boundaryElementSolver->setPointSources(solvScriptReader->getPointSources());
    boundaryElementSolver->setObservationFields(obsScriptReader->getObservationFields());
    boundaryElementSolver->setObservationPoints(obsScriptReader->getObservationPoints());
//    boundaryElementSolver->setFrequency(solvScriptReader->getFrequencies());
    boundaryElementSolver->prepareBoundaryElements();
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
    hasBeenSolved = false;
    setUpMenu();
    update();
}

void MainWindow::loadObservationScript()
{
    obsScriptReader->setElementSectionsNames(solvScriptReader->getElementSectionsNames());
    obsScriptReader->setMeshFilesAndAliases(projFileHandler->getMeshFileList(), projFileHandler->getMeshFileAliasList());
    obsScriptReader->setEdgeLength(solvScriptReader->getGlobalEdgelength());
    obsScriptReader->setWavespeed(solvScriptReader->getWavespeed());
    obsScriptReader->readObservationScript(projFileHandler->getObservScriptList());
    obsScriptReader->checkNodesIdentifiersAndSort();
    obsScriptReader->removeInvalidElements();
    obsScriptReader->setupObservationElements();

    boundaryElementSolver->setObservationFields(obsScriptReader->getObservationFields());
    boundaryElementSolver->setObservationPoints(obsScriptReader->getObservationPoints());
    openGlWidget->setObservationFields(obsScriptReader->getObservationFields());
    openGlWidget->setObservationPoints(obsScriptReader->getObservationPoints());
}

//void MainWindow::calculateSolution()
//{
//    if(hasBeenSolved)
//    {
//        QMessageBox::StandardButton reply;
//        reply = QMessageBox::question(this, "", "Project has already been solved.\r\n Solve again?", QMessageBox::Yes|QMessageBox::No);
//        if (reply == QMessageBox::No)
//        {
//            return;
//        }
//    }
//    setSolverParameters();
//    boundaryElementSolver->setBoundaryElements(solvScriptReader->getBoundaryElements());
//    boundaryElementSolver->setPointSources(solvScriptReader->getPointSources());
//    boundaryElementSolver->setObservationFields(obsScriptReader->getObservationFields());
//    boundaryElementSolver->setObservationPoints(obsScriptReader->getObservationPoints());
////    boundaryElementSolver->setWavenumber(std::complex<double>(0.37,0.0));
//    boundaryElementSolver->setFrequency(solvScriptReader->getFrequencies().first());
//    boundaryElementSolver->setWavespeed(solvScriptReader->getWavespeed());
//    boundaryElementSolver->setAirDensity(solvScriptReader->getAirDensity());
//    boundaryElementSolver->prepareBoundaryElements();
//    boundaryElementSolver->solve();
//    openGlWidget->setBoundaryElements(boundaryElementSolver->getBoundaryElements());
//    showCalculateFieldButton=true;
//    hasBeenSolved=true;
//    openGlWidget->setShowSolutionValue(true);
//    setUpMenu();
//    openGlWidget->update();
//}

void MainWindow::calculateSolution()
{
    if(hasBeenSolved)
    {
        QMessageBox::StandardButton reply;
        reply = QMessageBox::question(this, "", "Project has already been solved.\r\n Solve again?", QMessageBox::Yes|QMessageBox::No);
        if (reply == QMessageBox::No)
        {
            return;
        }
    }
    setSolverParameters();
    boundaryElementSolver->setBoundaryElements(solvScriptReader->getBoundaryElements());
    boundaryElementSolver->setPointSources(solvScriptReader->getPointSources());

    boundaryElementSolver->setWavespeed(solvScriptReader->getWavespeed());
    boundaryElementSolver->setAirDensity(solvScriptReader->getAirDensity());

    for(int i=0; i<freqWidget->getFrequencies().length(); i++)
    {
        boundaryElementSolver->setObservationFields(obsScriptReader->getObservationFields());
        boundaryElementSolver->setObservationPoints(obsScriptReader->getObservationPoints());
        boundaryElementSolver->setFrequency(freqWidget->getFrequencies().at(i));
//        boundaryElementSolver->setBoundaryElements(boundaryElementsFreq.at(i));
        boundaryElementSolver->prepareBoundaryElements();
        boundaryElementSolver->solve();
        freqWidget->setSolutions(i, boundaryElementSolver->getBoundaryElements().phiSolution, boundaryElementSolver->getBoundaryElements().dPhiSolution, boundaryElementSolver->getBoundaryElements().soundPressure);
    }
    std::pair<double, double> minMaxSoundPressure = freqWidget->getMinAndMaxSoundPressureOnBoundary();
    openGlWidget->setGlobalMinDbOnBoundary(minMaxSoundPressure.first);
    openGlWidget->setGlobalMaxDbOnBoundary(minMaxSoundPressure.second);

    openGlWidget->setBoundaryElements(boundaryElementSolver->getBoundaryElements());
//    openGlWidget->->getBoundaryElements().

    showCalculateFieldButton = true;
    hasBeenSolved = true;
    openGlWidget->setShowSolutionValue(true);
    setUpMenu();
    openGlWidget->update();
    freqWidget->show();
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
    openGlWidget->setShowSolutionValue(true);
    setUpMenu();
    openGlWidget->update();
}

void MainWindow::hideSolutionOnBoundaryElements()
{
    openGlWidget->setShowSolutionValue(false);
    setUpMenu();
    openGlWidget->update();
}

//void MainWindow::calculateSolutionOnField()
//{
//    setSolverParameters();
//    boundaryElementSolver->calculateFieldSolution();
//    openGlWidget->setObservationFields(boundaryElementSolver->getObservationFields());
//    openGlWidget->setDrawField(true);
//    setUpMenu();
//    openGlWidget->update();
//}

void MainWindow::calculateSolutionOnField()
{
    setSolverParameters();
//    boundaryElementSolver->setObservationFields(obsScriptReader->getObservationFields());

    for(int i=0; i<freqWidget->getFrequencies().length(); i++)
    {
        std::tuple<Eigen::VectorXcd, Eigen::VectorXcd, Eigen::VectorXcd> solution = freqWidget->getSolution(i);
        BoundaryElements elements = boundaryElementSolver->getBoundaryElements();
        elements.phiSolution = std::get<0>(solution);
        elements.dPhiSolution = std::get<1>(solution);
        elements.soundPressure = std::get<2>(solution);
        boundaryElementSolver->setBoundaryElements(elements);
        boundaryElementSolver->setFrequency(freqWidget->getFrequencies().at(i));
        boundaryElementSolver->prepareBoundaryElements();
        boundaryElementSolver->calculateFieldSolution();

        QVector<ObservationField> obsFields = boundaryElementSolver->getObservationFields();
        QVector<Eigen::VectorXcd> phiSolutionField(obsFields.size());
        QVector<Eigen::VectorXcd> soundPressureField(obsFields.size());
        for(int j=0; j<obsFields.size(); j++)
        {
            phiSolutionField[j] = obsFields.at(j).phiSolution;
            soundPressureField[j] = obsFields.at(j).soundPressure;
        }
        freqWidget->setSolutionsField(i, phiSolutionField, soundPressureField);
    }
    openGlWidget->setBoundaryElements(boundaryElementSolver->getBoundaryElements());
    openGlWidget->setObservationFields(boundaryElementSolver->getObservationFields());

    std::pair<double, double> minMaxSoundPressure = freqWidget->getMinAndMaxSoundPressureOnField();
    openGlWidget->setGlobalMinDbOnField(minMaxSoundPressure.first);
    openGlWidget->setGlobalMaxDbOnField(minMaxSoundPressure.second);
    openGlWidget->setDrawField(true);
    setUpMenu();
    openGlWidget->update();
    freqWidget->show();
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

void MainWindow::setHSolving()
{
    boundaryElementSolver->setHSolving(true);
    setUpMenu();
}

void MainWindow::setACARelativeError()
{
    bool ok;
    double acaRelError = boundaryElementSolver->getACArelativeError();
    acaRelError = QInputDialog::getDouble(this, tr("Set the relative ACA approximation error"),tr("relative Error:"), acaRelError, 0, 1, 6, &ok, Qt::WindowFlags(), 0.01);
    if(ok)
    {
        boundaryElementSolver->setACArelativeError(acaRelError);
    }
}

void MainWindow::setPreconditionerRank()
{
    bool ok;
    unsigned long rank = boundaryElementSolver->getPreconditionerRank();
    rank = QInputDialog::getInt(this, tr("Set the rank of the preconditioner. Rank zero is a direct solver."),tr("Order:"), rank, 0, std::numeric_limits<int>::max(), 1, &ok);
    if(ok)
    {
        boundaryElementSolver->setPreconditionerRank(rank);
    }
}

void MainWindow::setRegularSolving()
{
    boundaryElementSolver->setHSolving(false);
    setUpMenu();
}

void MainWindow::setHFieldSolving()
{
    boundaryElementSolver->setHFieldSolving(true);
    setUpMenu();
}

void MainWindow::setFieldACARelativeError()
{
    bool ok;
    double acaRelError = boundaryElementSolver->getFieldACARelError();
    acaRelError = QInputDialog::getDouble(this, tr("Set the relative ACA approximation error for the field calculatation"),tr("relative Error:"), acaRelError, 0, 1, 6, &ok, Qt::WindowFlags(), 0.01);
    if(ok)
    {
        boundaryElementSolver->setFieldACARelError(acaRelError);
    }
}

void MainWindow::setFieldACAMaxRank()
{
    bool ok;
    unsigned long rank = boundaryElementSolver->getFieldACAMaxRank();
    rank = QInputDialog::getInt(this, tr("Set the maximum local rank of the field ACA. Rank zero is umlimited."),tr("Order:"), rank, 0, std::numeric_limits<int>::max(), 1, &ok);
    if(ok)
    {
        boundaryElementSolver->setFieldACAMaxRank(rank);
    }
}

void MainWindow::setRegularFieldSolving()
{
    boundaryElementSolver->setHFieldSolving(false);
    setUpMenu();
}

void MainWindow::setRegularQuadrature()
{
    bool ok;
    regQuadOrder = QInputDialog::getInt(this, tr("Set the order of the regular quadrature method"),tr("Order:"), regQuadOrder, 1, 14, 1, &ok);
    if(ok)
    {
        boundaryElementSolver->setRegularOrderQuadratureRule(regQuadOrder);
    }
}

void MainWindow::setHighOrderQuadrature()
{
    bool ok;
    highQuadOrder = QInputDialog::getInt(this, tr("Set the order of the high order quadrature method"),tr("Order:"), highQuadOrder, 1, 14, 1, &ok);
    if(ok)
    {
        boundaryElementSolver->setHighOrderQuadratureRule(highQuadOrder);
    }
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
    solvinScript=registerScriptsTab->getSolvingScript();
    observationScripts=registerScriptsTab->getObservationScripts();
    meshFiles=registerScriptsTab->getMeshFiles();
    meshFileAlias=registerScriptsTab->getMeshFileAliases();
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

void MainWindow::reloadProjectQuery()
{
    if(hasBeenSolved)
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

void MainWindow::flipCouplingSign()
{
    if(couplingSignPos)
    {
        couplingSignPos = false;
    }
    else
    {
        couplingSignPos = true;
    }
    setUpMenu();
    boundaryElementSolver->flipCouplingSign();
}

void MainWindow::setSolverParameters()
{
    if(parameterDialog->getCoupling())
    {
        if(parameterDialog->getBurtonMillerCoupling())
        {
            boundaryElementSolver->setBurtonMillerCoupling();
        }
        else
        {
            boundaryElementSolver->setKirkupCoupling();
        }
    }
    else
    {
        boundaryElementSolver->setNoCoupling();
    }
    boundaryElementSolver->setHSolving(parameterDialog->getUseHSolver());
    boundaryElementSolver->setACArelativeError(parameterDialog->getACARelError());
    boundaryElementSolver->setACAMaxRank(parameterDialog->getACAMaxRank());
    boundaryElementSolver->setUsePreconditioner(parameterDialog->getUseHPreconditioning());

    boundaryElementSolver->setPreconditionerRelativeError(parameterDialog->getPreconditionerRelError());
    boundaryElementSolver->setPreconditionerRank(parameterDialog->getPreconditionerMaxRank());

    boundaryElementSolver->setCalculateNormAndConditionNumber(parameterDialog->getCalculeteNormAndConditionNumber());

    boundaryElementSolver->setHFieldSolving(parameterDialog->getUseHFieldSolver());
    boundaryElementSolver->setFieldACARelError(parameterDialog->getHFieldSolverRelError());
    boundaryElementSolver->setFieldACAMaxRank(parameterDialog->getHFieldSolverMaxRank());
}
