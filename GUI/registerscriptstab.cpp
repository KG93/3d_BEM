#include "registerscriptstab.h"

RegisterScriptsTab::RegisterScriptsTab(QWidget *parent) : QWidget(parent)
{
//    this->setUpdatesEnabled(true);
    QFileInfo argFileInfo;
    this->solvingScript=argFileInfo;
    QLabel *bemScript = new QLabel(tr("Boundary Element Script:"));
//    fileNameEdit = new QLineEdit(fileInfo.fileName());
    solvScriptLabel = new QLabel(solvingScript.absoluteFilePath());
    solvScriptLabel->setStyleSheet("QLabel { background-color : white; color : black; }");
    solvScriptLabel->setFrameStyle(QFrame::Panel | QFrame::Sunken);
    solvScriptLabel->setContextMenuPolicy(Qt::CustomContextMenu);
    solvScriptLabel->setWordWrap(true);
    connect(solvScriptLabel, SIGNAL(customContextMenuRequested(const QPoint&)), this, SLOT(showSolvContextMenu(const QPoint&)));

    QLabel *obsScript = new QLabel(tr("Observation Scripts:"));
    pathObsValueLabel = new QLabel(solvingScript.absoluteFilePath());
    pathObsValueLabel->setFrameStyle(QFrame::Panel | QFrame::Sunken);
    pathObsValueLabel->setWordWrap(true);

    observationListWidget = new QListWidget();
    observationListWidget->setContextMenuPolicy(Qt::CustomContextMenu);
    connect(observationListWidget, SIGNAL(customContextMenuRequested(const QPoint&)), this, SLOT(showObsContextMenu(const QPoint&)));

    QLabel *meshFiles = new QLabel(tr("MeshFiles:"));
    meshFilesLabel= new QLabel(solvingScript.absoluteFilePath());
    meshFilesLabel->setFrameStyle(QFrame::Panel | QFrame::Sunken);
    meshFilesLabel->setWordWrap(true);

    meshListWidget = new QListWidget();
    meshListWidget->setContextMenuPolicy(Qt::CustomContextMenu);
    connect(meshListWidget, SIGNAL(customContextMenuRequested(const QPoint&)), this, SLOT(showMeshContextMenu(const QPoint&)));

//    connect(listWidget, SIGNAL(customContextMenuRequested(QPoint)), this, SLOT(showContextMenu(QPoint)));
//    QLabel *sizeLabel = new QLabel(tr("Size:"));
//    qlonglong size = fileInfo.size()/1024;
//    QLabel *sizeValueLabel = new QLabel(tr("%1 K").arg(size));
//    sizeValueLabel->setFrameStyle(QFrame::Panel | QFrame::Sunken);

//    QLabel *lastReadLabel = new QLabel(tr("Last Read:"));
//    QLabel *lastReadValueLabel = new QLabel(fileInfo.lastRead().toString());
//    lastReadValueLabel->setFrameStyle(QFrame::Panel | QFrame::Sunken);

//    QLabel *lastModLabel = new QLabel(tr("Last Modified:"));
//    QLabel *lastModValueLabel = new QLabel(fileInfo.lastModified().toString());
//    lastModValueLabel->setFrameStyle(QFrame::Panel | QFrame::Sunken);

    QVBoxLayout *mainLayout = new QVBoxLayout;
    mainLayout->addWidget(bemScript);
//    mainLayout->addWidget(fileNameEdit);
//    mainLayout->addWidget(pathLabel);
    mainLayout->addWidget(solvScriptLabel);
    mainLayout->addWidget(obsScript);
    mainLayout->addWidget(pathObsValueLabel);
    mainLayout->addWidget(observationListWidget);
    mainLayout->addWidget(meshFiles);
    mainLayout->addWidget(meshFilesLabel);
    mainLayout->addWidget(meshListWidget);
//    mainLayout->addWidget(sizeLabel);
//    mainLayout->addWidget(sizeValueLabel);
//    mainLayout->addWidget(lastReadLabel);
//    mainLayout->addWidget(lastReadValueLabel);
//    mainLayout->addWidget(lastModLabel);
//    mainLayout->addWidget(lastModValueLabel);
//    mainLayout->addStretch(1);
    setLayout(mainLayout);
}

void RegisterScriptsTab::setSolvScript(const QFileInfo& fileInfo)
{
    this->solvingScript=fileInfo;
   updateSolvLabel();
}

void RegisterScriptsTab::setObsScripts(const QVector<QFileInfo>& obsFileInfos)
{
    this->observationScripts=obsFileInfos;
    updateObservationScripts();
}

void RegisterScriptsTab::setMeshFiles(const QVector<QFileInfo>& meshFiles, const QStringList& meshfileAliases)
{
    this->meshFiles=meshFiles;
    this->meshFileAliases=meshfileAliases;
    updateMeshFiles();
}

QStringList  RegisterScriptsTab::getObservationScriptPaths()
{
    QStringList obsScripts;
    for (int i=0;i<observationScripts.size();i++)
    {
        obsScripts.append(observationScripts.at(i).absoluteFilePath());
    }
    return obsScripts;
}

QStringList  RegisterScriptsTab::getMeshFilePaths()
{
    QStringList meshFilesPath;
    for (int i=0;i<meshFiles.size();i++)
    {
        meshFilesPath.append(meshFiles.at(i).absoluteFilePath());
    }
    return meshFilesPath;
}

void RegisterScriptsTab::clean()
{
    meshFiles.clear();
    meshFileAliases.clear();
    observationScripts.clear();
    observationListWidget->clear();
    meshListWidget->clear();
}

void RegisterScriptsTab::updateSolvLabel()
{
    solvScriptLabel->setText(solvingScript.absoluteFilePath());
    update();
}


void RegisterScriptsTab::updateMeshFiles()
{
    meshListWidget->clear();
    QString meshFileString;
    for (int i=0;i<meshFiles.size();i++)
    {
        meshFileString.append(meshFiles.at(i).absoluteFilePath());
        meshFileString.append("\r\n");
        QListWidgetItem* item =new QListWidgetItem(meshFiles.at(i).fileName()+ "," + meshFileAliases.at(i) );
        meshListWidget->addItem(item);
    }
    meshFilesLabel->setText(meshFileString);
    update();
}

void RegisterScriptsTab::updateObservationScripts()
{
    observationListWidget->clear();
    QString observationFiles;
    for (int i=0;i<observationScripts.size();i++)
    {
        observationFiles.append(observationScripts.at(i).absoluteFilePath());
        observationFiles.append("\r\n");
        QListWidgetItem* item =new QListWidgetItem(observationScripts.at(i).fileName());
        observationListWidget->addItem(item);
    }
    pathObsValueLabel->setText(observationFiles);
    update();
}

void RegisterScriptsTab::mousePressEvent(QMouseEvent *event)
{
    if(event->button() == Qt::RightButton)
    {
        emit customContextMenuRequested(event->pos());
    }

    if (event->buttons() != Qt::LeftButton)
    {
        return;
    }
}

void RegisterScriptsTab::showSolvContextMenu(const QPoint& pos)
{
    QPoint globalPos = solvScriptLabel->mapToGlobal(pos);

    QMenu* contextMenu=new QMenu(this);
    if(solvingScript.fileName().isEmpty())
    {
        contextMenu->addAction("Add SolvingScript", this, SLOT(addSolvScript()));
    }
    else
    {
        contextMenu->addAction("Replace SolvingScript",  this, SLOT(replaceSolvScript()));
        contextMenu->addAction("Open SolvingScript in Editor",  this, SLOT(solvScriptInEditor()));
    }
//    contextMenu->addAction("Erase SolvingScript",  this, SLOT(eraseMeshItem()));
    contextMenu->exec(globalPos);
}

void RegisterScriptsTab::showMeshContextMenu(const QPoint& pos)
{
    QPoint globalPos = meshListWidget->mapToGlobal(pos);

    QMenu* contextMenu=new QMenu(this);
    contextMenu->addAction("Add Meshfile", this, SLOT(addMeshItem()));
    if(meshListWidget->count()>=1)
    {
        contextMenu->addAction("Open MeshFile in Editor",  this, SLOT(meshFileInEditor()));
        contextMenu->addAction("Edit MeshFileAlias",  this, SLOT(editMeshItemAlias()));
        contextMenu->addAction("Remove MeshFile",  this, SLOT(eraseMeshItem()));
    }
    contextMenu->exec(globalPos);
}

void RegisterScriptsTab::showObsContextMenu(const QPoint& pos)
{
    QPoint globalPos = observationListWidget->mapToGlobal(pos);

    QMenu* contextMenu=new QMenu(this);
    contextMenu->addAction("Add ObservationScript", this, SLOT(addObservationItem()));
    if(observationListWidget->count()>=1)
    {
        contextMenu->addAction("Open ObservationScript in Editor",  this, SLOT(obsScriptInEditor()));
        contextMenu->addAction("Replace ObservationScript",  this, SLOT(replaceObservationItem()));
        contextMenu->addAction("Remove ObservationScript",  this, SLOT(eraseObservationItem()));
    }
    contextMenu->exec(globalPos);
}

void RegisterScriptsTab::addSolvScript()
{
    QString filename = QFileDialog::getOpenFileName(this, "Select observation script.","",tr("Text files (*.txt)"));
        // filter, show only txt files
    if ( filename.isEmpty() )
    {
        std::cout <<"Empty filename: "<<std::endl;
        return;
    }
    solvingScript=QFileInfo(filename);
    updateSolvLabel();
    emit scriptChanged();
}

void RegisterScriptsTab::replaceSolvScript()
{
    QString filename = QFileDialog::getOpenFileName(this, "Select observation script.","",tr("Text files (*.txt)"));
        // filter, show only txt files
    if ( filename.isEmpty() )
    {
        std::cout <<"Empty filename: "<<std::endl;
        return;
    }
    solvingScript=QFileInfo(filename);
    updateSolvLabel();
    emit scriptChanged();
}

void RegisterScriptsTab::addMeshItem()
{
    QString filename = QFileDialog::getOpenFileName(this, "Select mesh file.","",tr("Mesh files (*.msh)"));
        // filter, show only msh files
    if ( filename.isEmpty() )
    {
        std::cout <<"Empty filename: "<<std::endl;
        return;
    }
    bool ok;
    int meshAliasNumber=meshFiles.size()+1;
    QString defAlias=QString("Alias").append(QString::number(meshAliasNumber));
    QString alias =QInputDialog::getText(this, tr("Specify MeshFileAlias"), tr("MeshFileAlias: "), QLineEdit::Normal,defAlias, &ok);
    if (!ok | alias.isEmpty())
    {
        std::cout <<"Empty MeshFileAlias"<<std::endl;
        return;
    }
    if(meshFileAliases.contains(alias,Qt::CaseInsensitive))
    {
        std::cout <<"Non unique MeshFileAlias: "<<alias.toUtf8().constData()<<std::endl;
        alias =QInputDialog::getText(this, tr("Specify unique MeshFileAlias"), tr("MeshFileAlias: "), QLineEdit::Normal,defAlias, &ok);
            if (!ok | alias.isEmpty())
            {
                std::cout <<"Empty MeshFileAlias"<<std::endl;
                return;
            }
            if(meshFileAliases.contains(alias,Qt::CaseInsensitive))
            {
                std::cout <<"Non unique MeshFileAlias: "<<alias.toUtf8().constData()<<std::endl;
                return;
            }
    }
    meshFiles.append(QFileInfo(filename));
    meshFileAliases.append(alias);
    updateMeshFiles();
    emit scriptChanged();
}

void RegisterScriptsTab::editMeshItemAlias()
{
    if(meshListWidget->count()==0)
    {
        return;
    }
    QListWidgetItem* currentlistWidgetItem=meshListWidget->currentItem();
    int listIndex=meshListWidget->currentRow();
    QString filenameAndAlias = currentlistWidgetItem->text();
    int index=filenameAndAlias.indexOf(",");
    if(index<=1)
    {
        return;
    }
//    QString meshPath=filenameAndAlias.left(index);
    QString meshAlias=filenameAndAlias.mid(index+1);
    bool ok;
    QString newAlias =QInputDialog::getText(this, tr("Specify MeshFileAlias"), tr("MeshFileAlias: "), QLineEdit::Normal,meshAlias, &ok);
    if (!ok | newAlias.isEmpty())
    {
        std::cout <<"Empty MeshFileAlias"<<std::endl;
        return;
    }
    QStringList meshFileAliasesWithoutCurrent=meshFileAliases;
    if(listIndex>=0 && listIndex<meshFileAliases.size())
    {
        meshFileAliasesWithoutCurrent.removeAt(listIndex);
    }
    if(meshFileAliasesWithoutCurrent.contains(newAlias,Qt::CaseInsensitive))
    {
        std::cout <<"Non unique MeshFileAlias: "<<newAlias.toUtf8().constData()<<std::endl;
        newAlias =QInputDialog::getText(this, tr("Specify unique MeshFileAlias"), tr("MeshFileAlias: "), QLineEdit::Normal,meshAlias, &ok);
            if (!ok | newAlias.isEmpty())
            {
                std::cout <<"Empty MeshFileAlias"<<std::endl;
                return;
            }
    }
    else
    {
        if(listIndex>=0 && listIndex<meshFileAliases.size())
        {
            meshFileAliases[listIndex]=newAlias;
            updateMeshFiles();
            emit scriptChanged();
        }
    }
}

void RegisterScriptsTab::eraseMeshItem()
{
    if(meshListWidget->count()==0)
    {
        return;
    }
    int listIndex=meshListWidget->currentRow();
    if(listIndex>=0 && listIndex<meshFileAliases.size() && listIndex<meshFiles.size())
    {
        meshFiles.removeAt(listIndex);
        meshFileAliases.removeAt(listIndex);
        updateMeshFiles();
        emit scriptChanged();
    }
}

void RegisterScriptsTab::addObservationItem()
{
    QString filename = QFileDialog::getOpenFileName(this, "Select observation script.","",tr("Text files (*.txt)"));
        // filter, show only txt files
    if ( filename.isEmpty() )
    {
        std::cout <<"Empty filename: "<<std::endl;
        return;
    }
    if(observationScripts.contains(QFileInfo(filename)))
    {
        std::cout <<"ObservationScript is already registered."<<std::endl;
        return;
    }
    else
    {
        observationScripts.append(QFileInfo(filename));
        updateObservationScripts();
        emit obsScriptChanged();
    }
}

void RegisterScriptsTab::replaceObservationItem()
{
    if(observationListWidget->count()==0)
    {
        return;
    }
    QString filename = QFileDialog::getOpenFileName(this, "Select observation script.","",tr("Text files (*.txt)"));
        // filter, show only txt files
    if ( filename.isEmpty() )
    {
        std::cout <<"Empty filename: "<<std::endl;
        return;
    }
    if(observationScripts.contains(QFileInfo(filename)))
    {
        std::cout <<"ObservationScript is already registered."<<std::endl;
        return;
    }
    else
    {
        int listIndex=observationListWidget->currentRow();
        if(listIndex>=0 && listIndex<observationScripts.size())
        {
            observationScripts[listIndex]=QFileInfo(filename);
            updateObservationScripts();
            emit obsScriptChanged();
        }
    }
}

void RegisterScriptsTab::eraseObservationItem()
{
    if(observationListWidget->count()==0)
    {
        return;
    }
    int listIndex=observationListWidget->currentRow();
    if(listIndex>=0 && listIndex<observationScripts.size())
    {
        observationScripts.removeAt(listIndex);
        updateObservationScripts();
        emit obsScriptChanged();
    }
}

void RegisterScriptsTab::solvScriptInEditor()
{
    if(solvingScript.absoluteFilePath().isEmpty()==false)
    {
        QString solvScriptPath = QDir::toNativeSeparators(solvingScript.absoluteFilePath());
        QDesktopServices::openUrl(QUrl::fromLocalFile(solvScriptPath));
    }
}

void RegisterScriptsTab::obsScriptInEditor()
{
    int listIndex=observationListWidget->currentRow();
     if(listIndex>=0 && listIndex<observationScripts.size())
     {
         QString url=QDir::toNativeSeparators(observationScripts.at(listIndex).absoluteFilePath());
         if(url.isEmpty()==false)
         {
             QDesktopServices::openUrl(QUrl::fromLocalFile(url));
         }
     }
}

void RegisterScriptsTab::meshFileInEditor()
{
    int listIndex=meshListWidget->currentRow();
    if(listIndex>=0 && listIndex<meshFiles.size())
    {
        QString url=QDir::toNativeSeparators(meshFiles.at(listIndex).absoluteFilePath());
        if(url.isEmpty()==false)
        {
            QDesktopServices::openUrl(QUrl::fromLocalFile(url));
        }
    }
}
