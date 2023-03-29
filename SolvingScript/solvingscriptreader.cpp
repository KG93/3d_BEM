#include "solvingscriptreader.h"

SolvingScriptReader::SolvingScriptReader(QObject* parent): QObject(parent)
{

}

bool SolvingScriptReader::setMeshFilesAndAliases( const QStringList &meshFiles, const QStringList &meshFileAlias)
{
    this->meshFiles=meshFiles;
    this->meshFileAlias=meshFileAlias;
    meshfilesAndAliasesSet=true;
    if(meshFiles.length() != meshFileAlias.length())
    {
        meshfilesAndAliasesSet=false;
        emit errorMessage("The list of the meshfile links and aliases have different lengths!");
        logStrings::errorLogString.append("The list of the meshfile links and aliases have different lengths!\r\n");
        std::cout <<"The list of the meshfile links and aliases have different lengths!"<<std::endl;
    }
    clear();
    return meshfilesAndAliasesSet;
}

void SolvingScriptReader::readSolvingScript(const QString& fileName)
{
    if(!meshfilesAndAliasesSet){
        emit errorMessage("Meshfile links were not passed to solvingscript reader!");
        logStrings::errorLogString.append("Meshfile links were not passed to solvingscript reader!\r\n");
        std::cout <<"Meshfile links were not passed to solvingscript reader!"<<std::endl;
        return;
    }

    controlSolverSectionDeclared=false;
    clear();

    QFile projFile(fileName);
    if (!projFile.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        emit errorMessage("Couldn't read solving script : " + fileName);
        logStrings::errorLogString.append("Couldn't read solving script : " + fileName + "\r\n");
        std::cout <<"Couldn't read solving script : "<< fileName.toUtf8().constData()<<std::endl;
        setSolverParameters();
        return;
    }
    emit logMessage("\r\n** Reading Solving_Script **");
    logStrings::logString.append("\r\n** Reading Solving_Script **\r\n");
    logStrings::errorLogString.append("\r\n** Reading Solving_Script **\r\n");
    emit errorMessage("\r\nSolving_Script:\r\n");

    QTextStream in(&projFile);
    QString scriptFileString=in.readAll();

    QRegularExpression separator("(\r\n|\r|\n|;)"); //string to list of the substrings by ; and newline separator
    QStringList scriptFileStringList=scriptFileString.split(separator, Qt::SkipEmptyParts);
    QStringListIterator iterator(scriptFileStringList);
    while (iterator.hasNext())
    {
        QString currentLine=iterator.next();
        currentLine=currentLine.simplified();
        currentLine.remove(QRegularExpression("\\/\\/.*$")); //remove comments
        if(currentLine.isEmpty()){continue;}
//                cout << iterator.next().toLocal8Bit().constData() << endl;
        if(!newSectionDeclared(currentLine))
        {
             handleLineAccordingToCurrentSection(currentLine);
             std::cout <<"line: "<< currentLine.toUtf8().constData()<<std::endl;
        }
    }
    emit logMessage("\r\n");
    logStrings::logString.append("\r\n");
    logStrings::errorLogString.append("\r\n");
    setSolverParameters();
}

void SolvingScriptReader::setSolverParameters()
{
    controlSolverSection.allMandatoryParametersSet();
    frequencies=controlSolverSection.getFrequencies();
    globalEdgelength=controlSolverSection.getEdgelength();
    c=controlSolverSection.getWavespeed();
    airDensity=controlSolverSection.getDensity();
    mediumImpedance=frequencies.last()*airDensity;
    std::cout<<"Frequencies: ";
    global::printVector<QVector<double> >(frequencies);
    std::cout<<"c: "<<c<<std::endl;
    std::cout<<"Globlal max edgelength: "<<globalEdgelength<<std::endl;
    logStrings::logString.append("Globlal max edgelength: "+QString::number(globalEdgelength)+"\r\n");
}

bool SolvingScriptReader::newSectionDeclared(const QString& currentLine){
    QRegularExpressionMatch match;
    QRegularExpression identifier = sectionIdentifierRegEx("Control_Solver");

    if(currentLine.contains(identifier,&match))
    {
//        QString name=match.captured(1);
//        global::removeQuotesAndWhitespace(value);
        if(controlSolverSectionDeclared)
        {
            std::cout <<"Control_Solver section already declared!"<<std::endl;
            emit errorMessage("\r\nControl_Solver second time declared!");
            logStrings::errorLogString.append("\r\nControl_Solver second time declared! \r\n");
            currentSection = Control_Solver;
        }
        else{
            controlSolverSectionDeclared=true;
            std::cout <<"Control_Solver section read!"<<std::endl;
            emit logMessage("\r\nControl_Solver section read!");
            logStrings::logString.append("\r\nControl_Solver section read! \r\n");
            logStrings::errorLogString.append("\r\nControl_Solver section read! \r\n");
            currentSection = Control_Solver;
        }
        return true;
    }

    identifier=sectionIdentifierRegEx("Elements");
    if(currentLine.contains(identifier,&match))
    {

        QString name= match.captured(1);
        global::removeQuotes(name);
        elementsSections.push_back((ElementSection(name)));
        std::cout <<"Elements section read!"<<std::endl;
        emit logMessage("\r\nElements section " +name + " read!");
        logStrings::logString.append("\r\nElements section " +name + " read!\r\n");
        logStrings::errorLogString.append("\r\nElements section " +name + " read!\r\n");
        currentSection = Elements;
        return true;
    }

    identifier=sectionIdentifierRegEx("MeshFile_Properties");
    if(currentLine.contains(identifier,&match))
    {
        std::cout <<"MeshFile_Properties section read!"<<std::endl;
        QString name= match.captured(1);
        emit logMessage("\r\nMeshFile_Properties section " +name + " read!");
        logStrings::logString.append("\r\nMeshFile_Properties section " +name + " read!\r\n");
        logStrings::errorLogString.append("\r\nMeshFile_Properties section " +name + " read!\r\n");
        meshFilePropertiesSections.push_back(MeshFilePropertiesSection(name));
        currentSection = MeshFile_Properties;
        return true;
    }

    identifier=sectionIdentifierRegEx("Nodes");
    if(currentLine.contains(identifier,&match))
    {
        QString name= match.captured(1);
        global::removeQuotes(name);
        nodesSections.push_back(NodesSection(name));
        std::cout <<"Nodes section read! Name: "<<name.toUtf8().constData()<<std::endl;
        emit logMessage("\r\nNodes section " +name+ " read!");
        logStrings::logString.append("\r\nNodes section " +name+ " read!\r\n");
        logStrings::errorLogString.append("\r\nNodes section " +name+ " read!\r\n");
        currentSection = Nodes;
        return true;
    }

    identifier=sectionIdentifierRegEx("Subdomain_Properties");
    if(currentLine.contains(identifier,&match))
    {
        QString name= match.captured(1);
        subdomainPropertiesSections.push_back(SubdomainPropertiesSection(name));
        std::cout <<"Subdomain_Properties section read! Name: "<<name.toUtf8().constData()<<std::endl;
        emit logMessage("\r\nSubdomain_Properties section " +name+ " read!");
        logStrings::logString.append("\r\nSubdomain_Properties section " +name+ " read!\r\n");
        logStrings::errorLogString.append("\r\nSubdomain_Properties section " +name+ " read!\r\n");
        currentSection = Subdomain_Properties;
        return true;
    }

    identifier=sectionIdentifierRegEx("WallImpedance");
    if(currentLine.contains(identifier,&match))
    {
        QString name= match.captured(1);
        wallImpedanceSections.push_back(WallImpedanceSection(name));
        std::cout <<"WallImpedance section read! Name: "<<name.toUtf8().constData()<<std::endl;
        emit logMessage("\r\nWallImpedance section " +name+ " read!");
        logStrings::logString.append("\r\nWallImpedance section " +name+ " read!\r\n");
        logStrings::errorLogString.append("\r\nWallImpedance section " +name+ " read!\r\n");
        currentSection = WallImpedance;
        return true;
    }

    identifier=sectionIdentifierRegEx("Driving");
    if(currentLine.contains(identifier,&match))
    {
        QString name= match.captured(1);
        global::removeQuotes(name);
        drivingSections.push_back(DrivingSection(name));
        std::cout <<"Driving section read!"<<std::endl;
        emit logMessage("\r\nDriving section " +name+ " read!");
        logStrings::logString.append("\r\nDriving section " +name+ " read!\r\n");
        logStrings::errorLogString.append("\r\nDriving section " +name+ " read!\r\n");
        currentSection = Driving;
        return true;
    }

    identifier=sectionIdentifierRegEx("Pressure_Points");
    if(currentLine.contains(identifier,&match))
    {
        QString name= match.captured(1);
        global::removeQuotes(name);
        pressurePointsSections.push_back(PressurePointsSection(name));
        std::cout <<"Pressure_Points section read! Name: "<<name.toUtf8().constData()<<std::endl;
        emit logMessage("\r\nPressure_Points section " +name+ " read!");
        logStrings::logString.append("\r\nPressure_Points section " +name+ " read!\r\n");
        logStrings::errorLogString.append("\r\nPressure_Points section " +name+ " read!\r\n");
        currentSection = Pressure_Points;
        return true;
    }

    identifier=sectionIdentifierRegEx("Infinite_Baffle");
    if(currentLine.contains(identifier,&match))
    {
        QString name= match.captured(1);
        global::removeQuotes(name);
        infiniteBaffleSections.push_back(InfiniteBaffleSection(name));
        std::cout <<"Infinite_Baffle section read! Name: "<<name.toUtf8().constData()<<std::endl;
        emit logMessage("\r\nInfinite_Baffle section " +name+ " read!");
        logStrings::logString.append("\r\nInfinite_Baffle section " +name+ " read!\r\n");
        logStrings::errorLogString.append("\r\nInfinite_Baffle section " +name+ " read!\r\n");
        currentSection = Infinite_Baffle;
        return true;
    }

    identifier=sectionIdentifierRegEx("off");
    if(currentLine.contains(identifier,&match))
    {
        std::cout <<"off!"<<std::endl;
        currentSection = none;
        return true;
    }
    return false;
}

QRegularExpression SolvingScriptReader::sectionIdentifierRegEx(const QString& identifier)
{
//    return QRegularExpression("^\\s*"+identifier+"\\s*={1}\\s*([\"\'][\\w\\s]+[\"\']|\\w+)", QRegularExpression::CaseInsensitiveOption);
    return QRegularExpression("^\\s*"+identifier+"\\s*([\"\'][^\"\']+[\"\']|[\\S]+|$)", QRegularExpression::CaseInsensitiveOption);
}

void SolvingScriptReader::handleLineAccordingToCurrentSection(const QString& line){
    switch(currentSection)
    {
    case none:
        break;

    case Control_Solver:
        controlSolverSection.handleScriptLine(line);
        break;

    case Driving:
        drivingSections.last().handleScriptLine(line);
        break;

    case Subdomain_Properties:
        subdomainPropertiesSections.last().handleScriptLine(line);
        break;

    case MeshFile_Properties:
        meshFilePropertiesSections.last().handleScriptLine(line);
        break;

    case Elements:
        elementsSections.last().handleScriptLine(line);
        break;

    case Nodes:
        nodesSections.last().handleScriptLine(line);
        break;

    case WallImpedance:
        wallImpedanceSections.last().handleScriptLine(line);
        break;

    case Pressure_Points:
        pressurePointsSections.last().handleScriptLine(line);
        break;

    case Infinite_Baffle:
        infiniteBaffleSections.last().handleScriptLine(line);
        break;
    }
}


void SolvingScriptReader::clear()
{
    controlSolverSection.clear();
    elementsSections.clear();
    nodesSections.clear();
    subdomainPropertiesSections.clear();
    meshFilePropertiesSections.clear();
    wallImpedanceSections.clear();
    pressurePointsSections.clear();
    infiniteBaffleSections.clear();
    drivingSections.clear();

    pointSources.clear();
    boundaryElements.clear();

    currentSection=none;
    controlSolverSectionDeclared=false;
}

bool SolvingScriptReader::removeInvalidElements()
{
    bool noLinkProblem=true;
    for(int i=0;i<elementsSections.length();i++)
    {
        QListIterator<MeshFileElement> meshFileElementIterator(elementsSections.at(i).meshFileElements);
        quint64 iteratorPos=0;
        while (meshFileElementIterator.hasNext())
        {
            MeshFileElement currentMeshFileElement=meshFileElementIterator.next();
            if(!meshFileAlias.contains(currentMeshFileElement.MeshFileAlias,Qt::CaseInsensitive))
            {
                 noLinkProblem=false;
                 elementsSections[i].meshFileElements.removeAt(iteratorPos); //remove Element with invalid meshfile alias
                 std::cout <<"MeshFileAlias: "<<currentMeshFileElement.MeshFileAlias.toUtf8().constData()<<" from "<< elementsSections.at(i).name.toUtf8().constData()
                          <<" Element "<< currentMeshFileElement.elementIndex<<" not declared."<<std::endl;
                 QString error=QString("MeshFileAlias: "+currentMeshFileElement.MeshFileAlias+" from "+elementsSections.at(i).name+" Element "+ QString::number(currentMeshFileElement.elementIndex)+" not declared.");
                 emit errorMessage(error);
                 logStrings::errorLogString.append(error+"\r\n");
                 continue;
            }
            iteratorPos=iteratorPos+1;
        }

//        std::cout <<"Elemens name:"<<elementsSections.at(i).name.toUtf8().constData()<<std::endl;
//        std::cout <<"MeshFileAlias:"<<elementsSections.at(i).MeshFileAlias.toUtf8().constData()<<std::endl;
        if(!elementsSections.at(i).MeshFileAlias.isEmpty())
        {
            if(!(meshFileAlias.contains(elementsSections.at(i).MeshFileAlias,Qt::CaseInsensitive)))
            {
                 //noLinkProblem=false;
                 std::cout <<"MeshFileAlias: "<<elementsSections.at(i).MeshFileAlias.toUtf8().constData()<<" from "<< elementsSections.at(i).name.toUtf8().constData()<<" not declared."<<std::endl;
                 QString error=QString("MeshFileAlias: "+elementsSections.at(i).MeshFileAlias+" from "+ elementsSections.at(i).name+" not declared.");
                 emit errorMessage(error);
                 logStrings::errorLogString.append(error+"\r\n");
            }
        }
    }

    QStringList nodesSectionsNames=getNodesSectionsNames();
    for(int i=0;i<elementsSections.length();i++)
    {

        ElementSection currentElementSection=elementsSections[i];
        QListIterator<Triangle> triangleElementIterator(currentElementSection.triangles);
        quint64 triangleIteratorPos=0;
        while (triangleElementIterator.hasNext())
        {
            Triangle currentTriangleElement=triangleElementIterator.next();
            if(!nodesSectionsNames.contains(currentTriangleElement.nodesReference,Qt::CaseInsensitive))
            {
                 noLinkProblem=false;
                 elementsSections[i].triangles.removeAt(triangleIteratorPos);
                 std::cout <<"nodesReference: "<<currentTriangleElement.nodesReference.toUtf8().constData()<<" from "<< currentElementSection.name.toUtf8().constData()
                          <<" triangle element "<< currentTriangleElement.elementIndex<<" not declared."<<std::endl;
                 QString error=QString("nodesReference: "+currentTriangleElement.nodesReference+" from "+currentElementSection.name+" triangle element "+ QString::number(currentTriangleElement.elementIndex)+" not declared.");
                 emit errorMessage(error);
                 logStrings::errorLogString.append(error+"\r\n");
                 continue;
            }
            triangleIteratorPos=triangleIteratorPos+1;
        }
        QListIterator<Quadrilateral> quadrilateralElementIterator(currentElementSection.quadrilaterals);
        quint64 quadrilateralIteratorPos=0;
        while (quadrilateralElementIterator.hasNext())
        {
            Quadrilateral currentQuadrilateralElement=quadrilateralElementIterator.next();
            if(!nodesSectionsNames.contains(currentQuadrilateralElement.nodesReference,Qt::CaseInsensitive))
            {
                 noLinkProblem=false;
                 elementsSections[i].quadrilaterals.removeAt(quadrilateralIteratorPos);
                 std::cout <<"nodesReference: "<<currentQuadrilateralElement.nodesReference.toUtf8().constData()<<" from "<< currentElementSection.name.toUtf8().constData()
                          <<" quadrilateral element "<< currentQuadrilateralElement.elementIndex<<" not declared."<<std::endl;
                 QString error=QString("nodesReference: "+currentQuadrilateralElement.nodesReference+" from "+ currentElementSection.name+" quadrilateral element "+ QString::number(currentQuadrilateralElement.elementIndex)+" not declared.");
                 emit errorMessage(error);
                 logStrings::errorLogString.append(error+"\r\n");
                 continue;
            }
            quadrilateralIteratorPos=quadrilateralIteratorPos+1;
        }

        if(!currentElementSection.RefNodes.isEmpty())
        {
            if(!nodesSectionsNames.contains(currentElementSection.RefNodes,Qt::CaseInsensitive))
            {
                std::cout <<"nodesReference: "<<currentElementSection.RefNodes.toUtf8().constData()<<" from "<< currentElementSection.name.toUtf8().constData()<<" not declared."<<std::endl;
                QString error=QString("nodesReference: "+currentElementSection.RefNodes+" from "+ currentElementSection.name+" not declared.");
                emit errorMessage(error);
                logStrings::errorLogString.append(error+"\r\n");
            }
        }
    }
    return noLinkProblem;
}

void SolvingScriptReader::checkNodesIdentifiersAndSort()
{
    QStringList duplicates;


    for(int i =0;i<nodesSections.length();i++)
    {
        NodesSection currentNodesSection=nodesSections[i];
        nodesSections[i].scaleNodes(); //scale Nodes
//        if(!currentNodesSection.checkIdentifiersAllUnique(currentNodesSection.nodes))
            if(!currentNodesSection.checkIdentifiersAllUnique())
        {
            QVector<quint64> nonUniqueIdentifiers =currentNodesSection.getNonUniqueIdentifiers();
            QString errMess("Non unique nodes identifiers in NodesSection "+currentNodesSection.name+": ");
            std::cout <<"Non unique nodes identifiers in NodesSection "<<currentNodesSection.name.toUtf8().constData()<<": ";

            for (auto a : nonUniqueIdentifiers)
            {
                std::cout << a << " ";
                errMess.append(a+" ");
            }
            emit errorMessage(errMess);
            logStrings::errorLogString.append(errMess+"\r\n");
            std::cout <<std::endl;
        }
    }

    QStringList nodesSectionsNames=getNodesSectionsNames();
    for(int i=0;i<nodesSectionsNames.length()-1;i++)
    {
        for(int j=i+1;j<nodesSectionsNames.length();j++)
        {
            if(0==nodesSections.value(i).name.compare(nodesSections.value(j).name,Qt::CaseInsensitive))
            {
                duplicates.push_back(nodesSections.value(i).name);
            }
        }
    }
    duplicates.removeDuplicates();
    if(!duplicates.isEmpty())
    {
        std::cout <<"Non unique Nodes section names: ";
        QString errMess("Non unique Nodes section names: ");
        for(int i=0;i<duplicates.length();i++)
        {
            std::cout <<duplicates.value(i).toUtf8().constData() <<" ";
            errMess.append(duplicates.value(i)+" ");
        }
        emit errorMessage(errMess);
        logStrings::errorLogString.append(errMess+"\r\n");
        std::cout<<std::endl;
    }
}

void SolvingScriptReader::checkElementsIdentifiersAndSort()
{
    QStringList duplicates;
    QStringList elementSectionsNames=getElementSectionsNames();
    for(int i=0;i<elementSectionsNames.length()-1;i++)
    {
        for(int j=i+1;j<elementSectionsNames.length();j++)
        {
            if(0==nodesSections.value(i).name.compare(nodesSections.value(j).name,Qt::CaseInsensitive))
            {
                duplicates.push_back(nodesSections.value(i).name);
            }
        }
    }
    duplicates.removeDuplicates();
    if(!duplicates.isEmpty())
    {
        std::cout <<"Non unique Element section names: ";
        QString errMess("Non unique Element section names: ");

        for(int i=0;i<duplicates.length();i++)
        {
            std::cout <<duplicates.value(i).toUtf8().constData() <<" ";
            errMess.append(duplicates.value(i)+" ");
        }
        emit errorMessage(errMess);
        logStrings::errorLogString.append(errMess+"\r\n");
        std::cout<<std::endl;
    }
}

QStringList SolvingScriptReader::getNodesSectionsNames()
{
    QStringList nodesSectionsNames;
    QVectorIterator<NodesSection> nodesSectionsIterator(nodesSections);
    while (nodesSectionsIterator.hasNext())
    {
        NodesSection currentNodesSection=nodesSectionsIterator.next();
        nodesSectionsNames.push_back(currentNodesSection.name);
    }
    return nodesSectionsNames;
}

QStringList SolvingScriptReader::getElementSectionsNames()
{
    QStringList elementSectionsNames;
    QListIterator<ElementSection> elementSectionsIterator(elementsSections);
    while (elementSectionsIterator.hasNext())
    {
        ElementSection currentElementssSection=elementSectionsIterator.next();
        elementSectionsNames.push_back(currentElementssSection.name);
    }
    return elementSectionsNames;
}

QStringList SolvingScriptReader::getPressurePointsSectionsNames()
{
    QStringList pressurePointsSectionsNames;
    QVectorIterator<PressurePointsSection> pressurePointsSectionsIterator(pressurePointsSections);
    while (pressurePointsSectionsIterator.hasNext())
    {
        PressurePointsSection currentPressurePointsSection=pressurePointsSectionsIterator.next();
        pressurePointsSectionsNames.push_back(currentPressurePointsSection.name);
    }
    return pressurePointsSectionsNames;
}


bool SolvingScriptReader::setupConstBoundaryElements()
{
//    boundaryElements.clear();
    bool errors = false;
    QStringList nodesSectionsNames=getNodesSectionsNames();
    for(int i=0;i<elementsSections.length();i++)
    {
        ElementSection currentElementSection=elementsSections[i];
        if(currentElementSection.containsNoElements())
        {
            continue;
        }
//        boundaryElements.push_back(BoundaryElements(currentElementSection.name));
//        std::cout<<"aliases: "<<meshFileAlias.at(0).toUtf8().constData()<<std::endl;
        QListIterator<MeshFileElement> meshFileElementIterator(currentElementSection.meshFileElements);
        std::cout <<"meshelementlist length "<<currentElementSection.meshFileElements.length()<<std::endl;
        emit logMessage("Meshelementlist length "+QString::number(currentElementSection.meshFileElements.size()));
        while(meshFileElementIterator.hasNext())
        {
            MeshFileElement currentMeshFileElement = meshFileElementIterator.next();
            int index = meshFileAlias.indexOf(QRegularExpression(currentMeshFileElement.MeshFileAlias, QRegularExpression::CaseInsensitiveOption));
            QString meshFilePath = meshFiles.value(index);
//            std::cout<<"MeshFilePath: "<<meshFilePath.toUtf8().constData()<<std::endl;
//            std::cout<<"index: "<<index<<std::endl;
//            std::cout<<"aliases: "<<currentMeshFileElement.MeshFileAlias.toUtf8().constData()<<std::endl;
            meshFileReader.readMsh(meshFilePath, currentMeshFileElement.elementIndex, currentMeshFileElement.includeAll,currentMeshFileElement.include,currentMeshFileElement.exclude);
            QVector<VectorTriangle> elementsVectorTriangles = meshFileReader.getTriangles();
            QVector<VectorQuadrilateral> elementsVectorQuadrilaterals = meshFileReader.getQuadrilaterals();
            if(currentMeshFileElement.swapNormals ^ currentElementSection.swapNormals) // xor
            {
                std::cout<<"swapNormals"<<std::endl;
                emit logMessage("swapNormals");
                currentElementSection.swapTriangleNormals(elementsVectorTriangles);
                currentElementSection.swapQuadrilateralNormals(elementsVectorQuadrilaterals);
            }
            elementsSections[i].elementsVectorTriangles.append(elementsVectorTriangles);
            elementsSections[i].elementsVectorQuadrilaterals.append(elementsVectorQuadrilaterals);
        }
        elementsSections[i].meshFileElements.clear();

        QListIterator<Triangle> triangleElementIterator(currentElementSection.triangles);
        while (triangleElementIterator.hasNext())
        {
            Triangle currentTriangleElement=triangleElementIterator.next();
            QRegularExpression nodesReference("^"+currentTriangleElement.nodesReference+"$", QRegularExpression::CaseInsensitiveOption);
            int nodesSectionIndex=nodesSectionsNames.indexOf(nodesReference);
            if(nodesSectionIndex!=-1)
            {
                int firstNode=currentTriangleElement.node1;
                int secondNode=currentTriangleElement.node2;
                int thirdNode=currentTriangleElement.node3;

                bool allNodesDeclared=true;
                if(!nodesSections.at(nodesSectionIndex).nodes.contains(firstNode))
                {
                    std::cout <<"nodessection: "<<nodesSections.at(nodesSectionIndex).name.toUtf8().constData()<<" doesn't contain node "<<firstNode<<std::endl; allNodesDeclared=false;
                    QString error=QString("nodessection: "+nodesSections.at(nodesSectionIndex).name+" doesn't contain node "+QString::number(firstNode));
                    logStrings::errorLogString.append(error+"\r\n");
                    errors = true;
                }
                if(!nodesSections.at(nodesSectionIndex).nodes.contains(secondNode))
                {
                    std::cout <<"nodessection: "<<nodesSections.at(nodesSectionIndex).name.toUtf8().constData()<<" doesn't contain node "<<secondNode<<std::endl; allNodesDeclared=false;
                    QString error=QString("nodessection: "+nodesSections.at(nodesSectionIndex).name+" doesn't contain node "+QString::number(secondNode));
                    logStrings::errorLogString.append(error+"\r\n");
                    errors = true;
                }
                if(!nodesSections.at(nodesSectionIndex).nodes.contains(thirdNode))
                {
                    std::cout <<"nodessection: "<<nodesSections.at(nodesSectionIndex).name.toUtf8().constData()<<" doesn't contain node "<<thirdNode<<std::endl; allNodesDeclared=false;
                    QString error=QString("nodessection: "+nodesSections.at(nodesSectionIndex).name+" doesn't contain node "+QString::number(thirdNode));
                    logStrings::errorLogString.append(error+"\r\n");
                    errors = true;
                }

                if(allNodesDeclared==true)
                {
                    qint64 triangleIndex=currentTriangleElement.elementIndex;
                    Eigen::Vector3d node1=nodesSections.at(nodesSectionIndex).nodes.value(firstNode).coordinates;
                    Eigen::Vector3d node2=nodesSections.at(nodesSectionIndex).nodes.value(secondNode).coordinates;
                    Eigen::Vector3d node3=nodesSections.at(nodesSectionIndex).nodes.value(thirdNode).coordinates;
    //                boundaryElements.last().elementsVectorTriangles.push_back(VectorTriangle(node1,node2,node3));
                    if(currentElementSection.swapNormals)
                    {
                        elementsSections[i].elementsVectorTriangles.push_back(VectorTriangle(triangleIndex,node3,node2,node1));
                    }
                    else
                    {
                        elementsSections[i].elementsVectorTriangles.push_back(VectorTriangle(triangleIndex,node1,node2,node3));
                    }
                }
            }
            else
            {
                std::cout <<"nodesReference: "<<currentTriangleElement.nodesReference.toUtf8().constData()<<" from "<< currentElementSection.name.toUtf8().constData()
                          <<" triangle element "<< currentTriangleElement.elementIndex<<" not declared at time of boundary element setup."<<std::endl;
                QString error=QString("nodesReference: "+currentTriangleElement.nodesReference+currentElementSection.name+" triangle element "+QString::number(currentTriangleElement.elementIndex)+" not declared at time of boundary element setup.");
                emit errorMessage(error);
                logStrings::errorLogString.append(error+"\r\n");
                errors = true;
            }
        }
        elementsSections[i].triangles.clear();

        QListIterator<Quadrilateral> quadrilateralElementIterator(currentElementSection.quadrilaterals);
        while (quadrilateralElementIterator.hasNext())
        {
            Quadrilateral currentQuadrilateralElement=quadrilateralElementIterator.next();
            QRegularExpression nodesReference("^"+currentQuadrilateralElement.nodesReference+"$", QRegularExpression::CaseInsensitiveOption);
            int nodesSectionIndex=nodesSectionsNames.indexOf(nodesReference);
            if(nodesSectionIndex!=-1)
            {
                int firstNode=currentQuadrilateralElement.node1;
                int secondNode=currentQuadrilateralElement.node2;
                int thirdNode=currentQuadrilateralElement.node3;
                int fourthNode=currentQuadrilateralElement.node4;
                bool allNodesDeclared=true;
                if(!nodesSections.at(nodesSectionIndex).nodes.contains(firstNode))
                {
                    std::cout <<"nodessection: "<<nodesSections.at(nodesSectionIndex).name.toUtf8().constData()<<" doesn't contain node "<<firstNode<<std::endl; allNodesDeclared=false;
                    QString error=QString("nodessection: "+nodesSections.at(nodesSectionIndex).name+" doesn't contain node "+QString::number(firstNode));
                    logStrings::errorLogString.append(error+"\r\n");
                    errors = true;
                }
                if(!nodesSections.at(nodesSectionIndex).nodes.contains(secondNode))
                {
                    std::cout <<"nodessection: "<<nodesSections.at(nodesSectionIndex).name.toUtf8().constData()<<" doesn't contain node "<<secondNode<<std::endl; allNodesDeclared=false;
                    QString error=QString("nodessection: "+nodesSections.at(nodesSectionIndex).name+" doesn't contain node "+QString::number(secondNode));
                    logStrings::errorLogString.append(error+"\r\n");
                    errors = true;
                }
                if(!nodesSections.at(nodesSectionIndex).nodes.contains(thirdNode))
                {
                    std::cout <<"nodessection: "<<nodesSections.at(nodesSectionIndex).name.toUtf8().constData()<<" doesn't contain node "<<thirdNode<<std::endl; allNodesDeclared=false;
                    QString error=QString("nodessection: "+nodesSections.at(nodesSectionIndex).name+" doesn't contain node "+QString::number(thirdNode));
                    logStrings::errorLogString.append(error+"\r\n");
                    errors = true;
                }
                if(!nodesSections.at(nodesSectionIndex).nodes.contains(fourthNode))
                {
                    std::cout <<"nodessection: "<<nodesSections.at(nodesSectionIndex).name.toUtf8().constData()<<" doesn't contain node "<<fourthNode<<std::endl; allNodesDeclared=false;
                    QString error=QString("nodessection: "+nodesSections.at(nodesSectionIndex).name+" doesn't contain node "+QString::number(fourthNode));
                    logStrings::errorLogString.append(error+"\r\n");
                    errors = true;
                }

                if(allNodesDeclared==true)
                {
                    qint64 quadrilateralIndex=currentQuadrilateralElement.elementIndex;
                    Eigen::Vector3d node1=nodesSections.at(nodesSectionIndex).nodes.value(firstNode).coordinates;
                    Eigen::Vector3d node2=nodesSections.at(nodesSectionIndex).nodes.value(secondNode).coordinates;
                    Eigen::Vector3d node3=nodesSections.at(nodesSectionIndex).nodes.value(thirdNode).coordinates;
                    Eigen::Vector3d node4=nodesSections.at(nodesSectionIndex).nodes.value(fourthNode).coordinates;
//                    boundaryElements.last().elementsVectorQuadrilaterals.push_back(VectorQuadrilateral(node1,node2,node3,node4));
                    if(currentElementSection.swapNormals)
                    {
                        elementsSections[i].elementsVectorQuadrilaterals.push_back(VectorQuadrilateral(quadrilateralIndex,node4,node3,node2,node1));
                    }
                    else
                    {
                        elementsSections[i].elementsVectorQuadrilaterals.push_back(VectorQuadrilateral(quadrilateralIndex,node1,node2,node3,node4));
                    }
                }
            }
            else
            {
                std::cout <<"nodesReference: "<<currentQuadrilateralElement.nodesReference.toUtf8().constData()<<" from "<< currentElementSection.name.toUtf8().constData()
                               <<" of quadrilateral element "<< currentQuadrilateralElement.elementIndex<<" not found!"<<std::endl;
                QString error=QString("nodesReference: "+currentQuadrilateralElement.nodesReference+" from "+ currentElementSection.name+" of quadrilateral element "+ QString::number(currentQuadrilateralElement.elementIndex)+" not found!");
                logStrings::errorLogString.append(error+"\r\n");
                emit errorMessage(error);
                errors = true;
            }

        }
        elementsSections[i].quadrilaterals.clear();
    }
    for(int i=0;i<elementsSections.length();i++)
    {
        elementsSections[i].scaleAndTransformElements();
    }
    setupBoundaryConditions();
    setupDrivingElements();
    calculateCenterOfMass();
    return errors;
}

bool SolvingScriptReader::setupLinearBoundaryElements()
{
    LinearBoundaryElements boundaryElements;
    bool errors = false;
    QStringList nodesSectionsNames=getNodesSectionsNames();
    for(int i=0;i<elementsSections.length();i++)
    {
        ElementSection currentElementSection=elementsSections[i];
        if(currentElementSection.containsNoElements())
        {
            continue;
        }
//        boundaryElements.push_back(BoundaryElements(currentElementSection.name));
//        std::cout<<"aliases: "<<meshFileAlias.at(0).toUtf8().constData()<<std::endl;
        QListIterator<MeshFileElement> meshFileElementIterator(currentElementSection.meshFileElements);
        std::cout <<"meshelementlist length "<<currentElementSection.meshFileElements.length()<<std::endl;
        emit logMessage("Meshelementlist length "+QString::number(currentElementSection.meshFileElements.size()));
        while (meshFileElementIterator.hasNext())
        {
            MeshFileElement currentMeshFileElement=meshFileElementIterator.next();
            int index=meshFileAlias.indexOf(QRegularExpression(currentMeshFileElement.MeshFileAlias, QRegularExpression::CaseInsensitiveOption));
            QString meshFilePath=meshFiles.value(index);
//            std::cout<<"MeshFilePath: "<<meshFilePath.toUtf8().constData()<<std::endl;
//            std::cout<<"index: "<<index<<std::endl;
//            std::cout<<"aliases: "<<currentMeshFileElement.MeshFileAlias.toUtf8().constData()<<std::endl;
            meshFileReader.readMsh(meshFilePath, currentMeshFileElement.elementIndex, currentMeshFileElement.includeAll,currentMeshFileElement.include,currentMeshFileElement.exclude);
            QVector<VectorTriangle> elementsVectorTriangles=meshFileReader.getTriangles();
            QVector<VectorQuadrilateral> elementsVectorQuadrilaterals=meshFileReader.getQuadrilaterals();
            if(currentMeshFileElement.swapNormals ^ currentElementSection.swapNormals) // xor
            {
                std::cout<<"swapNormals"<<std::endl;
                emit logMessage("swapNormals");
                currentElementSection.swapTriangleNormals(elementsVectorTriangles);
                currentElementSection.swapQuadrilateralNormals(elementsVectorQuadrilaterals);
            }
            elementsSections[i].elementsVectorTriangles.append(elementsVectorTriangles);
            elementsSections[i].elementsVectorQuadrilaterals.append(elementsVectorQuadrilaterals);
        }
        elementsSections[i].meshFileElements.clear();

        QListIterator<Triangle> triangleElementIterator(currentElementSection.triangles);
        while (triangleElementIterator.hasNext())
        {
            Triangle currentTriangleElement=triangleElementIterator.next();
            QRegularExpression nodesReference("^"+currentTriangleElement.nodesReference+"$", QRegularExpression::CaseInsensitiveOption);
            int nodesSectionIndex=nodesSectionsNames.indexOf(nodesReference);
            if(nodesSectionIndex!=-1)
            {
                int firstNode=currentTriangleElement.node1;
                int secondNode=currentTriangleElement.node2;
                int thirdNode=currentTriangleElement.node3;

                bool allNodesDeclared=true;
                if(!nodesSections.at(nodesSectionIndex).nodes.contains(firstNode))
                {
                    std::cout <<"nodessection: "<<nodesSections.at(nodesSectionIndex).name.toUtf8().constData()<<" doesn't contain node "<<firstNode<<std::endl; allNodesDeclared=false;
                    QString error=QString("nodessection: "+nodesSections.at(nodesSectionIndex).name+" doesn't contain node "+QString::number(firstNode));
                    logStrings::errorLogString.append(error+"\r\n");
                    errors = true;
                }
                if(!nodesSections.at(nodesSectionIndex).nodes.contains(secondNode))
                {
                    std::cout <<"nodessection: "<<nodesSections.at(nodesSectionIndex).name.toUtf8().constData()<<" doesn't contain node "<<secondNode<<std::endl; allNodesDeclared=false;
                    QString error=QString("nodessection: "+nodesSections.at(nodesSectionIndex).name+" doesn't contain node "+QString::number(secondNode));
                    logStrings::errorLogString.append(error+"\r\n");
                    errors = true;
                }
                if(!nodesSections.at(nodesSectionIndex).nodes.contains(thirdNode))
                {
                    std::cout <<"nodessection: "<<nodesSections.at(nodesSectionIndex).name.toUtf8().constData()<<" doesn't contain node "<<thirdNode<<std::endl; allNodesDeclared=false;
                    QString error=QString("nodessection: "+nodesSections.at(nodesSectionIndex).name+" doesn't contain node "+QString::number(thirdNode));
                    logStrings::errorLogString.append(error+"\r\n");
                    errors = true;
                }

                if(allNodesDeclared==true)
                {
                    qint64 triangleIndex=currentTriangleElement.elementIndex;
                    Eigen::Vector3d node1=nodesSections.at(nodesSectionIndex).nodes.value(firstNode).coordinates;
                    Eigen::Vector3d node2=nodesSections.at(nodesSectionIndex).nodes.value(secondNode).coordinates;
                    Eigen::Vector3d node3=nodesSections.at(nodesSectionIndex).nodes.value(thirdNode).coordinates;
    //                boundaryElements.last().elementsVectorTriangles.push_back(VectorTriangle(node1,node2,node3));
                    if(currentElementSection.swapNormals)
                    {
                        elementsSections[i].elementsVectorTriangles.push_back(VectorTriangle(triangleIndex,node3,node2,node1));
                    }
                    else
                    {
                        elementsSections[i].elementsVectorTriangles.push_back(VectorTriangle(triangleIndex,node1,node2,node3));
                    }
                }
            }
            else
            {
                std::cout <<"nodesReference: "<<currentTriangleElement.nodesReference.toUtf8().constData()<<" from "<< currentElementSection.name.toUtf8().constData()
                          <<" triangle element "<< currentTriangleElement.elementIndex<<" not declared at time of boundary element setup."<<std::endl;
                QString error=QString("nodesReference: "+currentTriangleElement.nodesReference+currentElementSection.name+" triangle element "+QString::number(currentTriangleElement.elementIndex)+" not declared at time of boundary element setup.");
                emit errorMessage(error);
                logStrings::errorLogString.append(error+"\r\n");
                errors = true;
            }
        }
        elementsSections[i].triangles.clear();

        QListIterator<Quadrilateral> quadrilateralElementIterator(currentElementSection.quadrilaterals);
        while (quadrilateralElementIterator.hasNext())
        {
            Quadrilateral currentQuadrilateralElement=quadrilateralElementIterator.next();
            QRegularExpression nodesReference("^"+currentQuadrilateralElement.nodesReference+"$", QRegularExpression::CaseInsensitiveOption);
            int nodesSectionIndex=nodesSectionsNames.indexOf(nodesReference);
            if(nodesSectionIndex!=-1)
            {
                int firstNode=currentQuadrilateralElement.node1;
                int secondNode=currentQuadrilateralElement.node2;
                int thirdNode=currentQuadrilateralElement.node3;
                int fourthNode=currentQuadrilateralElement.node4;
                bool allNodesDeclared=true;
                if(!nodesSections.at(nodesSectionIndex).nodes.contains(firstNode))
                {
                    std::cout <<"nodessection: "<<nodesSections.at(nodesSectionIndex).name.toUtf8().constData()<<" doesn't contain node "<<firstNode<<std::endl; allNodesDeclared=false;
                    QString error=QString("nodessection: "+nodesSections.at(nodesSectionIndex).name+" doesn't contain node "+QString::number(firstNode));
                    logStrings::errorLogString.append(error+"\r\n");
                    errors = true;
                }
                if(!nodesSections.at(nodesSectionIndex).nodes.contains(secondNode))
                {
                    std::cout <<"nodessection: "<<nodesSections.at(nodesSectionIndex).name.toUtf8().constData()<<" doesn't contain node "<<secondNode<<std::endl; allNodesDeclared=false;
                    QString error=QString("nodessection: "+nodesSections.at(nodesSectionIndex).name+" doesn't contain node "+QString::number(secondNode));
                    logStrings::errorLogString.append(error+"\r\n");
                    errors = true;
                }
                if(!nodesSections.at(nodesSectionIndex).nodes.contains(thirdNode))
                {
                    std::cout <<"nodessection: "<<nodesSections.at(nodesSectionIndex).name.toUtf8().constData()<<" doesn't contain node "<<thirdNode<<std::endl; allNodesDeclared=false;
                    QString error=QString("nodessection: "+nodesSections.at(nodesSectionIndex).name+" doesn't contain node "+QString::number(thirdNode));
                    logStrings::errorLogString.append(error+"\r\n");
                    errors = true;
                }
                if(!nodesSections.at(nodesSectionIndex).nodes.contains(fourthNode))
                {
                    std::cout <<"nodessection: "<<nodesSections.at(nodesSectionIndex).name.toUtf8().constData()<<" doesn't contain node "<<fourthNode<<std::endl; allNodesDeclared=false;
                    QString error=QString("nodessection: "+nodesSections.at(nodesSectionIndex).name+" doesn't contain node "+QString::number(fourthNode));
                    logStrings::errorLogString.append(error+"\r\n");
                    errors = true;
                }

                if(allNodesDeclared==true)
                {
                    qint64 quadrilateralIndex=currentQuadrilateralElement.elementIndex;
                    Eigen::Vector3d node1=nodesSections.at(nodesSectionIndex).nodes.value(firstNode).coordinates;
                    Eigen::Vector3d node2=nodesSections.at(nodesSectionIndex).nodes.value(secondNode).coordinates;
                    Eigen::Vector3d node3=nodesSections.at(nodesSectionIndex).nodes.value(thirdNode).coordinates;
                    Eigen::Vector3d node4=nodesSections.at(nodesSectionIndex).nodes.value(fourthNode).coordinates;
//                    boundaryElements.last().elementsVectorQuadrilaterals.push_back(VectorQuadrilateral(node1,node2,node3,node4));
                    if(currentElementSection.swapNormals)
                    {
                        elementsSections[i].elementsVectorQuadrilaterals.push_back(VectorQuadrilateral(quadrilateralIndex,node4,node3,node2,node1));
                    }
                    else
                    {
                        elementsSections[i].elementsVectorQuadrilaterals.push_back(VectorQuadrilateral(quadrilateralIndex,node1,node2,node3,node4));
                    }
                }
            }
            else
            {
                std::cout <<"nodesReference: "<<currentQuadrilateralElement.nodesReference.toUtf8().constData()<<" from "<< currentElementSection.name.toUtf8().constData()
                               <<" of quadrilateral element "<< currentQuadrilateralElement.elementIndex<<" not found!"<<std::endl;
                QString error=QString("nodesReference: "+currentQuadrilateralElement.nodesReference+" from "+ currentElementSection.name+" of quadrilateral element "+ QString::number(currentQuadrilateralElement.elementIndex)+" not found!");
                logStrings::errorLogString.append(error+"\r\n");
                emit errorMessage(error);
                errors = true;
            }

        }
        elementsSections[i].quadrilaterals.clear();
    }
    for(int i=0;i<elementsSections.length();i++)
    {
        elementsSections[i].scaleAndTransformElements();
    }
    setupBoundaryConditions();
    setupDrivingElements();
    calculateCenterOfMass();
    return errors;
}

void SolvingScriptReader::setupBoundaryConditions()
{
    std::cout <<"Setting up wall impedances."<<std::endl;
    QString error=QString("Setting up wall impedances.");
    logStrings::errorLogString.append(error+"\r\n");
    emit errorMessage(error);

    QStringList elementSectionsNames = getElementSectionsNames();
    std::cout <<"Number of WallImpedance sections: "<<wallImpedanceSections.length()<<std::endl;
    for(int j=0;j<wallImpedanceSections.length();j++)
    {
        WallImpedanceSection& currentWallImpedanceSection=wallImpedanceSections[j];
        currentWallImpedanceSection.setWavenumer((2.0*frequencies.last() * global::PI)/c);
        currentWallImpedanceSection.setMediumImpedance(mediumImpedance);
        currentWallImpedanceSection.setUpValues();
        if(currentWallImpedanceSection.wallImpecanceElements.length()>=1)
        {
            for(int i=0;i<currentWallImpedanceSection.wallImpecanceElements.length();i++)
            {
                WallImpecanceElement currentWallImpecanceElement=currentWallImpedanceSection.wallImpecanceElements.at(i);
                QRegularExpression elementSectionReference("^"+currentWallImpecanceElement.RefElements+"$", QRegularExpression::CaseInsensitiveOption);
                int elementSectionIndex=elementSectionsNames.indexOf(elementSectionReference);
                if(elementSectionIndex!=-1 && elementSectionIndex<elementsSections.size())
                {
                    RobinBoundaryCondition currentBoundaryCondition=currentWallImpecanceElement.robinBoundaryCondition;
                    int elementIndex=currentWallImpecanceElement.elementIndex;
                    elementsSections[elementSectionIndex].setTrianglesWallimpedance(elementIndex,currentBoundaryCondition);
                    elementsSections[elementSectionIndex].setQuadrilateralsWallimpedance(elementIndex,currentBoundaryCondition);
                }
                else
                {
                    std::cout <<"ElementReference: "<<currentWallImpecanceElement.RefElements.toUtf8().constData()<<" from WallImpedance section "<< currentWallImpedanceSection.name.toUtf8().constData()<<" Element "<<currentWallImpecanceElement.index<<" not found!"<<std::endl;
                    QString error=QString("ElementReference: "+currentWallImpecanceElement.RefElements+" from WallImpedance section "+ currentWallImpedanceSection.name+" Element "+QString::number(currentWallImpecanceElement.index)+" not found!");
                    logStrings::errorLogString.append(error+"\r\n");
                    emit errorMessage(error);
                }
            }
        }
        else //set boundary condition for entire elements section
        {
            QRegularExpression elementSectionReference("^"+currentWallImpedanceSection.RefElements+"$", QRegularExpression::CaseInsensitiveOption);
            int elementSectionIndex=elementSectionsNames.indexOf(elementSectionReference);
            if(elementSectionIndex!=-1 && elementSectionIndex<elementsSections.size())
            {
                RobinBoundaryCondition currentBoundaryCondition = currentWallImpedanceSection.robinBoundaryCondition;
                elementsSections[elementSectionIndex].setLocalWallimpedance(currentBoundaryCondition);
            }
            else
            {
                std::cout <<"ElementReference: "<<currentWallImpedanceSection.RefElements.toUtf8().constData()<<" from "<< currentWallImpedanceSection.name.toUtf8().constData()<<" not found!"<<std::endl;
                QString error=QString("ElementReference: "+currentWallImpedanceSection.RefElements+" from "+ currentWallImpedanceSection.name+" not found!");
                logStrings::errorLogString.append(error+"\r\n");
                emit errorMessage(error);
            }
        }
    }
}

void SolvingScriptReader::setupDrivingElements()
{
    std::cout <<"Setting up driving elements."<<std::endl;
    QString error=QString("Setting up driving elements.");
    logStrings::errorLogString.append(error+"\r\n");
    emit errorMessage(error);

    QStringList elementSectionsNames= getElementSectionsNames();
    std::cout <<"Number of Driving sections: "<<drivingSections.length()<<std::endl;
    for(int j=0;j<drivingSections.length();j++)
    {
        DrivingSection& currentDrivingSection=drivingSections[j];
//        currentWallImpedanceSection.setWavenumer((2.0*frequencies.last()*global::PI)/c);
//        currentWallImpedanceSection.setMediumImpedance(mediumImpedance);
//        currentWallImpedanceSection.setUpValues();
        if(currentDrivingSection.drivingElements.length()>=1)
        {
            for(int i=0;i<currentDrivingSection.drivingElements.length();i++)
            {
                DrivingElement currentDrivingElement=currentDrivingSection.drivingElements.at(i);
                QRegularExpression elementSectionReference("^"+currentDrivingElement.RefElements+"$", QRegularExpression::CaseInsensitiveOption);
                int elementSectionIndex=elementSectionsNames.indexOf(elementSectionReference);
                if(elementSectionIndex!=-1 && elementSectionIndex<elementsSections.size())
                {
                    std::complex<double> drivingWeight=currentDrivingElement.Weight * exp(std::complex<double>(0,-1) * (currentDrivingElement.Delay * frequencies.last() * 2.0 * M_PI));
//                    std::cout<<"Delay influence: "<<exp(std::complex<double>(0,-1) * currentDrivingElement.Delay * frequencies.last()*2.0*M_PI)<<std::endl;
//                    std::cout<<"Delay in s: "<<currentDrivingSection.DrvDelay<<std::endl;
//                    std::cout<<"frequency in Hertz: "<<frequencies.last()<<std::endl;
                    int elementIndex=currentDrivingElement.elementIndex;
                    elementsSections[elementSectionIndex].setTrianglesDrivingWeight(elementIndex,drivingWeight);
                    elementsSections[elementSectionIndex].setQuadrilateralsDrivingWeight(elementIndex,drivingWeight);
                }
                else
                {
                    std::cout <<"ElementReference: "<<currentDrivingElement.RefElements.toUtf8().constData()<<" from Driving section "<< currentDrivingSection.name.toUtf8().constData()<<" Element "<<currentDrivingElement.index<<" not found!"<<std::endl;
                    QString error=QString("ElementReference: "+currentDrivingElement.RefElements+" from Driving section "+ currentDrivingSection.name+" Element "+QString::number(currentDrivingElement.index)+" not found!");
                    logStrings::errorLogString.append(error+"\r\n");
                    emit errorMessage(error);
                }
            }
        }
        else //set driving weights for entire elements section
        {
            QRegularExpression elementSectionReference("^"+currentDrivingSection.RefElements+"$", QRegularExpression::CaseInsensitiveOption);
            int elementSectionIndex=elementSectionsNames.indexOf(elementSectionReference);
            if(elementSectionIndex!=-1 && elementSectionIndex<elementsSections.size())
            {
                std::complex<double> drivingWeight=currentDrivingSection.DrvWeight * exp(std::complex<double>(0,-1)*(currentDrivingSection.DrvDelay * frequencies.last() * 2*M_PI));
                std::cout<<"Delay influence: "<<exp(std::complex<double>(0,-1)*(currentDrivingSection.DrvDelay * frequencies.last() * 2*M_PI))<<std::endl;
                elementsSections[elementSectionIndex].setLocalDrivingWeight(drivingWeight);
            }
            else
            {
                std::cout <<"ElementReference: "<<currentDrivingSection.RefElements.toUtf8().constData()<<" from "<< currentDrivingSection.name.toUtf8().constData()<<" not found!"<<std::endl;
                QString error=QString("ElementReference: "+currentDrivingSection.RefElements+" from "+ currentDrivingSection.name+" not found!");
                logStrings::errorLogString.append(error+"\r\n");
                emit errorMessage(error);
            }
        }
    }
}

void SolvingScriptReader::setupPressurePoints()
{
    QStringList nodesSectionsNames=getNodesSectionsNames();
    std::cout <<"Number of PressurePoints sections: "<<pressurePointsSections.length()<<std::endl;
    QString message=QString("Number of PressurePoints sections: "+QString::number(pressurePointsSections.length()));
    logStrings::logString.append(message+"\r\n");
    emit logMessage(message);

    for(int i=0;i<pressurePointsSections.length();i++)
    {
        PressurePointsSection currentPressurePointsSection=pressurePointsSections[i];
        for(int j=0;j<currentPressurePointsSection.pressurePoints.size();j++)
        {
            PressurePoint currentPressurePoint=currentPressurePointsSection.pressurePoints.at(j);
            QRegularExpression nodesReference("^"+currentPressurePoint.RefNodes+"$", QRegularExpression::CaseInsensitiveOption);
            int nodesSectionIndex=nodesSectionsNames.indexOf(nodesReference);
            if(nodesSectionIndex!=-1)
            {
                int node=currentPressurePoint.nodeIndex;

                bool nodeDeclared=true;
                if(!nodesSections.at(nodesSectionIndex).nodes.contains(node)){std::cout <<"nodessection: "<<nodesSections.at(nodesSectionIndex).name.toUtf8().constData()<<" doesn't contain node "<<node<<std::endl; nodeDeclared=false;}

                if(nodeDeclared==true)
                {
                    int pointSourceIndex=currentPressurePoint.index;
                    double weight=currentPressurePoint.weight;
                    double delay=currentPressurePoint.delay;
                    Eigen::Vector3d pointSourcePosition=nodesSections.at(nodesSectionIndex).nodes.value(node).coordinates;
                    pointSources.push_back(PointSource(pointSourceIndex,weight,delay,pointSourcePosition));
                }
            }
            else
            {
                std::cout <<"nodesReference: "<<currentPressurePoint.RefNodes.toUtf8().constData()<<" from Pressure_Points section "<< currentPressurePointsSection.name.toUtf8().constData()
                               <<"of pressure point "<< currentPressurePoint.index<<" not found."<<std::endl;
                QString error=QString("nodesReference: "+currentPressurePoint.RefNodes+" from Pressure_Points section "+ currentPressurePointsSection.name+"of pressure point "+ QString::number(currentPressurePoint.index)+" not found.");
                logStrings::errorLogString.append(error+"\r\n");
                emit errorMessage(error);
            }
        }
    }
}

BoundaryElements SolvingScriptReader::getBoundaryElements()
{
    BoundaryElements returnElements("Test");

    for(int i =0; i< elementsSections.length();i++)
    {
        returnElements.triangles.append(elementsSections.at(i).elementsVectorTriangles);
        returnElements.exteriorProblem=elementsSections.at(i).ElType;//not finished
    }
    if(infiniteBaffleSections.length() >= 1)
    {
        InfiniteBaffleSection infiniteBaffleSection = infiniteBaffleSections.at(0);
        ImpedancePlane impedancePlane;
        impedancePlane.halfSpacePlanePoint = infiniteBaffleSection.planePoint;
        impedancePlane.halfSpacePlaneNormal = infiniteBaffleSection.planeNormal.normalized();
        returnElements.impedancePlanes.append(impedancePlane);

        if(infiniteBaffleSection.distanceToParallelPlane > 0)
        {
            impedancePlane.halfSpacePlanePoint += infiniteBaffleSection.distanceToParallelPlane * infiniteBaffleSection.planeNormal.normalized();
            impedancePlane.halfSpacePlaneNormal *= -1;
            returnElements.impedancePlanes.append(impedancePlane);
            if(infiniteBaffleSections.length() >= 2)
            {
                std::cout << "Two parallel impedance planes have already been defined via the distanceToParallelPlane variable from Infinite_Baffle section " << infiniteBaffleSection.name.toStdString() << ". Further Infinite_Baffle sections will be ignored." << std::endl;
                QString error=QString("Two parallel impedance planes have already been defined via the distanceToParallelPlane variable from Infinite_Baffle section "+ infiniteBaffleSection.name+". Further Infinite_Baffle sections will be ignored.");
                emit errorMessage(error);
                logStrings::errorLogString.append(error+"\r\n");
            }
        }
        else if(infiniteBaffleSections.length() >= 2)
        {
            InfiniteBaffleSection infiniteBaffleSection = infiniteBaffleSections.at(1);

            impedancePlane.halfSpacePlanePoint = infiniteBaffleSection.planePoint;
            impedancePlane.halfSpacePlaneNormal *= -1;
            returnElements.impedancePlanes.append(impedancePlane);
        }

    }
    return returnElements;
}

void SolvingScriptReader::refineElements()
{
    for(int i=0; i<elementsSections.length(); i++)
    {
        ElementSection currentElementSection = elementsSections[i];
        if(currentElementSection.MeshFrequency==0 && currentElementSection.EdgeLength==0)
        {
            elementsSections[i].EdgeLength = globalEdgelength;
        }
        else if(currentElementSection.MeshFrequency != 0)
        {
            elementsSections[i].EdgeLength = (c/currentElementSection.MeshFrequency)/2;
        }
        QVector<VectorTriangle> triangles = elementsSections[i].refineQuadrilaterals(elementsSections[i].elementsVectorQuadrilaterals);
        elementsSections[i].elementsVectorQuadrilaterals.clear();
        elementsSections[i].elementsVectorTriangles.append(triangles);
        for(int j=0; j<elementsSections[i].elementsVectorTriangles.size(); j++)
        {
            elementsSections[i].elementsVectorTriangles[j].normal = global::normalVectorOfTriangle(elementsSections[i].elementsVectorTriangles[j]);
//            std::cout<<"normal; "<<elementsSections[i].elementsVectorTriangles[j].normal(0)<<" "<<elementsSections[i].elementsVectorTriangles[j].normal(1)<<" "<<elementsSections[i].elementsVectorTriangles[j].normal(2)<<std::endl;
        }
        elementsSections[i].elementsVectorTriangles=elementsSections[i].refineTriangles(elementsSections[i].elementsVectorTriangles);
//        elementsSections[i].elementsVectorTriangles=MeshFunctions::refineTriangles(elementsSections[i].elementsVectorTriangles,elementsSections[i].EdgeLength*0.9);
    }
}

void SolvingScriptReader::calculateCenterOfMass()
{
    QVector<VectorTriangle> triangles;
    for(int i=0; i<elementsSections.length(); i++)
    {
        ElementSection currentElementSection = elementsSections[i];
        QVector<VectorTriangle> tmpTriangles;
        tmpTriangles.append(currentElementSection.elementsVectorTriangles);
        tmpTriangles.append(global::quadrilateralsToTriangles(currentElementSection.elementsVectorQuadrilaterals));
        if(currentElementSection.ElType==false)//interior
        {
            currentElementSection.swapTriangleNormals(tmpTriangles);
        }
        triangles.append(tmpTriangles);
    }
    int numberOfTriangles = triangles.size();
    QPair<Eigen::Vector3d,double> centerOfMassAndVolume = global::centerAndVolumeOfMassOfTriangleObject(triangles);
    centerOfMass = centerOfMassAndVolume.first;
    objectVolume = centerOfMassAndVolume.second;

    containingSphereRadius = global::calculateMaxRelativeDistance(triangles, centerOfMass);

    if(numberOfTriangles == 0)
    {
        containingSphereRadius = 1;
        centerOfMass = Eigen::Vector3d(0,0,0);
    }
    std::cout << "Radius of containing sphere: " << containingSphereRadius << " m."<<std::endl;
    QString message=QString("Radius of containing sphere: "+QString::number(containingSphereRadius)+" m.");
    logStrings::logString.append(message+"\r\n");

    emit logMessage(message);

    triangles.clear();
}
