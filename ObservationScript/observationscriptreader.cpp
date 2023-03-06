#include "ObservationScript/observationscriptreader.h"

ObservationScriptReader::ObservationScriptReader(QObject* parent): QObject(parent)
{

}

bool ObservationScriptReader::setMeshFilesAndAliases( const QStringList& meshFiles, const QStringList& meshFileAlias)
{
    this->meshFiles=meshFiles;
    this->meshFileAlias=meshFileAlias;
    meshfilesAndAliasesSet=true;
    if(meshFiles.length() != meshFileAlias.length())
    {
        meshfilesAndAliasesSet=false;
        std::cout <<"The list of the meshfile links and aliases have different lengths!"<<std::endl;
        logStrings::errorLogString.append("The list of the meshfile links and aliases have different lengths!\r\n");
    }
    clear();
    return meshfilesAndAliasesSet;
}

void ObservationScriptReader::readObservationScript(const QStringList &fileName)
{
    if(!meshfilesAndAliasesSet){
        std::cout <<"Meshfile links were not passed to solvingscript reader!"<<std::endl;
        logStrings::errorLogString.append("Meshfile links were not passed to solvingscript reader!\r\n");
        emit errorMessage("Meshfile links were not passed to solvingscript reader!");

        return;
    }
    clear();
    controlSpectrumSectionDeclared=false;

    if(fileName.size()==0)
    {
        std::cout <<"No observation scripts declared!"<<std::endl;
        logStrings::errorLogString.append("No observation scripts declared!\r\n");
        emit logMessage("No observation scripts declared!");
        return;
    }
    QString scriptFileString;
    for(int i=0;i<fileName.size();i++)
    {
        QFile projFile(fileName.at(i));
        if (!projFile.open(QIODevice::ReadOnly | QIODevice::Text))
        {
            std::cout <<"Couldn't read observation script : "<< fileName.at(i).toUtf8().constData()<<std::endl;
            emit errorMessage("Couldn't read observation script : " + fileName.at(i));
            logStrings::errorLogString.append("Couldn't read observation script : " + fileName.at(i)+"\r\n");
            return;
        }
        QTextStream in(&projFile);
        QString tmpScriptFileString=in.readAll();
        scriptFileString.append(tmpScriptFileString);
    }

    emit logMessage("\r\n** Reading Observation_Scripts **");
    logStrings::logString.append("\r\n** Reading Observation_Scripts **\r\n");
    logStrings::errorLogString.append("\r\n** Reading Observation_Scripts **\r\n");
    emit errorMessage("\r\nObservation_Scripts:\r\n");

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
//    controlSolverSection.allMandatoryParametersSet();
}

bool ObservationScriptReader::newSectionDeclared(const QString& currentLine)
{
    QRegularExpressionMatch match;
    QRegularExpression identifier=sectionIdentifierRegEx("Control_Spectrum");

    if(currentLine.contains(identifier,&match))
    {
//        QString name=match.captured(1);
//        global::removeQuotesAndWhitespace(value);
        if(controlSpectrumSectionDeclared)
        {
            std::cout <<"Control_Spectrum section already declared!"<<std::endl;
            logStrings::errorLogString.append("\r\nControl_Spectrum section already declared!\r\n");
            emit errorMessage("\r\nControl_Spectrum section already declared!");
            currentSection=Control_Spectrum;
        }
        else{
            controlSpectrumSectionDeclared=true;
            std::cout <<"Control_Spectrum section read!"<<std::endl;
            currentSection=Control_Spectrum;
            emit logMessage("\r\nControl_Spectrum section read!");
            logStrings::logString.append("\r\nControl_Spectrum section read!\r\n");
            logStrings::errorLogString.append("\r\nControl_Spectrum section read!\r\n");
        }
        return true;
    }

    identifier=sectionIdentifierRegEx("Control_Field");
    if(currentLine.contains(identifier,&match))
    {
//        QString name=match.captured(1);
//        global::removeQuotesAndWhitespace(value);
        if(controlFieldSectionDeclared)
        {
            std::cout <<"Control_Field section already declared!"<<std::endl;
            emit errorMessage("\r\nControl_Field section already declared!");
            logStrings::errorLogString.append("\r\nControl_Field section already declared!\r\n");
            currentSection=Control_Field;
        }
        else
        {
            controlFieldSectionDeclared=true;
            std::cout <<"Control_Spectrum section read!"<<std::endl;
            logStrings::logString.append("Control_Spectrum section read!\r\n");
            logStrings::errorLogString.append("Control_Spectrum section read!\r\n");
            emit logMessage("\r\nControl_Spectrum section read!");
            currentSection=Control_Field;
        }
        return true;
    }

    identifier=sectionIdentifierRegEx("Field");
    if(currentLine.contains(identifier,&match))
    {

        QString name= match.captured(1);
        fieldSections.push_back((FieldSection(name)));
        std::cout <<"Field section read! Name: "<<name.toUtf8().constData()<<std::endl;
        logStrings::logString.append("\r\nField section " +name + " read!\r\n");
        logStrings::errorLogString.append("\r\nField section " +name + " read!\r\n");
        emit logMessage("\r\nField section " +name + " read!");
        currentSection=Field;
        return true;
    }

    identifier=sectionIdentifierRegEx("BE_Spectrum");
    if(currentLine.contains(identifier,&match))
    {
        QString name= match.captured(1);
        bESpectrumSections.push_back(BESpectrumSection(name));
        std::cout <<"BE_Spectrum section read! Name: "<<name.toUtf8().constData()<<std::endl;
        logStrings::logString.append("\r\nBE_Spectrum section " +name + " read!\r\n");
        logStrings::errorLogString.append("\r\nBE_Spectrum section " +name + " read!\r\n");
        emit logMessage("\r\nBE_Spectrum section " +name + " read!");
        currentSection=BE_Spectrum;
        return true;
    }

    identifier=sectionIdentifierRegEx("Nodes");
    if(currentLine.contains(identifier,&match))
    {
        QString name= match.captured(1);
        global::removeQuotes(name);
        nodesSections.push_back(NodesSection(name));
        std::cout <<"Nodes section read! Name: "<<name.toUtf8().constData()<<std::endl;
        logStrings::logString.append("\r\nNodes section " +name + " read!\r\n");
        logStrings::errorLogString.append("\r\nNodes section " +name + " read!\r\n");
        emit logMessage("\r\nNodes section " +name + " read!");
        currentSection=Nodes;
        return true;
    }

    identifier=sectionIdentifierRegEx("Driving_Values");
    if(currentLine.contains(identifier,&match))
    {
        QString name= match.captured(1);
        global::removeQuotes(name);
        std::cout <<"Driving_Values section read! Name: "<<name.toUtf8().constData()<<std::endl;
        logStrings::logString.append("\r\nDriving_Values section " +name + " read!\r\n");
        logStrings::errorLogString.append("\r\nDriving_Values section " +name + " read!\r\n");
        emit logMessage("\r\nDriving_Values section " +name + " read!");
        currentSection=Driving_Values;
        return true;
    }

    identifier=sectionIdentifierRegEx("MeshFile_Properties");
    if(currentLine.contains(identifier,&match))
    {
        QString name= match.captured(1);
        global::removeQuotes(name);
        meshFilePropertiesSections.push_back(MeshFilePropertiesSection(name));
        logStrings::logString.append("\r\nMeshFile_Properties " +name + " read!\r\n");
        logStrings::errorLogString.append("\r\nMeshFile_Properties " +name + " read!\r\n");
        emit logMessage("\r\nMeshFile_Properties " +name + " read!");
        currentSection=MeshFile_Properties;
        return true;
    }

    identifier=sectionIdentifierRegEx("Radiation_Impedance");
    if(currentLine.contains(identifier,&match))
    {
        QString name= match.captured(1);
        global::removeQuotes(name);
        std::cout <<"Radiation_Impedance section read!"<<std::endl;
        logStrings::logString.append("\r\nRadiation_Impedance section " +name + " read!\r\n");
        logStrings::errorLogString.append("\r\nRadiation_Impedance section " +name + " read!\r\n");
        emit logMessage("\r\nRadiation_Impedance section " +name + " read!");
        currentSection=Radiation_Impedance;
        return true;
    }

    identifier=sectionIdentifierRegEx("off");
    if(currentLine.contains(identifier,&match))
    {
        std::cout <<"off!"<<std::endl;
        currentSection=none;
        return true;
    }
    return false;
}

QRegularExpression ObservationScriptReader::sectionIdentifierRegEx(const QString& identifier)
{
//    return QRegularExpression("^\\s*"+identifier+"\\s*={1}\\s*([\"\'][\\w\\s]+[\"\']|\\w+)", QRegularExpression::CaseInsensitiveOption);
    return QRegularExpression("^\\s*"+identifier+"\\s*([\"\'][^\"\']+[\"\']|[\\S]+|$)", QRegularExpression::CaseInsensitiveOption);
}

void ObservationScriptReader::handleLineAccordingToCurrentSection(const QString& line)
{
            switch(currentSection)
            {
            case none:
                break;

            case Control_Spectrum:
                controlSpectrumSection.handleScriptLine(line);
                break;

            case Control_Field:
                controlFieldSection.handleScriptLine(line);
                break;

            case Field:
                fieldSections.last().handleScriptLine(line);
                break;

            case BE_Spectrum:
                bESpectrumSections.last().handleScriptLine(line);
                break;

            case Nodes:
                nodesSections.last().handleScriptLine(line);
                break;

            case Driving_Values:
//                nodesSections.last().handleScriptLine(line);
                break;

            case MeshFile_Properties:
                meshFilePropertiesSections.last().handleScriptLine(line);
                break;

            case Radiation_Impedance:
//                pressurePointsSections.last().handleScriptLine(line);
                break;

            }
}

void ObservationScriptReader::clear()
{

    controlFieldSection=ControlFieldSection();
    controlSpectrumSection=ControlSpectrumSection();
//    elementsSections.clear();
    nodesSections.clear();
//    subdomainPropertiesSections.clear();
    meshFilePropertiesSections.clear();
    bESpectrumSections.clear();
    fieldSections.clear();
//    wallImpedanceSections.clear();
//    pressurePointsSections.clear();

//    pointSources.clear();
//    boundaryElements.clear();

    currentSection=none;
    controlSpectrumSectionDeclared=false;
    controlFieldSectionDeclared=false;
}

QStringList ObservationScriptReader::getNodesSectionsNames()
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

void ObservationScriptReader::checkNodesIdentifiersAndSort()
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
            QString errMess("Non unique nodes identifiers in NodesSection "+currentNodesSection.name);
            std::cout <<"Non unique nodes identifiers in NodesSection "<<currentNodesSection.name.toUtf8().constData()<<": ";
            for (auto a : nonUniqueIdentifiers)
            {
                std::cout << a << " ";
                errMess.append(a+" ");
            }
            std::cout <<std::endl;
            emit errorMessage(errMess);
            logStrings::errorLogString.append(errMess+"\r\n");
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
        std::cout <<"Non unique NodesSection names: ";
        QString errMess("Non unique NodesSection names: ");

        for(int i=0;i<duplicates.length();i++)
        {
            std::cout <<duplicates.value(i).toUtf8().constData() <<" ";
            errMess.append(duplicates.value(i)+" ");
        }
        std::cout<<std::endl;
        emit errorMessage(errMess);
        logStrings::errorLogString.append(errMess+"\r\n");
    }
}

bool ObservationScriptReader::removeInvalidElements()
{
    std::cout <<"Removing possibly invalid/unlinkable observation elements."<<std::endl;
    emit errorMessage("Removing possibly invalid/unlinkable observation elements.");
    bool noLinkProblem=true;
    for(int i=0;i<fieldSections.length();i++)
    {

        QVectorIterator<MeshField> meshFieldIterator(fieldSections.at(i).meshFields);
        quint64 iteratorPos=0;
        while (meshFieldIterator.hasNext())
        {
            MeshField currentMeshField=meshFieldIterator.next();
            if(!meshFileAlias.contains(currentMeshField.MeshFileAlias, Qt::CaseInsensitive))
            {
                 noLinkProblem=false;
                 fieldSections[i].meshFields.removeAt(iteratorPos); //remove field with invalid meshfile alias
                 std::cout <<"MeshFileAlias: "<<currentMeshField.MeshFileAlias.toUtf8().constData()<<" from Field_Section "<< fieldSections.at(i).name.toUtf8().constData()
                          <<", field number "<< currentMeshField.index<<" not declared."<<std::endl;
                 emit errorMessage("MeshFileAlias: "+currentMeshField.MeshFileAlias+" from Field_Section "+ fieldSections.at(i).name+", field number "+QString::number(currentMeshField.index)+" not declared.");
                 logStrings::errorLogString.append("MeshFileAlias: "+currentMeshField.MeshFileAlias+" from Field_Section "+ fieldSections.at(i).name+", field number "+QString::number(currentMeshField.index)+" not declared.\r\n");
                 continue;
            }
            iteratorPos=iteratorPos+1;
        }


//        std::cout <<"Field name:"<<fieldSections.at(i).name.toUtf8().constData()<<std::endl;
//        std::cout <<"MeshFileAlias:"<<fieldSections.at(i).MeshFileAlias.toUtf8().constData()<<std::endl;
        if(!fieldSections.at(i).MeshFileAlias.isEmpty())
        {
            if(!(meshFileAlias.contains(fieldSections.at(i).MeshFileAlias, Qt::CaseInsensitive)))
            {
                 noLinkProblem=false;
                 std::cout <<"MeshFileAlias: "<<fieldSections.at(i).MeshFileAlias.toUtf8().constData()<<" from "<< fieldSections.at(i).name.toUtf8().constData()<<" not declared."<<std::endl;
                 emit errorMessage("MeshFileAlias: "+fieldSections.at(i).MeshFileAlias+" from section "+fieldSections.at(i).name+" not declared.");
                 logStrings::errorLogString.append("MeshFileAlias: "+fieldSections.at(i).MeshFileAlias+" from section "+fieldSections.at(i).name+" not declared.\r\n");
            }
        }
    }

    QStringList nodesSectionsNames=getNodesSectionsNames();
    for(int i=0;i<fieldSections.length();i++)
    {

        FieldSection currentFieldSection=fieldSections[i];
        QVectorIterator<NodesField> fieldIterator(currentFieldSection.nodesFields);
        quint64 fieldIteratorPos=0;
        while (fieldIterator.hasNext())
        {
            NodesField currentField=fieldIterator.next();
            if(!nodesSectionsNames.contains(currentField.RefNodes, Qt::CaseInsensitive))
            {
                 noLinkProblem=false;
                 std::cout <<"RefNodes: "<<currentField.RefNodes.toUtf8().constData()<<" from "<< currentFieldSection.name.toUtf8().constData()
                          <<" field element "<< currentField.index<<" not declared."<<std::endl;
                 emit errorMessage("RefNodes: "+currentField.RefNodes+" from section "+ currentFieldSection.name+" field element "+ QString::number(currentField.index)+" not declared.");
                 logStrings::errorLogString.append("RefNodes: "+currentField.RefNodes+" from section "+ currentFieldSection.name+" field element "+ QString::number(currentField.index)+" not declared.\r\n");
                 fieldSections[i].nodesFields.removeAt(fieldIteratorPos);
                 continue;
            }
            fieldIteratorPos=fieldIteratorPos+1;
        }

        if(!currentFieldSection.RefNodes.isEmpty())
        {
            if(!nodesSectionsNames.contains(currentFieldSection.RefNodes,Qt::CaseInsensitive))
            {
                std::cout <<"nodesReference: "<<currentFieldSection.RefNodes.toUtf8().constData()<<" from section "<< currentFieldSection.name.toUtf8().constData()<<" not declared."<<std::endl;
                emit errorMessage("nodesReference: "+currentFieldSection.RefNodes+" from section "+ currentFieldSection.name+" not declared.");
                logStrings::errorLogString.append("nodesReference: "+currentFieldSection.RefNodes+" from section "+ currentFieldSection.name+" not declared.\r\n");
            }
        }
//        std::cout <<"Field name:"<<currentFieldSection.name.toUtf8().constData()<<std::endl;
//        std::cout <<"nodesReference:"<<currentFieldSection.RefNodes.toUtf8().constData()<<std::endl;
    }
    for(int i=0;i<fieldSections.length();i++)
    {
        FieldSection currentFieldSection=fieldSections.at(i);
        QVectorIterator<BoundaryField> boundaryFieldIterator(currentFieldSection.boundaryFields);
        quint64 boundaryFieldIteratorPos=0;
        while (boundaryFieldIterator.hasNext())
        {
            BoundaryField currentBoundaryField=boundaryFieldIterator.next();
            if(!elementSectionsNames.contains(currentBoundaryField.Refelements,Qt::CaseInsensitive))
            {
                 noLinkProblem=false;
                 std::cout <<"Refelements: "<<currentBoundaryField.Refelements.toUtf8().constData()<<" from section "<< currentFieldSection.name.toUtf8().constData()
                          <<" boundary field "<< currentBoundaryField.index<<" not declared."<<std::endl;
                 emit errorMessage("Refelements: "+currentBoundaryField.Refelements+" from section "+ currentFieldSection.name+" boundary field "+QString::number(currentBoundaryField.index)+" not declared.");
                 logStrings::errorLogString.append("Refelements: "+currentBoundaryField.Refelements+" from section "+ currentFieldSection.name+" boundary field "+QString::number(currentBoundaryField.index)+" not declared.\r\n");
                 fieldSections[i].boundaryFields.removeAt(boundaryFieldIteratorPos);
                 continue;
            }
            boundaryFieldIteratorPos=boundaryFieldIteratorPos+1;
        }
    }
    for(int i=0;i<bESpectrumSections.length();i++)
    {

        BESpectrumSection currentBESpectrumSection=bESpectrumSections.at(i);
        QListIterator<BESpectrumItem> spectrumIterator(currentBESpectrumSection.bESpectrumItems);
        quint64 spectrumIteratorPos=0;
        while (spectrumIterator.hasNext())
        {
            BESpectrumItem currentSpectrum=spectrumIterator.next();
            if(!nodesSectionsNames.contains(currentSpectrum.RefNodes,Qt::CaseInsensitive))
            {
                noLinkProblem=false;
                bESpectrumSections[i].bESpectrumItems.removeAt(spectrumIteratorPos);
                continue;
            }
            spectrumIteratorPos=spectrumIteratorPos+1;
        }

        if(!currentBESpectrumSection.RefNodes.isEmpty())
        {
            if(!nodesSectionsNames.contains(currentBESpectrumSection.RefNodes,Qt::CaseInsensitive))
            {
                std::cout <<"nodesReference: "<<currentBESpectrumSection.RefNodes.toUtf8().constData()<<" from "<< currentBESpectrumSection.name.toUtf8().constData()<<" not declared."<<std::endl;
                emit errorMessage("nodesReference: "+currentBESpectrumSection.RefNodes+" from section "+ currentBESpectrumSection.name+" not declared.");
                logStrings::errorLogString.append("nodesReference: "+currentBESpectrumSection.RefNodes+" from section "+ currentBESpectrumSection.name+" not declared.\r\n");
            }
        }
//        std::cout <<"BE_Spectrum Section name:"<<currentBESpectrumSection.name.toUtf8().constData()<<std::endl;
//        std::cout <<"nodesReference:"<<currentBESpectrumSection.RefNodes.toUtf8().constData()<<std::endl;
    }
    return noLinkProblem;
}

void ObservationScriptReader::setupObservationElements(){
    std::cout <<"Setting up the observation elements."<<std::endl;
    QStringList nodesSectionsNames=getNodesSectionsNames();
//    double globalMinMeshFrequency=controlSolverSection.MeshFrequency;
//    QListIterator<ElementSection> elementIterator(elementsSections);
//    while (elementIterator.hasNext())
    for(int i=0;i<fieldSections.length();i++)
    {
        FieldSection currentfieldSection=fieldSections[i];
        if(currentfieldSection.containsNoElements())
        {
            continue;
        }
//        boundaryElements.push_back(BoundaryElements(currentElementSection.name));

        QVectorIterator<MeshField> meshFieldIterator(currentfieldSection.meshFields);
//        std::cout <<"meshfieldlist length "<<currentfieldSection.meshFields.length()<<std::endl;
        while (meshFieldIterator.hasNext())
        {
            MeshField currentMeshField=meshFieldIterator.next();
            int index=meshFileAlias.indexOf(currentMeshField.MeshFileAlias,Qt::CaseInsensitive);
            QString meshFilePath=meshFiles.value(index);
            meshFileReader.readMsh(meshFilePath, currentMeshField.index, currentMeshField.includeAll,currentMeshField.include,currentMeshField.exclude);
            QVector<VectorTriangle> elementsVectorTriangles=meshFileReader.getTriangles();
            QVector<VectorQuadrilateral> elementsVectorQuadrilaterals=meshFileReader.getQuadrilaterals();
//            if(currentMeshField.swapNormals^currentfieldSection.swapNormals)
//            {
//                std::cout<<"swapNormals"<<std::endl;
//                currentfieldSection.swapTriangleNormals(elementsVectorTriangles);
//                currentfieldSection.swapQuadrilateralNormals(elementsVectorQuadrilaterals);
//            }
            fieldSections[i].fieldTriangles.append(elementsVectorTriangles);
            fieldSections[i].fieldQuadrilaterals.append(elementsVectorQuadrilaterals);
        }
        fieldSections[i].meshFields.clear();

        QVectorIterator<NodesField> nodesFieldIterator(currentfieldSection.nodesFields);
        while (nodesFieldIterator.hasNext())
        {
            NodesField currentNodesField=nodesFieldIterator.next();
            QRegularExpression nodesReference("^"+currentNodesField.RefNodes+"$", QRegularExpression::CaseInsensitiveOption);
            int nodesSectionIndex=nodesSectionsNames.indexOf(nodesReference);
            if(nodesSectionIndex!=-1)
            {
                int firstNode=currentNodesField.node1;
                int secondNode=currentNodesField.node2;
                int thirdNode=currentNodesField.node3;
                int fourthNode=currentNodesField.node4;

                bool allNodesDeclared=true;
                if(!nodesSections.at(nodesSectionIndex).nodes.contains(firstNode))
                {
                    std::cout <<"nodessection: "<<nodesSections.at(nodesSectionIndex).name.toUtf8().constData()<<" doesn't contain node "<<firstNode<<std::endl; allNodesDeclared=false;
                    logStrings::errorLogString.append("nodessection: "+nodesSections.at(nodesSectionIndex).name+" doesn't contain node "+QString::number(firstNode)+"\r\n");
                }
                if(!nodesSections.at(nodesSectionIndex).nodes.contains(secondNode))
                {
                    std::cout <<"nodessection: "<<nodesSections.at(nodesSectionIndex).name.toUtf8().constData()<<" doesn't contain node "<<secondNode<<std::endl; allNodesDeclared=false;
                    logStrings::errorLogString.append("nodessection: "+nodesSections.at(nodesSectionIndex).name+" doesn't contain node "+QString::number(secondNode)+"\r\n");
                }
                if(!nodesSections.at(nodesSectionIndex).nodes.contains(thirdNode))
                {
                    std::cout <<"nodessection: "<<nodesSections.at(nodesSectionIndex).name.toUtf8().constData()<<" doesn't contain node "<<thirdNode<<std::endl; allNodesDeclared=false;
                    logStrings::errorLogString.append("nodessection: "+nodesSections.at(nodesSectionIndex).name+" doesn't contain node "+QString::number(thirdNode)+"\r\n");
                }
                if(!nodesSections.at(nodesSectionIndex).nodes.contains(fourthNode))
                {
                    std::cout <<"nodessection: "<<nodesSections.at(nodesSectionIndex).name.toUtf8().constData()<<" doesn't contain node "<<fourthNode<<std::endl; allNodesDeclared=false;
                    logStrings::errorLogString.append("nodessection: "+nodesSections.at(nodesSectionIndex).name+" doesn't contain node "+QString::number(fourthNode)+"\r\n");
                }

                if(allNodesDeclared==true)
                {
                    qint64 fieldIndex=currentNodesField.index;
                    Eigen::Vector3d node1=nodesSections.at(nodesSectionIndex).nodes.value(firstNode).coordinates;
                    Eigen::Vector3d node2=nodesSections.at(nodesSectionIndex).nodes.value(secondNode).coordinates;
                    Eigen::Vector3d node3=nodesSections.at(nodesSectionIndex).nodes.value(thirdNode).coordinates;
                    Eigen::Vector3d node4=nodesSections.at(nodesSectionIndex).nodes.value(fourthNode).coordinates;

    //                boundaryElements.last().elementsVectorTriangles.push_back(VectorTriangle(node1,node2,node3));
    //                if(currentfieldSection.swapNormals)
    //                {
    //                    fieldSections[i].fieldQuadrilaterals.push_back(NodesField(fieldIndex,node4,node3,node2,node1));
    //                }
    //                else
    //                {
                        fieldSections[i].fieldQuadrilaterals.push_back(VectorQuadrilateral(fieldIndex,node1,node2,node3,node4));
    //                }
                }
            }
            else
            {
                std::cout <<"nodesReference: "<<currentNodesField.RefNodes.toUtf8().constData()<<" from "<< currentfieldSection.name.toUtf8().constData()
                          <<" field with index "<< currentNodesField.index<<" not declared at time of field setup."<<std::endl;
                emit errorMessage("nodesReference: "+currentNodesField.RefNodes+" from section "+ currentfieldSection.name+" field with index "+QString::number(currentNodesField.index)+" not declared at time of field setup.");
                logStrings::errorLogString.append("nodesReference: "+currentNodesField.RefNodes+" from section "+ currentfieldSection.name+" field with index "+QString::number(currentNodesField.index)+" not declared at time of field setup.\r\n");
            }
        }
        fieldSections[i].nodesFields.clear();

//        QVectorIterator<BoundaryField> quadrilateralElementIterator(currentfieldSection.quadrilaterals);
//        while (quadrilateralElementIterator.hasNext())
//        {
//            Quadrilateral currentQuadrilateralElement=quadrilateralElementIterator.next();
//            QRegularExpression nodesReference("^"+currentQuadrilateralElement.nodesReference+"$", QRegularExpression::CaseInsensitiveOption);
//            int nodesSectionIndex=nodesSectionsNames.indexOf(nodesReference);
//            if(nodesSectionIndex!=-1)
//            {
//                int firstNode=currentQuadrilateralElement.node1;
//                int secondNode=currentQuadrilateralElement.node2;
//                int thirdNode=currentQuadrilateralElement.node3;
//                int fourthNode=currentQuadrilateralElement.node4;
//                bool allNodesDeclared=true;
//                if(!nodesSections.at(nodesSectionIndex).nodes.contains(firstNode)){std::cout <<"nodessection: "<<nodesSections.at(nodesSectionIndex).name.toUtf8().constData()<<" doesn't contain node "<<firstNode<<std::endl; allNodesDeclared=false;}
//                if(!nodesSections.at(nodesSectionIndex).nodes.contains(secondNode)){std::cout <<"nodessection: "<<nodesSections.at(nodesSectionIndex).name.toUtf8().constData()<<" doesn't contain node "<<secondNode<<std::endl; allNodesDeclared=false;}
//                if(!nodesSections.at(nodesSectionIndex).nodes.contains(thirdNode)){std::cout <<"nodessection: "<<nodesSections.at(nodesSectionIndex).name.toUtf8().constData()<<" doesn't contain node "<<thirdNode<<std::endl; allNodesDeclared=false;}
//                if(!nodesSections.at(nodesSectionIndex).nodes.contains(fourthNode)){std::cout <<"nodessection: "<<nodesSections.at(nodesSectionIndex).name.toUtf8().constData()<<" doesn't contain node "<<fourthNode<<std::endl; allNodesDeclared=false;}

//                if(allNodesDeclared==true)
//                {
//                    qint64 quadrilateralIndex=currentQuadrilateralElement.elementIndex;
//                    Eigen::Vector3d node1=nodesSections.at(nodesSectionIndex).nodes.value(firstNode).coordinates;
//                    Eigen::Vector3d node2=nodesSections.at(nodesSectionIndex).nodes.value(secondNode).coordinates;
//                    Eigen::Vector3d node3=nodesSections.at(nodesSectionIndex).nodes.value(thirdNode).coordinates;
//                    Eigen::Vector3d node4=nodesSections.at(nodesSectionIndex).nodes.value(fourthNode).coordinates;
////                    boundaryElements.last().elementsVectorQuadrilaterals.push_back(VectorQuadrilateral(node1,node2,node3,node4));
//                    if(currentfieldSection.swapNormals)
//                    {
//                        elementsSections[i].elementsVectorQuadrilaterals.push_back(VectorQuadrilateral(quadrilateralIndex,node4,node3,node2,node1));
//                    }
//                    else
//                    {
//                        elementsSections[i].elementsVectorQuadrilaterals.push_back(VectorQuadrilateral(quadrilateralIndex,node1,node2,node3,node4));
//                    }


//                }
//            }
//            else
//            {
//                std::cout <<"nodesReference: "<<currentQuadrilateralElement.nodesReference.toUtf8().constData()<<" from "<< currentfieldSection.name.toUtf8().constData()
//                               <<" of quadrilateral element "<< currentQuadrilateralElement.elementIndex<<" not found!"<<std::endl;
//            }

//        }
//        elementsSections[i].quadrilaterals.clear();

    }
    std::cout <<"Setting up the observation elements (spectrum part)."<<std::endl;

    for(int i=0;i<bESpectrumSections.length();i++)
    {
        BESpectrumSection currentBESpectrumSection=bESpectrumSections.at(i);
        if(currentBESpectrumSection.containsNoElements())
        {
            continue;
        }
        QListIterator<BESpectrumItem> spectrumItemsIterator(currentBESpectrumSection.bESpectrumItems);
        while(spectrumItemsIterator.hasNext())
        {
            BESpectrumItem spectrumItem=spectrumItemsIterator.next();
            QRegularExpression nodesReference("^"+spectrumItem.RefNodes+"$", QRegularExpression::CaseInsensitiveOption);
            int nodesSectionIndex=nodesSectionsNames.indexOf(nodesReference);
            if(nodesSectionIndex!=-1)
            {
                int nodeIndex=spectrumItem.nodeIndex;
                bool nodeDeclared=true;
                if(!nodesSections.at(nodesSectionIndex).nodes.contains(nodeIndex))
                {
                    std::cout <<"nodessection: "<<nodesSections.at(nodesSectionIndex).name.toUtf8().constData()<<" doesn't contain node "<<nodeIndex<<std::endl; nodeDeclared=false;
                    logStrings::errorLogString.append("nodessection: "+nodesSections.at(nodesSectionIndex).name+" doesn't contain node "+QString::number(nodeIndex)+"\r\n");
                }
                if(nodeDeclared==true)
                {
                    qint64 spectrumIndex=spectrumItem.index;
                    Eigen::Vector3d node=nodesSections.at(nodesSectionIndex).nodes.value(nodeIndex).coordinates;
                    bESpectrumSections[i].bESpectrumPoints.push_back(BESpectrumPoint(spectrumIndex,node));
                }
            }
            else
            {
                logStrings::errorLogString.append("nodesReference: "+spectrumItem.RefNodes+" from section "+ currentBESpectrumSection.name+" spectrum item with index "+QString::number(spectrumItem.index)+" not declared at time of field setup.\r\n");
            }
        }
    }
    if(c<global::tiny)
    {
        c=343.2;
    }
    if(EdgeLength<0.00000001)
    {
        EdgeLength=1;
    }
    for(int i=0;i<fieldSections.length();i++)
    {
        FieldSection currentFieldSection=fieldSections.at(i);
        if(currentFieldSection.MeshFrequency==0 && currentFieldSection.EdgeLength==0)
        {
            fieldSections[i].EdgeLength=EdgeLength;
        }
        else if(currentFieldSection.MeshFrequency!=0)
        {
            fieldSections[i].EdgeLength=(c/currentFieldSection.MeshFrequency)/2;
        }
        fieldSections[i].fieldTriangles.append(MeshFunctions::refineQuadrilaterals(fieldSections.at(i).fieldQuadrilaterals,fieldSections[i].EdgeLength));
        fieldSections[i].fieldTriangles=MeshFunctions::refineTriangles(fieldSections.at(i).fieldTriangles,fieldSections[i].EdgeLength);
        fieldSections[i].fieldQuadrilaterals.clear();
    }
}

QVector<ObservationPoint> ObservationScriptReader::getObservationPoints()
{
    QVector<ObservationPoint> observationPoints;
    for(int i=0;i<bESpectrumSections.size();i++)
    {
        BESpectrumSection currentSpectrumSection=bESpectrumSections.at(i);
        for(int j=0;j<currentSpectrumSection.bESpectrumPoints.size();j++)
        {
           BESpectrumPoint currentSpectrumPoint=currentSpectrumSection.bESpectrumPoints.at(j);
            observationPoints.push_back(ObservationPoint(currentSpectrumSection.name,currentSpectrumPoint.index,currentSpectrumPoint.coordinates));
        }
    }
    return observationPoints;
}

QVector<ObservationField> ObservationScriptReader::getObservationFields()
{
    QVector<ObservationField> observationFields;
    for(int i=0;i<fieldSections.size();i++)
    {
        FieldSection currentFieldSection=fieldSections.at(i);
        observationFields.push_back(ObservationField(currentFieldSection.name,currentFieldSection.fieldTriangles));
    }
    return observationFields;
}
