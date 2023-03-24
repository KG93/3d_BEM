#include "mshreader.h"

//mshReader::mshReader(QString fileName, QObject *parent) : QObject(parent)
//{
//    mshFilename = fileName;
//}

mshReader::mshReader(QObject* parent) : QObject(parent){}

void mshReader::readMsh(const QString &filename, qint64 elementIndex, bool includeAll, const QStringList &include, const QStringList &exclude)
{
   this->elementIndex = elementIndex;
    mshFilename = filename;
    this->includeAll = includeAll;
    includeList = include;
    excludeList = exclude;
    physicalGroupNumbersFilterList.clear();
    geometricalGroupNumbersFilterList.clear();
    readSections();
}

void mshReader::sortOutIncludesAndExcludes()
{

    quint64 maxValue = std::numeric_limits<quint64>::max();
    QVector<QPair<quint64,quint64>> geometricalGroupNumbersToBeIncluded;
    QVector<QPair<quint64,quint64>> geometricalGroupNumbersToBeExcluded;
    QVector<quint64> physicalGroupNumbersToBeIncluded;
    QVector<quint64> physicalGroupNumbersToBeExcluded;
    int index = includeList.indexOf("all", Qt::CaseInsensitive);

    if(index>=0)
    {
//        geometricalGroupNumbersToBeIncluded.push_back(qMakePair(1,maxValue));
        includeList.removeAt(index);
        includeAll = true;
    }
    if(includeList.isEmpty())
    {
        includeAll = true;
    }
    else
    {
        includeAll = false;
    }

    QRegularExpressionMatch match;
    for(int i=0; i<includeList.length(); i++)
    {
        quint64 start;
        quint64 end;
        QRegularExpression identifier("^\\s*(\\d*)\\.\\.\\.?(\\d*)\\s*$", QRegularExpression::CaseInsensitiveOption);

        if(includeList.value(i).contains(identifier,&match))
        {
            QString rangeStart = match.captured(1);
            QString rangeEnd = match.captured(2);
            if(rangeStart.isEmpty() && rangeEnd.isEmpty())
            {
                continue;
            }
            if(rangeStart.isEmpty())
            {
                start = 1;
            }
            else{
                start = rangeStart.toUInt();
            }
            if(rangeEnd.isEmpty())
            {
                end = maxValue;
            }
            else{
                end = rangeEnd.toUInt();

            }
            if(start <= end){
                geometricalGroupNumbersToBeIncluded.push_back(qMakePair(start, end));
            }
            else{
                geometricalGroupNumbersToBeIncluded.push_back(qMakePair(end, start));

            }
            continue;
        }
        identifier = QRegularExpression ("^\\s*(\\d+)\\s*$", QRegularExpression::CaseInsensitiveOption);
        if(includeList.value(i).contains(identifier, &match))
        {
//            std::cout<<"Number found: "<<std::endl;

            QString value = match.captured(1); //first numbers
            if(!value.isEmpty())
            {
                start = value.toUInt();
//                std::cout<<"includeNumber: "<<start<<std::endl;

                geometricalGroupNumbersToBeIncluded.push_back(qMakePair(start, start));
            }
            continue;
        }

        identifier = QRegularExpression ("^\\s*"+includeList.value(i)+"\\s*$", QRegularExpression::CaseInsensitiveOption);
        int indexOfPhysicalGroup = physicalGroupsNames.indexOf(identifier);
        if(indexOfPhysicalGroup >= 0 && indexOfPhysicalGroup < physicalGroupsNumbers.length())
        {
            int physGroupNumber = physicalGroupsNumbers.value(indexOfPhysicalGroup);
            physicalGroupNumbersToBeIncluded.push_back(physGroupNumber);
            continue;
        }


    }
//    std::cout <<"geometricalGroupNumbersToBeIncluded length : "<< geometricalGroupNumbersToBeIncluded.length()<<std::endl;

    for(int i=0; i<excludeList.length(); i++)
    {
        quint64 start;
        quint64 end;
        QRegularExpression identifier("^\\s*(\\d*)\\.\\.\\.?(\\d*)\\s*$", QRegularExpression::CaseInsensitiveOption);

        if(excludeList.value(i).contains(identifier,&match))
        {
            QString rangeStart = match.captured(1);
            QString rangeEnd = match.captured(2);
            if(rangeStart.isEmpty() && rangeEnd.isEmpty())
            {
                continue;
            }
            if(rangeStart.isEmpty())
            {
                start = 1;
            }
            else{
                start = rangeStart.toUInt();
            }
            if(rangeEnd.isEmpty())
            {
                end = std::numeric_limits<quint64>::max();
            }
            else{
                end = rangeEnd.toUInt();

            }
            if(start <= end){
                geometricalGroupNumbersToBeExcluded.push_back(qMakePair(start,end));
            }
            else{
                geometricalGroupNumbersToBeExcluded.push_back(qMakePair(end,start));

            }
            continue;
        }

        identifier = QRegularExpression ("^\\s*(\\d+)\\s*$", QRegularExpression::CaseInsensitiveOption);
        if(excludeList.value(i).contains(identifier,&match))
        {
//            std::cout<<"Number found: "<<std::endl;

            QString value = match.captured(1); //first numbers
            if(!value.isEmpty())
            {
                start = value.toUInt();
//                std::cout<<"excludeNumber: "<<start<<std::endl;
                geometricalGroupNumbersToBeExcluded.push_back(qMakePair(start,start));
            }
            continue;
        }

        identifier = QRegularExpression ("^\\s*"+excludeList.value(i)+"\\s*$", QRegularExpression::CaseInsensitiveOption);
        int indexOfPhysicalGroup = physicalGroupsNames.indexOf(identifier);
        int indexOfPhysicalGroupIfAlsoInInclude = includeList.indexOf(identifier); //if it's also on the include list, it will be ignored
        if(includeAll && indexOfPhysicalGroup >= 0 && indexOfPhysicalGroup < physicalGroupsNumbers.length() && indexOfPhysicalGroupIfAlsoInInclude == -1)
        {
            int physGroupNumber = physicalGroupsNumbers.value(indexOfPhysicalGroup);
            physicalGroupNumbersToBeExcluded.push_back(physGroupNumber);
            continue;
        }

    }
//    std::cout <<"geometricalGroupNumbersToBeExcluded length : "<< geometricalGroupNumbersToBeExcluded.length()<<std::endl;

//    geometricalGroupNumbersToBeIncluded = global::diffIntervals(geometricalGroupNumbersToBeIncluded, geometricalGroupNumbersToBeExcluded);
    if(includeAll)
    {
        physicalGroupNumbersFilterList = physicalGroupNumbersToBeExcluded;
        global::mergeIntervals(geometricalGroupNumbersToBeExcluded);
        geometricalGroupNumbersFilterList = geometricalGroupNumbersToBeExcluded;

    }
    else
    {
        physicalGroupNumbersFilterList = physicalGroupNumbersToBeIncluded;
        geometricalGroupNumbersFilterList = global::diffIntervals(geometricalGroupNumbersToBeIncluded, geometricalGroupNumbersToBeExcluded);

    }
    std::cout << "geometricalGroupNumbersFilterList length : " << geometricalGroupNumbersFilterList.length() << std::endl;
    std::cout << "physicalGroupNumbersFilterList length : " << physicalGroupNumbersFilterList.length() << std::endl;
}

void  mshReader::readSections()
{
    QFile meshFile(mshFilename);
    if(!meshFile.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        std::cout << "Incorrect file : " << mshFilename.toUtf8().constData() << std::endl;
        return;
        //    emit nachricht("Konnte die Datei nicht einlesen: "+mshFilename);
    }

    QTextStream in(&meshFile);
    QString line = in.readLine();

    if(line.startsWith("$MeshFormat",Qt::CaseInsensitive))
    {
        std::cout <<"correct .msh file : "<< mshFilename.toUtf8().constData()<<std::endl;
    }
    while(!line.startsWith("$EndMeshFormat",Qt::CaseInsensitive))
    {
        line = in.readLine();
        std::cout << line.toUtf8().constData() << std::endl;
    }

    while(!line.isNull())
    {
        if(line.startsWith("$Nodes",Qt::CaseInsensitive)) // Nodes section
        {
            QString tmp = in.readLine();
            numberOfNodes = tmp.toUInt(&validLine);

//            nodelist.clear();
            nodes.clear();
//            nodelist.squeeze();
//            nodelist.resize(numberOfNodes);
            nodes.reserve(numberOfNodes);
            quint64 nodenumber;
            double x,y,z;
            for(quint64 i=0; i<numberOfNodes; i++ )
            {
//                QString node = in.readLine();
                in >> nodenumber;
//                std::cout <<"nodenumber "<< nodenumber<<std::endl;

                in >> x;
                in >> y;
                in >> z;

//                std::array<double,3> array= {x,y,z};
//                nodelist.at(i) = array;
                nodes.insert(nodenumber,Node({x,y,z}));
            }
            line = in.readLine();// skip to next line after last >> operator run
            line = in.readLine();

            if(line.startsWith("$EndNodes",Qt::CaseInsensitive))
            {
                std::cout << "All "<< nodenumber << " nodes succesfully read." << std::endl;
            }
            else
            {
                std::cout <<"Error while reading the nodes. Number of nodes:  "<< nodenumber<<". Last read line: "<<line.toUtf8().constData() <<std::endl<<"next line"<<std::endl;
            }
        }

        if(line.startsWith("$Elements",Qt::CaseInsensitive)) // Elements section
        {
            sortOutIncludesAndExcludes();
            QString tmp = in.readLine();
            numberOfElements = tmp.toUInt(&validLine);
            triangles.clear();
            quadrilaterals.clear();
            quint64 elementNumber;
            quint64 elementType;
            quint64 numberOfTags;
            quint64 numberOfPhysicalGroup;
            quint64 numberOfGeometricalGroup;
            quint64 node1;
            quint64 node2;
            quint64 node3;
            quint64 node4;

            for(quint64 i=0; i<numberOfElements; i++ )
            {
                line = in.readLine(); //finish line
                QStringList lineAsList = line.split(' ', Qt::SkipEmptyParts);
                elementNumber = lineAsList.value(0).toUInt();
                elementType = lineAsList.value(1).toUInt();
                numberOfTags= lineAsList.value(2).toUInt();

                if(numberOfTags == 0)
                {
                    numberOfPhysicalGroup = 0;
                    numberOfGeometricalGroup = 0;//not sure
                }
                else if(numberOfTags == 1)
                {
                    numberOfPhysicalGroup = lineAsList.value(3).toUInt();
                    numberOfGeometricalGroup = 0;//not sure
                }
                else
                {
                    numberOfPhysicalGroup = lineAsList.value(3).toUInt();
                    numberOfGeometricalGroup = lineAsList.value(4).toUInt();
                }

                if(! elementIsIncluded(numberOfPhysicalGroup, numberOfGeometricalGroup)) //check wether element is filtered
                {
                    continue;
                }

//                std::cout<<"Element "<<elementNumber <<" with geometrical number "<<numberOfGeometricalGroup<<" and physical Number "<<numberOfPhysicalGroup<<" will be included."<<std::endl;
                if(elementType == 2) //Triangle
                {
                    node1 = lineAsList.value(2+numberOfTags+1).toUInt();
                    node2 = lineAsList.value(2+numberOfTags+2).toUInt();
                    node3 = lineAsList.value(2+numberOfTags+3).toUInt();
                    if(nodes.contains(node1) && nodes.contains(node2) && nodes.contains(node3))
                    {
                        triangles.push_back(VectorTriangle(elementIndex, nodes.value(node1).coordinates, nodes.value(node2).coordinates, nodes.value(node3).coordinates));
                    }
                    else
                    {
                        std::cout << "Element " << elementNumber << "from file " << mshFilename.toUtf8().constData() << " contains unspecified nodes." << std::endl;
                    }
                }
                if(elementType == 3) //Quadrilateral
                {
                    node1 = lineAsList.value(2+numberOfTags+1).toUInt();
                    node2 = lineAsList.value(2+numberOfTags+2).toUInt();
                    node3 = lineAsList.value(2+numberOfTags+3).toUInt();
                    node4 = lineAsList.value(2+numberOfTags+4).toUInt();
                    if(nodes.contains(node1) && nodes.contains(node2) && nodes.contains(node3) && nodes.contains(node4))
                    {
                        quadrilaterals.push_back(VectorQuadrilateral(elementIndex, nodes.value(node1).coordinates, nodes.value(node2).coordinates, nodes.value(node3).coordinates, nodes.value(node4).coordinates));
                    }
                    else
                    {
                        std::cout << "Element " << elementNumber << "from file " << mshFilename.toUtf8().constData() << " contains unspecified nodes." << std::endl;
                    }
                }
////                line=in.readLine(); //finish line
            }

//            line = in.readLine();
            line = in.readLine();
            if(line.startsWith("$EndElements", Qt::CaseInsensitive))
            {
                std::cout << "All " << numberOfElements << " elements succesfully read." << std::endl;
            }
            else
            {
                std::cout << "Error while reading the elements. Number of elements:  " << numberOfElements << ". Last read line: " << line.toUtf8().constData() << std::endl << "next line" << std::endl;

            }
        }

        if(line.startsWith("$PhysicalNames", Qt::CaseInsensitive)) // Physical names section
        {
            line = in.readLine();
            numberOfPhysicalNames = line.toUInt(&validLine);
            physicalGroupsNumbers.clear();
            physicalGroupsNames.clear();
            physicalGroupsNumbers.reserve(numberOfPhysicalNames);
            quint64 physicalGroupnumber;
            QString physicalGroupName;
            for(quint64 i=0; i<numberOfPhysicalNames; i++ )
            {
                line = in.readLine();
                QStringList stringList = global::stringToStringListQuotesIntact(line);
                physicalGroupnumber = stringList.value(1).toUInt();
                physicalGroupName = stringList.value(2);
                global::removeQuotes(physicalGroupName);
                physicalGroupsNumbers.push_back(physicalGroupnumber);
                physicalGroupsNames.push_back(physicalGroupName);
            }
            line = in.readLine();

            if(line.startsWith("$EndPhysicalNames", Qt::CaseInsensitive))
            {                
                std::cout << "All "<< numberOfPhysicalNames<<" physical names succesfully read." << std::endl;
            }
            else
            {
                std::cout << "Error while reading the physical names. Number of physical names:  " << numberOfPhysicalNames << ". Last read line: " << line.toUtf8().constData() << std::endl << "next line" << std::endl;

            }
            global::printQStringList(physicalGroupsNames);
            global::printVector(physicalGroupsNumbers);
            global::printQStringList(includeList);
            std::cout << "physicalGroupNumbersFilterList ";
            global::printVector(physicalGroupNumbersFilterList);
            std::cout << "geometricalGroupNumbersFilterList ";
            global::printIntervalList(geometricalGroupNumbersFilterList);
        }
        line = in.readLine();
    }
}

bool mshReader::elementIsIncluded(quint64 numberOfPhysicalGroup, quint64 numberOfGeometricalGroup)
{
    if(includeAll) //filter lists are used as exclude lists
    {
        if(physicalGroupNumbersFilterList.contains(numberOfPhysicalGroup))//element is on exclude physical group list
        {
            return false;
        }

        if(global::valueIsInInterval(numberOfGeometricalGroup,geometricalGroupNumbersFilterList))//element is on exclude geometrical group list
        {
            return false;
        }
        return true; //Element is not on exclude filter

    }
    else //filter lists are used as include lists
    {
        if(physicalGroupNumbersFilterList.contains(numberOfPhysicalGroup))//element is on include physical group list
        {
            return true;
        }
        if(global::valueIsInInterval(numberOfGeometricalGroup,geometricalGroupNumbersFilterList))
        {
//            std::cout<<"numberOfGeometricalGroup is in interval "<<numberOfGeometricalGroup<<std::endl;
            return true;
        }
        return false; //Element is not on include filter
    }
//    return false;
}
