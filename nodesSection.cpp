#include "nodesSection.h"

NodesSection::NodesSection()
{

}

bool NodesSection::handleScriptLine(const QString &currentLine)
{
    validLine = false;
    if(!readingNodesList)
    {
        validLine = handleParameterSectionLine(currentLine);
    }
    else
    {
        validLine = handleNodesListLine(currentLine);
    }
    return validLine;
}

bool NodesSection::handleNodesListLine(const QString &currentLine)
{

    QRegularExpression separator( "[ \\t,]+"); //identify any number of white spaces and tabs intermixed or commas
    QStringList list = currentLine.split(separator, Qt::SkipEmptyParts);
    if(list.length() == 3) //only dim3d is implemented, therefore 4 entries are required
    {
        validLine = false;
        return validLine;
    }
    if(list.length() == 4)  //syntax:     <Index>   <Coordinate1>   <Coordinate2>   <Coordinate3>
    {
        bool isInteger;
        int nodeName = list.first().toInt(&isInteger);
        if(isInteger) // current line starts with an element identifier
        {
            readingNodesList=true; //programm reads second part of Nodes section (declaration of list of nodes)
        }
        else{
            validLine=false; // current line doesn't start with a node identifier
            return validLine;

        }

        bool xisDouble;
        bool yisDouble;
        bool zisDouble;
//        double xCoordinate=list.at(1).toDouble(&xisDouble);
        double xCoordinate = global::stringWithUnitsToDouble(list.at(1),QString("m"),xisDouble);
//        double yCoordinate=list.at(2).toDouble(&yisDouble);
        double yCoordinate = global::stringWithUnitsToDouble(list.at(2),QString("m"),yisDouble);
//        double zCoordinate=list.at(3).toDouble(&zisDouble);
        double zCoordinate = global::stringWithUnitsToDouble(list.at(3),QString("m"),zisDouble);

        bool validCoordinateEntries = xisDouble && yisDouble && zisDouble;
        validLine = validCoordinateEntries;

        if(validCoordinateEntries)
        {
//            nodeNames.push_back(nodeName);
//            nodeCoordinates.push_back({xCoordinate, yCoordinate, zCoordinate});
//            nodes.push_back(Node(nodeName,{xCoordinate, yCoordinate, zCoordinate}));
            nodes.insert(nodeName,Node({xCoordinate, yCoordinate, zCoordinate}));
            return validLine;
        }
        else
        {
            return validLine;

        }


    }
    return validLine;
}

bool NodesSection::handleParameterSectionLine(const QString &currentLine)
{
    validLine = false;
    QRegularExpressionMatch match;
    QString value;
    QRegularExpression identifier;
    identifier = global::regExWithOneValue("Scale");
    if(currentLine.contains(identifier,&match))
    {
        value = match.captured(1);
        global::removeQuotes(value);
        bool alreadyCaptured = false;
        QRegularExpression separator( "[ \\t,]+");
        QStringList list = value.split(separator, Qt::SkipEmptyParts);
        if(list.length() == 1)
        {
            QString captured1 = list.at(0);
            double val = global::stringWithUnitsToDouble(captured1,QString("m"),validLine);
            Scale = {val,val,val};
            alreadyCaptured = true;
        }
        if(list.length()==2)
        {
            bool valid1,valid2;
            QString captured1 = list.at(0);
            QString captured2 = list.at(1);
            double val1 = global::stringWithUnitsToDouble(captured1,QString("m"),valid1);
            double val2 = global::stringWithUnitsToDouble(captured2,QString("m"),valid2);
            Scale = {val1,val2,1};
            validLine = (valid1 && valid2);
            alreadyCaptured = true;
        }
        if(list.length() == 3)
        {
            bool valid1,valid2,valid3;
            QString captured1 = list.at(0);
            QString captured2 = list.at(1);
            QString captured3 = list.at(2);
            double val1 = global::stringWithUnitsToDouble(captured1,QString("m"),valid1);
            double val2 = global::stringWithUnitsToDouble(captured2,QString("m"),valid2);
            double val3 = global::stringWithUnitsToDouble(captured3,QString("m"),valid3);
            Scale = {val1,val2,val3};
            validLine = (valid1 && valid2 && valid3);
            alreadyCaptured = true;
        }
        if(!alreadyCaptured)
        {
            validLine = false;
        }
        if(Scale(0) == 0) {Scale(0) = 1;}
        if(Scale(1) == 0) {Scale(1) = 1;}
        if(Scale(2) == 0) {Scale(2) = 1;}
//        std::cout<<"scale: "<<Scale(0)<<", "<<Scale(1)<<", "<<Scale(2)<<std::endl;
        return validLine;
    }

    identifier = global::regExWithOneValue("Shift");
    if(currentLine.contains(identifier,&match))
    {
        value = match.captured(1);
        global::removeQuotes(value);
        bool alreadyCaptured = false;
        QRegularExpression separator( "[ \\t,]+");
        QStringList list = value.split(separator, Qt::SkipEmptyParts);
        if(list.length() == 1)
        {
            QString captured1 = list.at(0);
            double val = global::stringWithUnitsToDouble(captured1,QString("m"),validLine);
            Shift = {val,0,0};
            alreadyCaptured = true;
        }
        if(list.length() == 2)
        {
            bool valid1,valid2;
            QString captured1 = list.at(0);
            QString captured2 = list.at(1);
            double val1 = global::stringWithUnitsToDouble(captured1,QString("m"),valid1);
            double val2 = global::stringWithUnitsToDouble(captured2,QString("m"),valid2);
            Shift = {val1,val2,0};
            validLine = (valid1 && valid2);
            alreadyCaptured = true;
        }
        if(list.length() == 3)
        {
            bool valid1,valid2,valid3;
            QString captured1 = list.at(0);
            QString captured2 = list.at(1);
            QString captured3 = list.at(2);
            double val1 = global::stringWithUnitsToDouble(captured1,QString("m"),valid1);
            double val2 = global::stringWithUnitsToDouble(captured2,QString("m"),valid2);
            double val3 = global::stringWithUnitsToDouble(captured3,QString("m"),valid3);
            Shift = {val1,val2,val3};
            validLine = (valid1 && valid2 && valid3);
            alreadyCaptured = true;
        }
        if(!alreadyCaptured)
        {
            validLine = false;
        }
        else if(validLine)
        {
            transformationMatrix.translate(Shift(0),Shift(1),Shift(2));
        }
//        std::cout<<"Shift: "<<Shift(0)<<", "<<Shift(1)<<", "<<Shift(2)<<std::endl;
        return validLine;
    }

    identifier = global::regExWithOneValue("Rotate");
    if(currentLine.contains(identifier,&match))
    {
        value = match.captured(1);
        global::removeQuotes(value);
        bool alreadyCaptured = false;
        QRegularExpression separator( "[ \\t,]+");
        QStringList list = value.split(separator, Qt::SkipEmptyParts);
        if(list.length() == 1)
        {
            QString captured1=list.at(0);
//            double val =global::stringWithUnitsToDouble(captured1,QString("m"),validLine);
            double val = captured1.toDouble(&validLine);
            Rotate = {val,0,0};
            alreadyCaptured=true;
        }
        if(list.length() == 2)
        {
            bool valid1,valid2;
            QString captured1 = list.at(0);
            QString captured2 = list.at(1);
            double val1 = captured1.toDouble(&valid1);
            double val2 = captured2.toDouble(&valid2);
            Rotate = {val1,val2,0};
            validLine = (valid1 && valid2);
            alreadyCaptured = true;
        }
        if(list.length() == 3)
        {
            bool valid1,valid2,valid3;
            QString captured1 = list.at(0);
            QString captured2 = list.at(1);
            QString captured3 = list.at(2);
            double val1 = captured1.toDouble(&valid1);
            double val2 = captured2.toDouble(&valid2);
            double val3 = captured3.toDouble(&valid3);
            Rotate = {val1,val2,val3};
            validLine = (valid1 && valid2 && valid3);
            alreadyCaptured = true;
        }
        if(!alreadyCaptured)
        {
            validLine = false;
        }
        else if(validLine)
        {
            transformationMatrix.rotate(Rotate(0),1,0,0);
            transformationMatrix.rotate(Rotate(1),0,1,0);
            transformationMatrix.rotate(Rotate(2),0,0,1);
        }
//        std::cout<<"Rotate: "<<Rotate(0)<<", "<<Rotate(1)<<", "<<Rotate(2)<<std::endl;
        return validLine;
    }

    QRegularExpression separator( "[ \\t,]+"); //identify any number of white spaces and tabs intermixed or commas
    QStringList list = currentLine.split(separator, Qt::SkipEmptyParts);
    if(!list.isEmpty())
    {
        bool isNumber;
        list.first().toInt(&isNumber);
        if(isNumber)
        {
            readingNodesList = true;
            validLine = handleNodesListLine(currentLine);
            return validLine;
        }
    }
    QString error = QString("Line in Nodes_Section section not recognized! "+ currentLine);
    logStrings::errorLogString.append(error+"\r\n");
    return validLine;
}

bool NodesSection::checkIdentifiersAllUnique()
{
    nonUniqueIdentifiers.clear();
    bool allIdentifiersUnique = true;
    nodeNames = nodes.keys();
    sortNodeNames();
    if(nodeNames.length()<2)
    {
        return true;
    }
    for(int i=0; i<=nodeNames.length()-2; i++)
    {
        if(nodeNames.at(i) == nodeNames.at(i+1))
        {
            nonUniqueIdentifiers.push_back(nodeNames.at(i));
            allIdentifiersUnique = false;
        }

    }
    return allIdentifiersUnique;
}

void NodesSection::scaleNodes()
{
    double xScale = Scale(0);
    double yScale = Scale(1);
    double zScale = Scale(2);

    bool reallyScale = true;
    if(xScale == 1 && yScale == 1 && zScale == 1)
    {
        reallyScale = false;
    }
    QMultiHash<quint64,Node>::iterator nodesIterator;
    if(reallyScale)
    {
        for (nodesIterator = nodes.begin(); nodesIterator != nodes.end(); ++nodesIterator)
        {
            double xCoordinate = nodesIterator.value().coordinates(0) * xScale;
            double yCoordinate = nodesIterator.value().coordinates(1) * yScale;
            double zCoordinate = nodesIterator.value().coordinates(2) * zScale;
            nodesIterator.value().coordinates = {xCoordinate, yCoordinate, zCoordinate};
        }
    }
    if(transformationMatrix.isIdentity() == false)
    {
        Eigen::Matrix4d transformationMatrixEigen = MeshFunctions::qtToEigenMatrix(transformationMatrix);
        Eigen::Vector3d node;
        Eigen::Vector4d homogNode;
        for (nodesIterator = nodes.begin(); nodesIterator != nodes.end(); ++nodesIterator)
        {
            node = nodesIterator.value().coordinates;
            homogNode = MeshFunctions::homogVec(node);
            homogNode = transformationMatrixEigen*homogNode;
            nodesIterator.value().coordinates = MeshFunctions::eigenVec(homogNode);
        }
    }
}
