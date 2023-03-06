#include "elementSection.h"

ElementSection::ElementSection()
{
}

bool ElementSection::handleScriptLine(const QString& currentLine)
{

    validLine=false;
    if(!readingElementList) //an elements section contains two parts: first the parameters and then the list of elements
    {
        validLine=handleParameterSectionLine(currentLine);
    }
    else {
        validLine=handleElementListLine(currentLine);
    }
    return validLine;
}

bool ElementSection::handleParameterSectionLine(const QString& currentLine)
{
    validLine=false;
    QRegularExpressionMatch match;
    QString value;
    QRegularExpression identifier;

    identifier=global::regExWithOneValue("Subdomain");
    if(currentLine.contains(identifier,&match))
    {
        value=match.captured(1);
        global::removeQuotesAndWhitespace(value);
        Subdomain=value.toUInt(&validLine);
//        if(Subdomain<0){validLine=false;}
        return validLine;
    }

    identifier=global::regExWithOneValue("ElType");
    if(currentLine.contains(identifier,&match))
    {
        value=match.captured(1);
        global::removeQuotesAndWhitespace(value);
        if(0==QString::compare("interior",value, Qt::CaseInsensitive))
        {
              validLine=true;
              ElType=interior;
              return validLine;
        }
        if(0==QString::compare("exterior",value, Qt::CaseInsensitive))
        {
            validLine=true;
              ElType=exterior;
              return validLine;
        }
        validLine=false;
        return validLine;

    }

    identifier=global::regExWithOneValue("MeshFrequency");
    if(currentLine.contains(identifier,&match))
    {
        value=match.captured(1);
        global::removeQuotesAndWhitespace(value);
        MeshFrequency=global::stringWithUnitsToDouble(value,QString("Hz"),validLine);
        if(MeshFrequency<=0){MeshFrequency=0; validLine=false;}
        return validLine;
    }

    identifier=global::regExWithOneValue("EdgeLength");
    if(currentLine.contains(identifier,&match))
    {
        value=match.captured(1);
        global::removeQuotesAndWhitespace(value);
//        EdgeLength=value.toDouble(&validLine);
        EdgeLength=global::stringWithUnitsToDouble(value,QString("m"),validLine);
//        std::cout<<"Edgelength: "<<EdgeLength<<std::endl;
        if(EdgeLength<=0){EdgeLength=0; validLine=false;}
        return validLine;
    }

    identifier=regExWithOptionalValue("swapNormals");
    if(currentLine.contains(identifier,&match))
    {
//        std::cout<<"normals Swapped!!!"<<std::endl;
          value=match.captured(2);
          std::cout<<value.toUtf8().constData()<<std::endl;
          global::removeQuotesAndWhitespace(value);
          if(!value.isEmpty())
          {
              if(value.compare("true",Qt::CaseInsensitive)==0||value.toInt()==1)
              {
                  swapNormals=true;
              }
          }
          else
          {
              swapNormals=true;
          }
          validLine=true;
          return validLine;
    }

    identifier=global::regExWithOneValue("Scale");
    if(currentLine.contains(identifier,&match))
    {
        value=match.captured(1);
        global::removeQuotes(value);
        bool alreadyCaptured=false;
        QRegularExpression separator( "[ \\t,]+");
        QStringList list =value.split(separator, Qt::SkipEmptyParts);
        if(list.length()==1)
        {
            QString captured1=list.at(0);
            double val =global::stringWithUnitsToDouble(captured1,QString("m"),validLine);
            Scale={val,val,val};
            alreadyCaptured=true;
        }
        if(list.length()==2)
        {
            bool valid1,valid2;
            QString captured1=list.at(0);
            QString captured2=list.at(1);
            double val1 =global::stringWithUnitsToDouble(captured1,QString("m"),valid1);
            double val2 =global::stringWithUnitsToDouble(captured2,QString("m"),valid2);
            Scale={val1,val2,1};
            validLine=(valid1 && valid2);
            alreadyCaptured=true;
        }
        if(list.length()==3)
        {
            bool valid1,valid2,valid3;
            QString captured1=list.at(0);
            QString captured2=list.at(1);
            QString captured3=list.at(2);
            double val1 =global::stringWithUnitsToDouble(captured1,QString("m"),valid1);
            double val2 =global::stringWithUnitsToDouble(captured2,QString("m"),valid2);
            double val3 =global::stringWithUnitsToDouble(captured3,QString("m"),valid3);
            Scale={val1,val2,val3};
            validLine=(valid1 && valid2 && valid3);
            alreadyCaptured=true;
        }
        if(!alreadyCaptured)
        {
            validLine=false;
        }
        if(Scale(0)==0){Scale(0)=1;}
        if(Scale(1)==0){Scale(1)=1;}
        if(Scale(2)==0){Scale(2)=1;}
        std::cout<<"scale: "<<Scale(0)<<", "<<Scale(1)<<", "<<Scale(2)<<std::endl;
        return validLine;
    }

    identifier=global::regExWithOneValue("Shift");
    if(currentLine.contains(identifier,&match))
    {
        value=match.captured(1);
        global::removeQuotes(value);
        bool alreadyCaptured=false;
        QRegularExpression separator( "[ \\t,]+");
        QStringList list =value.split(separator, Qt::SkipEmptyParts);
        if(list.length()==1)
        {
            QString captured1=list.at(0);
            double val =global::stringWithUnitsToDouble(captured1,QString("m"),validLine);
            Shift={val,0,0};
            alreadyCaptured=true;
        }
        if(list.length()==2)
        {
            bool valid1,valid2;
            QString captured1=list.at(0);
            QString captured2=list.at(1);
            double val1 =global::stringWithUnitsToDouble(captured1,QString("m"),valid1);
            double val2 =global::stringWithUnitsToDouble(captured2,QString("m"),valid2);
            Shift={val1,val2,0};
            validLine=(valid1 && valid2);
            alreadyCaptured=true;
        }
        if(list.length()==3)
        {
            bool valid1,valid2,valid3;
            QString captured1=list.at(0);
            QString captured2=list.at(1);
            QString captured3=list.at(2);
            double val1 =global::stringWithUnitsToDouble(captured1,QString("m"),valid1);
            double val2 =global::stringWithUnitsToDouble(captured2,QString("m"),valid2);
            double val3 =global::stringWithUnitsToDouble(captured3,QString("m"),valid3);
            Shift={val1,val2,val3};
            validLine=(valid1 && valid2 && valid3);
            alreadyCaptured=true;
        }
        if(!alreadyCaptured)
        {
            validLine=false;
        }
        else if(validLine)
        {
            transformationMatrix.translate(Shift(0),Shift(1),Shift(2));
        }
        std::cout<<"Shift: "<<Shift(0)<<", "<<Shift(1)<<", "<<Shift(2)<<std::endl;
        return validLine;
    }

    identifier=global::regExWithOneValue("Rotate");
    if(currentLine.contains(identifier,&match))
    {
        value=match.captured(1);
        global::removeQuotes(value);
        bool alreadyCaptured=false;
        QRegularExpression separator( "[ \\t,]+");
        QStringList list =value.split(separator, Qt::SkipEmptyParts);
        if(list.length()==1)
        {
            QString captured1=list.at(0);
//            double val =global::stringWithUnitsToDouble(captured1,QString("m"),validLine);
            double val=captured1.toDouble(&validLine);
            Rotate={val,0,0};
            alreadyCaptured=true;
        }
        if(list.length()==2)
        {
            bool valid1,valid2;
            QString captured1=list.at(0);
            QString captured2=list.at(1);
            double val1=captured1.toDouble(&valid1);
            double val2=captured2.toDouble(&valid2);
            Rotate={val1,val2,0};
            validLine=(valid1 && valid2);
            alreadyCaptured=true;
        }
        if(list.length()==3)
        {
            bool valid1,valid2,valid3;
            QString captured1=list.at(0);
            QString captured2=list.at(1);
            QString captured3=list.at(2);
            double val1=captured1.toDouble(&valid1);
            double val2=captured2.toDouble(&valid2);
            double val3=captured3.toDouble(&valid3);
            Rotate={val1,val2,val3};
            validLine=(valid1 && valid2 && valid3);
            alreadyCaptured=true;
        }
        if(!alreadyCaptured)
        {
            validLine=false;
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

    identifier=global::regExWithOneValue("RefNodes");
    if(currentLine.contains(identifier,&match))
    {
        value=match.captured(1);
        global::removeQuotesAndWhitespace(value);
        RefNodes=value;
        global::removeQuotes(RefNodes);
        if(RefNodes.isEmpty()){validLine=false;}
        else{validLine=true;}
        return validLine;
    }

    identifier=global::regExWithOneValue("MeshFileAlias");
    if(currentLine.contains(identifier,&match))
    {
        value=match.captured(1);
        global::removeQuotesAndWhitespace(value);
        MeshFileAlias=value;
        global::removeQuotes(MeshFileAlias);
        if(MeshFileAlias.isEmpty()){validLine=false;}
        else{validLine=true;}
        return validLine;
    }


    QRegularExpression separator( "[ \\t,]+"); //identify any number of white spaces and tabs intermixed or commas
    QStringList list =currentLine.split(separator, Qt::SkipEmptyParts);
    if(!list.isEmpty())
    {
        bool isNumber;
        /*int firstValue=*/list.first().toInt(&isNumber); // a first entry, which is a number, is an elements index
        if(isNumber)
        {
            readingElementList=true;
            validLine=handleElementListLine(currentLine);
            return validLine;
        }
    }
    QString error=QString("Line in Element_Section section not recognized! "+ currentLine);
    logStrings::errorLogString.append(error+"\r\n");
    return validLine;
}

bool ElementSection::handleElementListLine(const QString& currentLine)
{
    QStringList list =stringToStringList(currentLine);
    if(list.length()>=2) //minimum to entries for valid line, i. e. <identifier> <Mesh>
    {
        bool isNumber;
        /*int firstValue=*/list.first().toInt(&isNumber);
        if(isNumber) // current line starts with an element identifier
        {
            readingElementList=true; //programm reads second part of Elements section (declaration of list of elements)
//            index.push_back(firstValue);
        }
        else{
            validLine=false; // current line doesn't start with an element identifier
            return validLine;

        }
        /*int secondValue=*/list.at(1).toInt(&isNumber);
        if(isNumber) //quadrilateral or triangle element directly declared
                    //syntax:     <Index>    <Node-index1>   <Node-index2>   <Node-index3>   <Node-index4>  RefNodes=
        {
            validLine=handleDirectElementDeclaration(currentLine);
            return validLine;
        }
//        if(list.contains("Mesh", Qt::CaseInsensitive)|list.contains("include", Qt::CaseInsensitive)
//                |list.contains("exclude", Qt::CaseInsensitive)||list.contains("MeshFileAlias", Qt::CaseInsensitive))
//        {

//        }
        else{ //list of elements declared via mesh-file
            validLine=handleMeshFileDeclaration(currentLine);
            return validLine;
        }



    }
    else{
        validLine=false;
    }


    return validLine;

}

bool ElementSection::handleDirectElementDeclaration(const QString& currentLine)
{
    QString currentLineCopy=currentLine;
    QRegularExpressionMatch match;
//    QRegularExpression identifier("\\s+RefNodes{1}\\s*={1}\\s*([\"\'][\\w\\s]+[\"\']|\\w+)", QRegularExpression::CaseInsensitiveOption);
    QRegularExpression identifier("\\s+RefNodes{1}\\s*={1}\\s*([\"\'][^\"\']+[\"\']|\\S+)", QRegularExpression::CaseInsensitiveOption);
    QString nodesReference=RefNodes;

    int matchPosition=currentLine.indexOf(identifier, 0, &match);
    if(matchPosition>=0)
    {
//        nodesReference=(currentLine.mid(matchPosition+match.capturedLength())).simplified();
        nodesReference=match.captured(1);
        global::removeQuotes(nodesReference);
//        std::cout <<"Nodes reference found! "<< nodesReference.toUtf8().constData()<<std::endl;
        currentLineCopy.resize(matchPosition); //remove the refnodes part of the current line
//        std::cout <<"currentLineCopy: "<< currentLineCopy.toUtf8().constData()<<std::endl;

    }
    if(nodesReference.isEmpty())
    {
        validLine=false;
        std::cout <<"No nodes reference!: "<< currentLine.toUtf8().constData()<<std::endl;
        QString error=QString("No nodes reference in element declaration in line: "+ currentLine);
        logStrings::errorLogString.append(error+"\r\n");
        return validLine;
    }



    QStringList list =stringToStringList(currentLineCopy);

    bool isNumber;
    bool isNumber2;
    bool isNumber3;
    bool isNumber4;
    bool isNumber5;

    if(list.length()==3) //line only in 2d mode
    {
        //not implemented, therefore
        validLine=false;
        return validLine;
    }
    else{
        if(list.length()==4) //triangle
        {
            int triangleIndex=list.at(0).toInt(&isNumber);
            int firstNode=list.at(1).toInt(&isNumber2);
            int secondNode=list.at(2).toInt(&isNumber3);
            int thirdNode=list.at(3).toInt(&isNumber4);
            if(!(isNumber&&isNumber2&&isNumber3&&isNumber4))
            {
                validLine=false; std::cout <<"Invalid triangle declaration "<< currentLine.toUtf8().constData()<<std::endl;
                QString error=QString("Invalid triangle declaration in line: "+ currentLine);
                logStrings::errorLogString.append(error+"\r\n");
                return validLine;
            }
            if(firstNode==secondNode || firstNode==thirdNode || secondNode==thirdNode)
            {
                validLine=false; std::cout <<"Invalid triangle declaration. Node Indexes should be pairwise distinct: "<< currentLine.toUtf8().constData()<<std::endl;
                QString error=QString("Invalid triangle declaration. Node Indexes should be pairwise distinct. Line: "+ currentLine);
                logStrings::errorLogString.append(error+"\r\n");
                return validLine;
            }
            triangles.push_back(Triangle(triangleIndex,firstNode,secondNode,thirdNode,nodesReference));
            validLine=true;
            return validLine;
        }
        else
        {
            if(list.length()==5) // quadrilateral element
            {
                int quadrilateralIndex=list.at(0).toInt(&isNumber);
                int firstNode=list.at(1).toInt(&isNumber2);
                int secondNode=list.at(2).toInt(&isNumber3);
                int thirdNode=list.at(3).toInt(&isNumber4);
                int fourthNode=list.at(4).toInt(&isNumber5);
                if(!(isNumber&&isNumber2&&isNumber3&&isNumber4&&isNumber5))
                {
                    validLine=false; std::cout <<"Invalid quadrilateral declaration "<< currentLine.toUtf8().constData()<<std::endl;
                    QString error=QString("Invalid quadrilateral declaration in line: "+ currentLine);
                    logStrings::errorLogString.append(error+"\r\n");
                    return validLine;
                }
//                QList<int> nodeIndexes={firstNode,secondNode,thirdNode,fourthNode};
//                QSet<int> set = QSet<int>::fromList(nodeIndexes);
//                QStringList set = QStringList(nodeIndexes);
                if(!global::allIntsUnique(firstNode,secondNode,thirdNode,fourthNode))
                {
                    validLine = false; std::cout <<"Invalid quadrilateral declaration. Node Indexes should pairwise distinct: "<< currentLine.toUtf8().constData()<<std::endl;
                    QString error = QString("Invalid quadrilateral declaration. Node Indexes should be pairwise distinct. Line: "+ currentLine);
                    logStrings::errorLogString.append(error+"\r\n");
                    return validLine;
                }

                quadrilaterals.push_back(Quadrilateral(quadrilateralIndex,firstNode,secondNode,thirdNode,fourthNode,nodesReference));
                validLine = true;
                return validLine;
            }
            else
            { //invalid line
                validLine = false;
            }
        }
    }
    return validLine;
}

bool ElementSection::handleMeshFileDeclaration(const QString& currentLine)
{
    // line syntax:     <Index>    Mesh    MeshFileAlias=  Include <tags>   Exclude <tags>   SwapNormals

    QString currentLineCopy=currentLine;
    QRegularExpressionMatch match;
    QRegularExpression identifier("\\s+MeshFileAlias{1}\\s*={1}\\s*([\"\'][^\"\']+[\"\']|\\S+)", QRegularExpression::CaseInsensitiveOption);
    QString meshFileAlias=MeshFileAlias;

    int matchPosition=currentLine.indexOf(identifier, 0, &match);

     if(matchPosition>=0) //Meshfile alias declaration in element declaration line found
    {
        meshFileAlias=match.captured(1);
        global::removeQuotes(meshFileAlias);
//        std::cout <<"MeshFileAlias: "<< meshFileAlias.toUtf8().constData()<<std::endl;
        currentLineCopy.remove(matchPosition,match.capturedLength() );
//        std::cout <<"Current line with meshfilealias removed: "<< currentLineCopy.toUtf8().constData()<<std::endl;
    }

    if(meshFileAlias.isEmpty())
    {
        validLine=false;
        std::cout <<"No MeshFileAlias!: "<< currentLine.toUtf8().constData()<<std::endl;
        QString error=QString("No MeshFileAlias in mesh-element declaration in line: "+ currentLine);
        logStrings::errorLogString.append(error+"\r\n");
        return validLine;
    }


    // valid line syntax:     <Index>    Mesh  Include <tags>   Exclude <tags>   SwapNormals

    QStringList list =stringToStringListQuotesIntact(currentLineCopy);



    if(!list.isEmpty())
    {
        bool isNumber;
        int elementIndex=list.first().toInt(&isNumber);
        bool localSwapNormals=false;/*=swapNormals;*/
        if(!isNumber)
        {
            validLine=false;
//            readingElementList=false;
            return validLine;
        }
        if(!(currentLine.contains(QRegularExpression("\\s+include([\\s]|$)",QRegularExpression::CaseInsensitiveOption)) || currentLine.contains(QRegularExpression("\\s+Mesh([\\s]|$)",QRegularExpression::CaseInsensitiveOption)) || currentLine.contains(QRegularExpression("\\s+meshfileAlias\\s*=",QRegularExpression::CaseInsensitiveOption))))
        {
            validLine=false;
//            readingElementList=false;
            return validLine;
        }
        int includePos=list.indexOf(QRegularExpression("^\\s*include$", QRegularExpression::CaseInsensitiveOption));
        int excludePos=list.indexOf(QRegularExpression("^\\s*exclude$", QRegularExpression::CaseInsensitiveOption));
        for (int i=0; i<=list.length();i++){
            std::cout << list.value(i).toUtf8().constData()<<std::endl;

        }
//        std::cout<<std::endl <<"include Position: "<< includePos<<std::endl;
//        std::cout <<"exclude Position: "<< excludePos<<std::endl;

        bool includeAll=false;
        QStringList include;
        QStringList exclude;

        if( (includePos==-1) | (includePos==list.length()-1) )
        {
            includeAll=true;
        }
        else
        {
            for(int i=includePos+1;i<=list.length()-1;i++)
            {
                QString currentListElement=list.value(i);
                if(currentListElement.contains(QRegularExpression("^\\s*MeshFileAlias",QRegularExpression::CaseInsensitiveOption)))
                {
                    break;
                }
                if(currentListElement.contains(QRegularExpression("^\\s*exclude\\s*$",QRegularExpression::CaseInsensitiveOption)))
                {
                    break;
                }
                if(currentListElement.contains(QRegularExpression("^\\s*swapNormals\\s*$",QRegularExpression::CaseInsensitiveOption)))
                {
                    localSwapNormals=true;
                    break;
                }
                if(currentListElement.contains(QRegularExpression("^\\s*all\\s*$",QRegularExpression::CaseInsensitiveOption)))
                {
                    includeAll=true;
                    continue;
                }
//                if(0==currentListElement.compare("Mesh",Qt::CaseInsensitive))
//                {
//                    continue;
//                }
                include.push_back(currentListElement);
//                std::cout <<"include: "<< currentListElement.toUtf8().constData()<<std::endl;

            }
        }
        if(includeAll)
        {
//            std::cout <<"includeAll=true"<<std::endl;
        }

        if(excludePos>=1 && (excludePos<=list.length()-2))
        {
            for(int i=excludePos+1;i<=list.length()-1;i++)
            {
                QString currentListElement=list.value(i);
                if(currentListElement.contains(QRegularExpression("^\\s*MeshFileAlias",QRegularExpression::CaseInsensitiveOption)))
                {
                    break;
                }
                if(currentListElement.contains(QRegularExpression("^\\s*include\\s*$",QRegularExpression::CaseInsensitiveOption)))
                {
                    break;
                }
                if(currentListElement.contains(QRegularExpression("^\\s*swapNormals\\s*$",QRegularExpression::CaseInsensitiveOption)))
                {
                    localSwapNormals=true;
                    break;
                }
//                if(0==currentListElement.compare("Mesh",Qt::CaseInsensitive))
//                {
//                    continue;
//                }

                exclude.push_back(currentListElement);
//                std::cout <<"exclude: "<< currentListElement.toUtf8().constData()<<std::endl;
            }
        }
//        if(include.isEmpty())
//        {
//            includeAll=true;
//        }
        meshFileElements.push_back(MeshFileElement(elementIndex,meshFileAlias,include,exclude,localSwapNormals,includeAll));
        validLine=true;
    }
    else
    {
        validLine=false;
    }
    return validLine;
}

//QRegularExpression ElementSection::regExWithOneValue(const QString& identifier)
//{
////    return QRegularExpression("^\\s*"+identifier+"\\s*={1}\\s*([\"\'][\\w\\s]+[\"\']|\\w+)", QRegularExpression::CaseInsensitiveOption);
//    return QRegularExpression("^\\s*"+identifier+"\\s*={1}\\s*([\"\'][^\"\']+[\"\']|.+$)", QRegularExpression::CaseInsensitiveOption);

//}

QRegularExpression ElementSection::regExWithOptionalValue(const QString& identifier)
{
    return QRegularExpression("^\\s*"+identifier+"\\s*(={1}\\s*([\"\'][^\"\']+[\"\']|.+$))?", QRegularExpression::CaseInsensitiveOption);
}

QStringList ElementSection::stringToStringList(const QString& currentLine)
{
    QRegularExpression separator( "[ \\t,]+"); //identify any number of white spaces and tabs intermixed or commas
    QStringList list =currentLine.split(separator, Qt::SkipEmptyParts);
    return list;

}

QStringList ElementSection::stringToStringListQuotesIntact(const QString& currentLine)
{
    QStringList list;
    QRegularExpression strings( "[\"\']([^\"\']*)[\"\']|(\\S+)"); //identify any string between "or'-quotes or any string which isn't whitespace
    QRegularExpressionMatchIterator iterator = strings.globalMatch(currentLine);
    while (iterator.hasNext()) {
        QRegularExpressionMatch match = iterator.next();
        list.push_back(match.captured(0));
    }
    return list;
}

//void ElementSection::removeQuotes(QString& line)
//{
//    QRegularExpression removeChars("[\"\']");
//    line.remove(removeChars);
//}

//bool ElementSection::uniqueIdentifiers(){
//    sortElements();

//}

bool ElementSection::containsNoElements()
{
if(triangles.isEmpty() && quadrilaterals.isEmpty() && meshFileElements.isEmpty())
{
    return true;

}
return false;
}

bool ElementSection::allMandatoryParametersSet()
{
    if(Subdomain==0)
    {
        Subdomain=1;
        std::cout <<"Subdomain not declared in Element "<< name.toUtf8().constData()<<std::endl;

//        return false;
    }
    return true;

}
void ElementSection::swapTriangleNormals(QVector<VectorTriangle>& elementsVectorTriangles)
{
    for(int i=0;i<elementsVectorTriangles.length();i++)
    {
        VectorTriangle tmp=elementsVectorTriangles.at(i);
        elementsVectorTriangles[i] = VectorTriangle(tmp.elementIndex,tmp.node3,tmp.node2,tmp.node1,-tmp.normal);
    }
}

void ElementSection::swapQuadrilateralNormals(QVector<VectorQuadrilateral>& elementsVectorQuadrilaterals)
{
    for(int i=0;i<elementsVectorQuadrilaterals.length();i++)
    {
        VectorQuadrilateral tmp=elementsVectorQuadrilaterals.at(i);
        elementsVectorQuadrilaterals[i]=VectorQuadrilateral(tmp.elementIndex,tmp.node4,tmp.node3,tmp.node2,tmp.node1);
    }
}

QVector<VectorTriangle> ElementSection::refineTriangles(const QVector<VectorTriangle>& elementsVectorTriangles)
{
    QVector<VectorTriangle> returnVector;
    int vectorLength=elementsVectorTriangles.length();
    for(int i=0;i<vectorLength;i++)
    {
        VectorTriangle tmp=elementsVectorTriangles.at(i);
        double edgeLength1=lineLength(tmp.node1,tmp.node2);
        double edgeLength2=lineLength(tmp.node2,tmp.node3);
        double edgeLength3=lineLength(tmp.node3,tmp.node1);
        double maxEdge = 0;
//        Eigen::Vector3d maxEdgeMidpoint;
        int max=indexOfMaxOfThree(edgeLength1, edgeLength2, edgeLength3);
        switch(max) {
            case 1 : maxEdge=edgeLength1;
                     break;
            case 2 : maxEdge=edgeLength2;
                     break;
            case 3 : maxEdge=edgeLength3;
                     break;
        }
        if(maxEdge<=EdgeLength)
        {
            returnVector.push_back(tmp);
            continue;
        }
        else
        {
            recursiveRefineTriangle(tmp,max,returnVector);
        }
//        switch(max) {
//            case 1 : maxEdgeMidpoint=middleOfLine(tmp.node1,tmp.node2);
//                    elementsVectorTriangles[i]=VectorTriangle(tmp.elementIndex,tmp.node1,maxEdgeMidpoint,tmp.node3);
//                    elementsVectorTriangles.push_back(VectorTriangle(tmp.elementIndex,maxEdgeMidpoint,tmp.node2,tmp.node3));
//                     break;
//            case 2 : maxEdgeMidpoint=middleOfLine(tmp.node2,tmp.node3);
//                    elementsVectorTriangles[i]=VectorTriangle(tmp.elementIndex,tmp.node1,tmp.node2,maxEdgeMidpoint);
//                    elementsVectorTriangles.push_back(VectorTriangle(tmp.elementIndex,maxEdgeMidpoint,tmp.node3,tmp.node1));
//                     break;
//            case 3 : maxEdgeMidpoint=middleOfLine(tmp.node3,tmp.node1);
//                    elementsVectorTriangles[i]=VectorTriangle(tmp.elementIndex,tmp.node1,tmp.node2,maxEdgeMidpoint);
//                    elementsVectorTriangles.push_back(VectorTriangle(tmp.elementIndex,maxEdgeMidpoint,tmp.node2,tmp.node3));
//                     break;
//        }
//        if(max==1){

//            elementsVectorTriangles[i]=VectorTriangle(tmp.elementIndex,tmp.node1,tmp.node2,tmp.node3);
//        }
//        else if(max==2){
//            maxEdge=edgeLength2;
//        }
//        else if(max==3){
//            maxEdge=edgeLength3;
//        }
    }
    return returnVector;
}
void ElementSection::recursiveRefineTriangle(VectorTriangle parentTriangle, int longestEdgeIndex, QVector<VectorTriangle>& vectorToStoreTriangles)
{
    double maxEdge = 0;
    switch(longestEdgeIndex) {
        case 1 : maxEdge=lineLength(parentTriangle.node1,parentTriangle.node2);
                 break;
        case 2 : maxEdge=lineLength(parentTriangle.node2,parentTriangle.node3);
                 break;
        case 3 : maxEdge=lineLength(parentTriangle.node3,parentTriangle.node1);
                 break;
    }
    if(maxEdge<=EdgeLength)
    {
        vectorToStoreTriangles.push_back(parentTriangle);
        return;
    }
    else
    {
        VectorTriangle vectorTriangle1;
        VectorTriangle vectorTriangle2;
        Eigen::Vector3d maxEdgeMidpoint;

        switch(longestEdgeIndex) {
            case 1 :  maxEdgeMidpoint=middleOfLine(parentTriangle.node1,parentTriangle.node2);
                      vectorTriangle1=VectorTriangle(parentTriangle.elementIndex,parentTriangle.node1,maxEdgeMidpoint,parentTriangle.node3,parentTriangle.normal);
                      vectorTriangle2=VectorTriangle(parentTriangle.elementIndex,maxEdgeMidpoint,parentTriangle.node2,parentTriangle.node3,parentTriangle.normal);
                     break;
            case 2 : maxEdgeMidpoint=middleOfLine(parentTriangle.node2,parentTriangle.node3);
                     vectorTriangle1=VectorTriangle(parentTriangle.elementIndex,parentTriangle.node1,parentTriangle.node2,maxEdgeMidpoint,parentTriangle.normal);
                     vectorTriangle2=VectorTriangle(parentTriangle.elementIndex,maxEdgeMidpoint,parentTriangle.node3,parentTriangle.node1,parentTriangle.normal);
                     break;
            case 3 : maxEdgeMidpoint=middleOfLine(parentTriangle.node3,parentTriangle.node1);
                     vectorTriangle1=VectorTriangle(parentTriangle.elementIndex,parentTriangle.node1,parentTriangle.node2,maxEdgeMidpoint,parentTriangle.normal);
                     vectorTriangle2=VectorTriangle(parentTriangle.elementIndex,maxEdgeMidpoint,parentTriangle.node2,parentTriangle.node3,parentTriangle.normal);
                     break;
        }
//        vectorTriangle1.normal=parentTriangle.normal;
//        vectorTriangle2.normal=parentTriangle.normal;
        vectorTriangle1.robinBoundaryCondition=parentTriangle.robinBoundaryCondition;
        vectorTriangle2.robinBoundaryCondition=parentTriangle.robinBoundaryCondition;
        vectorTriangle1.color=parentTriangle.color;
        vectorTriangle2.color=parentTriangle.color;
        double edgeLength1=lineLength(vectorTriangle1.node1,vectorTriangle1.node2);
        double edgeLength2=lineLength(vectorTriangle1.node2,vectorTriangle1.node3);
        double edgeLength3=lineLength(vectorTriangle1.node3,vectorTriangle1.node1);
        int max=indexOfMaxOfThree(edgeLength1, edgeLength2, edgeLength3);
        recursiveRefineTriangle(vectorTriangle1,max,vectorToStoreTriangles);

         edgeLength1=lineLength(vectorTriangle2.node1,vectorTriangle2.node2);
         edgeLength2=lineLength(vectorTriangle2.node2,vectorTriangle2.node3);
         edgeLength3=lineLength(vectorTriangle2.node3,vectorTriangle2.node1);
         max=indexOfMaxOfThree(edgeLength1, edgeLength2, edgeLength3);
         recursiveRefineTriangle(vectorTriangle2,max,vectorToStoreTriangles);
    }
}


QVector<VectorTriangle> ElementSection::refineQuadrilaterals(const QVector<VectorQuadrilateral>& elementsVectorQuadrilaterals)
{
//    QVector<VectorTriangle> triangles;
    QVector<VectorQuadrilateral> refinedQuadrilaterals;
    for(int i=0;i<elementsVectorQuadrilaterals.length();i++) //split quadrilaterals with too high aspect ratio recursively in two quadrilaterals
    {
        VectorQuadrilateral parentQuadrilateral=elementsVectorQuadrilaterals.at(i);

        double edgeLength1=lineLength(parentQuadrilateral.node1,parentQuadrilateral.node2);
        double edgeLength2=lineLength(parentQuadrilateral.node2,parentQuadrilateral.node3);
        double edgeLength3=lineLength(parentQuadrilateral.node3,parentQuadrilateral.node4);
        double edgeLength4=lineLength(parentQuadrilateral.node4,parentQuadrilateral.node1);
        QPair<double,int> longestEdgeWithIndex=maxOfFour(edgeLength1, edgeLength2, edgeLength3, edgeLength4);
        double aspectRatio=global::quadrilateralAspectRatio(parentQuadrilateral);
        if(aspectRatio>3 && longestEdgeWithIndex.first>EdgeLength) //refine Quadrilateral if aspect ratio is high end the longer edge is longer than the specified max edge length
        {
            recursiveRefineQuadrilateral(parentQuadrilateral, longestEdgeWithIndex.second, refinedQuadrilaterals);

        }
        else
        {
            refinedQuadrilaterals.push_back(parentQuadrilateral);
        }
    }
    QVector<VectorTriangle> triangles;
    triangles=global::quadrilateralsToTriangles(refinedQuadrilaterals);
    refinedQuadrilaterals.clear();
//    triangles=refineTriangles(triangles);
    return triangles;
}

void ElementSection::recursiveRefineQuadrilateral(VectorQuadrilateral parentQuadrilateral, int longestEdgeIndex, QVector<VectorQuadrilateral>& refinedQuadrilateralsStorage)
{
    VectorQuadrilateral vectorQuadrilateral1;
    VectorQuadrilateral vectorQuadrilateral2;
    Eigen::Vector3d maxEdgeMidpoint;
    Eigen::Vector3d oppositeEdgeMidpoint;

    switch(longestEdgeIndex) {
        case 1 :  maxEdgeMidpoint=middleOfLine(parentQuadrilateral.node1,parentQuadrilateral.node2);
                  oppositeEdgeMidpoint=middleOfLine(parentQuadrilateral.node3,parentQuadrilateral.node4);
                  vectorQuadrilateral1=VectorQuadrilateral(parentQuadrilateral.elementIndex,parentQuadrilateral.node1,maxEdgeMidpoint,oppositeEdgeMidpoint,parentQuadrilateral.node4);
                  vectorQuadrilateral2=VectorQuadrilateral(parentQuadrilateral.elementIndex,parentQuadrilateral.node2,parentQuadrilateral.node3,oppositeEdgeMidpoint,maxEdgeMidpoint);
                 break;
        case 2 : maxEdgeMidpoint=middleOfLine(parentQuadrilateral.node2,parentQuadrilateral.node3);
                 oppositeEdgeMidpoint=middleOfLine(parentQuadrilateral.node4,parentQuadrilateral.node1);
                 vectorQuadrilateral1=VectorQuadrilateral(parentQuadrilateral.elementIndex,parentQuadrilateral.node2,maxEdgeMidpoint,oppositeEdgeMidpoint,parentQuadrilateral.node1);
                 vectorQuadrilateral2=VectorQuadrilateral(parentQuadrilateral.elementIndex,parentQuadrilateral.node3,parentQuadrilateral.node4,oppositeEdgeMidpoint,maxEdgeMidpoint);
                 break;
        case 3 : maxEdgeMidpoint=middleOfLine(parentQuadrilateral.node3,parentQuadrilateral.node4);
                 oppositeEdgeMidpoint=middleOfLine(parentQuadrilateral.node1,parentQuadrilateral.node2);
                 vectorQuadrilateral1=VectorQuadrilateral(parentQuadrilateral.elementIndex,parentQuadrilateral.node3,maxEdgeMidpoint,oppositeEdgeMidpoint,parentQuadrilateral.node2);
                 vectorQuadrilateral2=VectorQuadrilateral(parentQuadrilateral.elementIndex,parentQuadrilateral.node4,parentQuadrilateral.node1,oppositeEdgeMidpoint,maxEdgeMidpoint);
                 break;
        case 4 : maxEdgeMidpoint=middleOfLine(parentQuadrilateral.node4,parentQuadrilateral.node1);
                 oppositeEdgeMidpoint=middleOfLine(parentQuadrilateral.node2,parentQuadrilateral.node3);
                 vectorQuadrilateral1=VectorQuadrilateral(parentQuadrilateral.elementIndex,parentQuadrilateral.node4,maxEdgeMidpoint,oppositeEdgeMidpoint,parentQuadrilateral.node3);
                 vectorQuadrilateral2=VectorQuadrilateral(parentQuadrilateral.elementIndex,parentQuadrilateral.node1,parentQuadrilateral.node2,oppositeEdgeMidpoint,maxEdgeMidpoint);
                 break;
    }
    vectorQuadrilateral1.robinBoundaryCondition=parentQuadrilateral.robinBoundaryCondition;
    vectorQuadrilateral2.robinBoundaryCondition=parentQuadrilateral.robinBoundaryCondition;
    vectorQuadrilateral1.color=parentQuadrilateral.color;
    vectorQuadrilateral2.color=parentQuadrilateral.color;

    double aspectRatio=global::quadrilateralAspectRatio(vectorQuadrilateral1);
    if(aspectRatio>3 ) //refine Quadrilateral if aspect ratio is high end the longer edge is longer than the specified max edge length
    {
        double edgeLength1=lineLength(vectorQuadrilateral1.node1,vectorQuadrilateral1.node2);
        double edgeLength2=lineLength(vectorQuadrilateral1.node2,vectorQuadrilateral1.node3);
        double edgeLength3=lineLength(vectorQuadrilateral1.node3,vectorQuadrilateral1.node4);
        double edgeLength4=lineLength(vectorQuadrilateral1.node4,vectorQuadrilateral1.node1);
        QPair<double,int> longestEdgeWithIndex=maxOfFour(edgeLength1, edgeLength2, edgeLength3, edgeLength4);
        if(longestEdgeWithIndex.first>EdgeLength)
        {
            recursiveRefineQuadrilateral(vectorQuadrilateral1, longestEdgeWithIndex.second, refinedQuadrilateralsStorage);
        }
        else
        {
            refinedQuadrilateralsStorage.push_back(vectorQuadrilateral1);
        }
    }
    else
    {
        refinedQuadrilateralsStorage.push_back(vectorQuadrilateral1);
    }
    double aspectRatio2= global::quadrilateralAspectRatio(vectorQuadrilateral2);
    if(aspectRatio2>3)
    {
        double edgeLength1=lineLength(vectorQuadrilateral2.node1,vectorQuadrilateral2.node2);
        double edgeLength2=lineLength(vectorQuadrilateral2.node2,vectorQuadrilateral2.node3);
        double edgeLength3=lineLength(vectorQuadrilateral2.node3,vectorQuadrilateral2.node4);
        double edgeLength4=lineLength(vectorQuadrilateral2.node4,vectorQuadrilateral2.node1);
        QPair<double,int> longestEdgeWithIndex=maxOfFour(edgeLength1, edgeLength2, edgeLength3, edgeLength4);
        if(longestEdgeWithIndex.first>EdgeLength)
        {
        recursiveRefineQuadrilateral(vectorQuadrilateral2, longestEdgeWithIndex.second, refinedQuadrilateralsStorage);
        }
        else
        {
            refinedQuadrilateralsStorage.push_back(vectorQuadrilateral2);
        }

    }
    else
    {
        refinedQuadrilateralsStorage.push_back(vectorQuadrilateral2);
    }
}

quint64 ElementSection::indexOfMaxOfThree(double a, double b, double c)
{
    if(a>b)
    {
        if(a>c)
        {
            return 1;
        }
        else
        {
            return 3;
        }
    }
    else
    {
        if(b>c)
        {
            return 2;
        }
        else
        {
            return 3;
        }

    }
}
QPair<double,quint64> ElementSection::maxOfFour(double a, double b, double c,double d)
{
    if(a>b)
    {
        if(a>c)
        {
            if(a>d)
            {
                return qMakePair(a,1);
            }
            else
            {
                return qMakePair(d,4);
            }
        }
        else
        {
            if(c>d)
            {
                return qMakePair(c,3);
            }
            else
            {
                return qMakePair(d,4);
            }
        }
    }
    else
    {
        if(b>c)
        {
            if(b>d)
            {
                return qMakePair(b,2);
            }
            else
            {
                return qMakePair(d,4);
            }
        }
        else
        {
            if(c>d)
            {
                return qMakePair(c,3);
            }
            else
            {
                return qMakePair(d,4);
            }
        }

    }
}

Eigen::Vector3d ElementSection::middleOfLine(Eigen::Vector3d point1,Eigen::Vector3d point2)
{
    double x=(point1(0)+point2(0))/2;
    double y=(point1(1)+point2(1))/2;
    double z=(point1(2)+point2(2))/2;
    return Eigen::Vector3d {x,y,z};
}

double ElementSection::lineLength(Eigen::Vector3d point1,Eigen::Vector3d point2)
{
    double lineLength=std::sqrt(std::pow(point1(0)-point2(0),2) + std::pow(point1(1)-point2(1),2) + std::pow(point1(2)-point2(2),2));
    return lineLength;
}

void ElementSection::scaleAndTransformElements()
{
    if(Scale(0) == 0){Scale(0) = 1;}
    if(Scale(1) == 0){Scale(1) = 1;}
    if(Scale(2) == 0){Scale(2) = 1;}
    bool reallyScale=true;
    if(Scale(0) == 1 && Scale(1) == 1 && Scale(2) == 1)
    {
        reallyScale=false;
    }
    if(reallyScale)
    {
        for(int i=0;i<elementsVectorQuadrilaterals.length();i++)
        {
            elementsVectorQuadrilaterals[i].node1=elementsVectorQuadrilaterals[i].node1.array()*Scale.array();
            elementsVectorQuadrilaterals[i].node2=elementsVectorQuadrilaterals[i].node2.array()*Scale.array();
            elementsVectorQuadrilaterals[i].node3=elementsVectorQuadrilaterals[i].node3.array()*Scale.array();
            elementsVectorQuadrilaterals[i].node4=elementsVectorQuadrilaterals[i].node4.array()*Scale.array();
        }
        for(int i=0;i<elementsVectorTriangles.length();i++)
        {
            elementsVectorTriangles[i].node1=elementsVectorTriangles[i].node1.array()*Scale.array();
            elementsVectorTriangles[i].node2=elementsVectorTriangles[i].node2.array()*Scale.array();
            elementsVectorTriangles[i].node3=elementsVectorTriangles[i].node3.array()*Scale.array();
        }
    }
    if(transformationMatrix.isIdentity()==false)
    {
        Eigen::Matrix4d transformationMatrixEigen=MeshFunctions::qtToEigenMatrix(transformationMatrix);
        MeshFunctions::transformTriangles(elementsVectorTriangles,transformationMatrixEigen);
        MeshFunctions::transformQuadrilaterals(elementsVectorQuadrilaterals,transformationMatrixEigen);
    }
}

void ElementSection::setTrianglesWallimpedance(int index, RobinBoundaryCondition robinBoundaryCondition)
{
    for(int i=0;i<elementsVectorTriangles.size();i++)
    {
        if(elementsVectorTriangles.at(i).elementIndex==index)
        {
            elementsVectorTriangles[i].robinBoundaryCondition=robinBoundaryCondition;
            elementsVectorTriangles[i].color = {0.5, 1, 0.5};
//            std::cout<<"wallimpedance changed."<<std::endl;
        }
    }
}

void ElementSection::setQuadrilateralsWallimpedance(int index, RobinBoundaryCondition robinBoundaryCondition)
{
    for(int i=0;i<elementsVectorQuadrilaterals.size();i++)
    {
        if(elementsVectorQuadrilaterals.at(i).elementIndex==index)
        {
            elementsVectorQuadrilaterals[i].robinBoundaryCondition=robinBoundaryCondition;
            elementsVectorQuadrilaterals[i].color = {0.5, 1, 0.5};
        }
    }
}

void ElementSection::setLocalWallimpedance(RobinBoundaryCondition robinBoundaryCondition)
{
    for(int i=0; i<elementsVectorQuadrilaterals.size(); i++)
    {
//        std::cout << "robinBoundaryCondition: " << robinBoundaryCondition.a << " " << robinBoundaryCondition.b << " " << robinBoundaryCondition.g << std::endl;
        elementsVectorQuadrilaterals[i].robinBoundaryCondition = robinBoundaryCondition;
        elementsVectorQuadrilaterals[i].color = {0.5, 1, 0.5};
    }
    for(int i=0;i<elementsVectorTriangles.size();i++)
    {
//        std::cout << "robinBoundaryCondition: " << robinBoundaryCondition.a << " " << robinBoundaryCondition.b << " " << robinBoundaryCondition.g << std::endl;
        elementsVectorTriangles[i].robinBoundaryCondition=robinBoundaryCondition;
        elementsVectorTriangles[i].color = {0.5, 1, 0.5};
    }
}

void ElementSection::setTrianglesDrivingWeight(int index, std::complex<double> weight)
{
    for(int i=0;i<elementsVectorTriangles.size();i++)
    {
        if(elementsVectorTriangles.at(i).elementIndex==index)
        {
            elementsVectorTriangles[i].robinBoundaryCondition.g=weight;
            elementsVectorTriangles[i].robinBoundaryCondition.b=1;
            elementsVectorTriangles[i].color={0.85,0,0};
//            std::cout<<"driving value changed."<<std::endl;
        }
    }
}

void ElementSection::setQuadrilateralsDrivingWeight(int index, std::complex<double> weight)
{
    for(int i=0;i<elementsVectorQuadrilaterals.size();i++)
    {
        if(elementsVectorQuadrilaterals.at(i).elementIndex==index)
        {
            elementsVectorQuadrilaterals[i].robinBoundaryCondition.g=weight;
            elementsVectorQuadrilaterals[i].robinBoundaryCondition.b=1;
            elementsVectorQuadrilaterals[i].color={0.85,0,0};
        }
    }
}

void ElementSection::setLocalDrivingWeight(std::complex<double> weight)
{
    for(int i=0;i<elementsVectorQuadrilaterals.size();i++)
    {
        elementsVectorQuadrilaterals[i].robinBoundaryCondition.g = weight;
        elementsVectorQuadrilaterals[i].robinBoundaryCondition.a = 0;
        elementsVectorQuadrilaterals[i].robinBoundaryCondition.b = 1;
        elementsVectorQuadrilaterals[i].color={0.85,0,0};
    }
    for(int i=0;i<elementsVectorTriangles.size();i++)
    {
        elementsVectorTriangles[i].robinBoundaryCondition.g = weight;
        elementsVectorTriangles[i].robinBoundaryCondition.a = 0;
        elementsVectorTriangles[i].robinBoundaryCondition.b = 1;
        elementsVectorTriangles[i].color={0.85,0,0};
    }
}
