#include "pressurepointssection.h"

PressurePointsSection::PressurePointsSection()
{

}

bool PressurePointsSection::handleScriptLine(const QString& currentLine)
{

    validLine=false;
    if(!readingElementList) //an elements section contains two parts: first the parameters and then the list of elements
    {
        validLine=handleParameterSectionLine(currentLine);
    }
    else {
        validLine=handlePressurePointsListLine(currentLine);
    }
    return validLine;
}

bool PressurePointsSection::handleParameterSectionLine(const QString& currentLine)
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

    identifier=global::regExWithOneValue("DrvGroup");
    if(currentLine.contains(identifier,&match))
    {
        value=match.captured(1);
        global::removeQuotesAndWhitespace(value);
        DrvGroup=value.toUInt(&validLine);
        if(!validLine)
        {
            std::cout<<"Invalid DrvGroup in Line \""<<currentLine.toUtf8().constData()<<"\" of PressurePoints section "<< name.toUtf8().constData()<<"."<<std::endl;
            QString error=QString("Invalid DrvGroup in Line \""+currentLine+"\" of PressurePoints section "+ name+".");
            logStrings::errorLogString.append(error+"\r\n");
        }
        return validLine;
    }
//    identifier="swapNormals";
//    identifier=global::regExWithOneValue("swapNormals");
//    identifier=regExWithOptionalValue("swapNormals");

//    if(currentLine.contains(identifier,&match))
//    {
//        std::cout<<"normals Swapped!!!"<<std::endl;
//          value=match.captured(2);
//          std::cout<<value.toUtf8().constData()<<std::endl;
//          global::removeQuotesAndWhitespace(value);
//          if(!value.isEmpty())
//          {
//              if(value.compare("true",Qt::CaseInsensitive)==0||value.toInt()==1)
//              {
//                  swapNormals=true;
//              }
//          }
//          else
//          {
//              swapNormals=true;
//          }
//          validLine=true;
//          return validLine;
//    }

//    identifier="RefNodes=";
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

//    identifier="MeshFileAlias=";
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

    identifier=global::regExWithOneValue("MeshMode");
    if(currentLine.contains(identifier,&match))
    {
        value=match.captured(1);
        global::removeQuotesAndWhitespace(value);
        if(0==QString::compare("Center",value, Qt::CaseInsensitive))
        {
              validLine=true;
              pointsOnCenterOrOnElements=true;
              return validLine;
        }
        if(0==QString::compare("Elements",value, Qt::CaseInsensitive))
        {
            validLine=true;
              pointsOnCenterOrOnElements=false;
              return validLine;
        }
        validLine=false;
        return validLine;
    }

    identifier=global::regExWithOneValue("MeshMinDistance");
    if(currentLine.contains(identifier,&match))
    {
        value=match.captured(1);
        global::removeQuotesAndWhitespace(value);
        MeshMinDistance=global::stringWithUnitsToDouble(value,QString("m"),validLine);
        if(MeshMinDistance<0){validLine=false; std::cout<<"MeshMinDistance <0 in Line \""<<currentLine.toUtf8().constData()<<"\" of Pressure_Points section "<< name.toUtf8().constData()<<"."<<std::endl;}
        return validLine;
    }

    identifier=global::regExWithOneValue("DrvWeight");
    if(currentLine.contains(identifier,&match))
    {
        value=match.captured(1);
        global::removeQuotesAndWhitespace(value);
        DrvWeight=value.toDouble(&validLine);
        if(validLine)
        {
            if(DrvWeight<=0)
            {
                std::cout<<"DrvWeight<=0 in Line \""<<currentLine.toUtf8().constData()<<"\" of Pressure_Points section "<< name.toUtf8().constData()<<"."<<std::endl;
                QString error=QString("DrvWeight<=0 in Line \""+currentLine+"\" of Pressure_Points section "+name+".");
                logStrings::errorLogString.append(error+"\r\n");
            }
        }
        return validLine;
    }

    identifier=global::regExWithOneValue("DrvDelay");
    if(currentLine.contains(identifier,&match))
    {
        value=match.captured(1);
        global::removeQuotesAndWhitespace(value);
        DrvDelay=global::stringWithUnitsToDouble(value,QString("s"),validLine);
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
            validLine=handlePressurePointsListLine(currentLine);
            return validLine;
        }
    }
    return validLine;
}

bool PressurePointsSection::handlePressurePointsListLine(const QString& currentLine)
{
    QStringList list =global::stringToStringListQuotesIntact(currentLine);
    if(list.length()>=2) //minimum to entries for valid line, i. e. <identifier> <Mesh>
    {
        bool isNumber;
        list.first().toInt(&isNumber);
        if(isNumber) // current line starts with an element identifier
        {
            readingElementList=true; //programm reads second part of Pressure_Points section (declaration of list of elements)
//            index.push_back(firstValue);
        }
        else{
            validLine=false; // current line doesn't start with an element identifier
            return validLine;

        }

        list.at(1).toInt(&isNumber);
        if(isNumber) //pressure point directly declared
                    //syntax:     <Index>    <Node-index1>   DrvGroup=   Weight=   Delay=   RefNodes=
        {
            validLine=handleDirectPressurePointDeclaration(currentLine);
            return validLine;
        }

        else{ //list of pressure points declared via mesh-file
            validLine=handleMeshFileDeclaration(currentLine);
            return validLine;
        }



    }
    else{
        validLine=false;
    }


    return validLine;

}

bool PressurePointsSection::handleDirectPressurePointDeclaration(const QString& currentLine)
{
    QString currentLineCopy=currentLine;
    QRegularExpressionMatch match;
    //find RefElements
    QRegularExpression identifier("\\s+RefNodes{1}\\s*={1}\\s*([\"\'][^\"\']+[\"\']|\\S+)", QRegularExpression::CaseInsensitiveOption);
    QString localRefNodes=RefNodes;

    int matchPosition=currentLineCopy.indexOf(identifier, 0, &match);
    if(matchPosition>=0) //RefNodes declaration in pressure point declaration line found
    {
        localRefNodes=match.captured(1);
        global::removeQuotes(localRefNodes);
//        std::cout <<"RefElements: "<< localRefNodes.toUtf8().constData()<<std::endl;
        currentLineCopy.remove(matchPosition,match.capturedLength() );
//        std::cout <<"Current line with RefNodes removed: "<< currentLineCopy.toUtf8().constData()<<std::endl;
    }

    if(localRefNodes.isEmpty())
    {
        validLine=false;
        std::cout <<"No RefElements!: "<< currentLine.toUtf8().constData()<<std::endl;
        QString error=QString("No RefElements in Pressure_Points declaration. Line: "+currentLine);
        logStrings::errorLogString.append(error+"\r\n");
        return validLine;
    }

    //find Delay
    identifier=QRegularExpression("\\s+Delay{1}\\s*={1}\\s*([\"\'][^\"\']+[\"\']|\\S+)", QRegularExpression::CaseInsensitiveOption);
    double localDelay=DrvDelay;
    bool isNumber=true;
    matchPosition=currentLineCopy.indexOf(identifier, 0, &match);
    if(matchPosition>=0) //DrvDelay declaration in pressure point declaration line found
    {
        QString value=match.captured(1);
        global::removeQuotesAndWhitespace(value);
        localDelay=global::stringWithUnitsToDouble(value,QString("s"),isNumber);
//        std::cout <<"DrvDelay: "<< localDelay<<std::endl;
        currentLineCopy.remove(matchPosition,match.capturedLength() );
//        std::cout <<"Current line with DrvDelay removed: "<< currentLineCopy.toUtf8().constData()<<std::endl;
    }

    if(!isNumber)
    {
        validLine=false;
        std::cout <<"Invalid DrvDelay declaration: "<< currentLine.toUtf8().constData()<<std::endl;
        QString error=QString("Invalid DrvDelay declaration in line: "+currentLine);
        logStrings::errorLogString.append(error+"\r\n");
        return validLine;
    }

    //find Weight
    identifier=QRegularExpression("\\s+Weight{1}\\s*={1}\\s*([\"\'][^\"\']+[\"\']|\\S+)", QRegularExpression::CaseInsensitiveOption);
    double localWeight=DrvWeight;
    isNumber=true;
    matchPosition=currentLineCopy.indexOf(identifier, 0, &match);
    if(matchPosition>=0) //DrvWeight declaration in pressure point declaration line found
    {
        QString value=match.captured(1);
        global::removeQuotesAndWhitespace(value);
        localWeight=value.toDouble(&isNumber);
        std::cout <<"DrvWeight: "<< localWeight<<std::endl;
        if(DrvWeight<=0){std::cout<<"DrvWeight<=0 in Line \""<<currentLine.toUtf8().constData()<<"\" of Pressure_Points section "<< name.toUtf8().constData()<<"."<<std::endl;}
        currentLineCopy.remove(matchPosition,match.capturedLength() );
//        std::cout <<"Current line with DrvWeight removed: "<< currentLineCopy.toUtf8().constData()<<std::endl;
    }

    if(!isNumber)
    {
        validLine=false;
        std::cout <<"Invalid DrvWeight declaration: "<< currentLine.toUtf8().constData()<<std::endl;
        QString error=QString("Invalid DrvWeight declaration in line: "+currentLine);
        logStrings::errorLogString.append(error+"\r\n");
        return validLine;
    }

    //find DrvGroup
    identifier=QRegularExpression("\\s+DrvGroup{1}\\s*={1}\\s*([\"\'][^\"\']+[\"\']|\\S+)", QRegularExpression::CaseInsensitiveOption);
    int localDrvGroup=DrvGroup;
    isNumber=true;
    matchPosition=currentLineCopy.indexOf(identifier, 0, &match);
    if(matchPosition>=0) //DrvGroup declaration in pressure point declaration line found
    {
        QString value=match.captured(1);
        global::removeQuotesAndWhitespace(value);
        localDrvGroup=value.toUInt(&isNumber);
        std::cout <<"DrvGroup: "<< localDrvGroup<<std::endl;
        currentLineCopy.remove(matchPosition,match.capturedLength() );
//        std::cout <<"Current line with DrvGroup removed: "<< currentLineCopy.toUtf8().constData()<<std::endl;
    }

    if(!isNumber)
    {
        validLine=false;
        std::cout <<"Invalid DrvGroup declaration: "<< currentLine.toUtf8().constData()<<std::endl;
        QString error=QString("Invalid DrvGroup declaration in line: "+currentLine);
        logStrings::errorLogString.append(error+"\r\n");
        return validLine;
    }

    QStringList list =global::stringToStringListQuotesIntact(currentLineCopy);

    if(list.length()<2)
    {
        std::cout<<"No Node-index declaration. "<< currentLine.toUtf8().constData()<<std::endl;
        QString error=QString("No Node-index declaration in line: "+currentLine);
        logStrings::errorLogString.append(error+"\r\n");
        validLine=false;
        return validLine;
    }
    else
    {
        int index=list.at(0).toInt(&isNumber);
        if(!isNumber)
        {
            validLine=false;
            std::cout <<"Invalid index declaration: "<< currentLine.toUtf8().constData()<<std::endl;
            QString error=QString("Invalid index declaration in line: "+currentLine);
            logStrings::errorLogString.append(error+"\r\n");
            return validLine;
        }
        int nodeIndex=list.at(1).toInt(&isNumber);
        if(!isNumber)
        {
            validLine=false;
            std::cout <<"Invalid node-index declaration: "<< currentLine.toUtf8().constData()<<std::endl;
            QString error=QString("Invalid node-index declaration in line: "+currentLine);
            logStrings::errorLogString.append(error+"\r\n");
            return validLine;
        }

//        pressurePoints.insertMulti(index,PressurePoint(index,nodeIndex,localDrvGroup,localWeight,localDelay,localRefNodes));
        pressurePoints.push_back(PressurePoint(index,nodeIndex,localDrvGroup,localWeight,localDelay,localRefNodes));
    }
    return validLine;
}

bool PressurePointsSection::handleMeshFileDeclaration(const QString& currentLine){
    // line syntax:     <Index>    Mesh    MeshFileAlias=  Include <tags>   Exclude <tags>   SwapNormals

    std::cout<<"Mesh-file declaration of pressure points not yet implemented! Line: "<< currentLine.toUtf8().constData()<<std::endl;
//    QString currentLineCopy=currentLine;
//    QRegularExpressionMatch match;
//    QRegularExpression identifier("\\s+MeshFileAlias{1}\\s*={1}\\s*([\"\'][^\"\']+[\"\']|\\S+)", QRegularExpression::CaseInsensitiveOption);
//    QString meshFileAlias=MeshFileAlias;

//    int matchPosition=currentLine.indexOf(identifier, 0, &match);

//     if(matchPosition>=0) //Meshfile alias declaration in element declaration line found
//    {
//        meshFileAlias=match.captured(1);
//        global::removeQuotes(meshFileAlias);
//        std::cout <<"MeshFileAlias: "<< meshFileAlias.toUtf8().constData()<<std::endl;
//        currentLineCopy.remove(matchPosition,match.capturedLength() );
////        std::cout <<"Current line with meshfilealias removed: "<< currentLineCopy.toUtf8().constData()<<std::endl;


//    }

//    if(meshFileAlias.isEmpty())
//    {
//        validLine=false;
//        std::cout <<"No MeshFileAlias!: "<< currentLine.toUtf8().constData()<<std::endl;
//        return validLine;
//    }


//    // valid line syntax:     <Index>    Mesh  Include <tags>   Exclude <tags>   SwapNormals

//    QStringList list =global::stringToStringListQuotesIntact(currentLineCopy);



//    if(!list.isEmpty())
//    {
//        bool isNumber;
//        int elementIndex=list.first().toInt(&isNumber);
//        bool localSwapNormals=false;/*=swapNormals;*/
//        if(!isNumber)
//        {
//            validLine=false;
////            readingElementList=false;
//            return validLine;
//        }
//        if(!(currentLine.contains(QRegularExpression("\\s+include([\\s]|$)",Qt::CaseInsensitive)) || currentLine.contains(QRegularExpression("\\s+Mesh([\\s]|$)",Qt::CaseInsensitive)) || currentLine.contains(QRegularExpression("\\s+meshfileAlias\\s*=",Qt::CaseInsensitive))))
//        {
//            validLine=false;
////            readingElementList=false;
//            return validLine;
//        }
//        int includePos=list.indexOf(QRegularExpression("^\\s*include$", Qt::CaseInsensitive));
//        int excludePos=list.indexOf(QRegularExpression("^\\s*exclude$", Qt::CaseInsensitive));
//        for (int i=0; i<=list.length();i++){
//            std::cout << list.value(i).toUtf8().constData()<<std::endl;

//        }
////        std::cout<<std::endl <<"include Position: "<< includePos<<std::endl;
////        std::cout <<"exclude Position: "<< excludePos<<std::endl;

//        bool includeAll=false;
//        QStringList include;
//        QStringList exclude;

//        if( (includePos==-1) | (includePos==list.length()-1) )
//        {
//            includeAll=true;
//        }
//        else
//        {
//            for(int i=includePos+1;i<=list.length()-1;i++)
//            {
//                QString currentListElement=list.value(i);
//                if(currentListElement.contains(QRegularExpression("^\\s*MeshFileAlias",Qt::CaseInsensitive)))
//                {
//                    break;
//                }
//                if(currentListElement.contains(QRegularExpression("^\\s*exclude\\s*$",Qt::CaseInsensitive)))
//                {
//                    break;
//                }
//                if(currentListElement.contains(QRegularExpression("^\\s*swapNormals\\s*$",Qt::CaseInsensitive)))
//                {
//                    localSwapNormals=true;
//                    break;
//                }
//                if(currentListElement.contains(QRegularExpression("^\\s*all\\s*$",Qt::CaseInsensitive)))
//                {
//                    includeAll=true;
//                    continue;
//                }
////                if(0==currentListElement.compare("Mesh",Qt::CaseInsensitive))
////                {
////                    continue;
////                }
//                include.push_back(currentListElement);
//                std::cout <<"include: "<< currentListElement.toUtf8().constData()<<std::endl;

//            }
//        }
//        if(includeAll){std::cout <<"includeAll=true"<<std::endl;}

//        if(excludePos>=1 && (excludePos<=list.length()-2))
//        {
//            for(int i=excludePos+1;i<=list.length()-1;i++)
//            {
//                QString currentListElement=list.value(i);
//                if(currentListElement.contains(QRegularExpression("^\\s*MeshFileAlias",Qt::CaseInsensitive)))
//                {
//                    break;
//                }
//                if(currentListElement.contains(QRegularExpression("^\\s*include\\s*$",Qt::CaseInsensitive)))
//                {
//                    break;
//                }
//                if(currentListElement.contains(QRegularExpression("^\\s*swapNormals\\s*$",Qt::CaseInsensitive)))
//                {
//                    localSwapNormals=true;
//                    break;
//                }
////                if(0==currentListElement.compare("Mesh",Qt::CaseInsensitive))
////                {
////                    continue;
////                }

//                exclude.push_back(currentListElement);
//                std::cout <<"exclude: "<< currentListElement.toUtf8().constData()<<std::endl;

//            }

//        }
////        if(include.isEmpty())
////        {
////            includeAll=true;
////        }
////        meshFileElements.push_back(MeshFileElement(elementIndex,meshFileAlias,include,exclude,localSwapNormals,includeAll));
//        validLine=true;

//    }
//    else{
//        validLine=false;
//    }

    return validLine=true;

}
