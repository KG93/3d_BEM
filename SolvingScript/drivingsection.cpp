#include "drivingsection.h"

DrivingSection::DrivingSection()
{

}

bool DrivingSection::handleScriptLine(const QString& currentLine)
{

    validLine=false;
    if(!readingDrivingElementsList)
    {
        validLine=handleParameterSectionLine(currentLine);
    }
    else
    {
        validLine=handleDrivingElementsListLine(currentLine);
    }
    return validLine;
}

bool DrivingSection::handleDrivingElementsListLine(const QString& currentLine)
{
    // valid line syntax:      <Index>    <Element-index>   RefElements=  DataFileAlias=  Value=

    QString currentLineCopy=currentLine;
    QRegularExpressionMatch match;
    //find RefElements
    QRegularExpression identifier("\\s+RefElements{1}\\s*={1}\\s*([\"\'][^\"\']+[\"\']|\\S+)", QRegularExpression::CaseInsensitiveOption);
    QString localRefElements=RefElements;

    int matchPosition=currentLineCopy.indexOf(identifier, 0, &match);

     if(matchPosition>=0) //RefElements alias declaration in driving declaration line found
    {
        localRefElements=match.captured(1);
        global::removeQuotes(localRefElements);
//        std::cout <<"RefElements: "<< localRefElements.toUtf8().constData()<<std::endl;
        currentLineCopy.remove(matchPosition,match.capturedLength() );
//        std::cout <<"Current line with RefElements removed: "<< currentLineCopy.toUtf8().constData()<<std::endl;
    }

    if(localRefElements.isEmpty())
    {
        validLine=false;
        std::cout <<"No RefElements!: "<< currentLine.toUtf8().constData()<<std::endl;
        QString error=QString("No RefElements declaration in line: "+currentLine);
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
        localDelay += global::stringWithUnitsToDouble(value,QString("s"),isNumber); //add extra delay for element
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
    std::complex<double> localWeight=DrvWeight;
    isNumber=true;
    matchPosition=currentLineCopy.indexOf(identifier, 0, &match);
    if(matchPosition>=0) //DrvWeight declaration in pressure point declaration line found
    {
        QString value=match.captured(1);
        global::removeQuotesAndWhitespace(value);
        localWeight=value.toDouble(&isNumber);
        std::cout <<"DrvWeight: "<< localWeight<<std::endl;
//        if(DrvWeight<=0){std::cout<<"DrvWeight<=0 in Line \""<<currentLine.toUtf8().constData()<<"\" of Pressure_Points section "<< name.toUtf8().constData()<<"."<<std::endl;}
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

    if(!list.isEmpty())
    {
        bool isNumber;
        int index=list.first().toInt(&isNumber);
        if(!isNumber)
        {
            validLine=false;
            std::cout <<"Invalid Index: "<< currentLine.toUtf8().constData()<<std::endl;
            QString error=QString("Invalid index in line: "+currentLine);
            logStrings::errorLogString.append(error+"\r\n");
            return validLine;
        }
        if(list.length()>=2)
        {
            int elementIndex=list.value(1).toInt(&isNumber);
            if(!isNumber)
            {
                validLine=false;
                std::cout <<"Invalid Element Index: "<< currentLine.toUtf8().constData()<<std::endl;
                QString error=QString("Invalid element index in line: "+currentLine);
                logStrings::errorLogString.append(error+"\r\n");
                return validLine;
            }
            validLine=true;
            drivingElements.push_back(DrivingElement(index,elementIndex,localRefElements,localDrvGroup,localWeight,localDelay));
        }
        else
        {
            validLine=false;
            std::cout <<"No Element Index: "<< currentLine.toUtf8().constData()<<std::endl;
            QString error=QString("No element index in line: "+currentLine);
            logStrings::errorLogString.append(error+"\r\n");
            return validLine;
        }

    }
    return validLine;

}

bool DrivingSection::handleParameterSectionLine(const QString& currentLine)
{
    validLine=false;
    QRegularExpressionMatch match;
    QString capturedValue;
    QRegularExpression identifier;

//    identifier=global::regExWithOneValue("ImpType");
//    if(currentLine.contains(identifier,&match))
//    {
//        capturedValue=match.captured(1);
//        global::removeQuotesAndWhitespace(capturedValue);
//        if(0==QString::compare("Reflection",capturedValue, Qt::CaseInsensitive))
//        {
//              validLine=true;
//              ImpType=Reflection;
//              return validLine;
//        }
//        else if(0==QString::compare("Damping",capturedValue, Qt::CaseInsensitive))
//        {
//            validLine=true;
//              ImpType=Damping;
//              return validLine;
//        }
//        else if(0==QString::compare("Impedance",capturedValue, Qt::CaseInsensitive))
//        {
//              validLine=true;
//              ImpType=Impedance;
//              return validLine;
//        }
//        else if(0==QString::compare("Admittance",capturedValue, Qt::CaseInsensitive))
//        {
//            validLine=true;
//              ImpType=Admittance;
//              return validLine;
//        }
//        else
//        {
//            QString message=QString("Unrecognizable ImpType in line: "+currentLine);
//            logStrings::logString.append(message+"\r\n");
//            logStrings::errorLogString.append(message+"\r\n");
//            std::cout<<message.toUtf8().constData()<<std::endl;
//            validLine=false;
//        }
//        return validLine;
//    }

    identifier=global::regExWithOptionalValue("DirIsAxial");
    if(currentLine.contains(identifier,&match))
    {
//        std::cout<<"DirIsAxial!!!"<<std::endl;
          capturedValue=match.captured(2);
          std::cout<<capturedValue.toUtf8().constData()<<std::endl;
          global::removeQuotesAndWhitespace(capturedValue);
          if(!capturedValue.isEmpty())
          {
              if(capturedValue.compare("true",Qt::CaseInsensitive)==0||capturedValue.toInt()==1)
              {
                  DirIsAxial=true;
              }
              if(capturedValue.compare("false",Qt::CaseInsensitive)==0||capturedValue.toInt()==0)
              {
                  DirIsAxial=false;
              }
          }
          else
          {
              DirIsAxial=true;
          }
          validLine=true;
          return validLine;
    }

    identifier=global::regExWithOneValue("RefElements");
    if(currentLine.contains(identifier,&match))
    {
        capturedValue=match.captured(1);
        global::removeQuotesAndWhitespace(capturedValue);
        RefElements=capturedValue;
        global::removeQuotes(RefElements);
        if(RefElements.isEmpty())
        {
            std::cout<<"Empty RefElements in Line \""<<currentLine.toUtf8().constData()<<"\" of Driving section "<< name.toUtf8().constData()<<"."<<std::endl;
            QString error=QString("Empty RefElements in Line \""+currentLine+"\" of Driving section "+ name+".");
            logStrings::errorLogString.append(error+"\r\n");
            validLine=false;
        }
        else{validLine=true;}
        return validLine;
    }

    identifier=global::regExWithOneValue("RefNodes");
    if(currentLine.contains(identifier,&match))
    {
        capturedValue=match.captured(1);
        global::removeQuotesAndWhitespace(capturedValue);
        RefNodes=capturedValue;
        global::removeQuotes(RefNodes);
        if(RefNodes.isEmpty())
        {
            std::cout<<"Empty RefNodes in Line \""<<currentLine.toUtf8().constData()<<"\" of Driving section "<< name.toUtf8().constData()<<"."<<std::endl;
            QString error=QString("Empty RefNodes in Line \""+currentLine+"\" of Driving section "+ name+".");
            logStrings::errorLogString.append(error+"\r\n");
            validLine=false;
        }
        else{validLine=true;}
        return validLine;
    }

    identifier=global::regExWithOneValue("MeshFileAlias");
    if(currentLine.contains(identifier,&match))
    {
        capturedValue=match.captured(1);
        global::removeQuotesAndWhitespace(capturedValue);
        MeshFileAlias=capturedValue;
        global::removeQuotes(MeshFileAlias);
        if(MeshFileAlias.isEmpty())
        {
            std::cout<<"Empty MeshFileAlias in Line \""<<currentLine.toUtf8().constData()<<"\" of Driving section "<< name.toUtf8().constData()<<"."<<std::endl;
            QString error=QString("Empty MeshFileAlias in Line \""+currentLine+"\" of Driving section "+ name+".");
            logStrings::errorLogString.append(error+"\r\n");
            validLine=false;
        }
        else{validLine=true;}
        return validLine;
    }

    identifier=global::regExWithOneValue("DrvGroup");
    if(currentLine.contains(identifier,&match))
    {
        capturedValue=match.captured(1);
        global::removeQuotesAndWhitespace(capturedValue);
        DrvGroup=capturedValue.toUInt(&validLine);
        if(!validLine)
        {
            std::cout<<"Invalid DrvGroup in Line \""<<currentLine.toUtf8().constData()<<"\" of Driving section "<< name.toUtf8().constData()<<"."<<std::endl;
            QString error=QString("Invalid DrvGroup in Line \""+currentLine+"\" of Driving section "+ name+".");
            logStrings::errorLogString.append(error+"\r\n");
        }
        return validLine;
    }

    identifier=global::regExWithOneValue("DrvWeight");
    if(currentLine.contains(identifier,&match))
    {
        capturedValue=match.captured(1);
        global::removeQuotesAndWhitespace(capturedValue);
        DrvWeight=capturedValue.toDouble(&validLine);
        if(validLine)
        {
//            if(DrvWeight<=0)
//            {
//                std::cout<<"DrvWeight<=0 in Line \""<<currentLine.toUtf8().constData()<<"\" of Driving section "<< name.toUtf8().constData()<<"."<<std::endl;
//                QString error=QString("DrvWeight<=0 in Line \""+currentLine+"\" of Driving section "+ name+".");
//                logStrings::errorLogString.append(error+"\r\n");
//            }
        }
        else
        {
            std::cout<<"Invalid DrvWeight in Line \""<<currentLine.toUtf8().constData()<<"\" of Driving section "<< name.toUtf8().constData()<<"."<<std::endl;
            QString error=QString("Invalid DrvWeight in Line \""+currentLine+"\" of Driving section "+ name+".");
            logStrings::errorLogString.append(error+"\r\n");
        }
        return validLine;
    }

    identifier=global::regExWithOneValue("DrvDelay");
    if(currentLine.contains(identifier,&match))
    {
        capturedValue=match.captured(1);
        global::removeQuotesAndWhitespace(capturedValue);
        DrvDelay=global::stringWithUnitsToDouble(capturedValue,QString("s"),validLine);
        return validLine;
    }

//    identifier=global::regExWithOneValue("Value");
//    if(currentLine.contains(identifier,&match))
//    {
//        capturedValue=match.captured(1);
//        global::removeQuotesAndWhitespace(capturedValue);
//        QRegularExpression separator( "[ \\t,]+");
//        QStringList list =capturedValue.split(separator, Qt::SkipEmptyParts);
//        bool valid1=true;
//        bool valid2=true;
//        double realPart=global::stringWithUnitsToDouble(list.at(0),QString("Hz"),valid1);
//        double imaginaryPart;

//        if(list.length()==2)
//        {
//            imaginaryPart=global::stringWithUnitsToDouble(list.at(1),QString("Hz"),valid2);
//        }
//        validLine=valid1&&valid2;
//        Value=std::complex<double>(realPart,imaginaryPart);
//        if(Value==std::complex<double>(0,0)){validLine=false;}
//        return validLine;
//    }

//    identifier=global::regExWithOneValue("ft");
//    if(currentLine.contains(identifier,&match))
//    {
//        capturedValue=match.captured(1);
//        global::removeQuotesAndWhitespace(capturedValue);
////        MeshFrequency=value.toDouble(&validLine);
//        ft=global::stringWithUnitsToDouble(capturedValue,QString("Hz"),validLine);
//        if(ft<=0){validLine=false;}
//        return validLine;
//    }



    identifier=global::regExWithOneValue("Color");
    if(currentLine.contains(identifier,&match))
    {
        validLine=true;
        return validLine;
    }


    QRegularExpression separator( "[ \\t,]+"); //identify any number of white spaces and tabs intermixed or commas
    QStringList list =currentLine.split(separator, Qt::SkipEmptyParts);
    if(!list.isEmpty())
    {
        bool isNumber;
        list.first().toInt(&isNumber);
        if(isNumber)
        {
            readingDrivingElementsList=true;
            validLine=handleDrivingElementsListLine(currentLine);
            return validLine;
        }
    }
    QString error=QString("Line in Driving section section not recognized! "+ currentLine);
    logStrings::errorLogString.append(error+"\r\n");
    return validLine;
}

//bool DrivingSection::setUpValues()
//{
//    if(mediumImpedance==0.0)
//    {
//        std::cout<<"Impedance of the medium not set in WallimpedanceSection."<<std::endl;
//        QString error=QString("Impedance of the medium not set in WallimpedanceSection.");
//        logStrings::errorLogString.append(error+"\r\n");
//    }
//    if(ImpType==Reflection)
//    {
//        double alpha=(1-std::real(Value))/(1+std::real(Value));
//        robinBoundaryCondition.a=imaginaryUnit*wavenumber*alpha;
//        robinBoundaryCondition.b=1;
//        robinBoundaryCondition.g=0;

//        for(int i=0;i<wallImpecanceElements.size();i++)
//        {
//            double currentVal=std::real(wallImpecanceElements.at(i).Value);
//            double alpha=(1-currentVal)/(1+currentVal);
//            wallImpecanceElements[i].robinBoundaryCondition.a=imaginaryUnit*wavenumber*alpha;
//            wallImpecanceElements[i].robinBoundaryCondition.b=1;
//            wallImpecanceElements[i].robinBoundaryCondition.g=0;
//        }
//    }
//    else if(ImpType==Damping)
//    {
//        double d=std::real(Value);
//        double alpha=d/(2-d);
//        robinBoundaryCondition.a=imaginaryUnit*wavenumber*alpha;
//        robinBoundaryCondition.b=1;
//        robinBoundaryCondition.g=0;

//        for(int i=0;i<wallImpecanceElements.size();i++)
//        {
//            double currentVal=std::real(wallImpecanceElements.at(i).Value);
//            double d=std::real(currentVal);
//            double alpha=d/(2-d);
//            wallImpecanceElements[i].robinBoundaryCondition.a=imaginaryUnit*wavenumber*alpha;
//            wallImpecanceElements[i].robinBoundaryCondition.b=1;
//            wallImpecanceElements[i].robinBoundaryCondition.g=0;
//        }
//    }
//    else if(ImpType==Impedance)
//    {
//        if(Normalized==true)
//        {
//            robinBoundaryCondition.a=imaginaryUnit*wavenumber*(1.0/Value);
//            robinBoundaryCondition.b=1;
//            robinBoundaryCondition.g=0;
//        }
//        else
//        {
//            std::complex<double> alpha=mediumImpedance/(Value);
//            robinBoundaryCondition.a=imaginaryUnit*wavenumber*(1.0/alpha);
//            robinBoundaryCondition.b=1;
//            robinBoundaryCondition.g=0;
//        }
//        for(int i=0;i<wallImpecanceElements.size();i++)
//        {
//            if(Normalized==true)
//            {
//                wallImpecanceElements[i].robinBoundaryCondition.a=imaginaryUnit*wavenumber*(1.0/wallImpecanceElements[i].Value);
//                wallImpecanceElements[i].robinBoundaryCondition.b=1;
//                wallImpecanceElements[i].robinBoundaryCondition.g=0;
//            }
//            else
//            {
//                std::complex<double> alpha=mediumImpedance/wallImpecanceElements[i].Value;
//                wallImpecanceElements[i].robinBoundaryCondition.a=imaginaryUnit*wavenumber*(1.0/alpha);
//                wallImpecanceElements[i].robinBoundaryCondition.b=1;
//                wallImpecanceElements[i].robinBoundaryCondition.g=0;
//            }
//        }
//    }
//    else if(ImpType==Admittance)
//    {
//        if(Normalized==true)
//        {
//            robinBoundaryCondition.a=imaginaryUnit*wavenumber*(Value);
//            robinBoundaryCondition.b=1;
//            robinBoundaryCondition.g=0;
//        }
//        else
//        {
//            std::complex<double> alpha=mediumImpedance/Value;
//            robinBoundaryCondition.a=imaginaryUnit*wavenumber*(alpha);
//            robinBoundaryCondition.b=1;
//            robinBoundaryCondition.g=0;
//        }
//        for(int i=0;i<wallImpecanceElements.size();i++)
//        {
//            if(Normalized==true)
//            {
//                wallImpecanceElements[i].robinBoundaryCondition.a=imaginaryUnit*wavenumber*(wallImpecanceElements[i].Value);
//                wallImpecanceElements[i].robinBoundaryCondition.b=1;
//                wallImpecanceElements[i].robinBoundaryCondition.g=0;
//            }
//            else
//            {
//                std::complex<double> alpha=mediumImpedance/wallImpecanceElements[i].Value;
//                wallImpecanceElements[i].robinBoundaryCondition.a=imaginaryUnit*wavenumber*(alpha);
//                wallImpecanceElements[i].robinBoundaryCondition.b=1;
//                wallImpecanceElements[i].robinBoundaryCondition.g=0;
//            }
//        }
//    }
//}
