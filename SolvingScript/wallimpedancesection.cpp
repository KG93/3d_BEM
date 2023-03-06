#include "wallimpedancesection.h"

bool WallImpedanceSection::handleScriptLine(const QString& currentLine)
{
    validLine=false;
    if(!readingWallImpedancesList)
    {
        validLine=handleParameterSectionLine(currentLine);
    }
    else
    {
        validLine=handleWallImpedancesListLine(currentLine);
    }
    return validLine;
}

bool WallImpedanceSection::handleWallImpedancesListLine(const QString& currentLine)
{
    // valid line syntax:      <Index>    <Element-index>   RefElements=  DataFileAlias=  Value=

    QString currentLineCopy = currentLine;
    QRegularExpressionMatch match;
    //find RefElements
    QRegularExpression identifier("\\s+RefElements{1}\\s*={1}\\s*([\"\'][^\"\']+[\"\']|\\S+)", QRegularExpression::CaseInsensitiveOption);
    QString localRefElements = RefElements;

    int matchPosition = currentLineCopy.indexOf(identifier, 0, &match);

     if(matchPosition >= 0) //RefElements alias declaration in wallimpedance line found
    {
        localRefElements=match.captured(1);
        global::removeQuotes(localRefElements);
//        std::cout <<"RefElements: "<< localRefElements.toUtf8().constData()<<std::endl;
        currentLineCopy.remove(matchPosition,match.capturedLength() );
//        std::cout <<"Current line with RefElements removed: "<< currentLineCopy.toUtf8().constData()<<std::endl;
    }

    if(localRefElements.isEmpty())
    {
        validLine = false;
        std::cout << "No RefElements!: " << currentLine.toUtf8().constData() << std::endl;
        QString error = QString("No RefElements declaration in line: "+currentLine);
        logStrings::errorLogString.append(error+"\r\n");
        return validLine;
    }

    //find DataFileAlias
    identifier = QRegularExpression("\\s+DataFileAlias{1}\\s*={1}\\s*([\"\'][^\"\']+[\"\']|\\S+)", QRegularExpression::CaseInsensitiveOption);
    QString localDataFileAlias = DataFileAlias;

    matchPosition = currentLineCopy.indexOf(identifier, 0, &match);

     if(matchPosition >= 0) //DataFileAlias alias declaration in wallimpedance line found
    {
        localDataFileAlias = match.captured(1);
        global::removeQuotes(localDataFileAlias);
        std::cout <<"DataFileAlias: " << localDataFileAlias.toUtf8().constData() << std::endl;
        currentLineCopy.remove(matchPosition,match.capturedLength() );
//        std::cout <<"Current line with DataFileAlias removed: "<< currentLineCopy.toUtf8().constData()<<std::endl;
    }

    if(localDataFileAlias.isEmpty())
    {
        validLine = false;
        std::cout <<"No DataFileAlias!: "<< currentLine.toUtf8().constData()<<std::endl;
    }

    //find Value
//    identifier=QRegularExpression("\\s+Value{1}\\s*={1}\\s*([\"\'][^\"\']+[\"\']|\\S+)", QRegularExpression::CaseInsensitiveOption);
    identifier = QRegularExpression("\\s+Value{1}\\s*={1}\\s*([\"\'][^\"\']+[\"\']|.+$)", QRegularExpression::CaseInsensitiveOption);
    QString captureForValue;
    std::complex<double> localValue = Value;
    matchPosition = currentLineCopy.indexOf(identifier, 0, &match);

     if(matchPosition >= 0) //Value declaration in wallimpedance line found
    {
        captureForValue = match.captured(1);
        global::removeQuotes(captureForValue);
        QRegularExpression separator( "[ \\t,]+");
        QStringList list = captureForValue.split(separator, Qt::SkipEmptyParts);
        bool valid1 = true;
        bool valid2 = true;

//        double realPart=global::stringWithUnitsToDouble(list.at(0),QString("Hz"),valid1);
        double realPart = list.at(0).toDouble(&valid1);
        double imaginaryPart = 0;
        if(list.length() == 2)
        {
//            imaginaryPart=global::stringWithUnitsToDouble(list.at(1),QString("Hz"),valid2);
            imaginaryPart = list.at(1).toDouble(&valid2);
        }
        validLine = valid1&&valid2;
        localValue = std::complex<double>(realPart,imaginaryPart);
        if(localValue == std::complex<double>(0,0)){validLine = false;}
        currentLineCopy.remove(matchPosition,match.capturedLength() );
    }

//    if(localValue==std::complex<double>(0,0))
//    {
//        validLine=false;
//        std::cout <<"Value=(0,0)!: "<< currentLine.toUtf8().constData()<<std::endl;
//        QString error=QString("Value=(0,0) in line: "+currentLine);
//        logStrings::errorLogString.append(error+"\r\n");
//        return validLine;
//    }

    QStringList list =global::stringToStringListQuotesIntact(currentLineCopy);

    if(!list.isEmpty())
    {
        bool isNumber;
        int index = list.first().toInt(&isNumber);
        if(!isNumber)
        {
            validLine = false;
            std::cout << "Invalid Index: " << currentLine.toUtf8().constData() << std::endl;
            QString error = QString("Invalid index in line: "+currentLine);
            logStrings::errorLogString.append(error+"\r\n");
            return validLine;
        }
        if(list.length() >= 2)
        {
            int elementIndex = list.value(1).toInt(&isNumber);
            if(!isNumber)
            {
                validLine = false;
                std::cout <<"Invalid Element Index: "<< currentLine.toUtf8().constData()<<std::endl;
                QString error = QString("Invalid element index in line: "+currentLine);
                logStrings::errorLogString.append(error+"\r\n");
                return validLine;
            }
            validLine = true;
            wallImpecanceElements.push_back(WallImpecanceElement(index, elementIndex, localRefElements, localDataFileAlias, localValue));
        }
        else
        {
            validLine = false;
            std::cout << "No Element Index: " << currentLine.toUtf8().constData() << std::endl;
            QString error = QString("No element index in line: "+currentLine);
            logStrings::errorLogString.append(error+"\r\n");
            return validLine;
        }

    }
    return validLine;
}

bool WallImpedanceSection::handleParameterSectionLine(const QString& currentLine)
{
    validLine = false;
    QRegularExpressionMatch match;
    QString capturedValue;
    QRegularExpression identifier;

    identifier = global::regExWithOneValue("ImpType");
    if(currentLine.contains(identifier,&match))
    {
        capturedValue = match.captured(1);
        global::removeQuotesAndWhitespace(capturedValue);
        if(0==QString::compare("Reflection",capturedValue, Qt::CaseInsensitive))
        {
              validLine = true;
              ImpType = Reflection;
              return validLine;
        }
        else if(0==QString::compare("Damping",capturedValue, Qt::CaseInsensitive))
        {
            validLine = true;
              ImpType = Damping;
              return validLine;
        }
        else if(0==QString::compare("Impedance",capturedValue, Qt::CaseInsensitive))
        {
              validLine = true;
              ImpType = Impedance;
              return validLine;
        }
        else if(0==QString::compare("Admittance",capturedValue, Qt::CaseInsensitive))
        {
            validLine = true;
              ImpType = Admittance;
              return validLine;
        }
        else
        {
            QString message = QString("Unrecognizable ImpType in line: "+currentLine);
            logStrings::logString.append(message+"\r\n");
            logStrings::errorLogString.append(message+"\r\n");
            std::cout << message.toUtf8().constData() << std::endl;
            validLine = false;
        }
        return validLine;
    }

    identifier = global::regExWithOptionalValue("Normalized");
    if(currentLine.contains(identifier,&match))
    {
//        std::cout<<"Normalized!!!"<<std::endl;
          capturedValue = match.captured(2);
          std::cout << capturedValue.toUtf8().constData() << std::endl;
          global::removeQuotesAndWhitespace(capturedValue);
          if(!capturedValue.isEmpty())
          {
              if(capturedValue.compare("true",Qt::CaseInsensitive)==0 || capturedValue.toInt()==1)
              {
                  Normalized = true;
              }
              if(capturedValue.compare("false",Qt::CaseInsensitive)==0 || capturedValue.toInt()==0)
              {
                  Normalized = false;
              }
          }
          else
          {
              Normalized = true;
          }
          validLine = true;
          return validLine;
    }

    identifier = global::regExWithOneValue("RefElements");
    if(currentLine.contains(identifier,&match))
    {
        capturedValue = match.captured(1);
        global::removeQuotesAndWhitespace(capturedValue);
        RefElements = capturedValue;
        global::removeQuotes(RefElements);
        if(RefElements.isEmpty()){validLine = false;}
        else{validLine = true;}
        return validLine;
    }

    identifier = global::regExWithOneValue("Value");
    if(currentLine.contains(identifier,&match))
    {
        capturedValue = match.captured(1);
//        global::removeQuotesAndWhitespace(capturedValue);
        global::removeQuotes(capturedValue);

        QRegularExpression separator( "[ \\t,]+");
        QStringList list = capturedValue.split(separator, Qt::SkipEmptyParts);
        bool valid1 = true;
        bool valid2 = true;
//        double realPart=global::stringWithUnitsToDouble(list.at(0),QString("Hz"),valid1);
        double realPart = list.at(0).toDouble(&valid1);
        double imaginaryPart = 0;
        if(list.length()==2)
        {
//            imaginaryPart=global::stringWithUnitsToDouble(list.at(1),QString("Hz"),valid2);
            imaginaryPart = list.at(1).toDouble(&valid2);
        }
        validLine = valid1&&valid2;
        Value = std::complex<double>(realPart,imaginaryPart);

//        if(Value==std::complex<double>(0,0)){validLine=false;}
        return validLine;
    }

    identifier = global::regExWithOneValue("ft");
    if(currentLine.contains(identifier,&match))
    {
        capturedValue = match.captured(1);
        global::removeQuotesAndWhitespace(capturedValue);
//        MeshFrequency=value.toDouble(&validLine);
        ft = global::stringWithUnitsToDouble(capturedValue,QString("Hz"),validLine);
        if(ft<=0){validLine=false;}
        return validLine;
    }

    identifier = global::regExWithOneValue("DataFileAlias");
    if(currentLine.contains(identifier,&match))
    {
        capturedValue = match.captured(1);
        global::removeQuotesAndWhitespace(capturedValue);
        DataFileAlias = capturedValue;
        global::removeQuotes(DataFileAlias);
        if(DataFileAlias.isEmpty()){validLine = false;}
        else{validLine = true;}
        return validLine;
    }

    identifier = global::regExWithOneValue("Color");
    if(currentLine.contains(identifier,&match))
    {
        validLine = true;
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
            readingWallImpedancesList=true;
            validLine = handleWallImpedancesListLine(currentLine);
            return validLine;
        }
    }
    QString error = QString("Line in Wallimpedance section section not recognized! "+ currentLine);
    logStrings::errorLogString.append(error+"\r\n");
    return validLine;
}

 void  WallImpedanceSection::setUpValues()
{    
    if(mediumImpedance==0.0)
    {
        std::cout<<"Impedance of the medium not set in Wallimpedance section."<<std::endl;
        QString error=QString("Impedance of the medium not set in Wallimpedance section.");
        logStrings::errorLogString.append(error+"\r\n");
    }
    if(ImpType==Reflection)
    {
        double alpha = (1-std::real(Value))/(1+std::real(Value));
        robinBoundaryCondition.a = imaginaryUnit*wavenumber*alpha;
        robinBoundaryCondition.b = 1;
        robinBoundaryCondition.g = 0;

        for(int i=0; i<wallImpecanceElements.size(); i++)
        {
            double currentVal = std::real(wallImpecanceElements.at(i).Value);
            double alpha = (1-currentVal)/(1+currentVal);
            wallImpecanceElements[i].robinBoundaryCondition.a = imaginaryUnit*wavenumber*alpha;
            wallImpecanceElements[i].robinBoundaryCondition.b = 1;
            wallImpecanceElements[i].robinBoundaryCondition.g = 0;
        }
    }
    else if(ImpType==Damping)
    {
        double d = std::real(Value);
        double alpha = d/(2-d);
        robinBoundaryCondition.a = imaginaryUnit*wavenumber*alpha;
        robinBoundaryCondition.b = 1;
        robinBoundaryCondition.g = 0;

        for(int i=0; i<wallImpecanceElements.size(); i++)
        {
            double currentVal = std::real(wallImpecanceElements.at(i).Value);
            double d = std::real(currentVal);
            double alpha = d/(2-d);
            wallImpecanceElements[i].robinBoundaryCondition.a = imaginaryUnit*wavenumber*alpha;
            wallImpecanceElements[i].robinBoundaryCondition.b = 1;
            wallImpecanceElements[i].robinBoundaryCondition.g = 0;
        }
    }
    else if(ImpType==Impedance)
    {
        if(Normalized==true)
        {
            robinBoundaryCondition.a = imaginaryUnit*wavenumber*(1.0/Value);
            robinBoundaryCondition.b = 1;
            robinBoundaryCondition.g = 0;
        }
        else
        {
            std::complex<double> alpha = mediumImpedance/(Value);
            robinBoundaryCondition.a = imaginaryUnit*wavenumber*(/*1.0/*/alpha);
            robinBoundaryCondition.b = 1;
            robinBoundaryCondition.g = 0;
        }
        for(int i=0; i<wallImpecanceElements.size(); i++)
        {
            if(Normalized==true)
            {
                wallImpecanceElements[i].robinBoundaryCondition.a = imaginaryUnit*wavenumber*(1.0/wallImpecanceElements[i].Value);
                wallImpecanceElements[i].robinBoundaryCondition.b = 1;
                wallImpecanceElements[i].robinBoundaryCondition.g = 0;
            }
            else
            {
                std::complex<double> alpha = mediumImpedance/wallImpecanceElements[i].Value;
                wallImpecanceElements[i].robinBoundaryCondition.a = imaginaryUnit*wavenumber*(/*1.0/*/alpha);
                wallImpecanceElements[i].robinBoundaryCondition.b = 1;
                wallImpecanceElements[i].robinBoundaryCondition.g = 0;
            }
        }
    }
    else if(ImpType==Admittance)
    {
        if(Normalized==true)
        {
            robinBoundaryCondition.a = imaginaryUnit*wavenumber*(Value);
            robinBoundaryCondition.b = 1;
            robinBoundaryCondition.g = 0;
        }
        else
        {
            std::complex<double> alpha = mediumImpedance*Value;
            robinBoundaryCondition.a = imaginaryUnit*wavenumber*(alpha);
            robinBoundaryCondition.b = 1;
            robinBoundaryCondition.g = 0;
        }
        for(int i=0; i<wallImpecanceElements.size(); i++)
        {
            if(Normalized==true)
            {
                wallImpecanceElements[i].robinBoundaryCondition.a = imaginaryUnit*wavenumber*(wallImpecanceElements[i].Value);
                wallImpecanceElements[i].robinBoundaryCondition.b = 1;
                wallImpecanceElements[i].robinBoundaryCondition.g = 0;
            }
            else
            {
                std::complex<double> alpha = mediumImpedance*wallImpecanceElements[i].Value;
                wallImpecanceElements[i].robinBoundaryCondition.a = imaginaryUnit*wavenumber*(alpha);
                wallImpecanceElements[i].robinBoundaryCondition.b = 1;
                wallImpecanceElements[i].robinBoundaryCondition.g = 0;
            }
        }
    }
//    std::cerr << "robinBoundaryCondition: " << robinBoundaryCondition.a << " " << robinBoundaryCondition.b << " " << robinBoundaryCondition.g << std::endl;
}
