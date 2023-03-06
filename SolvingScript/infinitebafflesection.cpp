#include "infinitebafflesection.h"

InfiniteBaffleSection::InfiniteBaffleSection()
{

}

bool InfiniteBaffleSection::handleScriptLine(const QString& currentLine)
{
    validLine = false;
    QRegularExpressionMatch match;
    QString value;
    QRegularExpression identifier;

    identifier=global::regExWithOneValue("Subdomain");
    if(currentLine.contains(identifier,&match))
    {
        value = match.captured(1);
        global::removeQuotesAndWhitespace(value);
        Subdomain = value.toUInt(&validLine);
        return validLine;
    }

    identifier = global::regExWithOneValue("Scale");
    if(currentLine.contains(identifier,&match))
    {
        value = match.captured(1);
        global::removeQuotes(value);
        bool alreadyCaptured = false;
        QRegularExpression separator( "[ \\t,]+");
        QStringList list = value.split(separator, Qt::SkipEmptyParts);
        if(list.length()==1)
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
        if(list.length()==3)
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
        if(Scale(0) == 0){Scale(0) = 1;}
        if(Scale(1) == 0){Scale(1) = 1;}
        if(Scale(2) == 0){Scale(2) = 1;}
        std::cout<<"scale: "<<Scale(0)<<", "<<Scale(1)<<", "<<Scale(2)<<std::endl;
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
        std::cout<<"Shift: "<<Shift(0)<<", "<<Shift(1)<<", "<<Shift(2)<<std::endl;
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
            QString captured1 = list.at(0);
//            double val =global::stringWithUnitsToDouble(captured1,QString("m"),validLine);
            double val = captured1.toDouble(&validLine);
            Rotate = {val,0,0};
            alreadyCaptured = true;
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

    identifier = global::regExWithOneValue("Position");
    if(currentLine.contains(identifier,&match))
    {
        value = match.captured(1);
        global::removeQuotes(value);
        QRegularExpression separator( "[ \\t,]+");
        QStringList list = value.split(separator, Qt::SkipEmptyParts);

        if(list.length() == 3)
        {
            bool valid1,valid2,valid3;
            QString captured1 = list.at(0);
            QString captured2 = list.at(1);
            QString captured3 = list.at(2);
            double val1 = global::stringWithUnitsToDouble(captured1,QString("m"),valid1);
            double val2 = global::stringWithUnitsToDouble(captured2,QString("m"),valid2);
            double val3 = global::stringWithUnitsToDouble(captured3,QString("m"),valid3);
            planePoint = {val1,val2,val3};
            validLine = (valid1 && valid2 && valid3);
        }
        else
        {
            validLine = false;
        }
        return validLine;
    }

    identifier = global::regExWithOneValue("Normal");
    if(currentLine.contains(identifier,&match))
    {
        value = match.captured(1);
        global::removeQuotes(value);
        QRegularExpression separator( "[ \\t,]+");
        QStringList list = value.split(separator, Qt::SkipEmptyParts);

        if(list.length() == 3)
        {
            bool valid1,valid2,valid3;
            QString captured1 = list.at(0);
            QString captured2 = list.at(1);
            QString captured3 = list.at(2);
            double val1 = captured1.toDouble(&valid1);
            double val2 = captured2.toDouble(&valid2);
            double val3 = captured3.toDouble(&valid3);
            planeNormal = {val1,val2,val3};
            validLine = (valid1 && valid2 && valid3);
        }
        else
        {
            validLine = false;
        }
        return validLine;
    }



//    identifier="RefNodes=";
    identifier = global::regExWithOneValue("RefNodes");
    if(currentLine.contains(identifier,&match))
    {
        value = match.captured(1);
        global::removeQuotesAndWhitespace(value);
        RefNodes = value;
        global::removeQuotes(RefNodes);
        if(RefNodes.isEmpty()){validLine = false;}
        else{validLine = true;}
        return validLine;
    }

//    identifier = "MeshFileAlias = ";
    identifier = global::regExWithOneValue("MeshFileAlias");
    if(currentLine.contains(identifier,&match))
    {
        value = match.captured(1);
        global::removeQuotesAndWhitespace(value);
        MeshFileAlias = value;
        global::removeQuotes(MeshFileAlias);
        if(MeshFileAlias.isEmpty()){validLine = false;}
        else{validLine = true;}
        return validLine;
    }

    identifier = global::regExWithOneValue("distanceToParallelPlane");
    if(currentLine.contains(identifier,&match))
    {
        value = match.captured(1);
        global::removeQuotesAndWhitespace(value);
        distanceToParallelPlane = value.toDouble(&validLine);
        if(distanceToParallelPlane <= 0)
        {
            validLine = false;
            std::cout <<"distanceToParallelPlane <= 0 in Infinite_Baffle section: "<< currentLine.toUtf8().constData()<<std::endl;
            QString error=QString("distanceToParallelPlane <= 0 in Infinite_Baffle declaration in line: "+ currentLine);
            logStrings::errorLogString.append(error+"\r\n");
        }
        return validLine;
    }
    return validLine;
}
