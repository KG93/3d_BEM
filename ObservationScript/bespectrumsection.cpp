#include "bespectrumsection.h"

BESpectrumSection::BESpectrumSection()
{

}

bool BESpectrumSection::handleScriptLine(const QString& currentLine)
{

    validLine=false;
    if(!readingSpectrumItemList)
    {
        validLine=handleParameterSectionLine(currentLine);
    }
    else
    {
        validLine=handleSpectrumItemListLine(currentLine);
    }
    return validLine;
}

bool BESpectrumSection::handleSpectrumItemListLine(const QString& currentLine)
{
    validLine=false;
    QString currentLineCopy=currentLine;
    QRegularExpressionMatch match;
    //find RefElements
    QRegularExpression identifier("\\s+RefNodes{1}\\s*={1}\\s*([\"\'][^\"\']+[\"\']|\\S+)", QRegularExpression::CaseInsensitiveOption);
    QString localRefNodes=RefNodes;

    int matchPosition=currentLineCopy.indexOf(identifier, 0, &match);
    if(matchPosition>=0) //RefNodes declaration in BE_Spectrum item declaration line found
    {
        localRefNodes=match.captured(1);
        global::removeQuotes(localRefNodes);
        std::cout <<"RefElements: "<< localRefNodes.toUtf8().constData()<<std::endl;
        currentLineCopy.remove(matchPosition,match.capturedLength() );
        std::cout <<"Current line with RefNodes removed: "<< currentLineCopy.toUtf8().constData()<<std::endl;

    }

    if(localRefNodes.isEmpty())
    {
        validLine=false;
        std::cout <<"No RefElements!: "<< currentLine.toUtf8().constData()<<std::endl;
        logStrings::errorLogString.append("No RefElements in line: "+ currentLine+".\r\n");
        return validLine;
    }

    //find ID
    identifier=QRegularExpression("\\s+ID{1}\\s*={1}\\s*([\"\'][^\"\']+[\"\']|\\S+)", QRegularExpression::CaseInsensitiveOption);
    QString ID;
    matchPosition=currentLineCopy.indexOf(identifier, 0, &match);
    if(matchPosition>=0) //DrvGroups declaration in BE_Spectrum item declaration line found
    {
        ID=match.captured(1);
        global::removeQuotes(localRefNodes);
        std::cout <<"ID: "<< ID.toUtf8().constData()<<std::endl;
        currentLineCopy.remove(matchPosition,match.capturedLength() );
        std::cout <<"Current line with DrvGroups removed: "<< currentLineCopy.toUtf8().constData()<<std::endl;

    }
    QStringList list =global::stringToStringListQuotesIntact(currentLineCopy);
    if(!list.isEmpty())
    {
        bool isNumber;
        int itemIndex=list.first().toInt(&isNumber);
        if(!isNumber)
        {
            validLine=false;
            return validLine;
        }
        if(list.size()>=2)
        {
            int nodeIndex=list.at(1).toInt(&isNumber);
            if(!isNumber)
            {
                validLine=false;
                std::cout<<"Invalid NodeIndex in BE_Spectrum section "<<name.toUtf8().constData()<<". Line: "<<currentLineCopy.toUtf8().constData()<<std::endl;
                logStrings::errorLogString.append("Invalid NodeIndex in BE_Spectrum section "+name+". Line: "+ currentLine+".\r\n");
                return validLine;
            }
            validLine=true;
            bESpectrumItems.push_back(BESpectrumItem(itemIndex,nodeIndex,localRefNodes,ID));
        }
        else
        {
            validLine=false;
            std::cout<<"No Node declared in BE_Spectrum section "<<name.toUtf8().constData()<<". Line: "<<currentLineCopy.toUtf8().constData()<<std::endl;
            logStrings::errorLogString.append("No Node declared in BE_Spectrum section "+name+". Line: "+ currentLine+".\r\n");
            return validLine;
        }
    }
    return validLine;
}

bool BESpectrumSection::handleParameterSectionLine(const QString& currentLine)
{
    validLine=false;
    QRegularExpressionMatch match;
    QString value;
    QRegularExpression identifier;

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

//    identifier=global::regExWithOneValue("scale");
//    if(currentLine.contains(identifier,&match))
//    {
//        value=match.captured(1);
//        global::removeQuotes(value);
//        bool alreadyCaptured=false;
//        QRegularExpression separator( "[ \\t,]+");
//        QStringList list =value.split(separator, Qt::SkipEmptyParts);
//        if(list.length()==1)
//        {
//            QString captured1=list.at(0);
//            double val =global::stringWithUnitsToDouble(captured1,QString("m"),validLine);
//            Scale={val,val,val};
//            alreadyCaptured=true;
//        }
//        if(list.length()==2)
//        {
//            bool valid1,valid2;
//            QString captured1=list.at(0);
//            QString captured2=list.at(1);
//            double val1 =global::stringWithUnitsToDouble(captured1,QString("m"),valid1);
//            double val2 =global::stringWithUnitsToDouble(captured2,QString("m"),valid2);
//            Scale={val1,val2,1};
//            validLine=(valid1 && valid2);
//            alreadyCaptured=true;
//        }
//        if(list.length()==3)
//        {
//            bool valid1,valid2,valid3;
//            QString captured1=list.at(0);
//            QString captured2=list.at(1);
//            QString captured3=list.at(2);
//            double val1 =global::stringWithUnitsToDouble(captured1,QString("m"),valid1);
//            double val2 =global::stringWithUnitsToDouble(captured2,QString("m"),valid2);
//            double val3 =global::stringWithUnitsToDouble(captured3,QString("m"),valid3);
//            Scale={val1,val2,val3};
//            validLine=(valid1 && valid2 && valid3);
//            alreadyCaptured=true;
//        }
//        if(!alreadyCaptured)
//        {
//            validLine=false;
//        }
//        if(Scale(0)==0){Scale(0)=1;}
//        if(Scale(1)==0){Scale(1)=1;}
//        if(Scale(2)==0){Scale(2)=1;}
//        std::cout<<"scale: "<<Scale(0)<<", "<<Scale(1)<<", "<<Scale(2)<<std::endl;
//        return validLine;
//    }
    QRegularExpression separator( "[ \\t,]+"); //identify any number of white spaces and tabs intermixed or commas
    QStringList list =currentLine.split(separator, Qt::SkipEmptyParts);
    if(!list.isEmpty())
    {
        bool isNumber;
        list.first().toInt(&isNumber);
        if(isNumber)
        {
            readingSpectrumItemList=true;
            validLine=handleSpectrumItemListLine(currentLine);
            return validLine;
        }
    }
    return validLine;
}

bool BESpectrumSection::containsNoElements()
{
    return (bESpectrumItems.isEmpty());
}
