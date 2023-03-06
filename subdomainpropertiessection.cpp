#include "subdomainpropertiessection.h"

SubdomainPropertiesSection::SubdomainPropertiesSection()
{

}

bool SubdomainPropertiesSection::handleScriptLine(const QString currentLine)
{
    validLine=false;
    QRegularExpressionMatch match;
    QString value;
    QRegularExpression identifier;

    identifier=regExWithOneValue("Subdomain");
    if(currentLine.contains(identifier,&match))
    {
        value=match.captured(1);
        global::removeQuotesAndWhitespace(value);
        Subdomain=value.toInt(&validLine);
        std::cout <<"Subdomain: "<< value.toUtf8().constData()<<std::endl;
        return validLine;

    }
    identifier=regExWithOneValue("ElType");
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
    identifier=regExWithOneValue("Color");
    if(currentLine.contains(identifier,&match))
    {
        validLine=true;
        return validLine;

    }

    {
//        std::cout <<"Line in Meshfile_Properties section not recognized!"<< currentLine.toUtf8().constData()<<std::endl;

    }
    return validLine;

}

QRegularExpression SubdomainPropertiesSection::regExWithOneValue(const QString& identifier)
{
//    return QRegularExpression("^\\s*"+identifier+"\\s*={1}\\s*([\"\'][\\w\\s]+[\"\']|\\w+)", QRegularExpression::CaseInsensitiveOption);
    return QRegularExpression("^\\s*"+identifier+"\\s*={1}\\s*([\"\'][^\"\']+[\"\']|.+$)", QRegularExpression::CaseInsensitiveOption);

}
