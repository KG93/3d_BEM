#include "controlfieldsection.h"

ControlFieldSection::ControlFieldSection()
{

}

bool ControlFieldSection::handleScriptLine(const QString &currentLine)
{
        validLine=false;
        QRegularExpressionMatch match;
        QString value;
        QRegularExpression identifier;

        identifier=global::regExWithOneValue("f1");
        if(currentLine.contains(identifier,&match))
        {
            value=match.captured(1);
            global::global::removeQuotesAndWhitespace(value);
            std::cout <<"f1: "<< value.toUtf8().constData()<<std::endl;
//            if(MeshFileAlias.isEmpty()){validLine=false;}
//            else{validLine=true;}
//            f1=value.toDouble(&validLine);
            f1=global::stringWithUnitsToDouble(value,QString("Hz"),validLine);
            std::cout <<"f1: "<<f1<<std::endl;

            if(f1<=0){validLine=false;}
            return validLine;
        }

        identifier=global::regExWithOneValue("f2");
         if(currentLine.contains(identifier,&match))
        {
            value=match.captured(1);
            global::removeQuotesAndWhitespace(value);
//            f2=value.toDouble(&validLine);
            f2=global::stringWithUnitsToDouble(value,QString("Hz"),validLine);

            std::cout <<"f2: "<< value.toUtf8().constData()<<std::endl;

            if(f2<=0){validLine=false;}
            return validLine;

        }

         identifier=global::regExWithOneValue("NumFrequencies");
         if(currentLine.contains(identifier,&match))
        {
            value=match.captured(1);
            global::removeQuotesAndWhitespace(value);
            NumFrequencies=value.toUInt(&validLine);
//            NumFrequencies=global::stringWithUnitsToUInt(value,QString(""),validLine);

            std::cout <<"NumFrequencies: "<< value.toUtf8().constData()<<std::endl;
            if(NumFrequencies==0){NumFrequencies=1; validLine=false;}
            return validLine;
        }

        identifier=global::regExWithOneValue("Abscissa");
        if(currentLine.contains(identifier,&match))
        {
            value=match.captured(1);
            global::removeQuotesAndWhitespace(value);
            value=value.simplified();
            std::cout <<"Abscissa: "<< value.toUtf8().constData()<<std::endl;

            if(0==QString::compare("log",value, Qt::CaseInsensitive))
            {
                  validLine=true;
                  abscissa=log;
                  return validLine;
            }
            if(0==QString::compare("lin",value, Qt::CaseInsensitive))
            {
                validLine=true;
                  abscissa=lin;
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
            std::cout <<"MeshFrequency: "<< value.toUtf8().constData()<<std::endl;
//            MeshFrequency=value.toDouble(&validLine);
            MeshFrequency=global::stringWithUnitsToDouble(value,QString("Hz"),validLine);
            std::cout <<"MeshFrequency: "<< value.toUtf8().constData()<<std::endl;

            if(MeshFrequency<=0){validLine=false; MeshFrequency=0;}
            return validLine;
        }

        identifier=global::regExWithOneValue("EdgeLength");
        if(currentLine.contains(identifier,&match))
        {
            value=match.captured(1);
            global::global::removeQuotesAndWhitespace(value);
    //        EdgeLength=value.toDouble(&validLine);
            EdgeLength=global::stringWithUnitsToDouble(value,QString("m"),validLine);
            std::cout<<"Edgelength: "<<EdgeLength<<std::endl;
            if(EdgeLength<=0){EdgeLength=1; validLine=false;}
            return validLine;
        }

        identifier=global::regExWithOneValue("Meshing");
        if(currentLine.contains(identifier,&match))
        {

        }

        identifier=global::global::regExWithOneValue("rho");
        if(currentLine.contains(identifier,&match))
        {
            value=match.captured(1);
            global::removeQuotesAndWhitespace(value);
            if(value.contains("kg/m3",Qt::CaseSensitive)){}
//            c=value.toDouble(&validLine);
            rho=global::stringWithUnitsToDouble(value,QString("kg/m3"),validLine);

            if(rho<=0){validLine=false;}
            return validLine;


        }

        identifier=global::regExWithOneValue("c");
        if(currentLine.contains(identifier,&match))
        {
            value=match.captured(1);
            global::removeQuotesAndWhitespace(value);
            if(value.contains("m/s",Qt::CaseSensitive)){}
//            c=value.toDouble(&validLine);
            c=global::stringWithUnitsToDouble(value,QString("m/s"),validLine);

            if(c<=0){validLine=false;}
            return validLine;

        }

        identifier=global::regExWithOneValue("Altitude");
        if(currentLine.contains(identifier,&match))
        {
            value=match.captured(1);
            global::removeQuotesAndWhitespace(value);
//            Altitude=value.toDouble(&validLine);
            Altitude=global::stringWithUnitsToDouble(value,QString("m"),validLine);

            return validLine;

        }

        identifier=global::regExWithOneValue("StaticPressure");
        if(currentLine.contains(identifier,&match))
        {
            value=match.captured(1);
            global::removeQuotesAndWhitespace(value);
            StaticPressure=global::stringWithUnitsToDouble(value,QString("Pa"),validLine);
            if(StaticPressure<0){validLine=false;}
            return validLine;

        }

        identifier=global::regExWithOneValue("Temperature");
        if(currentLine.contains(identifier,&match))
        {
            value=match.captured(1);
            global::removeQuotesAndWhitespace(value);
            Temperature=global::stringWithUnitsToDouble(value,QString("C"),validLine);

            if(Temperature<-273.15){validLine=false;}
            return validLine;


        }

        std::cout <<"Line in Control_Field section not recognized! "<< currentLine.toUtf8().constData()<<std::endl;

      return validLine;
}
