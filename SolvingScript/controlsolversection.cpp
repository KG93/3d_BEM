#include "controlsolversection.h"

//ControlSolverSection::ControlSolverSection(QObject *parent) : QObject(parent)
//{

//}

ControlSolverSection::ControlSolverSection()
{

}

bool ControlSolverSection::handleScriptLine(const QString currentLine)
{
        validLine=false;
        QRegularExpressionMatch match;
        QString value;
        QRegularExpression identifier;

        identifier=global::regExWithOneValue("f1");
        if(currentLine.contains(identifier,&match))
        {
            value=match.captured(1);
            global::removeQuotesAndWhitespace(value);
            std::cout <<"f1: "<< value.toUtf8().constData()<<std::endl;
            f1=global::stringWithUnitsToDouble(value,QString("Hz"),validLine);
            std::cout <<"f1: "<<f1<<std::endl;
            QString message=QString("f1: "+ QString::number(f1)+"Hz");
            logStrings::logString.append(message+"\r\n");
            if(f1<=0){validLine=false;}
            return validLine;
        }

        identifier=global::regExWithOneValue("f2");
         if(currentLine.contains(identifier,&match))
        {
            value=match.captured(1);
            global::removeQuotesAndWhitespace(value);
            f2=global::stringWithUnitsToDouble(value,QString("Hz"),validLine);

            std::cout <<"f2: "<< value.toUtf8().constData()<<std::endl;
            QString message=QString("f2: "+ QString::number(f2)+"Hz");
            logStrings::logString.append(message+"\r\n");
            if(f2<=0){validLine=false;}
            return validLine;
        }

         identifier=global::regExWithOneValue("NumFrequencies");
         if(currentLine.contains(identifier,&match))
        {
            value=match.captured(1);
            global::removeQuotesAndWhitespace(value);
            NumFrequencies=value.toUInt(&validLine);
            std::cout <<"NumFrequencies: "<< value.toUtf8().constData()<<std::endl;
            QString message=QString("NumFrequencies: "+ value);
            logStrings::logString.append(message+"\r\n");
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
            QString message=QString("Abscissa: "+ value);
            logStrings::logString.append(message+"\r\n");
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
            MeshFrequency=global::stringWithUnitsToDouble(value,QString("Hz"),validLine);
            QString message=QString("Global MeshFrequency: "+QString::number(MeshFrequency));
            logStrings::logString.append(message+"\r\n");
            std::cout <<"Global MeshFrequency: "<< MeshFrequency /*value.toUtf8().constData()*/<<std::endl;
            if(MeshFrequency<=0){validLine=false; MeshFrequency=0;}
            return validLine;
        }

        identifier=global::regExWithOneValue("Meshing");
        if(currentLine.contains(identifier,&match))
        {
            return true;
        }

//        identifier= QRegularExpression("^\\s*Sym\\s*={1}\\s*([\"\'][\\w\\s]+[\"\']|\\w+)", QRegularExpression::CaseInsensitiveOption);
        identifier=global::regExWithOneValue("Sym");
        if(currentLine.contains(identifier,&match))
        {
            return true;
        }

//        identifier= QRegularExpression("^\\s*Dim\\s*={1}\\s*([\"\'][\\w\\s]+[\"\']|\\w+)", QRegularExpression::CaseInsensitiveOption);
        identifier=global::regExWithOneValue("Dim");
        if(currentLine.contains(identifier,&match))
        {
            return true;
        }

//        identifier= QRegularExpression("^\\s*NUC\\s*={1}\\s*([\"\'][\\w\\s]+[\"\']|\\w+)", QRegularExpression::CaseInsensitiveOption);
        identifier=global::regExWithOneValue("NUC");
        if(currentLine.contains(identifier,&match))
        {
            return true;
        }

        identifier=global::regExWithOneValue("ImpedanceNumFrequencies");
        if(currentLine.contains(identifier,&match))
        {
            return true;
        }

        identifier=global::regExWithOneValue("rho");
        if(currentLine.contains(identifier,&match))
        {
            value=match.captured(1);
            global::removeQuotesAndWhitespace(value);
            rho=global::stringWithUnitsToDouble(value,QString("kg/m3"),validLine);

            if(rho<=0){validLine=false;}
            return validLine;
        }

        identifier=global::regExWithOneValue("c");
        if(currentLine.contains(identifier,&match))
        {
            value=match.captured(1);
            global::removeQuotesAndWhitespace(value);
            c=global::stringWithUnitsToDouble(value,QString("m/s"),validLine);

            if(c<=0){validLine=false;}
            return validLine;

        }

        identifier=global::regExWithOneValue("Altitude");
        if(currentLine.contains(identifier,&match))
        {
            value=match.captured(1);
            global::removeQuotesAndWhitespace(value);
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

        std::cout <<"Line in Control_Solver section not recognized: "<< currentLine.toUtf8().constData()<<std::endl;
        QString error=QString("Line in Control_Solver section not recognized: "+ currentLine);
        logStrings::errorLogString.append(error+"\r\n");
        return validLine;
}

QVector<double> ControlSolverSection::getFrequencies()
{
    if(f1 != 0 && f2 != 0)
    {
        double tmp1=std::min(f1,f2);
        double tmp2=std::max(f1,f2);
        f1=tmp1;
        f2=tmp2;
    }

    if(NumFrequencies<=1 || f2 == 0)
    {
        double frequency=std::max(f1,f2);
        QVector<double> vector={frequency};
        return vector;
    }
    else
    {
        if(abscissa==log)//log
        {
            double logF1=std::log(f1);
            double logF2=std::log(f2);
            QVector<double> vector;
            vector.resize(NumFrequencies);
            double diff=(logF2-logF1)/(NumFrequencies-1);
            for(quint64 i=0;i<NumFrequencies;i++)
            {
                double tmp=logF1+i*diff;
                vector[i]=std::exp(tmp);
            }
            return vector;

        }
        else // lin
        {
            QVector<double> vector;
            vector.resize(NumFrequencies);
            double diff=(f2-f1)/(NumFrequencies-1);
            for(quint64 i=0;i<NumFrequencies;i++)
            {
                vector[i]=f1+i*diff;
            }
            return vector;
        }
    }
}

bool ControlSolverSection::allMandatoryParametersSet()
{
    bool val = true;
    if(f1==0)
    {
        std::cout <<"f1 in controlsolver section not declared."<<std::endl;
        QString error=QString("f1 in controlsolver section not declared.");
        logStrings::errorLogString.append(error+"\r\n");
        val = false;
    }
    if(f2==0)
    {
//        std::cout <<"f2 in controlsolver section not declared. Defaulted to 100."<<std::endl;
//        QString error=QString("f2 in controlsolver section not declared. Defaulted to 100.");
        std::cout <<"f2 in controlsolver section not declared. Defaulted to 0."<<std::endl;
        QString error=QString("f2 in controlsolver section not declared. Defaulted to 0.");
        logStrings::errorLogString.append(error+"\r\n");
        f2=0;
        val = false;
    }
    if(MeshFrequency==0)
    {
        std::cout <<"MeshFrequency in controlsolver section not declared."<<std::endl;
        QString error=QString("MeshFrequency in controlsolver section not declared.");
        logStrings::errorLogString.append(error+"\r\n");
        val = false;
    }
    if(NumFrequencies==0)
    {
        std::cout <<"NumFrequencies in controlsolver section not declared. Defaulted to 1."<<std::endl;
        QString error=QString("NumFrequencies in controlsolver section not declared. Defaulted to 1.");
        logStrings::errorLogString.append(error+"\r\n");
        NumFrequencies=1;
        val = false;
    }
    if(c==0)
    {
        std::cout <<"Wavespeed in controlsolver section not declared. Defaulted to 343.2m/s."<<std::endl;
        QString error=QString("Wavespeed in controlsolver section not declared. Defaulted to 343.2m/s.");
        logStrings::errorLogString.append(error+"\r\n");
        c=343.2;
        val = false;
    }
    return val;

//            double f2;
//    double MeshFrequency;
//    double c;
//    double Altitude;
//    double StaticPressure;
//    double Temperature;

//    quint64 NumFrequencies;
}

double ControlSolverSection::getWavespeed()
{
    if(c==0)
    {
        c=343.21;
    }
    return c;
}

double ControlSolverSection::getEdgelength()
{
    if(c==0)
    {
        c=343.21;
    }
    double maxFrequency = std::max(f1,f2);
    double edgelength = (c/maxFrequency)/2;
    if(MeshFrequency!=0)
    {
        edgelength=(c/MeshFrequency)/2;
    }
    return edgelength;
}

double ControlSolverSection::getTemperatureInC()
{
    if(Temperature==0)
    {
        Temperature=20;
    }
    return Temperature;
}

double ControlSolverSection::getDensity()
{
    if(rho==0)
    {
        rho=1.2041;
    }
    return rho;
//    double specificGasConstant=287.058;
//    double airDensity=staticPressure/((temperature+273.15)*specificGasConstant);
}

void ControlSolverSection::clear()
{
    f1=0;
    f2=0;
    MeshFrequency=0;
    c=0;
    Altitude=0;
    StaticPressure=0;
    Temperature=0;
    rho=0;
}
