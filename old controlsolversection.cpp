#include "controlsolversection.h"

ControlSolverSection::ControlSolverSection()
{

}

bool ControlSolverSection::handleLine(const QString& currentLine){

    validLine=false;
    QRegularExpressionMatch match;
    QRegularExpression identifier("^\\s*f1\\s*={1}\\s*([\"\'][\\w\\s]+[\"\']|\\w+)", QRegularExpression::CaseInsensitiveOption);

    if(currentLine.contains(identifier))
    {
//        nodesReference=(currentLine.mid(matchPosition+match.capturedLength())).simplified();
//        MeshFileAlias=match.captured(1);
//        removeQuotes(MeshFileAlias);
//        std::cout <<"MeshFileAlias: "<< MeshFileAlias.toUtf8().constData()<<std::endl;
//        if(MeshFileAlias.isEmpty()){validLine=false;}
//        else{validLine=true;}
//        return validLine;
    }


    if(currentLine.startsWith("f1=",Qt::CaseInsensitive))
    {

    }
    else if(currentLine.startsWith("f2=",Qt::CaseInsensitive))
    {

    }
    else if(currentLine.startsWith("NumFrequencies=",Qt::CaseInsensitive))
    {

    }
    else if(currentLine.startsWith("Abscissa=",Qt::CaseInsensitive))
    {

    }
    else if(currentLine.startsWith("MeshFrequency=",Qt::CaseInsensitive))
    {

    }
    else if(currentLine.startsWith("Meshing=",Qt::CaseInsensitive))
    {

    }
    else if(currentLine.startsWith("Sym=",Qt::CaseInsensitive))
    {

    }
    else if(currentLine.startsWith("Dim=",Qt::CaseInsensitive))
    {

    }
    else if(currentLine.startsWith("NUC=",Qt::CaseInsensitive))
    {

    }
    else if(currentLine.startsWith("ImpedanceNumFrequencies=",Qt::CaseInsensitive))
    {

    }
    else if(currentLine.startsWith("rho=",Qt::CaseInsensitive))
    {

    }
    else if(currentLine.startsWith("c=",Qt::CaseInsensitive))
    {

    }
    else if(currentLine.startsWith("Altitude=",Qt::CaseInsensitive))
    {

    }
    else if(currentLine.startsWith("StaticPressure=",Qt::CaseInsensitive))
    {

    }
    else if(currentLine.startsWith("Temperature=",Qt::CaseInsensitive))
    {

    }
    else{
        std::cout <<"Line in Control_Solver section not recognized!"<< currentLine.toUtf8().constData()<<std::endl;

    }
    return validLine;

}

void ControlSolverSection::removeQuotes(QString& line)
{
    QRegularExpression removeChars("[\"\']");
    line.remove(removeChars);
}

