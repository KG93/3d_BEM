#include "meshfilepropertiessection.h"

MeshFilePropertiesSection::MeshFilePropertiesSection()
{

}
bool MeshFilePropertiesSection::handleScriptLine(const QString& currentLine){
    validLine=false;
    QRegularExpressionMatch match;
    QRegularExpression identifier("^\\s*MeshFileAlias\\s*={1}\\s*([\"\'][\\w\\s]+[\"\']|\\w+)", QRegularExpression::CaseInsensitiveOption);

    if(currentLine.contains(identifier,&match))
    {
//        nodesReference=(currentLine.mid(matchPosition+match.capturedLength())).simplified();
        MeshFileAlias=match.captured(1);
        removeQuotes(MeshFileAlias);
        std::cout <<"MeshFileAlias: "<< MeshFileAlias.toUtf8().constData()<<std::endl;
        if(MeshFileAlias.isEmpty()){validLine=false;}
        else{validLine=true;}
        return validLine;
    }

//    identifier=QRegularExpression("^\\s*MeshFileAlias{1}\\s*={1}\\s*([\"\'][\\w\\s]+[\"\']|\\w+)", QRegularExpression::CaseInsensitiveOption);
//    if(currentLine.contains(identifier))
//    {
////        nodesReference=(currentLine.mid(matchPosition+match.capturedLength())).simplified();
//        MeshFileAlias=match.captured(1);
//        removeQuotes(MeshFileAlias);
//        std::cout <<"MeshFileAlias: "<< MeshFileAlias.toUtf8().constData()<<std::endl;
//        if(MeshFileAlias.isEmpty()){validLine=false;}
//        else{validLine=true;}
//        return validLine;
//    }
    return validLine;

}

void MeshFilePropertiesSection::removeQuotes(QString& line)
{
    QRegularExpression removeChars("[\"\']");
    line.remove(removeChars);
}
