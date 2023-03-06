#ifndef MESHFILEPROPERTIESSECTION_H
#define MESHFILEPROPERTIESSECTION_H

#include <iostream>
#include <QString>
#include <QRegularExpressionMatch>
#include <QRegularExpression>


class MeshFilePropertiesSection
{
public:
    MeshFilePropertiesSection();
    MeshFilePropertiesSection(QString name){this->name = name;}

    QString name;
    bool handleScriptLine(const QString& line);
    QString MeshFileAlias;
//        Tags =
//        MinAngle =
//        MinArea =
//        Scale =
//        Shift =
//        Rotate =
//        Mirror =


private:
    bool validLine = false;
    void removeQuotes(QString& line);


};

#endif // MESHFILEPROPERTIESSECTION_H
