#ifndef SUBDOMAINPROPERTIESSECTION_H
#define SUBDOMAINPROPERTIESSECTION_H

#include"global.h"
#include <iostream>
#include <QObject>
#include <QRegularExpressionMatch>
#include <QRegularExpression>


class SubdomainPropertiesSection
{
public:
    SubdomainPropertiesSection();
    SubdomainPropertiesSection(QString name){this->name = name;}
    bool handleScriptLine(const QString line);

    QString name;
    int Subdomain;
    int ElType = 0; //default interior
    enum ElTypeEnum {interior = 0, exterior};

private:
    bool validLine;
    QRegularExpression regExWithOneValue(const QString& identifier);
};

#endif // SUBDOMAINPROPERTIESSECTION_H
