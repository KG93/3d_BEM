#ifndef CONTROLSOLVERSECTION_H
#include <iostream>
#include <QString>
#include <QRegularExpressionMatch>
#include <QRegularExpression>
#define CONTROLSOLVERSECTION_H


class ControlSolverSection
{
public:
    ControlSolverSection();
    bool handleLine(const QString& currentLine);

private:

    bool validLine=false;
    void removeQuotes(QString& line);
};

#endif // CONTROLSOLVERSECTION_H
