#ifndef LOGWIDGET_H
#define LOGWIDGET_H
#include "globallogstrings.h"
#include <QPlainTextEdit>
#include <QFile>


class LogWidget: public QPlainTextEdit
{
    Q_OBJECT
public:
    explicit LogWidget(QWidget *parent = nullptr);

private:
    QFile m_logFile;

public slots:
    void writeOutLogString();
    void writeOutErrorLogString();
};

#endif // LOGWIDGET_H
