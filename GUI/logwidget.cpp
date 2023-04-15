#include "logwidget.h"

/**
* \class LogWidget
* \brief Provides a read-only plain text editor for displaying log strings.
*/
LogWidget::LogWidget(QWidget *parent) : QPlainTextEdit(parent)
{
    this->setReadOnly(true);
}

void LogWidget::writeOutLogString()
{
    this->clear();
    this->appendPlainText(logStrings::logString);
//    this->appendPlainText(text);
}

void LogWidget::writeOutErrorLogString( )
{
    this->clear();
    this->appendPlainText(logStrings::errorLogString);
//    this->appendPlainText(text);
}
