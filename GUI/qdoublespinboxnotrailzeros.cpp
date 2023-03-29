#include "qdoublespinboxnotrailzeros.h"

QDoubleSpinBoxNoTrailZeros::QDoubleSpinBoxNoTrailZeros(QWidget* parent) : QDoubleSpinBox(parent),
  validator(new QRegularExpressionValidator(QRegularExpression("-?\\d{1,}(?:[,.]{1})\\d*"), this))
{

}

QDoubleSpinBoxNoTrailZeros::~QDoubleSpinBoxNoTrailZeros()
{
    delete validator;
}

QValidator::State QDoubleSpinBoxNoTrailZeros::validate(QString &text, int &pos) const
{
    return validator->validate(text,pos);
}

double QDoubleSpinBoxNoTrailZeros::valueFromText(const QString &text) const
{
    QString tmp = text;
    tmp.replace(",", ".");
    return tmp.toDouble();
}

QString QDoubleSpinBoxNoTrailZeros::textFromValue(double value) const
{
    QString returnStr = QString::number(value, 'g', 15);
    return returnStr;
}
