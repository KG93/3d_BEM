#ifndef QDOUBLESPINBOXNOTRAILZEROS_H
#define QDOUBLESPINBOXNOTRAILZEROS_H

#include <QDoubleSpinBox>
#include <QValidator>
#include <QRegularExpressionValidator>

/**
* \class QDoubleSpinBoxNoTrailZeros
* \brief The QDoubleSpinBoxNoTrailZeros class is a reimplementation of QDoubleSpinBox, but with the trailing zeros removed.
*/
class QDoubleSpinBoxNoTrailZeros : public QDoubleSpinBox
{
    Q_OBJECT
public:
    QDoubleSpinBoxNoTrailZeros(QWidget* parent = nullptr);
    ~QDoubleSpinBoxNoTrailZeros();

    QValidator::State validate(QString &text, int &pos) const override;

    double valueFromText(const QString &text) const override;

    QString textFromValue(double value) const override;

private:
    QRegularExpressionValidator* validator;
};

#endif // QDOUBLESPINBOXNOTRAILZEROS_H
