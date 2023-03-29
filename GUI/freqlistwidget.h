#ifndef FREQLISTWIDGET_H
#define FREQLISTWIDGET_H

#include <QListWidget>
#include <QMenu>
#include <iostream>
#include <eigen3/Eigen/Geometry>

/**
* \class FreqListWidget
* \brief This class provides a frequency selector widget and also stores the frequency-dependent BEM solutions.
*/
class FreqListWidget : public QListWidget
{
    Q_OBJECT
public:
    FreqListWidget(QWidget* parent = nullptr);
    ~FreqListWidget();

    void setFrequencies(const QVector<double> frequenciesArg);
    QVector<double> getFrequencies(){return frequencies;}
    unsigned int getCurrentSelectedFreqIndex() {currentSelected = currentRow(); return currentSelected;}

private:
    int currentSelected = 0;
    QVector<double> frequencies;

public slots:
    void colorSolvedFreq(const int freqIndex);
};

#endif // FREQLISTWIDGET_H
