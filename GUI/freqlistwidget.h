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

    void setFrequencies(QVector<double> frequenciesArg);
    void setSolutions(long i, Eigen::VectorXcd phiSolution, Eigen::VectorXcd dPhiSolution, Eigen::VectorXcd soundPressure){
        if(phiSolutionFreq.size() != frequencies.size() || dPhiSolution.size() != frequencies.size() || soundPressure.size() != frequencies.size())
        {
            phiSolutionFreq.resize(frequencies.size());
            dPhiSolutionFreq.resize(frequencies.size());
            soundPressureFreq.resize(frequencies.size());
            minSoundPressureOnBoundary.resize(0);
            maxSoundPressureOnBoundary.resize(0);
        }
        if(i >= 0 && i < phiSolutionFreq.length() && i < dPhiSolutionFreq.length() && i < soundPressureFreq.length())
        {
            this->phiSolutionFreq[i] = phiSolution;
            this->dPhiSolutionFreq[i] = dPhiSolution;
            this->soundPressureFreq[i] = soundPressure;
        }
        fieldSolved = false;
    }

    void setSolutions(Eigen::VectorXcd phiSolution, Eigen::VectorXcd dPhiSolution, Eigen::VectorXcd soundPressure){
        if(phiSolutionFreq.size() != frequencies.size() || dPhiSolution.size() != frequencies.size() || soundPressure.size() != frequencies.size())
        {
            phiSolutionFreq.resize(frequencies.size());
            dPhiSolutionFreq.resize(frequencies.size());
            soundPressureFreq.resize(frequencies.size());
            minSoundPressureOnBoundary.resize(0);
            maxSoundPressureOnBoundary.resize(0);
        }
        if(currentSelected >= 0 && currentSelected < phiSolutionFreq.length() && currentSelected < dPhiSolutionFreq.length() && currentSelected < soundPressureFreq.length())
        {
            this->phiSolutionFreq[currentSelected] = phiSolution;
            this->dPhiSolutionFreq[currentSelected] = dPhiSolution;
            this->soundPressureFreq[currentSelected] = soundPressure;
        }
        fieldSolved = false;
    }

    void setSolutionsField(long i, QVector<Eigen::VectorXcd> phiSolutionFieldArg, QVector<Eigen::VectorXcd> soundPressureFieldArg){
        if(phiSolutionField.size() != frequencies.size() || soundPressureField.size() != frequencies.size())
        {
            phiSolutionField.resize(frequencies.size());
            soundPressureField.resize(frequencies.size());
            minSoundPressureOnField.resize(0);
            maxSoundPressureOnField.resize(0);
        }
        if(i >= 0 && i < phiSolutionField.length() && i < soundPressureField.length())
        {
            this->phiSolutionField[i] = phiSolutionFieldArg;
            this->soundPressureField[i] = soundPressureFieldArg;
        }
        fieldSolved = true;
    }

    QVector<double> getFrequencies(){return frequencies;}
    std::tuple<Eigen::VectorXcd, Eigen::VectorXcd, Eigen::VectorXcd> getSolution(int freqIndex);
    std::pair<double, double> getMinAndMaxSoundPressure();
    std::pair<double, double> getMinAndMaxSoundPressureOnBoundary(long index = -1);
    std::pair<double, double> getMinAndMaxSoundPressureOnField(long index = -1);

private:
    bool fieldSolved = false;
    int currentSelected = -1;
    QVector<double> frequencies;

//    QVector<QVector<Eigen::Vector4d>> triangleColors;
    QVector<Eigen::VectorXcd> phiSolutionFreq;
    QVector<Eigen::VectorXcd> dPhiSolutionFreq;
    QVector<Eigen::VectorXcd> soundPressureFreq;
    Eigen::VectorXd maxSoundPressureOnBoundary;
    double globalMaxSoundPressureOnBoundary;
    Eigen::VectorXd minSoundPressureOnBoundary;
    double globalMinSoundPressureOnBoundary;

    QVector<QVector<Eigen::VectorXcd> > phiSolutionField;
    QVector<QVector<Eigen::VectorXcd> > soundPressureField;
    Eigen::VectorXd maxSoundPressureOnField;
    double globalMaxSoundPressureOnField;
    Eigen::VectorXd minSoundPressureOnField;
    double globalMinSoundPressureOnField;

signals:
//    void signalNewSelection(double frequency);
    void signalNewSelection(Eigen::VectorXcd phiSolution, Eigen::VectorXcd dPhiSolution, Eigen::VectorXcd soundPressure);
    void signalNewSelection(Eigen::VectorXcd phiSolution, Eigen::VectorXcd dPhiSolution, Eigen::VectorXcd soundPressure, QVector<Eigen::VectorXcd> phiSolutionField, QVector<Eigen::VectorXcd> soundPressureField);

private slots:
    void emitNewSelectionSignal(int row);
//    void showCustomContextMenu(const QPoint& pos);
    void calcSolutionForFreq();
};

#endif // FREQLISTWIDGET_H
