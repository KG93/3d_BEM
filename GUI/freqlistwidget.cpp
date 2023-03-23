#include "freqlistwidget.h"

FreqListWidget::FreqListWidget(QWidget *parent) : QListWidget(parent)
{
    setWindowTitle(tr("Frequencies in Hz"));
    setContextMenuPolicy(Qt::CustomContextMenu);
//     QObject::connect(this, SIGNAL(customContextMenuRequested(const QPoint&)), this, SLOT(showCustomContextMenu(const QPoint&)));
     QObject::connect(this, SIGNAL(currentRowChanged(int)), this, SLOT(emitNewSelectionSignal(int)));
     setAttribute( Qt::WA_QuitOnClose, false);
}

FreqListWidget::~FreqListWidget()
{

}

void FreqListWidget::setFrequencies(QVector<double> frequenciesARg)
{
    this->clear();
    this->frequencies = frequenciesARg;

    for(int i=0; i<frequencies.length(); i++)
    {
        QListWidgetItem* frequency = new QListWidgetItem(QString::number(frequencies.at(i)/*, 'g', 15*/));
        addItem(frequency);
    }
    currentSelected = currentRow();
    frequencySolved = QVector<bool>(frequencies.length(), false);
    phiSolutionFreq.clear();
    dPhiSolutionFreq.clear();
    soundPressureFreq.clear();
    maxSoundPressureOnBoundary.resize(0);
    minSoundPressureOnBoundary.resize(0);

    phiSolutionField.clear();
    soundPressureField.clear();
    maxSoundPressureOnField.resize(0);
    minSoundPressureOnField.resize(0);
    fieldSolved = false;
}

std::tuple<Eigen::VectorXcd, Eigen::VectorXcd, Eigen::VectorXcd> FreqListWidget::getSolution(int freqIndex)
{
    if(freqIndex >= 0 && freqIndex < phiSolutionFreq.length() && freqIndex < dPhiSolutionFreq.length() && freqIndex < soundPressureFreq.length())
    {
        return {phiSolutionFreq.at(freqIndex), dPhiSolutionFreq.at(freqIndex), soundPressureFreq.at(freqIndex)};
    }
    else
    {
        std::cerr << "Invalid freqIndex parameter in getSolution() call! Returning empty vectors." << std::endl;
        Eigen::VectorXcd empty;
        return {empty, empty, empty};
    }
}

std::pair<double, double> FreqListWidget::getMinAndMaxSoundPressure()
{
    double min = 0, max = 0;
    if(frequencies.size() > 0 && soundPressureFreq.size() == frequencies.size())
    {
        Eigen::VectorXd minPerFreq(frequencies.size());
        Eigen::VectorXd maxPerFreq(frequencies.size());
        for(int i=0; i<frequencies.size(); i++)
        {
            minPerFreq(i) = soundPressureFreq.at(i).cwiseAbs().minCoeff();
            maxPerFreq(i) = soundPressureFreq.at(i).cwiseAbs().maxCoeff();
        }
        if(fieldSolved)
        {
            for(int i=0; i<soundPressureField.size(); i++)
            {
                for(int j=0; j<soundPressureField.at(i).size(); j++)
                {
                    double tmpMin = soundPressureField.at(i).at(j).cwiseAbs().minCoeff();
                    double tmpMax = soundPressureField.at(i).at(j).cwiseAbs().maxCoeff();
                    if(tmpMin < minPerFreq(i))
                    {
                        minPerFreq(i) = tmpMin;
                    }
                    if(tmpMax > maxPerFreq(i))
                    {
                        maxPerFreq(i) = tmpMax;
                    }
                }
            }
        }
        min = minPerFreq.minCoeff();
        max = maxPerFreq.maxCoeff();
    }
    return {min,max};
}

std::pair<double, double> FreqListWidget::getMinAndMaxSoundPressureOnBoundary(long index)
{
    long n = frequencies.size();
    if(n > 0 && soundPressureFreq.size() == n)
    {
        if(minSoundPressureOnBoundary.size() != n || maxSoundPressureOnBoundary.size() != n)
        {
            minSoundPressureOnBoundary.resize(n);
            maxSoundPressureOnBoundary.resize(n);

            #pragma omp parallel for
            for(int i=0; i<frequencies.size(); i++)
            {

                minSoundPressureOnBoundary(i) = soundPressureFreq.at(i).cwiseAbs().minCoeff();
                maxSoundPressureOnBoundary(i) = soundPressureFreq.at(i).cwiseAbs().maxCoeff();
            }
            globalMinSoundPressureOnBoundary = minSoundPressureOnBoundary.minCoeff();
            globalMaxSoundPressureOnBoundary = maxSoundPressureOnBoundary.maxCoeff();
        }
        if(index >= 0 && index < n)
        {
            return {minSoundPressureOnBoundary(index), maxSoundPressureOnBoundary(index)};
        }
        else
        {
            return {globalMinSoundPressureOnBoundary, globalMaxSoundPressureOnBoundary};
        }
    }
    else
    {
        return {0,0};
    }
}

std::pair<double, double> FreqListWidget::getMinAndMaxSoundPressureOnField(long index)
{
    long n = frequencies.size();
    if(fieldSolved && n > 0 && soundPressureField.size() == n)
    {
        if(minSoundPressureOnField.size() != n || maxSoundPressureOnField.size() != n)
        {
            minSoundPressureOnField.resize(n);
            maxSoundPressureOnField.resize(n);

            #pragma omp parallel for
            for(int i=0; i<soundPressureField.size(); i++)
            {
                if(soundPressureField.at(i).size() >= 1)
                {
                    minSoundPressureOnField(i) = soundPressureField.at(i).at(0).cwiseAbs().minCoeff();
                    maxSoundPressureOnField(i) = soundPressureField.at(i).at(0).cwiseAbs().maxCoeff();
                }
                for(int j=1; j<soundPressureField.at(i).size(); j++)
                {
                    double tmpMin = soundPressureField.at(i).at(j).cwiseAbs().minCoeff();
                    double tmpMax = soundPressureField.at(i).at(j).cwiseAbs().maxCoeff();
                    if(tmpMin < minSoundPressureOnField(i))
                    {
                        minSoundPressureOnField(i) = tmpMin;
                    }
                    if(tmpMax > maxSoundPressureOnField(i))
                    {
                        maxSoundPressureOnField(i) = tmpMax;
                    }
                }
            }
            globalMinSoundPressureOnField = minSoundPressureOnField.minCoeff();
            globalMaxSoundPressureOnField = maxSoundPressureOnField.maxCoeff();
        }
        if(index >= 0 && index < n)
        {
            return {minSoundPressureOnField(index), maxSoundPressureOnField(index)};
        }
        else
        {
            return {globalMinSoundPressureOnField, globalMaxSoundPressureOnField};
        }
    }
    else
    {
        return {0,0};
    }
}

bool FreqListWidget::allContainersSameSize()
{
    return phiSolutionFreq.size() == frequencies.size() && dPhiSolutionFreq.size() == frequencies.size() && soundPressureFreq.size() == frequencies.size() && frequencySolved.size() == frequencies.size();
}

void FreqListWidget::emitNewSelectionSignal(int row)
{
    if(!fieldSolved)
    {
        if(row >= 0 && row < phiSolutionFreq.length() && row < dPhiSolutionFreq.length() && row < soundPressureFreq.length())
        {
            emit signalNewSelection(phiSolutionFreq.at(row), dPhiSolutionFreq.at(row), soundPressureFreq.at(row));
        }
    }
    else
    {
        if(row >= 0 && row < phiSolutionFreq.length() && row < dPhiSolutionFreq.length() && row < soundPressureFreq.length()&& row < phiSolutionField.length()&& row < soundPressureField.length())
        {
            emit signalNewSelection(phiSolutionFreq.at(row), dPhiSolutionFreq.at(row), soundPressureFreq.at(row), phiSolutionField.at(row), soundPressureField.at(row));
        }
    }
}

//void FreqListWidget::showCustomContextMenu(const QPoint& pos)
//{
////    QPoint globalPos = this->mapToGlobal(pos);

////    QMenu* contextMenu = new QMenu(this);
////    contextMenu->addAction("Calculate solution", this, SLOT(calcSolutionForFreq()));
////    contextMenu->exec(globalPos);
//}

void FreqListWidget::calcSolutionForFreq()
{
    std::cout << "calcSolutionForFreq() called." << std::endl;
}
