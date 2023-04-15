#include "freqlistwidget.h"

FreqListWidget::FreqListWidget(QWidget* parent) : QListWidget(parent)
{
    setWindowTitle(tr("Frequencies in Hz"));
    setContextMenuPolicy(Qt::CustomContextMenu);
//     QObject::connect(this, SIGNAL(customContextMenuRequested(const QPoint&)), this, SLOT(showCustomContextMenu(const QPoint&)));
//     QObject::connect(this, SIGNAL(currentRowChanged(int)), this, SLOT(emitNewSelectionSignal(int)));
     setAttribute( Qt::WA_QuitOnClose, false);
     setUniformItemSizes(true);
}

FreqListWidget::~FreqListWidget()
{

}

void FreqListWidget::setFrequencies(const QVector<double> frequenciesArg)
{
    this->clear();
    this->frequencies = frequenciesArg;

    for(int i=0; i<frequencies.size(); i++)
    {
        QListWidgetItem* frequency = new QListWidgetItem(QString::number(frequencies.at(i)));
        addItem(frequency);
    }
    if(currentSelected < frequencies.size())
    {
        setCurrentRow(currentSelected);
    }
    else if(frequencies.size() > 0)
    {
        setCurrentRow(0);
    }
    currentSelected = currentRow();
}

unsigned int FreqListWidget::getCurrentSelectedFreqIndex()
{
    currentSelected = currentRow();
    if(currentSelected < 0 || currentSelected >= count())
    {
        return 0;
    }
    return currentSelected;
}

void FreqListWidget::colorSolvedFreq(const int freqIndex)
{
    int numberOfFrequencies = this->count();
    if(0 <= freqIndex && freqIndex < numberOfFrequencies)
    {
        QListWidgetItem* freqItem = this->item(freqIndex);
        freqItem->setForeground(Qt::darkGreen);
    }
}
