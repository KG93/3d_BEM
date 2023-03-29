#include "parameterdialog.h"

ParameterDialog::ParameterDialog(QWidget* parent) : QDialog(parent)
{
////////////////////// CouplingParameters

    couplingCheckBox = new QCheckBox(tr("Coupling"));
    couplingCheckBox->setChecked(true);

    burtonMillerCheckBox = new QRadioButton(tr("Burton Miller"));
    burtonMillerCheckBox->setChecked(false);
    kirkupCheckBox = new QRadioButton(tr("Kirkup"));
    kirkupCheckBox->setChecked(true);

    QGridLayout* couplingParameterLayout = new QGridLayout;
    couplingParameterLayout->addWidget(burtonMillerCheckBox, 0, 1);
    couplingParameterLayout->addWidget(kirkupCheckBox, 0, 2);

    QGroupBox* couplingParameterGroupBox = new QGroupBox(tr("Coupling parameters"));
    couplingParameterGroupBox->setLayout(couplingParameterLayout);


////////////////////// HSolver Parameters
    hsolverCheckBox = new QCheckBox(tr("H-Matrix solver"));
    hsolverCheckBox->setChecked(true);

    acaErrorLabel = new QLabel(tr("relavive ACA error:"));

    acaErrorQSpinBox = new QDoubleSpinBoxNoTrailZeros;
    acaErrorQSpinBox->setValue(0.01);
    acaErrorQSpinBox->setSingleStep(0.01);
    acaErrorQSpinBox->setStepType(QAbstractSpinBox::AdaptiveDecimalStepType);
    acaErrorQSpinBox->setRange(0,1);
//    acaErrorQSpinBox->setDecimals(16);
//    acaErrorLabel->setBuddy(acaErrorQSpinBox);

    acaMaxRankLabel = new QLabel(tr("max ACA Rank:"));
    acaMaxRankQSpinBox= new QSpinBox;
    acaMaxRankQSpinBox->setMinimum(0);
    acaMaxRankQSpinBox->setMaximum(global::maxInt);
    acaMaxRankQSpinBox->setValue(0);
    acaMaxRankLabel->setBuddy(acaMaxRankQSpinBox);

    usePreconditionerCheckBox = new QCheckBox(tr("Use preconditioner"));
    usePreconditionerCheckBox->setChecked(false);

    precRankLabel = new QLabel(tr("preconditioner Rank:"));
    precRankQSpinBox = new QSpinBox;
    precRankQSpinBox->setMinimum(0);
    precRankQSpinBox->setMaximum(global::maxInt);
    precRankQSpinBox->setValue(2);
    precRankLabel->setBuddy(precRankQSpinBox);

    precErrorLabel = new QLabel(tr("Relative error of preconditioner"));
    precErrorQSpinBox = new QDoubleSpinBoxNoTrailZeros;
    precErrorQSpinBox->setValue(0.0);
    precErrorQSpinBox->setSingleStep(0.1);
    precErrorQSpinBox->setStepType(QAbstractSpinBox::AdaptiveDecimalStepType);
    precErrorQSpinBox->setRange(0,1);
//    precErrorQSpinBox->setDecimals(16);

    calculateNormCheckBox = new QCheckBox(tr("Calculate matrix norm and condition number"));
    calculateNormCheckBox->setChecked(false);

    QGridLayout* hParamLayout = new QGridLayout;
//    hParamLayout->setSizeConstraint(QLayout::SetFixedSize);
    QGroupBox* hParameterWidget = new QGroupBox;
    hParameterWidget->setLayout(hParamLayout);

    hParamLayout->addWidget(acaErrorLabel, 1, 1);
    hParamLayout->addWidget(acaErrorQSpinBox, 1, 2);

    hParamLayout->addWidget(acaMaxRankLabel, 2, 1);
    hParamLayout->addWidget(acaMaxRankQSpinBox, 2, 2);

    hParamLayout->addWidget(usePreconditionerCheckBox, 3, 1);

    hParamLayout->addWidget(precRankLabel, 4, 1);
    hParamLayout->addWidget(precRankQSpinBox, 4, 2);

    hParamLayout->addWidget(precErrorLabel, 5, 1);
    hParamLayout->addWidget(precErrorQSpinBox, 5, 2);

    hParamLayout->addWidget(calculateNormCheckBox, 7, 1,1,2);

////////////////////// Field Solving Parameters

    HFieldSolverCheckBox = new QCheckBox(tr("H-Matrix field solver"));
    HFieldSolverCheckBox->setChecked(true);

    QLabel* fieldACALabel = new QLabel(tr("Field ACA Rank:"));
    fieldACAMaxRankQSpinBox = new QSpinBox;
    fieldACAMaxRankQSpinBox->setMinimum(0);
    fieldACAMaxRankQSpinBox->setMaximum(global::maxInt);
    fieldACAMaxRankQSpinBox->setValue(0);
    fieldACALabel->setBuddy(fieldACAMaxRankQSpinBox);

    QLabel* fieldACAErrorLabel = new QLabel(tr("Relative field ACA error"));
    fieldACARelErrorQSpinBox = new QDoubleSpinBoxNoTrailZeros;
    fieldACARelErrorQSpinBox->setValue(0.01);
    fieldACARelErrorQSpinBox->setSingleStep(0.01);
    fieldACARelErrorQSpinBox->setStepType(QAbstractSpinBox::AdaptiveDecimalStepType);
    fieldACARelErrorQSpinBox->setRange(0,1);
//    fieldACARelErrorQSpinBox->setDecimals(16);
    fieldACAErrorLabel->setBuddy(fieldACARelErrorQSpinBox);

    QGridLayout* fieldParameterLayout = new QGridLayout;
    QGroupBox* fieldParameterGroupBox = new QGroupBox(tr("Field solver parameters"));
    fieldParameterGroupBox->setLayout(fieldParameterLayout);

    fieldParameterLayout->addWidget(fieldACALabel,1,1);
    fieldParameterLayout->addWidget(fieldACAMaxRankQSpinBox,1,2/*,Qt::AlignRight*/);
    fieldParameterLayout->addWidget(fieldACAErrorLabel,2,1);
    fieldParameterLayout->addWidget(fieldACARelErrorQSpinBox,2,2/*,Qt::AlignRight*/);
//        fieldACAMaxRankQSpinBox->resize(fieldACARelErrorQSpinBox->size());

////////////////////// MainLayout

    QGridLayout* mainLayout = new QGridLayout;
    mainLayout->setSizeConstraint(QLayout::SetFixedSize);
    mainLayout->addWidget(couplingCheckBox, 0, 1);
    mainLayout->addWidget(couplingParameterGroupBox,1,1);
    mainLayout->addWidget(hsolverCheckBox, 2, 1);
    mainLayout->addWidget(hParameterWidget,3,1);
    mainLayout->addWidget(HFieldSolverCheckBox,4,1);
    mainLayout->addWidget(fieldParameterGroupBox,5,1);

    connect(couplingCheckBox, &QAbstractButton::toggled, couplingParameterGroupBox, &QWidget::setVisible);
    connect(hsolverCheckBox, &QAbstractButton::toggled, hParameterWidget, &QWidget::setVisible);
    connect(usePreconditionerCheckBox, &QAbstractButton::toggled, precRankLabel, &QWidget::setVisible);
    connect(usePreconditionerCheckBox, &QAbstractButton::toggled, precRankQSpinBox, &QWidget::setVisible);
    connect(usePreconditionerCheckBox, &QAbstractButton::toggled, precErrorLabel, &QWidget::setVisible);
    connect(usePreconditionerCheckBox, &QAbstractButton::toggled, precErrorQSpinBox, &QWidget::setVisible);
    connect(HFieldSolverCheckBox, &QAbstractButton::toggled, fieldParameterGroupBox, &QWidget::setVisible);

    setLayout(mainLayout);
    setWindowTitle(tr("Solver parameters"));
//    this->show();
}

BemCoupling ParameterDialog::getCoupling()
{
    if(couplingCheckBox->isChecked())
    {
        if(burtonMillerCheckBox->isChecked())
        {
            return BurtonMillerCoupling;
        }
        else
        {
            return KirkupCoupling;
        }
    }
    else
    {
        return NoCoupling;
    }
}

bool ParameterDialog::getUseHSolver()
{
    return hsolverCheckBox->isChecked();
}

double ParameterDialog::getACARelError()
{
    return acaErrorQSpinBox->value();
}

double ParameterDialog::getACAMaxRank()
{
    return acaMaxRankQSpinBox->value();
}

bool ParameterDialog::getUseHPreconditioning()
{
    return usePreconditionerCheckBox->isChecked();
}

double ParameterDialog::getPreconditionerRelError()
{
    return precErrorQSpinBox->value();
}

double ParameterDialog::getPreconditionerMaxRank()
{
    return precRankQSpinBox->value();
}

double ParameterDialog::getCalculateNormAndConditionNumber()
{
    return calculateNormCheckBox->isChecked();
}

bool ParameterDialog::getUseHFieldSolver()
{
    return HFieldSolverCheckBox->isChecked();
}

int ParameterDialog::getHFieldSolverMaxRank()
{
    return fieldACAMaxRankQSpinBox->value();
}

double ParameterDialog::getHFieldSolverRelError()
{
    return fieldACARelErrorQSpinBox->value();
}

BemSolverParameters ParameterDialog::getBemSolverParameters()
{
    BemSolverParameters parameters;

    parameters.hMatSolving = getUseHSolver(); // Use the (fast) H-matrix solver.
    parameters.acaRelativeError = getACARelError(); // Controls the approximation error for the adaptive cross approximation.
    parameters.acaMaxRank = getACAMaxRank(); // Controls the maximum rank for the adaptive cross approximation.
    parameters.usePreconditioner = getUseHPreconditioning(); // Use HLU preconditioning for the H-matrix solver.
    parameters.preconditionerRank = getPreconditionerMaxRank(); // Sets the (maximum) local rank of the H-LU preconditioner for the GMRES. If set to zero, the acaRelativeError tolerance is used.
    parameters.preconditionerRelError = getPreconditionerRelError(); // Controls the accuracy during the preconditioner calculation.
    parameters.hMatFieldSolving = getUseHFieldSolver(); // Use the (fast) H-matrix solver for the field calculation.
    parameters.fieldACARelativeError = getHFieldSolverRelError(); // Controls the approximation error for the adaptive cross approximation in the field calculation.
    parameters.fieldACAMaxRank = getHFieldSolverMaxRank(); // Sets the (maximum) local rank of the ACA for the field calculation. If set to zero, only the fieldACARelativeError tolerance is used.
    parameters.calculateNormAndCond = getCalculateNormAndConditionNumber(); // Controls whether to estimate the norm and condition number of the dphi bem matrix.
    parameters.coupling = getCoupling(); // Holds the selected BEM coupling method to solve the non-uniqueness problem.

    return parameters;
}
