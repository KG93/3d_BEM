#ifndef PARAMETERDIALOG_H
#define PARAMETERDIALOG_H

#include "global.h"
#include <QButtonGroup>
#include <QCheckBox>
#include <QDialog>
#include <QGridLayout>
#include <QGroupBox>
#include <QLabel>
#include <QLineEdit>
#include <QRadioButton>
#include <QSpinBox>
#include <QValidator>

#include "qdoublespinboxnotrailzeros.h"
#include "../bemparameters.h"

/**
* \class ParameterDialog
* \brief The ParameterDialog class provides a dialog for setting the simulation parameters.
*/
class ParameterDialog : public QDialog {
    Q_OBJECT
public:
    ParameterDialog(QWidget* parent = nullptr);

    /**
    * \brief Returns the selected BEM coupling method to solve the non-uniqueness problem.
    * \return Returns one of NoCoupling, BurtonMillerCoupling and KirkupCoupling.
    */
    BemCoupling getCoupling();

    /**
    * \brief Returns whether the H-Matrix solver is enabled.
    * \return true if the H-Matrix solver is enabled, false otherwise
    */
    bool getUseHSolver();

    /**
    * \brief Returns the relative ACA error specified by the user.
    * \return The relative ACA error.
    */
    double getACARelError();

    /**
    * \brief Returns the maximum ACA rank specified by the user.
    * \return The maximum ACA rank.
    */
    double getACAMaxRank();

    /**
    * \brief Returns whether the preconditioner is enabled.
    * \return True if the preconditioner is enabled, false otherwise.
    */
    bool getUseHPreconditioning();

    /**
    * \brief Returns the relative error of the preconditioner specified by the user.
    * \return The relative error of the preconditioner.
    */
    double getPreconditionerRelError();

    /**
    * \brief Returns the maximum rank of the preconditioner specified by the user.
    * \return The maximum rank of the preconditioner.
    */
    double getPreconditionerMaxRank();

    /**
    * \brief Returns whether the user has chosen to calculate the matrix norm and condition number.
    * \return True if the matrix norm and condition number should be calculated, false otherwise.
    */
    double getCalculateNormAndConditionNumber();

    /**
    * \brief Returns whether the H-Field solver is enabled.
    * \return True if the H-Field solver is enabled, false otherwise.
    */
    bool getUseHFieldSolver();

    /**
    * \brief Returns the maximum rank of the H-Field solver specified by the user.
    * \return The maximum rank of the H-Field solver.
    */
    int getHFieldSolverMaxRank();

    /**
    * \brief Returns the maximum rank of the H-field solver.
    * \return The maximum rank of the H-field solver as an integer.
    */
    double getHFieldSolverRelError();

    BemSolverParameters getBemSolverParameters();

private:
  // checkboxes for the bem coupling parameter
  QCheckBox* couplingCheckBox;
  QRadioButton* burtonMillerCheckBox;
  QRadioButton* kirkupCheckBox;

  QLabel* couplingParamLabel;

  QCheckBox* hsolverCheckBox;
  QLabel* acaErrorLabel;
  QLineEdit* acaErrorLineEdit;
  QDoubleSpinBoxNoTrailZeros* acaErrorQSpinBox;
  QLabel* acaMaxRankLabel;
  QSpinBox* acaMaxRankQSpinBox;

  QCheckBox* usePreconditionerCheckBox;

  QLabel* precRankLabel;
  QSpinBox* precRankQSpinBox;
  QLabel* precErrorLabel;
  QDoubleSpinBoxNoTrailZeros* precErrorQSpinBox;

  QCheckBox* calculateNormCheckBox;

  QDoubleValidator doubleValidator = QDoubleValidator(0, 1, 16);
  QIntValidator intValidator = QIntValidator(0, global::maxInt);

  QCheckBox* HFieldSolverCheckBox;
  QSpinBox* fieldACAMaxRankQSpinBox;
  QDoubleSpinBoxNoTrailZeros* fieldACARelErrorQSpinBox;
};

#endif // PARAMETERDIALOG_H
