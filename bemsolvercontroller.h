#ifndef BEMSOLVERCONTROLLER_H
#define BEMSOLVERCONTROLLER_H

#include "boundaryelementsolver.h"

#include <QThread>

///**
//* \brief The solution to the BEM equations (the potential and its normal derivative) at the collocation points.
//*/
//struct BoundarySolution
//{
//    Eigen::VectorXcd phiSolution; /*!< The acoustic potential (is proportional to the sound pressure). */
//    Eigen::VectorXcd dPhiSolution; /*!< The normal derivative of the acoustic potential (is proportional to the particle velocity). */
//    Eigen::VectorXcd soundPressure; /*!< The sound pressure. */
//};

///**
//* \brief The solution to the BEM equations (the potential) at the field points.
//*/
//struct FieldSolutions
//{
//    QVector<Eigen::VectorXcd> phiSolutions; /*!< The acoustic potential (is proportional to the sound pressure). */
//    QVector<Eigen::VectorXcd> soundPressures; /*!< The sound pressure. */
//};

enum BemSolutionStates : int {NothingSolved, AtLeastOneSolved};
enum BemControllerStates : int {NotWorking, Working};
enum SolveObjective : int {SolveBoundary, SolveField, SolveBoundaryAndField};

class BemSolverController : public QThread
{
    Q_OBJECT
public:
    explicit BemSolverController(QObject* parent = nullptr);

    /**
    * \brief Set the boundary elements for the solver.
    * \param elements The boundary elements.
    */
    void setBoundaryElements(BoundaryElements elements) {this->boundaryElements = elements; solutionState = NothingSolved;}
    /**
    \brief Set the point sources for the solver.
    \param newPointSources A QVector of PointSource objects representing the point sources.
    */
    void setPointSources(QVector<PointSource> newPointSources) {this->pointSources = newPointSources; solutionState = NothingSolved;}
    /**
    * \brief Set observation fields.
    * \param observationFields The observation fields to set.
    */
    void setObservationFields(QVector<ObservationField> observationFields){this->observationFields = observationFields;}

    void setBemParameters(BemParameters parameters) {bemParameters = parameters; solutionState = NothingSolved;}
    void setBemSolverParameters(BemSolverParameters parameters) {bemSolverParameters = parameters;}
    void clear();
    void clearFields();
    void setFrequencies(const QVector<double> frequenciesArg);

    void setBoundarySolutions(const long i, const BoundarySolution solution)
    {
        if(boundarySolutions.size() != frequencies.size())
        {
            boundarySolutions.resize(frequencies.size());
            minSoundPressureOnBoundary.resize(frequencies.size());
            maxSoundPressureOnBoundary.resize(frequencies.size());
        }
        if(i >= 0 && i < boundarySolutions.length())
        {
            boundarySolutions[i] = solution;
        }
    }

    void setFieldSolutions(const long i, const FieldSolutions solutions)
    {
        if(fieldSolutions.size() != frequencies.size())
        {
            fieldSolutions.resize(frequencies.size());
            minSoundPressureOnFields.resize(frequencies.size());
            maxSoundPressureOnFields.resize(frequencies.size());
        }
        if(i >= 0 && i < fieldSolutions.length())
        {
            fieldSolutions[i] = solutions;
        }
    }

    void run() override;
    std::pair<double, double> getMinAndMaxSoundPressureOnBoundary(const long index = -1);
    BoundaryElements getBoundaryElements() {return boundaryElementSolver.getBoundaryElements();}
    QVector<ObservationField> getObservationFields() {return boundaryElementSolver.getObservationFields();}
    bool isBoundarySolvedAtFreq(const int freqIndex);
    bool isFieldSolvedAtFreq(const int freqIndex);
    BemSolutionStates getSolutionState() {return solutionState;}
    BemControllerStates getControllerState() {return controllerState;}
    BoundarySolution getBoundarySolutionAtFrequencyIndex(const int freqIndex);
    FieldSolutions getFieldSolutionAtFrequencyIndex(const int freqIndex);
    void setSolverObjective(SolveObjective objective) {solveObjective = objective;}

private:

    BemParameters bemParameters;
    BemSolverParameters bemSolverParameters;
    BoundaryElements boundaryElements;
    BoundaryElementSolver boundaryElementSolver;

    QVector<PointSource> pointSources;  /** \brief List of point sources for the simulation. */
    QVector<ObservationField> observationFields;

    QVector<double> frequencies;
    QVector<bool> boundariesSolvedAtFrequency;
    QVector<BoundarySolution> boundarySolutions;
    Eigen::VectorXd maxSoundPressureOnBoundary;
    double globalMaxSoundPressureOnBoundary;
    Eigen::VectorXd minSoundPressureOnBoundary;
    double globalMinSoundPressureOnBoundary;

    QVector<bool> fieldsSolvedAtFrequency;
    QVector<FieldSolutions> fieldSolutions;
    Eigen::VectorXd maxSoundPressureOnFields;
    double globalMaxSoundPressureOnField;
    Eigen::VectorXd minSoundPressureOnFields;
    double globalMinSoundPressureOnField;

    BemSolutionStates solutionState = NothingSolved;
    BemControllerStates controllerState = NotWorking;
    SolveObjective solveObjective = SolveBoundary;

public slots:
    void terminateThread();  /** Forcefully terminate the bem calculation thread. Not safe yet. Might crash. */
    void getUpdateLogSignalFromSolver(){emit updateLog();}
signals:
    void boundarySolvedForFreqIndex(int freqIndex);
    void fieldSolvedForFreqIndex(int freqIndex);
    void updateLog(); /*!< Update the GUI log. */
};

#endif // BEMSOLVERCONTROLLER_H
