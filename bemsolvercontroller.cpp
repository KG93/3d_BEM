#include "bemsolvercontroller.h"

BemSolverController::BemSolverController(QObject* parent)
    : QThread{parent}
{
    this->setTerminationEnabled(true); // not safe yet
}

void BemSolverController::clear()
{
    terminateThread();
    boundariesSolvedAtFrequency = QVector<bool> (frequencies.size(), false);
    boundarySolutions.clear();
    maxSoundPressureOnBoundary.resize(0);
    minSoundPressureOnBoundary.resize(0);

    clearFields();

    solutionState = NothingSolved;
    controllerState = NotWorking;
}

void BemSolverController::clearFields()
{
    fieldsSolvedAtFrequency = QVector<bool> (frequencies.size(), false);
    fieldSolutions.clear();
    maxSoundPressureOnFields.resize(0);
    minSoundPressureOnFields.resize(0);
}

void BemSolverController::setFrequencies(const QVector<double> frequenciesArg)
{
    terminateThread();
    frequencies = frequenciesArg;
    const int n = frequencies.size();
    boundarySolutions.resize(n);
    boundariesSolvedAtFrequency.resize(n);
    maxSoundPressureOnBoundary.resize(n);
    minSoundPressureOnBoundary.resize(n);

    fieldSolutions.resize(n);
    fieldsSolvedAtFrequency.resize(n);
    maxSoundPressureOnFields.resize(n);
    minSoundPressureOnFields.resize(n);
    solutionState = NothingSolved;
}

void BemSolverController::run()
{
    controllerState = Working;
    solutionState = NothingSolved;

    boundaryElementSolver.setBemParameters(bemParameters);
    boundaryElementSolver.setBemSolverParameters(bemSolverParameters);
    boundaryElementSolver.setBoundaryElements(boundaryElements);
    boundaryElementSolver.setPointSources(pointSources);

    if(solveObjective == SolveBoundary || solveObjective == SolveBoundaryAndField)
    {
        boundariesSolvedAtFrequency = QVector<bool>(frequencies.size(), false);

        for(int i=0; i<frequencies.size(); i++)
        {
            boundaryElementSolver.setFrequency(frequencies.at(i));
            boundaryElementSolver.calculateBoundarySolution();
            setBoundarySolutions(i, boundaryElementSolver.getBoundaryElements().getSolution());
            boundaryElements = boundaryElementSolver.getBoundaryElements();
            boundariesSolvedAtFrequency[i] = true;
            solutionState = AtLeastOneSolved;
            emit boundarySolvedForFreqIndex(i);
        }
        std::pair<double, double> minMaxSoundPressure = getMinAndMaxSoundPressureOnBoundary();
        globalMinSoundPressureOnBoundary = minMaxSoundPressure.first;
        globalMaxSoundPressureOnBoundary = minMaxSoundPressure.second;

        boundaryElements = boundaryElementSolver.getBoundaryElements();
    }
    if(solveObjective == SolveField || solveObjective == SolveBoundaryAndField)
    {
        fieldsSolvedAtFrequency = QVector<bool>(frequencies.size(), false);
        boundaryElementSolver.setObservationFields(observationFields);

        for(int i=0; i<frequencies.size(); i++)
        {
            if(!isBoundarySolvedAtFreq(i))
            {
                continue;
            }
            boundaryElements.setSolution(boundarySolutions.at(i)); // set the to the frequency corresponding solution to the boundary elements
            boundaryElementSolver.setBoundaryElements(boundaryElements); // pass boundary elements to the solver
            boundaryElementSolver.setFrequency(frequencies.at(i));
            boundaryElementSolver.calculateFieldSolution(); // start field solver
            // save boundary solution
            setBoundarySolutions(i, boundaryElementSolver.getBoundaryElements().getSolution());

            setFieldSolutions(i, ObservationField::getFieldSolutions(boundaryElementSolver.getObservationFields()));

            fieldsSolvedAtFrequency[i] = true;
            solutionState = AtLeastOneSolved;
            emit fieldSolvedForFreqIndex(i);
        }
//        std::pair<double, double> minMaxSoundPressure = getMinAndMaxSoundPressureOnBoundary();
//        globalMinSoundPressureOnField = minMaxSoundPressure.first;
//        globalMaxSoundPressureOnField = minMaxSoundPressure.second;
        observationFields = boundaryElementSolver.getObservationFields();
    }
    controllerState = NotWorking;
}

std::pair<double, double> BemSolverController::getMinAndMaxSoundPressureOnBoundary(const long index)
{
    long n = frequencies.size();
    if(n > 0 && boundarySolutions.size() == n)
    {
        if(minSoundPressureOnBoundary.size() != n || maxSoundPressureOnBoundary.size() != n)
        {
            minSoundPressureOnBoundary.resize(n);
            maxSoundPressureOnBoundary.resize(n);

            #pragma omp parallel for
            for(int i=0; i<frequencies.size(); i++)
            {

                minSoundPressureOnBoundary(i) = boundarySolutions.at(i).soundPressure.cwiseAbs().minCoeff();
                maxSoundPressureOnBoundary(i) = boundarySolutions.at(i).soundPressure.cwiseAbs().maxCoeff();
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

BoundarySolution BemSolverController::getBoundarySolutionAtFrequencyIndex(const int freqIndex)
{
    if(freqIndex >= 0 && freqIndex < boundarySolutions.size())
    {
        return boundarySolutions.at(freqIndex);
    }
    else
    {
        std::cerr << "Invalid freqIndex in getBoundarySolutionAtFrequencyIndex() call." << std::endl;
        return BoundarySolution();
    }
}

FieldSolutions BemSolverController::getFieldSolutionAtFrequencyIndex(const int freqIndex)
{
    if(freqIndex >= 0 && freqIndex < fieldSolutions.size())
    {
        return fieldSolutions.at(freqIndex);
    }
    else
    {
        std::cerr << "Invalid freqIndex in getFieldSolutionAtFrequencyIndex() call." << std::endl;
        return FieldSolutions();
    }
}

bool BemSolverController::isBoundarySolvedAtFreq(const int freqIndex)
{
    if(freqIndex >= 0 && freqIndex < frequencies.size() && freqIndex < boundariesSolvedAtFrequency.size())
    {
        return boundariesSolvedAtFrequency.at(freqIndex);
    }
    else
    {
        std::cerr << "Invalid freqIndex passed to isBoundarySolvedAtFreq() call." << std::endl;
        return false;
    }
}

bool BemSolverController::isFieldSolvedAtFreq(const int freqIndex)
{
    if(freqIndex >= 0 && freqIndex < frequencies.size() && freqIndex < fieldsSolvedAtFrequency.size())
    {
        return fieldsSolvedAtFrequency.at(freqIndex);
    }
    else
    {
        std::cerr << "Invalid freqIndex passed to isFieldSolvedAtFreq() call." << std::endl;
        return false;
    }
}

//FieldSolutions BemSolverController::getFieldSolutions(const QVector<ObservationField> &fields)
//{
//    FieldSolutions fieldSolutions;
//    fieldSolutions.phiSolutions.resize(fields.size());
//    fieldSolutions.soundPressures.resize(fields.size());
//    for(int i=0; i<fields.size(); i++)
//    {
//        fieldSolutions.phiSolutions[i] = fields.at(i).phiSolution;
//        fieldSolutions.soundPressures[i] = fields.at(i).soundPressure;
//    }
//    return fieldSolutions;
//}

void BemSolverController::terminateThread()
{
    if(controllerState != NotWorking)
    {
        terminate();
        controllerState = NotWorking;
        std::cout << "Terminating worker thread!" << std::endl;
    }
}
