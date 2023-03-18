//The  (regular LU-) boundary element solver is implemented after the text
//"the boundary element method in acoustics" by Stephen Kirkup
//1997/2007

#include "boundaryelementsolver.h"
//#include <valgrind/callgrind.h>
BoundaryElementSolver::BoundaryElementSolver()
{

}

void BoundaryElementSolver::solve()
{
    /////// for generating images /////////////
//    std::shuffle(boundaryElements.triangles.begin(), boundaryElements.triangles.end(), std::random_device());
    /////// for generating images /////////////
    global::trimMemory();
    std::cout<< "Number of OpenMP devices: "  << omp_get_num_devices() << std::endl;
//    //omp_set_default_device(1);
    std::cout<< "Number of default OpenMP device: "  << omp_get_default_device() << std::endl;
//    omp_set_max_active_levels(1);
//    std::cout<< "Number of active OpenMP levels: "  << omp_get_max_active_levels() << std::endl;

    QString message;
    message = QString("** Boundary element solver started. **");
    logStrings::logString.append("\r\n" + message);
    logStrings::errorLogString.append("\r\n" + message + "\r\n");
    std::cout << std::endl << message.toUtf8().constData() << std::endl;

    wavenumber = (2.0 * frequency * global::PI) / waveSpeed;
    wavenumberSquared = wavenumber * wavenumber;
    wavenumberSquaredHalf = wavenumberSquared / 2.0;
    iWavenumber = wavenumber * imaginaryUnit;
    setUpQuadratureRules();
    if(burtonMillerCoupling == true)
    {
        couplingParameter = imaginaryUnit / wavenumber;
        if(std::abs(wavenumber) < global::global::tiny)
        {
            couplingParameter = imaginaryUnit;
        }
        std::cout<<"Burton and Miller coupling with parameter: "<<couplingParameter<<std::endl;
        message = QString("Burton and Miller coupling with parameter: "+QString::number(std::real(couplingParameter))+" "+QString::number(std::imag(couplingParameter)));
    }
    else if(kirkupCoupling == true)
    {
        couplingParameter = imaginaryUnit / (wavenumber + 1.0);
        std::cout<<"Kirkup coupling with parameter: "<<couplingParameter<<std::endl;
        message = QString("Kirkup coupling with parameter: "+QString::number(std::real(couplingParameter))+" "+QString::number(std::imag(couplingParameter)));
    }
    else
    {
        couplingParameter = 0.0;
        message = QString("No coupling. Coupling parameter: "+QString::number(std::real(couplingParameter))+" "+QString::number(std::imag(couplingParameter)));
        std::cout<<"No coupling. Coupling parameter: "<<couplingParameter<<std::endl;
    }
    logStrings::logString.append("\r\n" + message + "\r\n");

    couplingParameterHalf = couplingParameter / 2.0;
    if(negativeCoupling)
    {
        message = QString("Coupling parameter sign set to negative.");
        std::cout << "Coupling parameter sign set to negative." << std::endl;
        logStrings::logString.append("\r\n" + message + "\r\n");
        setCouplingParameterNegative();
    }
    std::cout << "Number of point sources: "<<pointSources.size() << std::endl;
    message = QString("Number of point sources: "+QString::number(pointSources.size()));
    logStrings::logString.append(message + "\r\n");

    std::cout << "Wavenumber: " << wavenumber << std::endl;
    message = QString("Wavenumber: "+QString::number(std::real(wavenumber))+" "+QString::number(std::imag(wavenumber)));
    logStrings::logString.append(message+"\r\n");
    std::cout<<"Frequency: "<<frequency<<std::endl;
    message=QString("Frequency: "+QString::number(frequency));
    logStrings::logString.append(message+"\r\n");
//    std::cout<<"NumberOfRegularWeightsAndAbscissa: "<<numberOfRegularWeightsAndAbscissa<<std::endl;

    if(hMatSolving == true) // use H-matrix solver
    {
        Timer timer;


        timer.start();
        hMatrixSolve();
        timer.stop();
        std::cout << "Runtime of hMatrixSolve: " << timer.secs() << " seconds" << std::endl;

        if(global::activeProgramming) // some diagnostics for testing purposes
        {
            //        timer.start();
            //        hMatrixSolve(preconditionerRank, acaRelativeError);
            //        timer.stop();
            //        std::cerr << "Runtime of hMatrixSolve: " << timer.secs() << " seconds" << std::endl;
            //        timer.reset();

            //        regularSolveWithGMRES(0.01);  //compare with regular solver

            if(testACA)
            {
                regularSolveWithLU();
            }

            if(testDirect.size() == testPhiHMat.size() && testPhiHMat.size() >= 1) // for testing
            {
                long row;
                long col;
                double maxError = (testDPhiDirect - testDPhiHMat).cwiseAbs().maxCoeff(&row, &col);
                std::cerr <<"(test1 - test2).cwiseAbs().maxCoeff(): " << maxError << std::endl;
                std::cerr <<"(testDPhiDirect - testDPhiHMat).cwiseAbs()/ std::abs(testDPhiDirect(row, col)): " << maxError/ std::abs(testDPhiDirect(row, col))<< std::endl;
                std::cerr <<"row: " << row << std::endl;
                std::cerr <<"col: " << col << std::endl;

                std::cerr <<"phi norm accuracy: " << (testDirect - testPhiHMat).norm() / testDirect.norm() << std::endl;
                testDirect.resize(0,0);
                testPhiHMat.resize(0,0);
                std::cerr <<"dPhi norm accuracy: " << (testDPhiDirect - testDPhiHMat).norm() / testDPhiDirect.norm() << std::endl;
                testDPhiDirect.resize(0,0);
                testDPhiHMat.resize(0,0);
            }
            if(testSol1.size() == testSol2.size())
            {
                std::cerr <<"solution accuracy: " << (testSol1 - testSol2).norm() / testSol1.norm() << std::endl;
            }
            if(testSolLU.size() == testSol1.size() && testSol1.size() != 0)
            {
                std::cerr <<"testSolLU accuracy: " << (testSolLU - testSol1).norm() / testSol1.norm() << std::endl;
            }
            if(testSolRandGMRES.size() == testSol1.size() && testSol1.size() != 0)
            {
                std::cerr <<"testSolRandGMRES accuracy: " << (testSolRandGMRES - testSol1).norm() / testSol1.norm() << std::endl;
            }
            if(directRightHandSide.size() == HMatRightHandSide.size())
            {
                std::cerr <<"HMatRightHandSide accuracy: " << (HMatRightHandSide - directRightHandSide).norm() / directRightHandSide.norm() << std::endl;
            }
        }
    }
    else //use  regular matrix solver
    {
        regularSolveWithLU();
    }
    if(frequency != 0) // calculate sound pressure from acoustic potential
    {
        boundaryElements.soundPressure = (imaginaryUnit*airDensity * PI2 * frequency) * boundaryElements.phiSolution;
    }
    else // laplace equation
    {
        boundaryElements.soundPressure = (imaginaryUnit*airDensity * PI2) * boundaryElements.phiSolution;
    }
    emit updateLog();
}

void BoundaryElementSolver::regularSolveWithLU()
{
    int numberOfElements = boundaryElements.triangles.length();
    if(numberOfElements == 0)
    {
        std::cout<<"No elements!"<<std::endl;
        QString message=QString("No elements!");
        logStrings::logString.append(message+"\r\n");
        return;
    }

    prepareBoundaryElements();
    calculateBoundaryConditionAndSourceTermVector();

    if(polarIntegOfSingOps && frequency != 0)
    {
        std::cout<<"Number of elements in solver with polar integration: "<<numberOfElements<<std::endl;
        QString message=QString("Number of elements in solver with polar integration: "+QString::number(numberOfElements));
        logStrings::logString.append(message+"\r\n");
    }
    else
    {
        std::cout<<"Number of elements in solver with singularity substraction: "<<numberOfElements<<std::endl;
        QString message=QString("Number of elements in solver with singularity substraction: "+QString::number(numberOfElements));
        logStrings::logString.append(message+"\r\n");
    }
    Timer timer;
    timer.start();
    std::cout<<"Start of matrix initialization."<<std::endl;
    QString message=QString("Start of matrix initialization.");
    logStrings::logString.append(message+"\r\n");
    emit updateLog();
    Eigen::MatrixXcd matrixPhi(numberOfElements,numberOfElements);
    Eigen::MatrixXcd matrixDPhi(numberOfElements,numberOfElements);

    std::complex<double> Lk;
    std::complex<double> Mk;
    std::complex<double> Mtk;
    std::complex<double> Nk;
    bool onPanel = false;

    if(polarIntegOfSingOps && frequency != 0)
    {
        if(!avoidGLS)
        {
            std::cerr << "only for testing." << std::endl;
            #pragma omp parallel for private(Mk,Nk,Lk,Mtk,onPanel)// parallelizes matrix initialization
            for(int column = 0; column < numberOfElements; column++)
            {
                for(int row = 0; row < numberOfElements; row++)
                {
                    BemOperatorsConst(row,column,Lk,Mk,Mtk,Nk);

                    matrixPhi(row,column) = Mk + couplingParameter * Nk;
                    matrixDPhi(row,column) = Lk + couplingParameter * Mtk;
                }
            }
            matrixPhi.diagonal() -= 0.5 * Eigen::VectorXcd::Ones(numberOfElements);
            matrixDPhi.diagonal() += couplingParameterHalf * Eigen::VectorXcd::Ones(numberOfElements);
        }
        else
        {
            std::cout << "Implicit substitutions." << std::endl;
            #pragma omp parallel for private(Mk,Nk,Lk,Mtk,onPanel)// parallelizes matrix initialization
            for(int column = 0; column < numberOfElements; column++)
            {
                for(int row = 0; row < numberOfElements; row++)
                {
                    BemOperatorsConst(row,column,Lk,Mk,Mtk,Nk);
                    onPanel = (row == column);
                    if(substituteDPhiWithPhi(column))
                    {
                        matrixPhi(row,column) = - (Lk + couplingParameter * Mtk + ((double)onPanel) * couplingParameterHalf) / beta(column);
                        matrixDPhi(row,column) = - alpha(column) / beta(column) * (Lk + couplingParameter * Mtk + ((double)onPanel) * couplingParameterHalf) - (Mk + couplingParameter * Nk - ((double)onPanel) * 0.5);
                    }
                    else
                    {
                        matrixPhi(row,column) = (Mk + couplingParameter * Nk - ((double)onPanel) * 0.5) / alpha(column);
                        matrixDPhi(row,column) = Lk + couplingParameter * Mtk + ((double)onPanel) * couplingParameterHalf + beta(column)/alpha(column) * (Mk + couplingParameter * Nk - ((double)onPanel) * 0.5);
                    }
                }
            }
        }
    }
    else // treat the singularitias with singularity subtraction instead of polar integration (which isn't implemented for f=0)
    {
        #pragma omp parallel for private(Mk,Nk,Lk,Mtk,onPanel)// parallelizes matrix initialization
        for(int column=0;column<numberOfElements;column++)
        {
            for(int row=0;row<numberOfElements;row++)
            {
                if(row==column)
                {
                    BemOperatorsSingularitySubtraction(row,column,Lk,Mk,Mtk,Nk);
                }
                else
                {
                    BemOperatorsConst(row,column,Lk,Mk,Mtk,Nk);
                }

                matrixPhi(row,column) = Mk + couplingParameter*Nk;
                matrixDPhi(row,column) = Lk + couplingParameter*Mtk;
            }
        }
        matrixPhi.diagonal() -= 0.5 * Eigen::VectorXcd::Ones(numberOfElements);
        matrixDPhi.diagonal() += couplingParameterHalf * Eigen::VectorXcd::Ones(numberOfElements);
    }

    timer.stop();
    std::cout << "Runtime of matrix initialization: " << timer.secs() << " seconds" << std::endl;
    message=QString("Runtime of matrix initialization: " + QString::number(timer.secs())+" seconds.");
    logStrings::logString.append(message + "\r\n");
    emit updateLog();
    timer.reset();
    if(testACA)
    {
        testDirect = matrixPhi; //tmp
        matrixPhi.resize(0,0);
        testDPhiDirect = matrixDPhi;
        matrixDPhi.resize(0,0);
        return;
    }
    timer.start();

    /////////////////////////////////////// generate matrix imaga for presentation ////////////////////////////////
//    hMatrixVisuals::matrixPhaseImage(matrixPhi);
//    hMatrixVisuals::matrixToGreyscaleImageHigherContrast(matrixPhi);
//    hMatrixVisuals::matrixPhaseImage(matrixDPhi);
//    hMatrixVisuals::matrixToGreyscaleImageHigherContrast(matrixDPhi);
    //////////////////////////////////////////////// end //////////////////////////////////////////////////////////

    boundaryElements.phiSolution = Eigen::VectorXcd (numberOfElements);
    boundaryElements.dPhiSolution = Eigen::VectorXcd (numberOfElements);

//    debug = false;
//    if(debug)
//    {
//        std::cout << "Sources Vector: " <<sourceTermVector << std::endl;
//        std::cout << "Phi Matrix: " <<matrixPhi << std::endl;
//        std::cout << "dPhi Matrix: " <<matrixDPhi << std::endl;
//        std::cout << "alpha: " <<alpha << std::endl;
//        std::cout << "beta: " <<beta << std::endl;
//        std::cout << "f: " <<f << std::endl;
//    }
//    std::cerr << "matrix dphi frobenius condition number: " << matrixDPhi.norm() * matrixDPhi.inverse().norm() << std::endl;
//    std::cerr << "matrix phi frobenius condition number: " << matrixPhi.norm() * matrixPhi.inverse().norm() << std::endl;
    if(avoidGLS)
    {
        Eigen::PartialPivLU<Eigen::MatrixXcd> luPP(matrixDPhi);
        directRightHandSide = sourceTermVector + matrixPhi * f;
        boundaryElements.dPhiSolution = luPP.solve(directRightHandSide);
        boundaryElements.phiSolution = boundaryElements.dPhiSolution;
        for(long i = 0; i < substituteDPhiWithPhi.size(); i++)
        {
            if (substituteDPhiWithPhi(i))
            {
                boundaryElements.dPhiSolution(i) = (f(i) - alpha(i) * boundaryElements.phiSolution(i)) / beta(i);
            }
            else
            {
                boundaryElements.phiSolution(i) = (f(i) - beta(i) * boundaryElements.dPhiSolution(i)) / alpha(i);
            }
        }
    }
    else
    {
        GLS::calculateGLS(matrixDPhi, boundaryElements.dPhiSolution, matrixPhi, boundaryElements.phiSolution, sourceTermVector, beta, alpha, f);
    }
    timer.stop();
    std::cout << "Runtime of gls solver: " << timer.secs() << " seconds" << std::endl;
    message=QString("Runtime of gls solver: " + QString::number(timer.secs()) + " seconds.");
    logStrings::logString.append(message+"\r\n");
    emit updateLog();

    bool showPhiSolution = false;
    if(showPhiSolution)
    {
        std::cout << "Phi Matrix: " << matrixPhi(0,0) << std::endl;
        std::cout << "dPhi Matrix: " << matrixDPhi(0,0) << std::endl;
        std::cout << "sources: " << sourceTermVector(0) << std::endl;

        boundaryElements.soundPressure = (imaginaryUnit*airDensity*PI2*frequency) * boundaryElements.phiSolution;
        for(int i=0;i<numberOfElements;i++)
        {
           std::cout<<"index: "<< boundaryElements.triangles.at(i).elementIndex<<" phi: "<<boundaryElements.phiSolution(i)<<" velocity: "<<boundaryElements.dPhiSolution(i)<<"sound pressure: "<<boundaryElements.soundPressure(i)<<std::endl;
        }
    }
    // for testing
    testSol1 = boundaryElements.phiSolution;

    matrixPhi.resize(0, 0);
    matrixDPhi.resize(0, 0);
}

// for benchmarking; full matrix method withe iterative solver
void BoundaryElementSolver::regularSolveWithGMRES(double error)
{
    int numberOfElements = boundaryElements.triangles.length();
    if(numberOfElements == 0)
    {
        std::cout<<"No elements!"<<std::endl;
        QString message=QString("No elements!");
        logStrings::logString.append(message+"\r\n");
        return;
    }

    std::cout<<"Number of elements in regularSolveWithGMRES: "<<numberOfElements<<std::endl;
    QString message=QString("Number of elements in regularSolveWithGMRES: "+QString::number(numberOfElements));
    logStrings::logString.append(message+"\r\n");

    Timer timer;
    timer.start();
    std::cout<<"Start of matrix initialization."<<std::endl;
    message=QString("Start of matrix initialization.");
    logStrings::logString.append(message+"\r\n");
    emit updateLog();

    prepareBoundaryElements();
    calculateBoundaryConditionAndSourceTermVector();

    Eigen::MatrixXcd matrixPhi(numberOfElements,numberOfElements);
    Eigen::MatrixXcd matrixDPhi(numberOfElements,numberOfElements);
    std::complex<double> Mk;
    std::complex<double> Nk;
    std::complex<double> Lk;
    std::complex<double> Mtk;
    bool onPanel = false;

    // assamble system matrices
    #pragma omp parallel for private(Mk,Nk,Lk,Mtk,onPanel)// parallelizes matrix initialization
    for(int column = 0; column < numberOfElements; column++)
    {
        for(int row = 0; row < numberOfElements; row++)
        {
            BemOperatorsConst(row,column,Lk,Mk,Mtk,Nk);
            onPanel = (row == column);

            if(substituteDPhiWithPhi(column))
            {
                matrixPhi(row,column) = -(Lk + couplingParameter * Mtk + ((double)onPanel) * couplingParameterHalf) / beta(column);
                matrixDPhi(row,column) = -alpha(column) / beta(column) * (Lk + couplingParameter * Mtk + ((double)onPanel) * couplingParameterHalf) - (Mk + couplingParameter * Nk - ((double)onPanel) * 0.5);
            }
            else
            {
                matrixPhi(row,column) = (Mk + couplingParameter * Nk - ((double)onPanel) * 0.5) / alpha(column);
                matrixDPhi(row,column) = Lk + couplingParameter * Mtk + ((double)onPanel) * couplingParameterHalf + beta(column)/alpha(column) * (Mk + couplingParameter * Nk - ((double)onPanel) * 0.5);
            }
        }
    }

    timer.stop();
    std::cout << "Runtime of matrix initialization: " << timer.secs() << " seconds" << std::endl;
    message=QString("Runtime of matrix initialization: " + QString::number(timer.secs())+" seconds.");
    logStrings::logString.append(message + "\r\n");
    emit updateLog();
//    std::cerr << "matrix dphi frobenius condition number: " << matrixDPhi.norm() * matrixDPhi.inverse().norm() << std::endl;
//    std::cerr << "matrix phi frobenius condition number: " << matrixPhi.norm() * matrixPhi.inverse().norm() << std::endl;
    timer.reset();
    if(testACA)
    {
        testDirect = matrixPhi; //tmp
        matrixPhi.resize(0,0);
        testDPhiDirect = matrixDPhi;
        matrixDPhi.resize(0,0);
        return;
    }

    boundaryElements.phiSolution = Eigen::VectorXcd (numberOfElements);
    boundaryElements.dPhiSolution = Eigen::VectorXcd (numberOfElements);

    timer.start();
    directRightHandSide = sourceTermVector + matrixPhi * f;
    boundaryElements.dPhiSolution = GMRES::gmresSolve(matrixDPhi, Eigen::VectorXcd::Random(boundaryElements.triangles.length()), directRightHandSide, 0.001*error);
    boundaryElements.phiSolution = boundaryElements.dPhiSolution;
    for(long i = 0; i < substituteDPhiWithPhi.size(); i++)
    {
        if (substituteDPhiWithPhi(i))
        {
            boundaryElements.dPhiSolution(i) = (f(i) - alpha(i) * boundaryElements.phiSolution(i)) / beta(i);
        }
        else
        {
            boundaryElements.phiSolution(i) = (f(i) - beta(i) * boundaryElements.dPhiSolution(i)) / alpha(i);
        }
    }
    timer.stop();
    std::cout << "Runtime of full GMRES solver: " << timer.secs() << " seconds" << std::endl;
    message=QString("Runtime of full GMRES solver: " + QString::number(timer.secs()) + " seconds.");
    logStrings::logString.append(message+"\r\n");
    emit updateLog();

    bool showPhiSolution = false;
    if(showPhiSolution)
    {
        std::cout << "Phi Matrix: " << matrixPhi(0,0) << std::endl;
        std::cout << "dPhi Matrix: " << matrixDPhi(0,0) << std::endl;
        std::cout << "sources: " << sourceTermVector(0) << std::endl;

        boundaryElements.soundPressure = (imaginaryUnit*airDensity*PI2*frequency) * boundaryElements.phiSolution;
        for(int i=0;i<numberOfElements;i++)
        {
           std::cout<<"index: "<< boundaryElements.triangles.at(i).elementIndex<<" phi: "<<boundaryElements.phiSolution(i)<<" velocity: "<<boundaryElements.dPhiSolution(i)<<"sound pressure: "<<boundaryElements.soundPressure(i)<<std::endl;
        }
    }
    ////////// for testing
    ///
    testSol1 = boundaryElements.phiSolution;


    matrixPhi.resize(0, 0);
    matrixDPhi.resize(0, 0);
}

//void BoundaryElementSolver::regularSolveWithLULinear()
//{
//    LinearBoundaryElements boundaryElements;
//    int numberOfNodes = boundaryElements.nodes.length();
//    int numberOfTriangles = boundaryElements.triangles.length();
////    int numberOfCollocationPoints = numberOfTriangles;

//    if(numberOfNodes == 0)
//    {
//        std::cout<<"No elements!"<<std::endl;
//        QString message=QString("No elements!");
//        logStrings::logString.append(message+"\r\n");
//        return;
//    }

//    std::cout<<"Number of elements in regularSolveWithGMRES: "<<numberOfNodes<<std::endl;
//    QString message=QString("Number of elements in regularSolveWithGMRES: "+QString::number(numberOfNodes));
//    logStrings::logString.append(message+"\r\n");

//    Timer timer;
//    timer.start();
//    std::cout<<"Start of matrix initialization."<<std::endl;
//    message=QString("Start of matrix initialization.");
//    logStrings::logString.append(message+"\r\n");
//    emit updateLog();

////    prepareBoundaryElements();
////    calculateBoundaryConditionAndSourceTermVector();

////    boundaryElements.calculateAssociatedTriangles();

//    Eigen::MatrixXcd matrixPhi = Eigen::MatrixXcd::Zero(numberOfTriangles,numberOfNodes);
//    Eigen::MatrixXcd matrixDPhi = Eigen::MatrixXcd::Zero(numberOfTriangles,numberOfNodes);
//    std::complex<double> Mk;
//    std::complex<double> Nk;
//    std::complex<double> Lk;
//    std::complex<double> Mtk;


//    #pragma omp parallel for private(Mk,Nk,Lk,Mtk)// parallelizes matrix initialization
//    for(int i=0; i< numberOfTriangles; i++)
//    {
//        LinearTriangle colloPoiTriangle = boundaryElements.triangles.at(i);
//        matrixPhi(i, colloPoiTriangle.node1Index) -= 1/6.0;
//        matrixPhi(i, colloPoiTriangle.node2Index) -= 1/6.0;
//        matrixPhi(i, colloPoiTriangle.node3Index) -= 1/6.0;

//        matrixDPhi(i, colloPoiTriangle.node1Index) += couplingParameterHalf/3.0;
//        matrixDPhi(i, colloPoiTriangle.node2Index) += couplingParameterHalf/3.0;
//        matrixDPhi(i, colloPoiTriangle.node3Index) += couplingParameterHalf/3.0;

//        for(int j=0; j<numberOfTriangles; j++)
//        {
//            LinearTriangle tmpTriangle = boundaryElements.triangles.at(j);
//            long node1Index = tmpTriangle.node1Index;
//            long node2Index = tmpTriangle.node2Index;
//            long node3Index = tmpTriangle.node3Index;

//            bool collocationPointOnTriangle = (i==j);

//            BemOperatorsLinearNode1(i,j,Lk,Mk,Mtk,Nk,boundaryElements, collocationPointOnTriangle);
//            matrixPhi(i,node1Index) = Mk + couplingParameter*Nk;
//            matrixDPhi(i,node1Index) = Lk + couplingParameter*Mtk;

//            BemOperatorsLinearNode2(i,j,Lk,Mk,Mtk,Nk,boundaryElements, collocationPointOnTriangle);
//            matrixPhi(i,node2Index) = Mk + couplingParameter*Nk;
//            matrixDPhi(i,node2Index) = Lk + couplingParameter*Mtk;

//            BemOperatorsLinearNode3(i,j,Lk,Mk,Mtk,Nk,boundaryElements, collocationPointOnTriangle);
//            matrixPhi(i,node3Index) = Mk + couplingParameter*Nk;
//            matrixDPhi(i,node3Index) = Lk + couplingParameter*Mtk;
//        }
//    }
//}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//void BoundaryElementSolver::hMatrixSolve()
//{
//    randomGenerator.seed(QRandomGenerator::system()->generate());

////    boundaryElements.calculateTriangleMidPoints();
//    std::cout<<"Number of elements in H-matrix solve: " << boundaryElements.triangles.length() << std::endl;

//    ClusterTree clusterTree(&boundaryElements);

//    HMatrix phiHMat;
//    HMatrix dPhiHMat;

//    phiHMat.setClusterTrees(&clusterTree, &clusterTree);
//    phiHMat.populateBlockClusterTree();

//    dPhiHMat.setClusterTrees(&clusterTree, &clusterTree);
//    dPhiHMat.populateBlockClusterTree();

//    prepareBoundaryElements();
//    calculateBoundaryConditionAndSourceTermVector();

//    QVector<BlockCluster*> phiHMatPartition = phiHMat.getMinPartition();
//    QVector<BlockCluster*> dPhiHMatPartition = dPhiHMat.getMinPartition();
//    if(phiHMatPartition.size() != dPhiHMatPartition.size())
//    {
//        std::cerr<< "phiHMatPartition.size() != dPhiHMatPartition.size() in hMatrixSolve()" << std::endl;
//    }

//    Timer timer;
//    timer.start();
//    std::cout<< "Started the block assembly." << std::endl;
//    #pragma omp parallel master
//    for(long blockIndex = 0; blockIndex < phiHMatPartition.size(); blockIndex++)
//    {
//        #pragma omp task
//        hBlockAssembly(phiHMatPartition.at(blockIndex), dPhiHMatPartition.at(blockIndex), acaMaxRank, acaRelativeError);
//    }
//    timer.stop();
//    std::cout << "runtime of hBlockAssembly: " << timer.secs() << std::endl;
//    timer.reset();

//    double normPhiResiduum = 0;
//    if(calculateNormAndCond)
//    {
//        timer.reset();
//        timer.start();
//        normPhiResiduum = HArithm::spectralNormLeastSignificant(dPhiHMat); // Approximates the spectral norm of the residuum of dPhiHMat
//        timer.stop();
//        std::cout << "Spectral norm of the dphi residuum: " << normPhiResiduum << std::endl;
//        std::cout << "Runtime of norm calculations: " << timer.secs() << std::endl;
//    }

//    timer.start();
//    double compressionRatio = phiHMat.getCompressionRatio();
//    timer.stop();
//    std::cout << std::endl << "phi compressionRatio: " << compressionRatio <<std::endl << "calculated in " << timer.secs() << " seconds" << std::endl;
//    timer.reset();

//    timer.start();
//    double dPhiCompressionRatio = dPhiHMat.getCompressionRatio();
//    timer.stop();
//    std::cout << "dPhi compressionRatio: " << dPhiCompressionRatio <<std::endl << "calculated in " << timer.secs() << " seconds" << std::endl;
//    timer.reset();
//    std::cout << "maximum local rank of phi matrix: " << phiHMat.calculateMaxBlockRank() << std::endl;
//    std::cout << "maximum local rank of dPhi matrix: " << dPhiHMat.calculateMaxBlockRank() << std::endl;

//    ////////////////////////////////////// Copy hMat for accuracy comparison
//    if(testACA)
//    {
//        testDPhiHMat = HArithm::hMatToFullMat(dPhiHMat);
//        testPhiHMat = HArithm::hMatToFullMat(phiHMat);

////        phiHMat.clear(true);
////        dPhiHMat.clear(true);
////        return;
//    }
//    ////////////////////////////////////// end Copy hMat for accuracy comparison

//    timer.start();
//    //re-compress matrices with lower relative error
//    HArithm::compressHMat(phiHMat, acaMaxRank, 0.1 * acaRelativeError);
//    HArithm::compressHMat(dPhiHMat, acaMaxRank, 0.1 * acaRelativeError);
////    std::cout << "phi.norm(): " << phiHMat.norm() << std::endl;
//    std::cout << "frobenius norm of phi matrix: " << phiHMat.normInSVDForm() << std::endl;
//    timer.stop();
//    compressionRatio = phiHMat.getCompressionRatio();
//    std::cout << std::endl << "phi compressionRatio: " << compressionRatio << " after additional compression calculated in " << timer.secs() << " seconds." << std::endl;
//    dPhiCompressionRatio = dPhiHMat.getCompressionRatio();
//    std::cout << "dPhi compressionRatio: " << dPhiCompressionRatio << " after additional compression calculated in " << timer.secs() << " seconds." << std::endl;
//    timer.reset();
//    std::cout << std::endl << "maximum local rank of phi matrix: " << phiHMat.calculateMaxBlockRank() << std::endl;
//    std::cout << "maximum local rank of dPhi matrix: " << dPhiHMat.calculateMaxBlockRank() << std::endl;

//    ////////////////////////////////////// Calculate the righthand side for the linear system

//    Eigen::VectorXcd rightHandside = sourceTermVector;
//    HArithm::MVM(rightHandside, phiHMat, f);
//    Eigen::VectorXcd rightHandsideCopy = rightHandside;
//    phiHMat.clear(true); // phimat in not needed anymore

//    HMatRightHandSide = rightHandside;
//    //////////////////////////////////////
//    if(boundaryElements.impedancePlanes.length() == 1)
//    {
//        std::cout << std::endl << "Calculating reflection matrices for a single impedance planes." << std::endl;
//        BoundaryElements reflectedBoundaryElements = boundaryElements;
//        reflectedBoundaryElements.reflectGeometry(boundaryElements.impedancePlanes.at(0));

//        ClusterTree reflectedClusterTree(clusterTree);  //copy the original clustertree
//        reflectedClusterTree.updateMinCuboids(&reflectedBoundaryElements);

//        HMatrix halfSpaceReflMatPhi;
//        halfSpaceReflMatPhi.setClusterTrees(&clusterTree, &reflectedClusterTree);
//        halfSpaceReflMatPhi.populateBlockClusterTree();

//        HMatrix halfSpaceReflMatDPhi;
//        halfSpaceReflMatDPhi.setClusterTrees(&clusterTree, &reflectedClusterTree);
//        halfSpaceReflMatDPhi.populateBlockClusterTree();

//        QVector<BlockCluster*> reflectedPhiHMatPartition = halfSpaceReflMatPhi.getMinPartition();
//        QVector<BlockCluster*> reflectedDPhiHMatPartition = halfSpaceReflMatDPhi.getMinPartition();
//        if(reflectedPhiHMatPartition.size() != reflectedDPhiHMatPartition.size())
//        {
//            std::cerr<< "reflectedPhiHMatPartition.size() != reflectedDPhiHMatPartition.size() in hMatrixSolve()" << std::endl;
//        }
//        std::cout<< "Started the block assembly of the reflection Matrix" << std::endl;
//        timer.start();
//        std::cout<< "reflectedPhiHMatPartition.size(): " << reflectedPhiHMatPartition.size() << std::endl;

//        #pragma omp parallel master
//        for(long blockIndex = 0; blockIndex < reflectedPhiHMatPartition.size(); blockIndex++)
//        {
//            #pragma omp task
//            hBlockAssemblyReflMat(reflectedPhiHMatPartition.at(blockIndex), reflectedDPhiHMatPartition.at(blockIndex), acaMaxRank, acaRelativeError, boundaryElements, reflectedBoundaryElements);
//        }
//        std::cout<< "Finished the block assembly of the reflection Matrix" << std::endl;
//        timer.stop();
//        std::cout << "Runtime of hBlockAssemblyReflectedMatrix: " << timer.secs() << std::endl;
//        timer.reset();

//        HArithm::MVM(rightHandside, halfSpaceReflMatPhi, f); // for fully reflective plane
////        HArithm::MVM(rightHandside, halfSpaceReflMatPhi, -f); // for absorbtive plane
//        halfSpaceReflMatPhi.clear(true); // halfSpaceReflMatPhi in not needed anymore
//        rightHandsideCopy = rightHandside;
//        HMatRightHandSide = rightHandside;

//        timer.start();
//        #pragma omp parallel master
//        {
//            HArithm::recursiveHMatAddition(dPhiHMat.getRootBlock(), halfSpaceReflMatDPhi.getRootBlock(), acaMaxRank, 0.1 * acaRelativeError); // for fully reflective plane
//        }
//        timer.stop();
//        std::cout << "Runtime of recursiveHMatAddition: " << timer.secs() << std::endl;
//        timer.reset();

////        HArithm::recursiveHMatSubstraction(dPhiHMat.getRootBlock(), halfSpaceReflMatDPhi.getRootBlock(), 0, 0.1 * error); // for absorbtive plane
//        dPhiCompressionRatio = halfSpaceReflMatDPhi.getCompressionRatio();
//        std::cout << "halfSpaceReflMatDPhi compressionRatio: " << dPhiCompressionRatio << ". " << std::endl;
//        halfSpaceReflMatDPhi.clear(true); // halfSpaceReflMatDPhi in not needed anymore
//        reflectedClusterTree.clear();
//        dPhiCompressionRatio = dPhiHMat.getCompressionRatio();
//        std::cout << "dPhi compressionRatio: " << dPhiCompressionRatio << " after additional compression calculated in " << timer.secs() << " seconds." << std::endl;
//        std::cout << "maximum local rank of dPhi matrix: " << dPhiHMat.calculateMaxBlockRank() << std::endl;
//    }
//    else if(boundaryElements.impedancePlanes.length() >= 2)
//    {
//        std::cout << "Calculating reflection matrices for two parallel impedance planes." << std::endl;
//        Timer refMatTimer;
//        refMatTimer.start();
//        HMatrix halfSpaceReflMatPhi;
//        HMatrix halfSpaceReflMatDPhi;

//        calcReflectionMatrices(halfSpaceReflMatPhi, halfSpaceReflMatDPhi, acaMaxRank, acaRelativeError, boundaryElements, boundaryElements, clusterTree);
//        refMatTimer.stop();
//        std::cout << "Runtime of new reflection matrix calculation: " << refMatTimer.secs() << std::endl;
//        std::cout << "Compression ratio of the dPhi reflection Matrix: " << halfSpaceReflMatDPhi.getCompressionRatio() << std::endl;
//        HArithm::MVM(rightHandside, halfSpaceReflMatPhi, f);
//        rightHandsideCopy = rightHandside;
//        HMatRightHandSide = rightHandside;
//        #pragma omp parallel master
//        {HArithm::recursiveHMatAddition(dPhiHMat.getRootBlock(), halfSpaceReflMatDPhi.getRootBlock(), acaMaxRank, 0.1 * acaRelativeError);}
//        halfSpaceReflMatPhi.clear(); // halfSpaceReflMatDPhi in not needed anymore
//        halfSpaceReflMatDPhi.clear(true); // halfSpaceReflMatDPhi in not needed anymore
//    }
//    ////////////////////////////////////// Copy dPhiHMat for later GMRES usage, as LUSolver modifies th original dPhiHMat

//    HMatrix dPhiHMatCopy;
//    std::pair<HMatrix,HMatrix> LUPair;
//    if(usePreconditioner)
//    {
//        timer.start();
//        dPhiHMatCopy = HMatrix(dPhiHMat, preconditionerRank, preconditionerRelError, true);
//        timer.stop();
//        std::cout << std::endl << "dPhiHMatCopy(dPhiHMat) took " << timer.secs() << " seconds." << std::endl;
//        timer.reset();
//        timer.start();
//        compressionRatio = dPhiHMatCopy.getCompressionRatio();
//        timer.stop();
//        std::cout << std::endl << "dPhiHMatCopy compressionRatio: " << compressionRatio <<std::endl << "calculated in " << timer.secs() << " seconds" << std::endl;
//        timer.reset();

//        timer.start();
//        std::cout << "Calculating the HLU preconditioner." << std::endl;
//        LUPair = HArithm::LUDecomposition(dPhiHMatCopy, preconditionerRank, preconditionerRelError); // calculate the LU decomposition of dPhiHMatCopy
//        dPhiHMatCopy.clear(true); // delete low rank copy
//        timer.stop();
//        std::cout << "Runtime of LUDecomposition: " << timer.secs() << std::endl;
//        timer.reset();
//    }

//    ////////////////////////////////////// Testing section start ////////////////////////////////
//    if(testACA)
//    {
//        phiHMat.clear(true);
//        dPhiHMat.clear(true);
//        return;
//    }
//    ////////////////////////////////////// Testing section end ////////////////////////////////

//    if(LUPair.first.hasNan() || LUPair.second.hasNan() || usePreconditioner == false) // LU decomposition has NaN; -> don't use the preconditioner
//    {
//        if(usePreconditioner == true)
//        {
//            std::cerr << " There are NaNs in the LU decomposition. Continuing without preconditioner." << std::endl;
//        }
//        timer.start();
//        HMatrixWrapper dPhiHMatForGMRES(&dPhiHMat); // set up a wrapper class for the hmatrix to be used in the matrix free gmres routine
//        boundaryElements.dPhiSolution = GMRES::gmresSolve(dPhiHMatForGMRES, Eigen::VectorXcd::Random(boundaryElements.triangles.length()), rightHandsideCopy, 0.001*acaRelativeError); // solve dPhiHMat * x = righthandside
//        timer.stop();
//        std::cout << "Runtime of GMRES: " << timer.secs() << std::endl;
//        timer.reset();
//        if(calculateNormAndCond)
//        {
//            timer.start();
//            double normPhi = HArithm::spectralNorm(dPhiHMat); // Approximate the spectral norm of dPhiHMat via power iteration
//            timer.stop();
//            std::cout << "Spectral norm of the dphi matrix: " << normPhi << std::endl;
//            std::cout << "Relative spectral approximation error of the dphi matrix: " << normPhiResiduum/normPhi << std::endl;
//            std::cout << "Runtime of norm calculation: " << timer.secs() << std::endl;
//        }
//        Eigen::VectorXcd b = Eigen::VectorXcd::Zero(boundaryElements.dPhiSolution.size());
//        HArithm::MVM(b, dPhiHMat, boundaryElements.dPhiSolution);
//        std::cout << "Relative error of the residuum: " << (b - rightHandsideCopy).norm()/ rightHandsideCopy.norm() << std::endl;
//    }
//    else
//    {
//        Eigen::VectorXcd gmresGuess = HArithm::LUSubstitutionSolve(LUPair.first, LUPair.second, rightHandsideCopy);
//        timer.start();
//        LUPrecondidionedHMatrixWrapper dPhiHMatForLUPreconditionedGMRES(&dPhiHMat, &LUPair.first, &LUPair.second);
//        boundaryElements.dPhiSolution = GMRES::gmresLUPreconditionedSolve(dPhiHMatForLUPreconditionedGMRES, gmresGuess, rightHandsideCopy, 0.001*acaRelativeError);
//        timer.stop();
//        std::cout << "Runtime of preconditioned GMRES: " << timer.secs() << std::endl;
//        timer.reset();

//        Eigen::VectorXcd b = Eigen::VectorXcd::Zero(boundaryElements.dPhiSolution.size());
//        HArithm::MVM(b, dPhiHMat, boundaryElements.dPhiSolution);
//        std::cout << "Relative error of the residuum: " << (b - rightHandsideCopy).norm()/ rightHandsideCopy.norm() << std::endl;
//        if(calculateNormAndCond)
//        {
//            timer.start();
//            double normPhi = HArithm::spectralNorm(dPhiHMat); // Approximates the spectral norm of dPhiHMat via power iteration
//            double normPhiInverse = HArithm::spectralNormFromLU(LUPair.first, LUPair.second); // Approximates the spectral norm of the inverse of of dPhiHMat via power iteration with its LU decomposition.
//            timer.stop();
//            std::cout << "Spectral norm of the dphi matrix: " << normPhi << std::endl;
//            std::cout << "Relative spectral approximation error of the dphi matrix: " << normPhiResiduum/normPhi << std::endl;
//            std::cout << "Spectral norm of the inverse dphi matrix: " << normPhiInverse /*<< ". Calculation may be unreliable for large low rank preconditioning matrices."*/ << std::endl;
//            std::cout << "Spectral condition number of the dphi matrix: " << normPhi * normPhiInverse<< std::endl;
//            std::cout << "Runtime of norm calculations: " << timer.secs() << std::endl;
//        }
//    }

//    dPhiHMat.clear(false); // delete matrix
//    LUPair.first.clear(false); // delete preconditioner
//    LUPair.second.clear(false);

//    /////////////////////////////////////////////////////////// reconstruct phiSolution from dPhiSolution via the boundary conditions

//    boundaryElements.phiSolution = boundaryElements.dPhiSolution;
//    for(long i = 0; i < substituteDPhiWithPhi.size(); i++)
//    {
//        if (substituteDPhiWithPhi(i))
//        {
//            boundaryElements.dPhiSolution(i) = (f(i) - alpha(i) * boundaryElements.phiSolution(i)) / beta(i);
//        }
//        else
//        {
//            boundaryElements.phiSolution(i) = (f(i) - beta(i) * boundaryElements.dPhiSolution(i)) / alpha(i);
//        }
//    }
//    testSol2 = boundaryElements.phiSolution; // for testing
//    clusterTree.clear();
//    global::trimMemory();
//}

void BoundaryElementSolver::hMatrixSolve()
{
    Timer timer;
    randomGenerator.seed(QRandomGenerator::system()->generate());
    std::cout<<"Number of elements in H-matrix solve: " << boundaryElements.triangles.length() << std::endl;

    ClusterTree clusterTree(&boundaryElements);

    prepareBoundaryElements();
    calculateBoundaryConditionAndSourceTermVector();

    HMatrix phiHMat;
    phiHMat.setClusterTrees(&clusterTree, &clusterTree);
    phiHMat.populateBlockClusterTree();

    timer.start();
    if(boundaryElements.impedancePlanes.length() == 0) // no impedance plane
    {
        phiHMat.assembleBlocks(acaMaxRank, acaRelativeError, std::bind(&BoundaryElementSolver::implicitPhiMatrix2, this, std::placeholders::_1, std::placeholders::_2), false);
    }
    else if(boundaryElements.impedancePlanes.length() == 1)
    {
        BoundaryElements reflectedBoundaryElements = boundaryElements;
        reflectedBoundaryElements.reflectGeometry(boundaryElements.impedancePlanes.at(0));
        phiHMat.assembleBlocks(acaMaxRank, acaRelativeError, std::bind(&BoundaryElementSolver::implicitPhiMatrixPlusReflected, this, std::placeholders::_1, std::placeholders::_2, std::ref(boundaryElements), std::ref(reflectedBoundaryElements)), false);
    }
    timer.stop();
    std::cout << "phiHMat assembled in " << timer.secs() << " seconds" << std::endl;
    //    std::cout << "frobenius norm of phi matrix: " << phiHMat.normInSVDForm() << std::endl;
//        std::cout << "frobenius norm of phi matrix: " << phiHMat.norm() << std::endl;

    double normPhiResiduum = 0;
    if(calculateNormAndCond)
    {
        timer.reset();
        timer.start();
        normPhiResiduum = HArithm::spectralNormLeastSignificant(phiHMat); // Approximates the spectral norm of the residuum of dPhiHMat
        timer.stop();
        std::cout << "Spectral norm of the phi residuum: " << normPhiResiduum << std::endl;
        std::cout << "Runtime of norm calculations: " << timer.secs() << std::endl;
    }


    timer.reset();

    timer.start();
    double compressionRatio = phiHMat.getCompressionRatio();
    timer.stop();
    std::cout << std::endl << "phi compressionRatio: " << compressionRatio <<std::endl << "calculated in " << timer.secs() << " seconds" << std::endl;
    timer.reset();
    std::cout << "maximum local rank of phi matrix: " << phiHMat.calculateMaxBlockRank() << std::endl;

    Eigen::VectorXcd rightHandside = sourceTermVector;
    HArithm::MVM(rightHandside, phiHMat, f);
    if(testACA) // just for testing
    {
        testPhiHMat = phiHMat.toFullMat();
    }
    phiHMat.clear(true); // phimat in not needed anymore

    Eigen::VectorXcd rightHandsideCopy = rightHandside;
    HMatRightHandSide = rightHandside;

    HMatrix dPhiHMat;
    dPhiHMat.setClusterTrees(&clusterTree, &clusterTree);
    dPhiHMat.populateBlockClusterTree();

    timer.start();
    if(boundaryElements.impedancePlanes.length() == 0) // no impedance plane
    {
//         std::bind(&BoundaryElementSolver::getNormalFilteredPivotIndices, this, std::placeholders::_1, std::placeholders::_2), std::bind(&BoundaryElementSolver::implicitDPhiMatrix, this, std::placeholders::_1, std::placeholders::_2)
        dPhiHMat.assembleBlocksExtra(acaMaxRank, acaRelativeError, std::bind(&BoundaryElementSolver::getNormalFilteredPivotIndices, this, std::placeholders::_1, std::placeholders::_2), std::bind(&BoundaryElementSolver::implicitDPhiMatrix2, this, std::placeholders::_1, std::placeholders::_2), false);
    }
    else if(boundaryElements.impedancePlanes.length() == 1)
    {
        BoundaryElements reflectedBoundaryElements = boundaryElements;
        reflectedBoundaryElements.reflectGeometry(boundaryElements.impedancePlanes.at(0));
        dPhiHMat.assembleBlocksExtra(acaMaxRank, acaRelativeError, std::bind(&BoundaryElementSolver::getNormalFilteredPivotIndices, this, std::placeholders::_1, std::placeholders::_2), std::bind(&BoundaryElementSolver::implicitDPhiMatrixPlusReflected, this, std::placeholders::_1, std::placeholders::_2, std::ref(boundaryElements), std::ref(reflectedBoundaryElements)), false);
    }
    timer.stop();
    std::cout << "dPhiHMat assembled in " << timer.secs() << " seconds" << std::endl;
    timer.reset();

    if(testACA) // just for testing
    {
        testDPhiHMat = dPhiHMat.toFullMat();
    }

    double normDPhiResiduum = 0;
    if(calculateNormAndCond)
    {
        timer.reset();
        timer.start();
        normDPhiResiduum = HArithm::spectralNormLeastSignificant(dPhiHMat); // Approximates the spectral norm of the residuum of dPhiHMat
        timer.stop();
        std::cout << "Spectral norm of the dphi residuum: " << normDPhiResiduum << std::endl;
        std::cout << "Runtime of norm calculations: " << timer.secs() << std::endl;
    }

    timer.start();
    compressionRatio = dPhiHMat.getCompressionRatio();
    timer.stop();
    std::cout << std::endl << "dphi compressionRatio: " << compressionRatio <<std::endl << "calculated in " << timer.secs() << " seconds" << std::endl;
    timer.reset();
    std::cout << "maximum local rank of dphi matrix: " << dPhiHMat.calculateMaxBlockRank() << std::endl;

    ////////////////////////////////////// Copy dPhiHMat for later GMRES usage, as LUSolver modifies the original dPhiHMat

    HMatrix dPhiHMatCopy;
    std::pair<HMatrix,HMatrix> LUPair;
    if(usePreconditioner)
    {
        timer.start();
        dPhiHMatCopy = HMatrix(dPhiHMat, preconditionerRank, preconditionerRelError, false);
        timer.stop();
        std::cout << std::endl << "dPhiHMatCopy(dPhiHMat) took " << timer.secs() << " seconds." << std::endl;
        timer.reset();
        timer.start();
        compressionRatio = dPhiHMatCopy.getCompressionRatio();
        timer.stop();
        std::cout << std::endl << "dPhiHMatCopy compressionRatio: " << compressionRatio <<std::endl << "calculated in " << timer.secs() << " seconds" << std::endl;
        timer.reset();

        timer.start();
        std::cout << "Calculating the HLU preconditioner." << std::endl;
        bool zeroBlocksUninitalized = true; // keep the zero blocks in the triangular block matrices (L and U) uninitialized
        LUPair = HArithm::LUDecomposition(dPhiHMatCopy, preconditionerRank, preconditionerRelError, zeroBlocksUninitalized); // calculate the LU decomposition of dPhiHMatCopy
        dPhiHMatCopy.clear(true); // delete low rank copy
        timer.stop();
        std::cout << "Runtime of LUDecomposition: " << timer.secs() << std::endl;
        timer.reset();
    }

    ////////////////////////////////////// Testing section start ////////////////////////////////
    if(testACA)
    {
        phiHMat.clear(true);
        dPhiHMat.clear(true);
        return;
    }
    ////////////////////////////////////// Testing section end ////////////////////////////////

    if(LUPair.first.hasNan() || LUPair.second.hasNan() || usePreconditioner == false) // LU decomposition has NaN; -> don't use the preconditioner
    {
        std::cerr << "GMRES start" << std::endl;
        if(usePreconditioner == true)
        {
            std::cerr << " There are NaNs in the LU decomposition. Continuing without preconditioner." << std::endl;
        }
        timer.start();
        HMatrixWrapper dPhiHMatForGMRES(&dPhiHMat); // set up a wrapper class for the hmatrix to be used in the matrix free gmres routine
        boundaryElements.dPhiSolution = GMRES::gmresSolve(dPhiHMatForGMRES, Eigen::VectorXcd::Random(boundaryElements.triangles.length()), rightHandsideCopy, 0.001*acaRelativeError); // solve dPhiHMat * x = righthandside
        timer.stop();
        std::cout << "Runtime of GMRES: " << timer.secs() << std::endl;
        timer.reset();
        if(calculateNormAndCond)
        {
            timer.start();
            double normPhi = HArithm::spectralNorm(dPhiHMat); // Approximate the spectral norm of dPhiHMat via power iteration
            timer.stop();
            std::cout << "Spectral norm of the dphi matrix: " << normPhi << std::endl;
            std::cout << "Relative spectral approximation error of the dphi matrix: " << normPhiResiduum/normPhi << std::endl;
            std::cout << "Runtime of norm calculation: " << timer.secs() << std::endl;
        }
        Eigen::VectorXcd b = Eigen::VectorXcd::Zero(boundaryElements.dPhiSolution.size());
        HArithm::MVM(b, dPhiHMat, boundaryElements.dPhiSolution);
        std::cout << "Relative error of the residuum: " << (b - rightHandsideCopy).norm()/ rightHandsideCopy.norm() << std::endl;
    }
    else
    {
        Eigen::VectorXcd gmresGuess = HArithm::LUSubstitutionSolve(LUPair.first, LUPair.second, rightHandsideCopy);
        timer.start();
        LUPrecondidionedHMatrixWrapper dPhiHMatForLUPreconditionedGMRES(&dPhiHMat, &LUPair.first, &LUPair.second);
        boundaryElements.dPhiSolution = GMRES::gmresLUPreconditionedSolve(dPhiHMatForLUPreconditionedGMRES, gmresGuess, rightHandsideCopy, 0.001*acaRelativeError);
        timer.stop();
        std::cout << "Runtime of preconditioned GMRES: " << timer.secs() << std::endl;
        timer.reset();

        Eigen::VectorXcd b = Eigen::VectorXcd::Zero(boundaryElements.dPhiSolution.size());
        HArithm::MVM(b, dPhiHMat, boundaryElements.dPhiSolution);
        std::cout << "Relative error of the residuum: " << (b - rightHandsideCopy).norm()/ rightHandsideCopy.norm() << std::endl;
        if(calculateNormAndCond)
        {
            timer.start();
            double normPhi = HArithm::spectralNorm(dPhiHMat); // Approximates the spectral norm of dPhiHMat via power iteration
            double normPhiInverse = HArithm::spectralNormFromLU(LUPair.first, LUPair.second); // Approximates the spectral norm of the inverse of of dPhiHMat via power iteration with its LU decomposition.
            timer.stop();
            std::cout << "Spectral norm of the dphi matrix: " << normPhi << std::endl;
            std::cout << "Relative spectral approximation error of the dphi matrix: " << normPhiResiduum/normPhi << std::endl;
            std::cout << "Spectral norm of the inverse dphi matrix: " << normPhiInverse /*<< ". Calculation may be unreliable for large low rank preconditioning matrices."*/ << std::endl;
            std::cout << "Spectral condition number of the dphi matrix: " << normPhi * normPhiInverse<< std::endl;
            std::cout << "Runtime of norm calculations: " << timer.secs() << std::endl;

//            timer.reset();
//            timer.start();
//            normPhi = dPhiHMat.norm();
//            normPhiInverse = HArithm::frobeniusNormFromLU(LUPair.first, LUPair.second, 5); // Approximates the frobenius norm of the inverse of of dPhiHMat.
//            timer.stop();
//            std::cout << "Frobenius norm of the dphi matrix: " << normPhi << std::endl;
////            std::cout << "Relative spectral approximation error of the dphi matrix: " << normPhiResiduum/normPhi << std::endl;
//            std::cout << "Frobenius norm of the inverse dphi matrix: " << normPhiInverse /*<< ". Calculation may be unreliable for large low rank preconditioning matrices."*/ << std::endl;
//            std::cout << "Frobenius condition number of the dphi matrix: " << normPhi * normPhiInverse<< std::endl;
//            std::cout << "Runtime of norm calculations: " << timer.secs() << std::endl;
        }
    }

    dPhiHMat.clear(); // delete matrix
    LUPair.first.clear(); // delete preconditioner
    LUPair.second.clear();

    /////////////////////////////////////////////////////////// reconstruct phiSolution from dPhiSolution via the boundary conditions

    boundaryElements.phiSolution = boundaryElements.dPhiSolution;
    for(long i = 0; i < substituteDPhiWithPhi.size(); i++)
    {
        if (substituteDPhiWithPhi(i))
        {
            boundaryElements.dPhiSolution(i) = (f(i) - alpha(i) * boundaryElements.phiSolution(i)) / beta(i);
        }
        else
        {
            boundaryElements.phiSolution(i) = (f(i) - beta(i) * boundaryElements.dPhiSolution(i)) / alpha(i);
        }
    }
    testSol2 = boundaryElements.phiSolution; // for testing
    clusterTree.clear();
    global::trimMemory();
}

void BoundaryElementSolver::calcReflectionMatrices(HMatrix &reflMatPhi, HMatrix &reflMatDPhi, const long maxRank, const double relativeError, const BoundaryElements &elements, BoundaryElements reflElements, ClusterTree &clusterTree, int lastPlaneIndex)
{
    for(int planeIndex=0; planeIndex<elements.impedancePlanes.length() && planeIndex < 2; planeIndex++)
    {
        if(lastPlaneIndex == planeIndex)
        {
            continue;
        }
        BoundaryElements reflectedBoundaryElements = reflElements;
        reflectedBoundaryElements.reflectGeometry(elements.impedancePlanes.at(planeIndex));

        ClusterTree reflectedClusterTree(clusterTree);  //copy the original clustertree
        reflectedClusterTree.updateMinCuboids(&reflectedBoundaryElements);

        HMatrix newReflMatPhi;
        newReflMatPhi.setClusterTrees(&clusterTree, &reflectedClusterTree);
        newReflMatPhi.populateBlockClusterTree();

        HMatrix newReflMatDPhi;
        newReflMatDPhi.setClusterTrees(&clusterTree, &reflectedClusterTree);
        newReflMatDPhi.populateBlockClusterTree();

        QVector<BlockCluster*> newReflMatPhiPartition = newReflMatPhi.getMinPartition();
        QVector<BlockCluster*> newReflMatDPhiPartition = newReflMatDPhi.getMinPartition();

        Timer timer;
        timer.start();

        #pragma omp parallel master
        for(long blockIndex = 0; blockIndex < newReflMatPhiPartition.size(); blockIndex++)
        {
            #pragma omp task
            hBlockAssemblyReflMat(newReflMatPhiPartition.at(blockIndex), newReflMatDPhiPartition.at(blockIndex), 0, relativeError, elements, reflectedBoundaryElements);
        }
        timer.stop();
        std::cout << "Runtime of hBlockAssemblyReflectedMatrix: " << timer.secs() << std::endl;
        timer.reset();

        omp_set_max_active_levels(2); // raise the number of active parallelization levels for the current section
        #pragma omp parallel sections
        {
            #pragma omp section
            {
                HArithm::compressHMat(newReflMatPhi, maxRank, 0.1 * relativeError);
            }
            #pragma omp section
            {
                HArithm::compressHMat(newReflMatDPhi, maxRank, 0.1 * relativeError);
            }
        }
        omp_set_max_active_levels(1);


        double newPhiNorm = newReflMatPhi.normInSVDForm();
        double newDPhiNorm = newReflMatDPhi.normInSVDForm();

//        double newPhiNorm = newReflMatPhi.norm();
//        double newDPhiNorm = newReflMatDPhi.norm();
        std::cout << "newPhiNorm: " << newPhiNorm << std::endl;

        if(lastPlaneIndex == -1 && planeIndex == 0)
        {
            reflMatPhi = newReflMatPhi;
            reflMatDPhi = newReflMatDPhi;

            calcReflectionMatrices(reflMatPhi, reflMatDPhi, maxRank, relativeError, elements, reflectedBoundaryElements, clusterTree, planeIndex);
//            reflectedClusterTree.clear();  // can't be deleted here -> still a memory leak
            continue;
        }
        omp_set_max_active_levels(2); // raise the number of active parallelization levels for the current section
        #pragma omp parallel sections
        {
            #pragma omp section
            HArithm::recursiveHMatAddition(* reflMatPhi.getRootBlock(), * newReflMatPhi.getRootBlock(), maxRank, 0.1 * relativeError); // for reflective planes
            #pragma omp section
            HArithm::recursiveHMatAddition(* reflMatDPhi.getRootBlock(), * newReflMatDPhi.getRootBlock(), maxRank, 0.1 * relativeError);
        }
        omp_set_max_active_levels(1);

//        HArithm::recursiveHMatSubstraction(reflMatPhi.getRootBlock(), newReflMatPhi.getRootBlock(), maxRank, 0.1 * relativeError); // for absorptive planes
//        HArithm::recursiveHMatSubstraction(reflMatDPhi.getRootBlock(), newReflMatDPhi.getRootBlock(), maxRank, 0.1 * relativeError);

        newReflMatDPhi.clear(); // newReflMatDPhi in not needed anymore
        newReflMatPhi.clear(); // newReflMatPhi in not needed anymore
        reflectedClusterTree.clear();

        double reflMatPhiNorm = reflMatPhi.norm();
        double reflMatDPhiNorm = reflMatDPhi.norm();
        std::cout << "reflMatPhiNorm: " << reflMatPhiNorm << std::endl;

        std::cout << "newPhiNorm / reflMatPhiNorm: " << newPhiNorm / reflMatPhiNorm << std::endl;
        std::cout << "newDPhiNorm / reflMatDPhiNorm: " << newDPhiNorm / reflMatDPhiNorm << std::endl;

        if(newPhiNorm / reflMatPhiNorm < relativeError && newDPhiNorm / reflMatDPhiNorm < relativeError)
        {
            continue;
        }
        else
        {
            calcReflectionMatrices(reflMatPhi, reflMatDPhi, maxRank, relativeError, elements, reflectedBoundaryElements, clusterTree, planeIndex);
        }
    }
    if(lastPlaneIndex == -1)
    {
        global::trimMemory();
    }
}

void BoundaryElementSolver::hBlockAssembly(BlockCluster* phiBlock, BlockCluster* dPhiBlock, const long maxRank, const double relativeError)
{
    if(phiBlock -> isAdmissible) // Adaptive Cross Approximation for low rank block assembly
    {
            #pragma omp task
            partialPivotACA(phiBlock, maxRank, relativeError, std::bind(&BoundaryElementSolver::implicitPhiMatrix, this, std::placeholders::_1, std::placeholders::_2));
            #pragma omp task
            partialPivotACAextra(dPhiBlock, maxRank, relativeError, std::bind(&BoundaryElementSolver::getNormalFilteredPivotIndices, this, std::placeholders::_1, std::placeholders::_2), std::bind(&BoundaryElementSolver::implicitDPhiMatrix, this, std::placeholders::_1, std::placeholders::_2));
    }
    else // nearfield block full matrix assembly
    {
        long rowStartIndex = phiBlock->rowStartIndex();
        long columnStartIndex = phiBlock->colStartIndex();
        long blockRows = phiBlock->rows(); //assumes contiguous ascending indexes
        long blockColumns = phiBlock->cols();

        phiBlock->fullMat = Eigen::MatrixXcd(blockRows, blockColumns);
        dPhiBlock->fullMat = Eigen::MatrixXcd(blockRows, blockColumns);

        std::complex<double> Mk;
        std::complex<double> Nk;
        std::complex<double> Lk;
        std::complex<double> Mtk;
        bool onPanel = false;

//        #pragma omp parallel for private(Mk,Nk,Lk,Mtk,onPanel)// parallelizes matrix initialization
        for(long row = 0; row < blockRows; row++)
        {
            for(long column = 0; column < blockColumns; column++)
            {
                long globalRowIndex = rowStartIndex + row;
                long globalColumnIndex = columnStartIndex + column;
                onPanel = (globalRowIndex == globalColumnIndex);

                BemOperatorsConst(rowStartIndex + row, columnStartIndex + column, Lk, Mk, Mtk, Nk);

                if(substituteDPhiWithPhi(globalColumnIndex))
                {
                    phiBlock->fullMat(row,column) = - (Lk + couplingParameter * Mtk + ((double)onPanel) * couplingParameterHalf) / beta(globalColumnIndex);
                    dPhiBlock->fullMat(row,column) = -alpha(globalColumnIndex) / beta(globalColumnIndex) * (Lk + couplingParameter * Mtk + ((double)onPanel) * couplingParameterHalf) - (Mk + couplingParameter * Nk - ((double)onPanel) * 0.5);
                }
                else
                {
                    phiBlock->fullMat(row,column) = (Mk + couplingParameter * Nk - ((double)onPanel) * 0.5) / alpha(globalColumnIndex);
                    dPhiBlock->fullMat(row,column) = Lk + couplingParameter * Mtk + ((double)onPanel) * couplingParameterHalf + beta(globalColumnIndex)/alpha(globalColumnIndex) * (Mk + couplingParameter * Nk - ((double)onPanel) * 0.5);
                }
            }
        }
    }
}

void BoundaryElementSolver::hBlockAssemblyReflMat(BlockCluster* phiBlock, BlockCluster* dPhiBlock, const long maxRank, const double relativeError, const BoundaryElements &boundaryElements, const BoundaryElements &reflectedElements)
{
    if(phiBlock -> isAdmissible) // Adaptive Cross Approximation
    {
        #pragma omp task
        partialPivotACA(phiBlock, maxRank, relativeError, std::bind(&BoundaryElementSolver::implicitPhiMatrixReflected, this, std::placeholders::_1, std::placeholders::_2, std::ref(boundaryElements), std::ref(reflectedElements)));
        #pragma omp task
        partialPivotACA(dPhiBlock, maxRank, relativeError, std::bind(&BoundaryElementSolver::implicitDPhiMatrixReflected, this, std::placeholders::_1, std::placeholders::_2, std::ref(boundaryElements), std::ref(reflectedElements)));
    }
    else // nearfield block full matrix assembly
    {
        long rowStartIndex = phiBlock->rowStartIndex();
        long columnStartIndex = phiBlock->colStartIndex();
        long blockRows = phiBlock->rows(); //assumes contiguous ascending indexes
        long blockColumns = phiBlock->cols();

        phiBlock->fullMat = Eigen::MatrixXcd(blockRows, blockColumns);
        dPhiBlock->fullMat = Eigen::MatrixXcd(blockRows, blockColumns);

        std::complex<double> Mk;
        std::complex<double> Nk;
        std::complex<double> Lk;
        std::complex<double> Mtk;

//        #pragma omp parallel for private(Mk,Nk,Lk,Mtk)// parallelizes matrix initialization
        for(long row = 0; row < blockRows; row++)
        {
            for(long column = 0; column < blockColumns; column++)
            {
                long globalColumnIndex = columnStartIndex + column;
                BemOperatorsReflected(rowStartIndex + row, columnStartIndex + column, Lk, Mk, Mtk, Nk, boundaryElements, reflectedElements);

                if(substituteDPhiWithPhi(globalColumnIndex))
                {
                    phiBlock->fullMat(row,column) = - (Lk + couplingParameter * Mtk) / beta(globalColumnIndex);
                    dPhiBlock->fullMat(row,column) = -alpha(globalColumnIndex) / beta(globalColumnIndex) * (Lk + couplingParameter * Mtk) - (Mk + couplingParameter * Nk);
                }
                else
                {
                    phiBlock->fullMat(row,column) = (Mk + couplingParameter * Nk) / alpha(globalColumnIndex);
                    dPhiBlock->fullMat(row,column) = Lk + couplingParameter * Mtk + beta(globalColumnIndex)/alpha(globalColumnIndex) * (Mk + couplingParameter * Nk);
                }
            }
        }
    }
}

void BoundaryElementSolver::fullPivotACA(BlockCluster* block, const long rank, double relativeError, std::function<std::complex<double> (long, long)> implicitMatrix)
{
    block->VAdjMat.resize(0,0);
    block->singularValues.resize(0);
    block->UMat.resize(0,0);

    double norm = 0; // full rank matrix norm approximation by the R1 matrix norm multiplied by relative error
    double residuumNorm = 0;
    long rowStartIndex = block->rowStartIndex();
    long columnStartIndex = block->colStartIndex();
    long blockRows = block->rows(); //assumes contiguous ascending indexes
    long blockColumns = block->cols();

    long maxRank = std::min(blockRows, blockColumns);
    long initialRankReservation = 5;
    if(rank > 0)
    {
        maxRank = std::min(maxRank, rank);
        block->UMat.resize(blockRows, maxRank);
        block->singularValues.resize(maxRank);
        block->VAdjMat.resize(maxRank, blockColumns);
    }
    else
    {
        long reservationRank = std::min(initialRankReservation, maxRank);
        block->UMat.resize(blockRows, reservationRank);
        block->singularValues.resize(reservationRank);
        block->VAdjMat.resize(reservationRank, blockColumns);
    }
    long tmpRank = 0;
    Eigen::RowVectorXcd rowVector(blockColumns);
    Eigen::MatrixXcd Residuum(blockRows, blockColumns);
    for (int i = 0; i < blockRows ; i++)
    {
        calcHBlockRowVector(rowVector, rowStartIndex + i, columnStartIndex, implicitMatrix);
        Residuum.row(i) = rowVector;
    }
    norm = Residuum.norm();
    long rowIndex;
    long columnIndex;

    Eigen::VectorXcd columnVector(blockRows);
    while(tmpRank < maxRank)
    {
        Residuum.cwiseAbs().maxCoeff(&rowIndex, &columnIndex);

        calcHBlockRowVector(rowVector, rowStartIndex + rowIndex, columnStartIndex, implicitMatrix);

        if(tmpRank >= 1)
        {
            rowVector -= (block->UMat.leftCols(tmpRank).row(rowIndex) * block->singularValues.head(tmpRank).asDiagonal()) * block->VAdjMat.topRows(tmpRank);
        }

        calcHBlockColumnVector(columnVector, rowStartIndex, columnStartIndex + columnIndex,  implicitMatrix);

        std::complex<double> alpha = 1.0 / rowVector(columnIndex);

        if(tmpRank >= 1)
        {
            columnVector -= block->UMat.leftCols(tmpRank) * (block->singularValues.head(tmpRank).asDiagonal() * block->VAdjMat.topRows(tmpRank).col(columnIndex));
        }

        if(block->singularValues.size() <= tmpRank)
        {
            long newSize = std::min(2*tmpRank, maxRank);
            block->UMat.conservativeResize(Eigen::NoChange, newSize);
            block->singularValues.conservativeResize(newSize);
            block->VAdjMat.conservativeResize(newSize, Eigen::NoChange);
         }

        block->UMat.col(tmpRank) = columnVector;
        block->singularValues(tmpRank) = alpha;
        block->VAdjMat.row(tmpRank) = rowVector;
        Residuum -= block->UMat.col(tmpRank) * block->singularValues(tmpRank) * block->VAdjMat.row(tmpRank);
        tmpRank++; // increment rank counter to actual current rank of the low rank matrix


        residuumNorm = Residuum.norm();

        if(relativeError > 0 && residuumNorm / norm < relativeError)
        {
            break;
        }
    }
    block->UMat.conservativeResize(Eigen::NoChange, tmpRank);
    block->singularValues.conservativeResize(tmpRank);
    block->VAdjMat.conservativeResize(tmpRank, Eigen::NoChange);

    if(tmpRank == 0) // if nothing has been sampled
    {
         block->UMat = Eigen::MatrixXcd::Zero(blockRows, 1);
         block->singularValues = Eigen::VectorXcd::Zero(1);
         block->VAdjMat = Eigen::MatrixXcd::Zero(1, blockColumns);
    }
}

void BoundaryElementSolver::partialPivotACA(BlockCluster* block, const long rank, double relativeError, std::function<std::complex<double> (long, long)> implicitMatrix)
{
    block->VAdjMat.resize(0,0);
    block->singularValues.resize(0);
    block->UMat.resize(0,0);

    std::complex<double> normEstimate = 0; // full rank matrix norm approximation by the R1 matrix norm multiplied by relative error
    double relErrorTarget = 0; // full rank matrix norm approximation by the R1 matrix norm multiplied by relative error
    long rowStartIndex = block->rowStartIndex();
    long columnStartIndex = block->colStartIndex();
    long blockRows = block->rows(); //assumes contiguous ascending indexes
    long blockColumns = block->cols();

    long maxRank = std::min(blockRows, blockColumns);
    long initialRankReservation = 5;
    if(rank > 0)
    {
        maxRank = std::min(maxRank, rank);
        block->UMat.resize(blockRows, maxRank);
        block->singularValues.resize(maxRank);
        block->VAdjMat.resize(maxRank, blockColumns);
    }
    else
    {
        long reservationRank = std::min(initialRankReservation, maxRank);
        block->UMat.resize(blockRows, reservationRank);
        block->singularValues.resize(reservationRank);
        block->VAdjMat.resize(reservationRank, blockColumns);
    }
    long tmpRank = 0;

    Eigen::RowVectorXcd rowVector(blockColumns);
    Eigen::VectorXcd columnVector(blockRows);

    long rowIndex;
    long columnIndex;

    columnIndex = QRandomGenerator::system()->bounded((qint32) blockColumns); // find random row index, indices can not be reused
    calcHBlockColumnVector(columnVector, rowStartIndex, columnStartIndex + columnIndex,  implicitMatrix);

    while(tmpRank < maxRank /*&& error condition*/)
    {
        if(tmpRank >= 1)
        {
            columnVector[rowIndex] = 0;
        }
        else //generate new random row index
        {
//            rowIndex = QRandomGenerator::global()->bounded((qint32) blockRows); // find random row index, indices can not be reused
            rowIndex = QRandomGenerator::system()->bounded((qint32) blockRows); // find random row index, indices can not be reused
        }
        columnVector.cwiseAbs().maxCoeff(&rowIndex);

        if(std::abs(columnVector(rowIndex)) <= global::tiny && tmpRank >= 1) //the max absolute element is too small -> no suitable column index for rowIndex
        {
            std::cerr << "max(columnVector) <= global::Tiny" << std::endl;
            break;    //change later
        }

        // calculate the block row rowIndex
        calcHBlockRowVector(rowVector, rowStartIndex + rowIndex, columnStartIndex, implicitMatrix);

        if(tmpRank >= 1)
        {
            rowVector -= (block->UMat.leftCols(tmpRank).row(rowIndex) * block->singularValues.head(tmpRank).asDiagonal()) * block->VAdjMat.topRows(tmpRank);
        }

         // find column index for largest absolute element in row with above rowindex
        double maxRowVal = rowVector.cwiseAbs().maxCoeff(&columnIndex); //  columnIndex is now max column index of phiAbsVec

        if(maxRowVal <= global::tiny && tmpRank >= 1) //the max absolute element is too small -> no suitable column index for rowIndex
        {
            std::cerr << "max(rowAbsVec) <= Tiny" << std::endl;
            break;    //change later
        }

        // calculate the block column columnIndex
        calcHBlockColumnVector(columnVector, rowStartIndex, columnStartIndex + columnIndex, implicitMatrix);

        std::complex<double> alpha = 1.0 / rowVector(columnIndex);

        if(tmpRank >= 1)
        {
            columnVector -= block->UMat.leftCols(tmpRank) * (block->singularValues.head(tmpRank).asDiagonal() * block->VAdjMat.topRows(tmpRank).col(columnIndex));
        }

        if(block->singularValues.size() <= tmpRank)
        {
            long newSize = std::min(2*tmpRank, maxRank);
            block->UMat.conservativeResize(Eigen::NoChange, newSize);
            block->singularValues.conservativeResize(newSize);
            block->VAdjMat.conservativeResize(newSize, Eigen::NoChange);
         }

        block->UMat.col(tmpRank) = columnVector;
        block->VAdjMat.row(tmpRank) = rowVector;
        block->singularValues(tmpRank) = alpha;
        tmpRank++; // increment rank counter to actual current rank of the low rank matrix

        normEstimate += std::abs(std::pow(alpha, 2)) *  rowVector.squaredNorm() * columnVector.squaredNorm();
        for(long i = 0; i < tmpRank - 1; i++)
        {
            normEstimate += 2.0 * std::conj(alpha) * block->singularValues(i) * (columnVector.adjoint() * block->UMat.col(i))(0,0) * (( block->VAdjMat.row(i)) * ( rowVector).adjoint())(0,0);
        }
        relErrorTarget = std::real(std::pow(relativeError, 2) * normEstimate);
        if(relativeError > 0 && std::pow(std::abs(alpha), 2) * rowVector.squaredNorm() * columnVector.squaredNorm() < relErrorTarget)
        {
            break;
        }
    }
    block->UMat.conservativeResize(Eigen::NoChange, tmpRank);
    block->singularValues.conservativeResize(tmpRank);
    block->VAdjMat.conservativeResize(tmpRank, Eigen::NoChange);

    block->frobeniusNorm = std::sqrt(std::abs(normEstimate));

    if(tmpRank == 0) // if the first pivoting element was already do small -> nothing has been sampled
    {
         block->UMat = Eigen::MatrixXcd::Zero(blockRows, 1);
         block->singularValues = Eigen::VectorXcd::Zero(1);
         block->VAdjMat = Eigen::MatrixXcd::Zero(1, blockColumns);
    }
}

void BoundaryElementSolver::partialPivotACAextra(BlockCluster* block, long rank, double relativeError, std::function<QVector<std::pair<long,long>>(BlockCluster*,double)> getPivotIndices, std::function<std::complex<double> (long, long)> implicitMatrix) /*!< Low-rank assembly of an h-block by ACA with heuristic partial pivoting. The (guaranteed) accuracy is improved for dPhi-blocks. */
{
    block->VAdjMat.resize(0,0);
    block->singularValues.resize(0);
    block->UMat.resize(0,0);

    std::complex<double> normEstimate = 0; // full rank matrix norm approximation by the R1 matrix norm multiplied by relative error
    double relErrorTarget = 0; // full rank matrix norm approximation by the R1 matrix norm multiplied by relative error
    long rowStartIndex = block->rowStartIndex();
    long columnStartIndex = block->colStartIndex();
    long blockRows = block->rows(); //assumes contiguous ascending indexes
    long blockColumns = block->cols();
    long maxRank = std::min(blockRows, blockColumns);
    long initialRankReservation = 5;
    if(rank > 0)
    {
        maxRank = std::min(maxRank, rank);
        block->UMat.resize(blockRows, maxRank);
        block->singularValues.resize(maxRank);
        block->VAdjMat.resize(maxRank, blockColumns);
    }
    else
    {
        long reservationRank = std::min(initialRankReservation, maxRank);
        block->UMat.resize(blockRows, reservationRank);
        block->singularValues.resize(reservationRank);
        block->VAdjMat.resize(reservationRank, blockColumns);
    }

    QVector<std::pair<long,long>> pivotIndices = getPivotIndices(block, 0.44);

    Eigen::RowVectorXcd rowVector(blockColumns);
    Eigen::VectorXcd columnVector(blockRows);
    long tmpRank = 0;

    for(long i = 0; i < pivotIndices.size() && tmpRank < maxRank; i++)
    {

        long rowIndex = pivotIndices.at(i).first;
        long columnIndex = pivotIndices.at(i).second;

        std::complex<double> pivotElement = implicitMatrix(rowStartIndex + rowIndex, columnStartIndex + columnIndex);

        if(std::abs(pivotElement) < global::tiny)
        {
            continue; // pivot element is too small
        }
        if(std::abs((pivotElement - ((block->UMat.row(rowIndex).head(tmpRank).transpose().array() * block->singularValues.head(tmpRank).array()).array() * block->VAdjMat.col(columnIndex).head(tmpRank).array()).sum()) / pivotElement) < relativeError)
        {
            continue; // pivot element is already sufficiently approximated
        }
        calcHBlockColumnVector(columnVector, rowStartIndex, columnStartIndex + columnIndex, implicitMatrix);
        if(tmpRank >= 1)
        {
            columnVector -= block->UMat.leftCols(tmpRank) * (block->singularValues.head(tmpRank).asDiagonal() * block->VAdjMat.topRows(tmpRank).col(columnIndex));
        }
        bool firstIterationForCurrentPivotElement = true;
        while(tmpRank < maxRank)
        {
//            std::cerr<<"tmpRank " <<tmpRank <<std::endl;

            Eigen::VectorXd colAbsVec;
            colAbsVec  = columnVector.cwiseAbs();
            if(!firstIterationForCurrentPivotElement)
            {
                colAbsVec[rowIndex] = 0; //at this point rowIndex is still the number from the last iteration
            }
            colAbsVec.maxCoeff(&rowIndex); //  columnIndex is now max column index of phiAbsVec
            firstIterationForCurrentPivotElement = false;
            if(std::abs(columnVector(rowIndex)) <= global::tiny && tmpRank >= 1) //the max absolute element is too small -> no suitable column index for rowIndex
            {
                 std::cerr << "max(columnVector) <= global::Tiny" << std::endl;
                 break;    //change later
            }

            // calculate the block row
            calcHBlockRowVector(rowVector, rowStartIndex + rowIndex, columnStartIndex, implicitMatrix);

            if(tmpRank >= 1)
            {
                rowVector -= (block->UMat.leftCols(tmpRank).row(rowIndex) * block->singularValues.head(tmpRank).asDiagonal()) * block->VAdjMat.topRows(tmpRank);
            }

             // find column index for largest absolute element in row with above rowindex

            Eigen::RowVectorXd rowAbsVec = rowVector.cwiseAbs();

             double maxRowVal = rowAbsVec.maxCoeff(&columnIndex); //  columnIndex is now max column index of phiAbsVec

            if(maxRowVal <= global::tiny && tmpRank >= 1) //the max absolute element is too small -> no suitable column index for rowIndex
            {
                std::cerr << "max(rowAbsVec) <= Tiny" << std::endl;
            //             std::cerr << "maxRowVal: " << maxRowVal << std::endl;
                break;    //change later
            }

            // calculate the block column
            calcHBlockColumnVector(columnVector, rowStartIndex, columnStartIndex + columnIndex, implicitMatrix);

            std::complex<double> alpha = 1.0 / rowVector(columnIndex);

            if(tmpRank >= 1)
            {
                columnVector -= block->UMat.leftCols(tmpRank) * (block->singularValues.head(tmpRank).asDiagonal() * block->VAdjMat.topRows(tmpRank).col(columnIndex));
            }

            if(block->singularValues.size() <= tmpRank)
            {
                long newSize = std::min(2*tmpRank, maxRank);
                block->UMat.conservativeResize(Eigen::NoChange, newSize);
                block->singularValues.conservativeResize(newSize);
                block->VAdjMat.conservativeResize(newSize, Eigen::NoChange);
             }

            block->UMat.col(tmpRank) = columnVector;
            block->VAdjMat.row(tmpRank) = rowVector;
            block->singularValues(tmpRank) = alpha;
            tmpRank++; // increment rank counter to actual current rank of the low rank matrix

            std::complex<double> oldNorm = normEstimate;

            normEstimate += std::abs(std::pow(alpha, 2)) * rowVector.squaredNorm() * columnVector.squaredNorm();
            for(long i = 0; i < tmpRank - 1; i++)
            {
                normEstimate += 2.0 * std::conj(alpha) * block->singularValues(i) * (columnVector.adjoint() * block->UMat.col(i))(0,0) * (( block->VAdjMat.row(i)) * ( rowVector).adjoint())(0,0);
            }
            relErrorTarget = std::real(std::pow(relativeError, 2) * normEstimate);
            double normOfLastGain2 = std::pow(std::abs(alpha), 2) * rowVector.squaredNorm() * columnVector.squaredNorm();

            if(relativeError > 0 && normOfLastGain2 < relErrorTarget)
            {
                if(relativeError > 0 && normOfLastGain2 < 0.1 * relErrorTarget)
                {
                    tmpRank--;
                    normEstimate = oldNorm;
                }
                break;
            }
        }
    }
    block->UMat.conservativeResize(Eigen::NoChange, tmpRank);
    block->singularValues.conservativeResize(tmpRank);
    block->VAdjMat.conservativeResize(tmpRank, Eigen::NoChange);

    block->frobeniusNorm = std::sqrt(std::abs(normEstimate));

    if(tmpRank == 0)
    {
         block->UMat = Eigen::MatrixXcd::Zero(blockRows, 1);
         block->singularValues = Eigen::VectorXcd::Zero(1);
         block->VAdjMat = Eigen::MatrixXcd::Zero(1, blockColumns);
    }
}

void BoundaryElementSolver::BemOperatorsConst(const long row, const long column, std::complex<double>& Lk, std::complex<double>& Mk, std::complex<double>& Mtk, std::complex<double>& Nk)
{
    if(row == column) //observation Point p is on the same Panel as source points(q)
    {
        BemOperatorsSingularPolarInt(row,  Lk,  Mk,  Mtk,  Nk); // operators for singular panels
        return;
    }
    else if((boundaryElements.collocationPoints.at(row) - boundaryElements.collocationPoints.at(column)).norm() < 4 * boundaryElements.averageTriangleDim) // observation point lies close to singularity
    {
        BemOperatorsNearSing(row, column,  Lk,  Mk,  Mtk,  Nk); // operators for almost singular panels
        return;
    }
    else //observation Point p is on different Panel than the source points(q)
    {
        Lk = 0.0;
        Mk = 0.0;
        Mtk = 0.0;
        Nk = 0.0;

        Eigen::Vector3d np = boundaryElements.triangles.at(row).normal;
        Eigen::Vector3d nq = boundaryElements.triangles.at(column).normal;
        double upNq = np.dot(nq);

        for(int i = 0; i < numberOfRegularWeightsAndAbscissa; i++) //gauss quadrature of standard triangle
        {
            double quadratureWeight = regularWeightsAndAbscissa(i,2);
            double quadratureAbscissaNode1 = regularWeightsAndAbscissa(i,0);
            double quadratureAbscissaNode2 = regularWeightsAndAbscissa(i,1);
            double quadratureAbscissaNode3 = (1.0-quadratureAbscissaNode2-quadratureAbscissaNode1);
            Eigen::Vector3d q = (quadratureAbscissaNode1*boundaryElements.triangles.at(column).node1) + (quadratureAbscissaNode2*boundaryElements.triangles.at(column).node2) + (quadratureAbscissaNode3*boundaryElements.triangles.at(column).node3);
            Eigen::Vector3d rVector = boundaryElements.collocationPoints.at(row) - q; //p-q
            double rr = rVector.squaredNorm();
            double r = std::sqrt(rr);
            double rnq = - (rVector.dot(nq))/r;
            double rnp = (rVector.dot(np))/r;
            std::complex<double> ikr = iWavenumber*r;
            std::complex<double> greensFunction = (std::exp(ikr))/(PI4*r);
            Lk += quadratureWeight*greensFunction;
            std::complex<double> dGreensDR = (greensFunction/r)*(ikr-1.0);
            Mk += quadratureWeight*dGreensDR*rnq;
            Mtk += quadratureWeight*dGreensDR*rnp;
            std::complex<double> d2GreensDR2 = (greensFunction/rr)*(2.0-2.0*ikr-wavenumberSquared*rr);
            std::complex<double> d2Gkupnq = dGreensDR*(-(upNq+rnp*rnq)/r)+d2GreensDR2*rnp*rnq;
            Nk += quadratureWeight*d2Gkupnq;
        }
        Lk *= boundaryElements.trianglesArea(column);
        Mk *= boundaryElements.trianglesArea(column);
        Mtk *= boundaryElements.trianglesArea(column);
        Nk *= boundaryElements.trianglesArea(column);
    }
}

void BoundaryElementSolver::BemOperatorsSingularPolarInt(const long triangleIndex, std::complex<double>& Lk, std::complex<double>& Mk, std::complex<double>& Mtk, std::complex<double>& Nk)
{
    Lk = 0.0;
    Mk = 0.0;
    Mtk = 0.0;
    Nk = 0.0;

    //observation Point p is on the same Panel as source points(q)

    static double xAbszissa[16] = {0.0483076656877383162348126,0.1444719615827964934851864,0.2392873622521370745446032,0.3318686022821276497799168,0.4213512761306353453641194,0.5068999089322293900237475,0.5877157572407623290407455,0.6630442669302152009751152,0.7321821187402896803874267,0.7944837959679424069630973,0.8493676137325699701336930,0.8963211557660521239653072,0.9349060759377396891709191,0.9647622555875064307738119,0.9856115115452683354001750,0.9972638618494815635449811};
    static double xWeights[16] = {0.0965400885147278005667648,0.0956387200792748594190820,0.0938443990808045656391802,0.0911738786957638847128686,0.0876520930044038111427715,0.0833119242269467552221991,0.0781938957870703064717409,0.0723457941088485062253994,0.0658222227763618468376501,0.0586840934785355471452836,0.0509980592623761761961632,0.0428358980222266806568786,0.0342738629130214331026877,0.0253920653092620594557526,0.0162743947309056706051706,0.0070186100094700966004071};
    int numberOfXWeightsAndAbszissas = 16;
    VectorTriangle triangle1(boundaryElements.triangles.at(triangleIndex).node1,boundaryElements.triangles.at(triangleIndex).node2,boundaryElements.collocationPoints.at(triangleIndex));
    VectorTriangle triangle2(boundaryElements.triangles.at(triangleIndex).node2,boundaryElements.triangles.at(triangleIndex).node3,boundaryElements.collocationPoints.at(triangleIndex));
    VectorTriangle triangle3(boundaryElements.triangles.at(triangleIndex).node3,boundaryElements.triangles.at(triangleIndex).node1,boundaryElements.collocationPoints.at(triangleIndex));
    triangle1.normal=boundaryElements.triangles.at(triangleIndex).normal;
    triangle2.normal=triangle1.normal;
    triangle3.normal=triangle1.normal;

    std::complex<double> Ltest1,Ltest2,Ltest3;
    std::complex<double> Nktest1,Nktest2,Nktest3;

    Eigen::Vector3d triangleCSide=(triangle1.node3-triangle1.node1);
    Eigen::Vector3d triangleBSide=(triangle1.node3-triangle1.node2);
    Eigen::Vector3d triangleASide=(triangle1.node2-triangle1.node1);
    double bSquared=triangleBSide.dot(triangleBSide);
    double cSquared=triangleCSide.dot(triangleCSide);
    double aSquared=triangleASide.dot(triangleASide);
    double b=std::sqrt(bSquared);
    double c=std::sqrt(cSquared);
//        double a=std::sqrt(aSquared);
    if(c<b)
    {
        double tmp=c;
        c=b;
        b=tmp;
    }
    double A=std::acos((bSquared+cSquared-aSquared)/(2.0*b*c));
    double B=std::atan((b*std::sin(A))/(c-b*std::cos(A)));


    double lower = 0.5*(A-0);
    double upper = 0.5*(A+0) ;
    for(int i=0; i<numberOfXWeightsAndAbszissas; i++)
    {
        double xVal1 = xAbszissa[i]*lower+upper;
        double xVal2 = -xAbszissa[i]*lower+upper;
        double r1 = c*sin(B)/(sin(xVal1+B));
        double r2 = c*sin(B)/(sin(xVal2+B));
        std::complex<double> ikr1 = iWavenumber*r1;
        std::complex<double> ikr2 = iWavenumber*r2;
        Ltest1+=xWeights[i]*(std::exp(ikr1)+std::exp(ikr2));
        Nktest1+=xWeights[i]*(std::exp(ikr1)/r1+std::exp(ikr2)/r2);
    }
    Ltest1 = Ltest1*lower;
    Nktest1 = Nktest1*lower;

    triangleCSide = (triangle2.node3-triangle2.node1);
    triangleBSide = (triangle2.node3-triangle2.node2);
    triangleASide = (triangle2.node2-triangle2.node1);
    bSquared = triangleBSide.dot(triangleBSide);
    cSquared = triangleCSide.dot(triangleCSide);
    aSquared = triangleASide.dot(triangleASide);
    b = std::sqrt(bSquared);
    c = std::sqrt(cSquared);
    if(c<b)
    {
        double tmp = c;
        c = b;
        b = tmp;
    }
    double AA = A;
    A = std::acos((bSquared+cSquared-aSquared)/(2.0*b*c));
    B = std::atan((b*std::sin(A))/(c-b*std::cos(A)));

    lower  =  0.5*(A);
    upper  =  0.5*(A+AA+AA);
    for(int i = 0; i<numberOfXWeightsAndAbszissas; i++)
    {
        double xVal1 = xAbszissa[i]*lower+upper;
        double xVal2 = -xAbszissa[i]*lower+upper;
        double r1 = c*sin(B)/(sin(xVal1-AA+B));
        double r2 = c*sin(B)/(sin(xVal2-AA+B));
        std::complex<double> ikr1 = iWavenumber*r1;
        std::complex<double> ikr2 = iWavenumber*r2;
        Ltest2+=xWeights[i]*(std::exp(ikr1)+std::exp(ikr2));
        Nktest2+=xWeights[i]*(std::exp(ikr1)/r1+std::exp(ikr2)/r2);
    }
    Ltest2 = Ltest2*lower;
    Nktest2 = Nktest2*lower;

    triangleCSide = (triangle3.node3-triangle3.node1);
    triangleBSide = (triangle3.node3-triangle3.node2);
    triangleASide = (triangle3.node2-triangle3.node1);
    bSquared = triangleBSide.dot(triangleBSide);
    cSquared = triangleCSide.dot(triangleCSide);
    aSquared = triangleASide.dot(triangleASide);
    b = std::sqrt(bSquared);
    c = std::sqrt(cSquared);
    if(c<b)
    {
        double tmp = c;
        c = b;
        b = tmp;
    }
    double AAA = AA+A;
    A = std::acos((bSquared+cSquared-aSquared)/(2.0*b*c));
    B = std::atan((b*std::sin(A))/(c-b*std::cos(A)));

    lower  =  0.5*(A);
    upper  =  0.5*(A+AAA+AAA);
    for(int i = 0; i<numberOfXWeightsAndAbszissas; i++)
    {
        double xVal1 = xAbszissa[i]*lower+upper;
        double xVal2 = -xAbszissa[i]*lower+upper;
        double r1 = c*sin(B)/(sin(xVal1-AAA+B));
        double r2 = c*sin(B)/(sin(xVal2-AAA+B));
        std::complex<double> ikr1 = iWavenumber*r1;
        std::complex<double> ikr2 = iWavenumber*r2;
        Ltest3+=xWeights[i]*(std::exp(ikr1)+std::exp(ikr2));
        Nktest3+=xWeights[i]*(std::exp(ikr1)/r1+std::exp(ikr2)/r2);
    }

    Ltest3 = Ltest3*lower;
    Nktest3 = Nktest3*lower;
    std::complex<double> LTest = (Ltest1+Ltest2+Ltest3)/PI2;
    std::complex<double> NkTest = (Nktest1+Nktest2+Nktest3)/PI4;

    LTest = (imaginaryUnit/(2.0*wavenumber))*(1.0-LTest);
    NkTest = (iWavenumber/2.0-NkTest);

    Lk = LTest;
    Nk = NkTest;
}

void BoundaryElementSolver::BemOperatorsSingularitySubtraction(const long row, const long column, std::complex<double>& Lk, std::complex<double>& Mk, std::complex<double>& Mtk, std::complex<double>& Nk)
{        //observation Point p is on the same Panel as source points(q)
    Lk = 0.0;
    Mk = 0.0;
    Mtk = 0.0;
    Nk = 0.0;

    Eigen::Vector3d np = boundaryElements.triangles.at(row).normal;

    VectorTriangle triangle1(boundaryElements.triangles.at(column).node1,boundaryElements.triangles.at(column).node2,boundaryElements.collocationPoints.at(column));
    VectorTriangle triangle2(boundaryElements.triangles.at(column).node2,boundaryElements.triangles.at(column).node3,boundaryElements.collocationPoints.at(column));
    VectorTriangle triangle3(boundaryElements.triangles.at(column).node3,boundaryElements.triangles.at(column).node1,boundaryElements.collocationPoints.at(column));
    triangle1.normal = boundaryElements.triangles.at(column).normal;
    triangle2.normal = triangle1.normal;
    triangle3.normal = triangle1.normal;
    double triangle1Area = global::areaOfTriangle(triangle1);
    double triangle2Area = global::areaOfTriangle(triangle2);
    double triangle3Area = global::areaOfTriangle(triangle3);

    std::complex<double> Nk1,Nk2,Nk3;
    std::complex<double> Lk1,Lk2,Lk3;

    Eigen::Vector3d triangleCSide = (triangle1.node3-triangle1.node1);
    Eigen::Vector3d triangleBSide = (triangle1.node3-triangle1.node2);
    Eigen::Vector3d triangleASide = (triangle1.node2-triangle1.node1);
    double bSquared = triangleBSide.dot(triangleBSide);
    double cSquared = triangleCSide.dot(triangleCSide);
    double aSquared = triangleASide.dot(triangleASide);
    double b = std::sqrt(bSquared);
    double c = std::sqrt(cSquared);
    if(c<b)
    {
        double tmp = c;
        c = b;
        b = tmp;
    }
    double A = std::acos((bSquared+cSquared-aSquared)/(2.0*b*c));
    double B = std::atan((b*std::sin(A))/(c-b*std::cos(A)));
    double cSinB = c*std::sin(B);
    std::complex<double> L0e = cSinB*(std::log((std::complex<double>)std::tan((B+A)/2.0))-std::log((std::complex<double>)std::tan(B/2.0)));
    double N0e = (std::cos(B+A)-std::cos(B))/(cSinB);

    triangleCSide = (triangle2.node3-triangle2.node1);
    triangleBSide = (triangle2.node3-triangle2.node2);
    triangleASide = (triangle2.node2-triangle2.node1);
    bSquared = triangleBSide.dot(triangleBSide);
    cSquared = triangleCSide.dot(triangleCSide);
    aSquared = triangleASide.dot(triangleASide);
    b = std::sqrt(bSquared);
    c = std::sqrt(cSquared);
    if(c<b)
    {
        double tmp = c;
        c = b;
        b = tmp;
    }
    A = std::acos((bSquared+cSquared-aSquared)/(2.0*b*c));
    B = std::atan((b*std::sin(A))/(c-b*std::cos(A)));

    cSinB = c*std::sin(B);
    L0e += cSinB*(std::log((std::complex<double>)std::tan((B+A)/2.0))-std::log((std::complex<double>)std::tan(B/2.0)));
    N0e += (std::cos(B+A)-std::cos(B))/(cSinB);

    triangleCSide = (triangle3.node3-triangle3.node1);
    triangleBSide = (triangle3.node3-triangle3.node2);
    triangleASide = (triangle3.node2-triangle3.node1);
    bSquared = triangleBSide.dot(triangleBSide);
    cSquared = triangleCSide.dot(triangleCSide);
    aSquared = triangleASide.dot(triangleASide);
    b = std::sqrt(bSquared);
    c = std::sqrt(cSquared);
    if(c<b)
    {
        double tmp=c;
        c=b;
        b=tmp;
    }
    A = std::acos((bSquared+cSquared-aSquared)/(2.0*b*c));
    B = std::atan((b*std::sin(A))/(c-b*std::cos(A)));
    cSinB = c*std::sin(B);
    L0e += cSinB*(std::log((std::complex<double>)std::tan((B+A)/2.0))-std::log((std::complex<double>)std::tan(B/2.0)));
    N0e += (std::cos(B+A)-std::cos(B))/(cSinB);

    L0e = L0e/PI4;
    N0e = N0e/PI4;

    for(int i=0; i<numberOfHighOrderWeightsAndAbscissa; i++) //gauss quadrature of three standard triangles
    {
        double quadratureWeight = highOrderweightsAndAbscissa(i,2);
        double quadratureAbscissaNode1 = highOrderweightsAndAbscissa(i,0);
        double quadratureAbscissaNode2 = highOrderweightsAndAbscissa(i,1);
        double quadratureAbscissaNode3 = (1.0-quadratureAbscissaNode2-quadratureAbscissaNode1);

        Eigen::Vector3d q1 = (quadratureAbscissaNode1*triangle1.node1) + (quadratureAbscissaNode2*triangle1.node2) + (quadratureAbscissaNode3*triangle1.node3);
        Eigen::Vector3d q2 = (quadratureAbscissaNode1*triangle2.node1) + (quadratureAbscissaNode2*triangle2.node2) + (quadratureAbscissaNode3*triangle2.node3);
        Eigen::Vector3d q3 = (quadratureAbscissaNode1*triangle3.node1) + (quadratureAbscissaNode2*triangle3.node2) + (quadratureAbscissaNode3*triangle3.node3);
        Eigen::Vector3d rVector1 = boundaryElements.collocationPoints.at(row)-q1;
        Eigen::Vector3d rVector2 = boundaryElements.collocationPoints.at(row)-q2;
        Eigen::Vector3d rVector3 = boundaryElements.collocationPoints.at(row)-q3;
        double upNq = np.dot(triangle1.normal);

        double rr1 = rVector1.squaredNorm();
        double rr2 = rVector2.squaredNorm();
        double rr3 = rVector3.squaredNorm();
        double r1 = std::sqrt(rr1);
        double r2 = std::sqrt(rr2);
        double r3 = std::sqrt(rr3);
        double rrr1 = rr1*r1;
        double rrr2 = rr2*r2;
        double rrr3 = rr3*r3;
        double rnq1 = -(rVector1.dot(triangle1.normal))/r1;
        double rnq2 = -(rVector2.dot(triangle2.normal))/r2;
        double rnq3 = -(rVector3.dot(triangle3.normal))/r3;
        double rup1 = (rVector1.dot(np))/r1;
        double rup2 = (rVector2.dot(np))/r2;
        double rup3 = (rVector3.dot(np))/r3;

        std::complex<double> ikr1 = iWavenumber*r1;
        std::complex<double> ikr2 = iWavenumber*r2;
        std::complex<double> ikr3 = iWavenumber*r3;
        std::complex<double> greensFunction1 = (std::exp(ikr1))/(PI4*r1);
        std::complex<double> greensFunction2 = (std::exp(ikr2))/(PI4*r2);
        std::complex<double> greensFunction3 = (std::exp(ikr3))/(PI4*r3);
        std::complex<double> greensK0Function1 = (1.0)/(PI4*r1);
        std::complex<double> greensK0Function2 = (1.0)/(PI4*r2);
        std::complex<double> greensK0Function3 = (1.0)/(PI4*r3);
        Lk1 += quadratureWeight*(greensFunction1-greensK0Function1);
        Lk2 += quadratureWeight*(greensFunction2-greensK0Function2);
        Lk3 += quadratureWeight*(greensFunction3-greensK0Function3);


        std::complex<double> dGreensDR1 = (greensFunction1/r1)*(ikr1-1.0);
        std::complex<double> dGreensDR2 = (greensFunction2/r2)*(ikr2-1.0);
        std::complex<double> dGreensDR3 = (greensFunction3/r3)*(ikr3-1.0);

        std::complex<double> d2GreensDR21 = (greensFunction1/rr1)*(2.0-2.0*ikr1-wavenumberSquared*rr1);
        std::complex<double> d2GreensDR22 = (greensFunction2/rr2)*(2.0-2.0*ikr2-wavenumberSquared*rr2);
        std::complex<double> d2GreensDR23 = (greensFunction3/rr3)*(2.0-2.0*ikr3-wavenumberSquared*rr3);
        std::complex<double> d2GkG0upnq1 = (dGreensDR1+1.0/(PI4*rr1))*(-(upNq+rup1*rnq1)/r1) + (d2GreensDR21-1.0/(PI2*rrr1))*rup1*rnq1;
        std::complex<double> d2GkG0upnq2 = (dGreensDR2+1.0/(PI4*rr2))*(-(upNq+rup2*rnq2)/r2) + (d2GreensDR22-1.0/(PI2*rrr2))*rup2*rnq2;
        std::complex<double> d2GkG0upnq3 = (dGreensDR3+1.0/(PI4*rr3))*(-(upNq+rup3*rnq3)/r3) + (d2GreensDR23-1.0/(PI2*rrr3))*rup3*rnq3;
        Nk1 += quadratureWeight*(d2GkG0upnq1+greensK0Function1*wavenumberSquaredHalf);
        Nk2 += quadratureWeight*(d2GkG0upnq2+greensK0Function2*wavenumberSquaredHalf);
        Nk3 += quadratureWeight*(d2GkG0upnq3+greensK0Function3*wavenumberSquaredHalf);
    }

    Lk1 *= triangle1Area;
    Nk1 *= triangle1Area;

    Lk2 *= triangle2Area;
    Nk2 *= triangle2Area;

    Lk3 *= triangle3Area;
    Nk3 *= triangle3Area;


    Lk = Lk1 + Lk2 + Lk3 + L0e;
//        Mk=0; // -> follows from asymptotic properties of opereators when lim p->q
//        Mtk=0; // -> follows from asymptotic properties of opereators when lim p->q
    Nk = Nk1 + Nk2 + Nk3 + N0e - wavenumberSquaredHalf * L0e;
}

void BoundaryElementSolver::BemOperatorsNearSingSinh(const long row, const long column, std::complex<double>& Lk, std::complex<double>& Mk, std::complex<double>& Mtk, std::complex<double>& Nk)
{
    // routine to evaluate the nearly singular BEM operators found in
    // "A new method for the numerical evaluation of nearly singular integrals on triangular elements in the 3D boundary element method"

    Lk = 0.0;
    Mk = 0.0;
    Mtk = 0.0;
    Nk = 0.0;

    Eigen::Vector3d p = boundaryElements.collocationPoints.at(row);
    Eigen::Vector3d np = boundaryElements.triangles.at(row).normal;
    VectorTriangle triangleDomain = boundaryElements.triangles.at(column);
    Eigen::Vector3d nq = boundaryElements.triangles.at(column).normal;

    const auto& [pProjectedOnTriangleDomain, projectionRegion] = global::projectPointOnTriangleFaster(triangleDomain, p); // project the point p onto the triangle triangleDomain
    // pProjectedOnTriangleDomain is the projected point
    // projectionRegion specifies the region of the triangle in which the projection point lies
    double distToSingularity = (p - pProjectedOnTriangleDomain).norm();

    QVector<VectorTriangle> splitTriangles;

    switch (projectionRegion)
    {
        case global::inInterior:
            splitTriangles.resize(3);
            splitTriangles[0] = VectorTriangle(pProjectedOnTriangleDomain, triangleDomain.node1, triangleDomain.node2);
            splitTriangles[1] = VectorTriangle(pProjectedOnTriangleDomain, triangleDomain.node2, triangleDomain.node3);
            splitTriangles[2] = VectorTriangle(pProjectedOnTriangleDomain, triangleDomain.node3, triangleDomain.node1);
            break;
        case global::onNode1:
            splitTriangles.resize(1);
            splitTriangles[0] = triangleDomain;
            break;
        case global::onNode2:
            splitTriangles.resize(1);
            splitTriangles[0] = VectorTriangle(triangleDomain.node2, triangleDomain.node3, triangleDomain.node1);
            break;
        case global::onNode3:
            splitTriangles.resize(1);
            splitTriangles[0] = VectorTriangle(triangleDomain.node3, triangleDomain.node1, triangleDomain.node2);
            break;
        case global::onSide21:
            splitTriangles.resize(2);
            splitTriangles[0] = VectorTriangle(pProjectedOnTriangleDomain, triangleDomain.node3, triangleDomain.node1);
            splitTriangles[1] = VectorTriangle(pProjectedOnTriangleDomain, triangleDomain.node2, triangleDomain.node3);
            break;
        case global::onSide32:
            splitTriangles.resize(2);
            splitTriangles[0] = VectorTriangle(pProjectedOnTriangleDomain, triangleDomain.node3, triangleDomain.node1);
            splitTriangles[1] = VectorTriangle(pProjectedOnTriangleDomain, triangleDomain.node1, triangleDomain.node2);
            break;
        case global::onSide13:
            splitTriangles.resize(2);
            splitTriangles[0] = VectorTriangle(pProjectedOnTriangleDomain, triangleDomain.node1, triangleDomain.node2);
            splitTriangles[1] = VectorTriangle(pProjectedOnTriangleDomain, triangleDomain.node2, triangleDomain.node3);
            break;
    }

    for(int i = 0; i < splitTriangles.size(); i++) // quadrature over each triangle
    {
        std::complex<double> tmpMk = 0.0;
        std::complex<double> tmpNk = 0.0;
        std::complex<double> tmpLk = 0.0;
        std::complex<double> tmpMtk = 0.0;
        VectorTriangle currentTriangle = splitTriangles.at(i);
//        std::cerr << "triangleArea: " << triangleArea << std:: endl;

        if(false)
        {
            double triangleArea = global::areaOfTriangle(currentTriangle);
            double upNq = np.dot(nq);

            for(int i = 0; i < numberOfRegularWeightsAndAbscissa; i++)
            {
                double quadratureWeight = regularWeightsAndAbscissa(i,2);
                double quadratureAbscissaNode1 = regularWeightsAndAbscissa(i,0);
                double quadratureAbscissaNode2 = regularWeightsAndAbscissa(i,1);
                double quadratureAbscissaNode3 = (1.0-quadratureAbscissaNode2-quadratureAbscissaNode1);
                Eigen::Vector3d q = (quadratureAbscissaNode1*currentTriangle.node1) + (quadratureAbscissaNode2*currentTriangle.node2) + (quadratureAbscissaNode3*currentTriangle.node3);
                Eigen::Vector3d rVector = p - q; //p-q
                double rr = rVector.squaredNorm();
                double r = std::sqrt(rr);
                double rnq = - (rVector.dot(nq))/r;
                double rnp = (rVector.dot(np))/r;
                std::complex<double> ikr = iWavenumber*r;
                std::complex<double> greensFunction = (std::exp(ikr))/(PI4*r);
                tmpLk += quadratureWeight*greensFunction;
                std::complex<double> dGreensDR = (greensFunction/r)*(ikr-1.0);
                tmpMk += quadratureWeight*dGreensDR*rnq;
                tmpMtk += quadratureWeight*dGreensDR*rnp;
                std::complex<double> d2GreensDR2 = (greensFunction/rr)*(2.0-2.0*ikr-wavenumberSquared*rr);
                std::complex<double> d2Gkupnq = dGreensDR*(-(upNq+rnp*rnq)/r)+d2GreensDR2*rnp*rnq;
                tmpNk += quadratureWeight*d2Gkupnq;
            }
            Lk += tmpLk * triangleArea;
            Mk += tmpMk * triangleArea;
            Mtk += tmpMtk * triangleArea;
            Nk += tmpNk * triangleArea;
        }
        else
        {
            Eigen::Vector3d triangle13Side = (currentTriangle.node3-currentTriangle.node1);
            Eigen::Vector3d triangle23Side = (currentTriangle.node3-currentTriangle.node2);
            Eigen::Vector3d triangle12Side = (currentTriangle.node2-currentTriangle.node1);
            double side13Squared = triangle13Side.squaredNorm();
            double side23Squared = triangle23Side.squaredNorm();
            double side12Squared = triangle12Side.squaredNorm();
            double side13Length = std::sqrt(side13Squared);
            double side23Length = std::sqrt(side23Squared);
            double side12Length = std::sqrt(side12Squared);
            double angleAtProjectionPoint = std::acos((side13Squared+side12Squared-side23Squared)/(2.0*side13Length*side12Length));
            double angleAtNode2 = std::acos((side23Squared+side12Squared-side13Squared)/(2.0*side23Length*side12Length));

            ////////////// for the conversion to global coordinates
            ///
//            double triangleHeightAtNode3 = std::sin(angleAtProjectionPoint) * side13Length;
//            double node3ProjectedOnSide12Length = std::cos(angleAtProjectionPoint) * side13Length;

            Eigen::Matrix3Xd rotationMatrix(3,2); // rotation matrix to convert from coordinates in triangle plane to the global coordinates
            rotationMatrix.col(0) = (triangle12Side.normalized());
//            rotationMatrix.col(2) = rotationMatrix.col(0).cross(triangle13Side).normalized();
//            rotationMatrix.col(1) = rotationMatrix.col(2).cross(rotationMatrix.col(0)).normalized();
            rotationMatrix.col(1) = (triangle13Side - rotationMatrix.col(0).dot(triangle13Side) * rotationMatrix.col(0)).normalized();

            ///
            /// ///////////////////////////////////////////////
//            double angleAtProjectionPoint =
            for(int i = 0; i < numberOfWeightsAndAbscissaLineQuadrature; i++)
            {
                double quadratureAbszissaForAngle = weightsAndAbscissaLineQuadrature(i,0);
                double quadratureWeightForAngle = weightsAndAbscissaLineQuadrature(i,1);


                /////////////////////// standard polar coordinates
//                double angle1 = (angleAtProjectionPoint + quadratureAbszissaForAngle * angleAtProjectionPoint)/2.0;
//                double maxRForAngle1 = side12Length*std::sin(angleAtNode2)/(std::sin(angle1 + angleAtNode2));
//                double mu1 = std::asinh(maxRForAngle1/distToSingularity) / 2.0;
                //////////////////////
                ///
                double sigma = (quadratureAbszissaForAngle + 1) / 2.0;
                double gamma = std::pow(sigma, 2) / (std::pow(sigma, 2) + std::pow(1 - sigma, 2) );
                double dGamma = 2.0 * sigma * (1 - sigma) /  std::pow(std::pow(sigma, 2) + std::pow(1 - sigma, 2), 2);
                double angle1 = angleAtProjectionPoint * gamma;
                double maxRForAngle1 = side12Length*std::sin(angleAtNode2)/(std::sin(angle1 + angleAtNode2));
                double mu1 = std::asinh(maxRForAngle1/distToSingularity) / 2.0;

                for(int j = 0; j < numberOfWeightsAndAbscissaLineQuadrature; j++)
                {
                    double s = weightsAndAbscissaLineQuadrature(j,0);
                    double quadratureWeightForS = weightsAndAbscissaLineQuadrature(j,1);

                    double rs = distToSingularity * std::sinh(mu1 * (s + 1));
//                    double rs = (maxRForAngle1 + weightsAndAbscissaLineQuadrature(j,0) * maxRForAngle1) / 2.0;

                    double x = rs * std::cos(angle1);
                    double y = rs * std::sin(angle1);

                    Eigen::Vector3d q = currentTriangle.node1 + x * rotationMatrix.col(0) + y * rotationMatrix.col(1);

//                    double factorAndWeights = rs * quadratureWeightForAngle * quadratureWeightForS * distToSingularity * mu1 * std::cosh( mu1 * (s + 1));
                    double factorAndWeights = dGamma * rs * quadratureWeightForAngle * quadratureWeightForS * distToSingularity * mu1 * std::cosh( mu1 * (s + 1));

                    Eigen::Vector3d rVector = p - q; //p-q
                    double upNq = np.dot(nq);
//                    double rr = x*x + y*y + distToSingularity*distToSingularity;
                    double rr = rVector.squaredNorm();
                    double r = std::sqrt(rr);
                    double rnq = - (rVector.dot(nq))/r;
                    double rnp = (rVector.dot(np))/r;
                    std::complex<double> ikr = iWavenumber*r;
                    std::complex<double> greensFunction = (std::exp(ikr))/(PI4*r);
                    tmpLk += factorAndWeights*greensFunction;
                    std::complex<double> dGreensDR = (greensFunction/r)*(ikr-1.0);
                    tmpMk += factorAndWeights*dGreensDR*rnq;
                    tmpMtk += factorAndWeights*dGreensDR*rnp;
                    std::complex<double> d2GreensDR2 = (greensFunction/rr)*(2.0-2.0*ikr-wavenumberSquared*rr);
                    std::complex<double> d2Gkupnq = dGreensDR*(-(upNq+rnp*rnq)/r)+d2GreensDR2*rnp*rnq;
                    tmpNk += factorAndWeights*d2Gkupnq;
                }
            }
//            Lk += std::pow(distToSingularity, (2-l)) * angleAtProjectionPoint / 2.0 * tmpLk;
//            Mk += std::pow(distToSingularity, (2-l)) * angleAtProjectionPoint / 2.0 * tmpMk;
//            Mtk += std::pow(distToSingularity, (2-l)) * angleAtProjectionPoint / 2.0 * tmpMtk;
//            Nk += std::pow(distToSingularity, (2-l)) * angleAtProjectionPoint / 2.0 * tmpNk;

            Lk += angleAtProjectionPoint / 2.0 * tmpLk;
            Mk += angleAtProjectionPoint / 2.0 * tmpMk;
            Mtk += angleAtProjectionPoint / 2.0 * tmpMtk;
            Nk += angleAtProjectionPoint / 2.0 * tmpNk;
        }
    }
}

void BoundaryElementSolver::BemOperatorsNearSing(const long row, const long column, std::complex<double>& Lk, std::complex<double>& Mk, std::complex<double>& Mtk, std::complex<double>& Nk)
{
    // routine to evaluate the nearly singular BEM operators found in
    // "A new method for the numerical evaluation of nearly singular integrals on triangular elements in the 3D boundary element method"

    Lk = 0.0;
    Mk = 0.0;
    Mtk = 0.0;
    Nk = 0.0;

    Eigen::Vector3d p = boundaryElements.collocationPoints.at(row);
    Eigen::Vector3d np = boundaryElements.triangles.at(row).normal;
    Eigen::Vector3d nq = boundaryElements.triangles.at(column).normal;
    double upNq = np.dot(nq);
    VectorTriangle triangleDomain = boundaryElements.triangles.at(column);


    const auto& [pProjectedOnTriangleDomain, projectionRegion] = global::projectPointOnTriangleFaster(triangleDomain, p); // project the point p onto the triangle triangleDomain
    // pProjectedOnTriangleDomain is the projected point
    // projectionRegion specifies the region of the triangle in which the projection point lies

    QVector<VectorTriangle> splitTriangles;

    switch (projectionRegion)
    {
        case global::inInterior:
            splitTriangles.resize(3);
            splitTriangles[0] = VectorTriangle(pProjectedOnTriangleDomain, triangleDomain.node1, triangleDomain.node2);
            splitTriangles[1] = VectorTriangle(pProjectedOnTriangleDomain, triangleDomain.node2, triangleDomain.node3);
            splitTriangles[2] = VectorTriangle(pProjectedOnTriangleDomain, triangleDomain.node3, triangleDomain.node1);
            break;
        case global::onNode1:
            splitTriangles.resize(1);
            splitTriangles[0] = triangleDomain;
            break;
        case global::onNode2:
            splitTriangles.resize(1);
            splitTriangles[0] = VectorTriangle(triangleDomain.node2, triangleDomain.node3, triangleDomain.node1);
            break;
        case global::onNode3:
            splitTriangles.resize(1);
            splitTriangles[0] = VectorTriangle(triangleDomain.node3, triangleDomain.node1, triangleDomain.node2);
            break;
        case global::onSide21:
            splitTriangles.resize(2);
            splitTriangles[0] = VectorTriangle(pProjectedOnTriangleDomain, triangleDomain.node3, triangleDomain.node1);
            splitTriangles[1] = VectorTriangle(pProjectedOnTriangleDomain, triangleDomain.node2, triangleDomain.node3);
            break;
        case global::onSide32:
            splitTriangles.resize(2);
            splitTriangles[0] = VectorTriangle(pProjectedOnTriangleDomain, triangleDomain.node3, triangleDomain.node1);
            splitTriangles[1] = VectorTriangle(pProjectedOnTriangleDomain, triangleDomain.node1, triangleDomain.node2);
            break;
        case global::onSide13:
            splitTriangles.resize(2);
            splitTriangles[0] = VectorTriangle(pProjectedOnTriangleDomain, triangleDomain.node1, triangleDomain.node2);
            splitTriangles[1] = VectorTriangle(pProjectedOnTriangleDomain, triangleDomain.node2, triangleDomain.node3);
            break;
    }

    for(int i = 0; i < splitTriangles.size(); i++) // quadrature over each triangle
    {
        std::complex<double> tmpMk = 0.0;
        std::complex<double> tmpNk = 0.0;
        std::complex<double> tmpLk = 0.0;
        std::complex<double> tmpMtk = 0.0;
        VectorTriangle currentTriangle = splitTriangles.at(i);

        Eigen::Vector3d triangle13Side = (currentTriangle.node3-currentTriangle.node1);
        Eigen::Vector3d triangle23Side = (currentTriangle.node3-currentTriangle.node2);
        Eigen::Vector3d triangle12Side = (currentTriangle.node2-currentTriangle.node1);
        double side13Squared = triangle13Side.squaredNorm();
        double side23Squared = triangle23Side.squaredNorm();
        double side12Squared = triangle12Side.squaredNorm();
        double side13Length = std::sqrt(side13Squared);
        double side23Length = std::sqrt(side23Squared);
        double side12Length = std::sqrt(side12Squared);
        double angleAtProjectionPoint = std::acos((side13Squared+side12Squared-side23Squared)/(2.0*side13Length*side12Length));
        double angleAtNode2 = std::acos((side23Squared+side12Squared-side13Squared)/(2.0*side23Length*side12Length));

        ////////////// for the conversion to global coordinates
        Eigen::Matrix3Xd rotationMatrix(3,2); // rotation matrix to convert from coordinates in triangle plane to the global coordinates
        rotationMatrix.col(0) = (triangle12Side.normalized());
        rotationMatrix.col(1) = (triangle13Side - rotationMatrix.col(0).dot(triangle13Side) * rotationMatrix.col(0)).normalized();

        for(int i = 0; i < numberOfWeightsAndAbscissaLineQuadrature; i++)
        {
            double quadratureAbszissaForAngle = weightsAndAbscissaLineQuadrature(i,0);
            double quadratureWeightForAngle = weightsAndAbscissaLineQuadrature(i,1);

            double angle1 = (angleAtProjectionPoint + quadratureAbszissaForAngle * angleAtProjectionPoint)/2.0;
            double maxRForAngle1 = side12Length*std::sin(angleAtNode2)/(std::sin(angle1 + angleAtNode2));

            for(int j = 0; j < numberOfWeightsAndAbscissaLineQuadrature; j++)
            {
                double quadratureWeightForS = weightsAndAbscissaLineQuadrature(j,1);

                double rs = (maxRForAngle1 + weightsAndAbscissaLineQuadrature(j,0) * maxRForAngle1) / 2.0;

                double x = rs * std::cos(angle1);
                double y = rs * std::sin(angle1);

                Eigen::Vector3d q = currentTriangle.node1 + x * rotationMatrix.col(0) + y * rotationMatrix.col(1);

                double factorAndWeights = maxRForAngle1/2.0 * rs * quadratureWeightForAngle * quadratureWeightForS;

                Eigen::Vector3d rVector = p - q; //p-q
                double rr = rVector.squaredNorm();
                double r = std::sqrt(rr);
                double rnq = - (rVector.dot(nq))/r;
                double rnp = (rVector.dot(np))/r;
                std::complex<double> ikr = iWavenumber*r;
                std::complex<double> greensFunction = (std::exp(ikr))/(PI4*r);
                tmpLk += factorAndWeights*greensFunction;
                std::complex<double> dGreensDR = (greensFunction/r)*(ikr-1.0);
                tmpMk += factorAndWeights*dGreensDR*rnq;
                tmpMtk += factorAndWeights*dGreensDR*rnp;
                std::complex<double> d2GreensDR2 = (greensFunction/rr)*(2.0-2.0*ikr-wavenumberSquared*rr);
                std::complex<double> d2Gkupnq = dGreensDR*(-(upNq+rnp*rnq)/r)+d2GreensDR2*rnp*rnq;
                tmpNk += factorAndWeights*d2Gkupnq;
            }
        }
        Lk += angleAtProjectionPoint / 2.0 * tmpLk;
        Mk += angleAtProjectionPoint / 2.0 * tmpMk;
        Mtk += angleAtProjectionPoint / 2.0 * tmpMtk;
        Nk += angleAtProjectionPoint / 2.0 * tmpNk;
    }
}

void BoundaryElementSolver::BemOperatorsReflected(const long row, const long column, std::complex<double>& Lk, std::complex<double>& Mk, std::complex<double>& Mtk, std::complex<double>& Nk, const BoundaryElements &boundaryElements, const BoundaryElements &reflectedElements)
{
    Lk = 0.0;
    Mk = 0.0;
    Mtk = 0.0;
    Nk = 0.0;

    Eigen::Vector3d np = boundaryElements.triangles.at(row).normal;

    for(int i = 0; i < numberOfRegularWeightsAndAbscissa; i++) //gauss quadrature of standard triangle
    {
        Eigen::Vector3d nq = reflectedElements.triangles.at(column).normal;
        double quadratureWeight = regularWeightsAndAbscissa(i,2);
        double quadratureAbscissaNode1 = regularWeightsAndAbscissa(i,0);
        double quadratureAbscissaNode2 = regularWeightsAndAbscissa(i,1);
        double quadratureAbscissaNode3 = (1.0-quadratureAbscissaNode2-quadratureAbscissaNode1);
        Eigen::Vector3d q = (quadratureAbscissaNode1*reflectedElements.triangles.at(column).node1) + (quadratureAbscissaNode2*reflectedElements.triangles.at(column).node2) + (quadratureAbscissaNode3*reflectedElements.triangles.at(column).node3);
        Eigen::Vector3d rVector = boundaryElements.collocationPoints.at(row)-q; //p-q
        double upNq = np.dot(nq);
        double rr = rVector.squaredNorm();
        double r = std::sqrt(rr);
        double rnq = - (rVector.dot(nq))/r;
        double rnp = (rVector.dot(np))/r;
        std::complex<double> ikr = iWavenumber*r;
        std::complex<double> greensFunction = (std::exp(ikr))/(PI4*r);
        Lk += quadratureWeight*greensFunction;
        std::complex<double> dGreensDR = (greensFunction/r)*(ikr-1.0);
        Mk += quadratureWeight*dGreensDR*rnq;
        Mtk += quadratureWeight*dGreensDR*rnp;
        std::complex<double> d2GreensDR2 = (greensFunction/rr)*(2.0-2.0*ikr-wavenumberSquared*rr);
        std::complex<double> d2Gkupnq = dGreensDR*(-(upNq+rnp*rnq)/r)+d2GreensDR2*rnp*rnq;
        Nk += quadratureWeight*d2Gkupnq;
    }
    Lk = Lk*reflectedElements.trianglesArea(column);
    Mk = Mk*reflectedElements.trianglesArea(column);
    Mtk = Mtk*reflectedElements.trianglesArea(column);
    Nk = Nk*reflectedElements.trianglesArea(column);
}


std::complex<double> BoundaryElementSolver::implicitPhiMatrixReflected(const long rowIndex, const long columnIndex, const BoundaryElements &boundaryElements, const BoundaryElements &reflectedElements)
{
    std::complex<double> Lk = 0.0;
    std::complex<double> Mk = 0.0;
    std::complex<double> Mtk = 0.0;
    std::complex<double> Nk = 0.0;

    Eigen::Vector3d np = boundaryElements.triangles.at(rowIndex).normal;
    Eigen::Vector3d nq = reflectedElements.triangles.at(columnIndex).normal;
    double upNq = np.dot(nq);

    for(int i = 0; i < numberOfRegularWeightsAndAbscissa; i++) //gauss quadrature of standard triangle
    {
        double quadratureWeight = regularWeightsAndAbscissa(i,2);
        double quadratureAbscissaNode1 = regularWeightsAndAbscissa(i,0);
        double quadratureAbscissaNode2 = regularWeightsAndAbscissa(i,1);
        double quadratureAbscissaNode3 = (1.0-quadratureAbscissaNode2-quadratureAbscissaNode1);
        Eigen::Vector3d q = (quadratureAbscissaNode1*reflectedElements.triangles.at(columnIndex).node1) + (quadratureAbscissaNode2*reflectedElements.triangles.at(columnIndex).node2) + (quadratureAbscissaNode3*reflectedElements.triangles.at(columnIndex).node3);
        Eigen::Vector3d rVector = boundaryElements.collocationPoints.at(rowIndex) - q; //p-q
        double rr = rVector.squaredNorm();
        double r = std::sqrt(rr);
        double rnq = - (rVector.dot(nq))/r;
        double rnp = (rVector.dot(np))/r;
        std::complex<double> ikr = iWavenumber*r;
        std::complex<double> greensFunction = (std::exp(ikr))/(PI4*r);
        Lk += quadratureWeight*greensFunction;
        std::complex<double> dGreensDR = (greensFunction/r)*(ikr-1.0);
        Mtk += quadratureWeight*dGreensDR*rnp;
        if(!substituteDPhiWithPhi(columnIndex))
        {
            Mk += quadratureWeight*dGreensDR*rnq;
            std::complex<double> d2GreensDR2 = (greensFunction/rr)*(2.0-2.0*ikr-wavenumberSquared*rr);
            std::complex<double> d2Gkupnq = dGreensDR*(-(upNq+rnp*rnq)/r)+d2GreensDR2*rnp*rnq;
            Nk += quadratureWeight*d2Gkupnq;
        }
    }
    if(substituteDPhiWithPhi(columnIndex))
    {
        return - (Lk + couplingParameter * Mtk ) * reflectedElements.trianglesArea(columnIndex) / beta(columnIndex);
    }
    else
    {
        return (Mk + couplingParameter * Nk) * reflectedElements.trianglesArea(columnIndex) / alpha(columnIndex);
    }
}

std::complex<double> BoundaryElementSolver::implicitDPhiMatrixReflected(const long rowIndex, const long columnIndex, const BoundaryElements &boundaryElements, const BoundaryElements &reflectedElements)
{
    std::complex<double> Lk = 0.0;
    std::complex<double> Mk = 0.0;
    std::complex<double> Mtk = 0.0;
    std::complex<double> Nk = 0.0;

    Eigen::Vector3d np = boundaryElements.triangles.at(rowIndex).normal;
    Eigen::Vector3d nq = reflectedElements.triangles.at(columnIndex).normal;
    double upNq = np.dot(nq);

    for(int i = 0; i < numberOfRegularWeightsAndAbscissa; i++) //gauss quadrature of standard triangle
    {
        double quadratureWeight = regularWeightsAndAbscissa(i,2);
        double quadratureAbscissaNode1 = regularWeightsAndAbscissa(i,0);
        double quadratureAbscissaNode2 = regularWeightsAndAbscissa(i,1);
        double quadratureAbscissaNode3 = (1.0-quadratureAbscissaNode2-quadratureAbscissaNode1);
        Eigen::Vector3d q = (quadratureAbscissaNode1*reflectedElements.triangles.at(columnIndex).node1) + (quadratureAbscissaNode2*reflectedElements.triangles.at(columnIndex).node2) + (quadratureAbscissaNode3*reflectedElements.triangles.at(columnIndex).node3);
        Eigen::Vector3d rVector = boundaryElements.collocationPoints.at(rowIndex) - q; //p-q
        double rr = rVector.squaredNorm();
        double r = std::sqrt(rr);
        double rnq = - (rVector.dot(nq))/r;
        double rnp = (rVector.dot(np))/r;
        std::complex<double> ikr = iWavenumber*r;
        std::complex<double> greensFunction = (std::exp(ikr))/(PI4*r);
        Lk += quadratureWeight*greensFunction;
        std::complex<double> dGreensDR = (greensFunction/r)*(ikr-1.0);
        Mk += quadratureWeight*dGreensDR*rnq;
        Mtk += quadratureWeight*dGreensDR*rnp;
        std::complex<double> d2GreensDR2 = (greensFunction/rr)*(2.0-2.0*ikr-wavenumberSquared*rr);
        std::complex<double> d2Gkupnq = dGreensDR*(-(upNq+rnp*rnq)/r)+d2GreensDR2*rnp*rnq;
        Nk += quadratureWeight*d2Gkupnq;
    }
    if(substituteDPhiWithPhi(columnIndex))
    {
        return ( -alpha(columnIndex) / beta(columnIndex) * (Lk + couplingParameter * Mtk ) - (Mk + couplingParameter * Nk)) * reflectedElements.trianglesArea(columnIndex);
    }
    else
    {
        return (Lk + couplingParameter * Mtk + beta(columnIndex)/alpha(columnIndex) * (Mk + couplingParameter * Nk)) * reflectedElements.trianglesArea(columnIndex);
    }
}

std::complex<double> BoundaryElementSolver::implicitPhiMatrix(long rowIndex, long columnIndex)
{
    if(rowIndex == columnIndex)
    {
        static double xAbszissa[16] = {0.0483076656877383162348126,0.1444719615827964934851864,0.2392873622521370745446032,0.3318686022821276497799168,0.4213512761306353453641194,0.5068999089322293900237475,0.5877157572407623290407455,0.6630442669302152009751152,0.7321821187402896803874267,0.7944837959679424069630973,0.8493676137325699701336930,0.8963211557660521239653072,0.9349060759377396891709191,0.9647622555875064307738119,0.9856115115452683354001750,0.9972638618494815635449811};
        static double xWeights[16] = {0.0965400885147278005667648,0.0956387200792748594190820,0.0938443990808045656391802,0.0911738786957638847128686,0.0876520930044038111427715,0.0833119242269467552221991,0.0781938957870703064717409,0.0723457941088485062253994,0.0658222227763618468376501,0.0586840934785355471452836,0.0509980592623761761961632,0.0428358980222266806568786,0.0342738629130214331026877,0.0253920653092620594557526,0.0162743947309056706051706,0.0070186100094700966004071};
        int numberOfXWeightsAndAbszissas = 16;
        VectorTriangle triangle1(boundaryElements.triangles.at(columnIndex).node1,boundaryElements.triangles.at(columnIndex).node2,boundaryElements.collocationPoints.at(columnIndex));
        VectorTriangle triangle2(boundaryElements.triangles.at(columnIndex).node2,boundaryElements.triangles.at(columnIndex).node3,boundaryElements.collocationPoints.at(columnIndex));
        VectorTriangle triangle3(boundaryElements.triangles.at(columnIndex).node3,boundaryElements.triangles.at(columnIndex).node1,boundaryElements.collocationPoints.at(columnIndex));
        triangle1.normal=boundaryElements.triangles.at(columnIndex).normal;
        triangle2.normal=triangle1.normal;
        triangle3.normal=triangle1.normal;

        std::complex<double> Ltest1,Ltest2,Ltest3;
        std::complex<double> Nktest1,Nktest2,Nktest3;

        Eigen::Vector3d triangleCSide=(triangle1.node3-triangle1.node1);
        Eigen::Vector3d triangleBSide=(triangle1.node3-triangle1.node2);
        Eigen::Vector3d triangleASide=(triangle1.node2-triangle1.node1);
        double bSquared=triangleBSide.dot(triangleBSide);
        double cSquared=triangleCSide.dot(triangleCSide);
        double aSquared=triangleASide.dot(triangleASide);
        double b=std::sqrt(bSquared);
        double c=std::sqrt(cSquared);
//        double a=std::sqrt(aSquared);
        if(c<b)
        {
            double tmp=c;
            c=b;
            b=tmp;
        }
        double A=std::acos((bSquared+cSquared-aSquared)/(2.0*b*c));
        double B=std::atan((b*std::sin(A))/(c-b*std::cos(A)));

        double lower = 0.5*(A-0);
        double upper = 0.5*(A+0) ;
        for(int i=0; i<numberOfXWeightsAndAbszissas; i++)
        {
            double xVal1 = xAbszissa[i]*lower+upper;
            double xVal2 = -xAbszissa[i]*lower+upper;
            double r1 = c*sin(B)/(sin(xVal1+B));
            double r2 = c*sin(B)/(sin(xVal2+B));
            std::complex<double> ikr1 = iWavenumber*r1;
            std::complex<double> ikr2 = iWavenumber*r2;
            Ltest1+=xWeights[i]*(std::exp(ikr1)+std::exp(ikr2));
            Nktest1+=xWeights[i]*(std::exp(ikr1)/r1+std::exp(ikr2)/r2);
        }
        Ltest1 = Ltest1*lower;
        Nktest1 = Nktest1*lower;

        triangleCSide = (triangle2.node3-triangle2.node1);
        triangleBSide = (triangle2.node3-triangle2.node2);
        triangleASide = (triangle2.node2-triangle2.node1);
        bSquared = triangleBSide.dot(triangleBSide);
        cSquared = triangleCSide.dot(triangleCSide);
        aSquared = triangleASide.dot(triangleASide);
        b = std::sqrt(bSquared);
        c = std::sqrt(cSquared);
        if(c<b)
        {
            double tmp = c;
            c = b;
            b = tmp;
        }
        double AA = A;
        A = std::acos((bSquared+cSquared-aSquared)/(2.0*b*c));
        B = std::atan((b*std::sin(A))/(c-b*std::cos(A)));

        lower  =  0.5*(A);
        upper  =  0.5*(A+AA+AA);
        for(int i = 0; i<numberOfXWeightsAndAbszissas; i++)
        {
            double xVal1 = xAbszissa[i]*lower+upper;
            double xVal2 = -xAbszissa[i]*lower+upper;
            double r1 = c*sin(B)/(sin(xVal1-AA+B));
            double r2 = c*sin(B)/(sin(xVal2-AA+B));
            std::complex<double> ikr1 = iWavenumber*r1;
            std::complex<double> ikr2 = iWavenumber*r2;
            Ltest2+=xWeights[i]*(std::exp(ikr1)+std::exp(ikr2));
            Nktest2+=xWeights[i]*(std::exp(ikr1)/r1+std::exp(ikr2)/r2);
        }
        Ltest2 = Ltest2*lower;
        Nktest2 = Nktest2*lower;

        triangleCSide = (triangle3.node3-triangle3.node1);
        triangleBSide = (triangle3.node3-triangle3.node2);
        triangleASide = (triangle3.node2-triangle3.node1);
        bSquared = triangleBSide.dot(triangleBSide);
        cSquared = triangleCSide.dot(triangleCSide);
        aSquared = triangleASide.dot(triangleASide);
        b = std::sqrt(bSquared);
        c = std::sqrt(cSquared);
        if(c<b)
        {
            double tmp = c;
            c = b;
            b = tmp;
        }
        double AAA = AA+A;
        A = std::acos((bSquared+cSquared-aSquared)/(2.0*b*c));
        B = std::atan((b*std::sin(A))/(c-b*std::cos(A)));

        lower  =  0.5*(A);
        upper  =  0.5*(A+AAA+AAA);
        for(int i = 0; i<numberOfXWeightsAndAbszissas; i++)
        {
            double xVal1 = xAbszissa[i]*lower+upper;
            double xVal2 = -xAbszissa[i]*lower+upper;
            double r1 = c*sin(B)/(sin(xVal1-AAA+B));
            double r2 = c*sin(B)/(sin(xVal2-AAA+B));
            std::complex<double> ikr1 = iWavenumber*r1;
            std::complex<double> ikr2 = iWavenumber*r2;
            Ltest3+=xWeights[i]*(std::exp(ikr1)+std::exp(ikr2));
            Nktest3+=xWeights[i]*(std::exp(ikr1)/r1+std::exp(ikr2)/r2);
        }

        Ltest3 = Ltest3*lower;
        Nktest3 = Nktest3*lower;
        std::complex<double> Lk = (Ltest1+Ltest2+Ltest3)/PI2;
        std::complex<double> Nk = (Nktest1+Nktest2+Nktest3)/PI4;

        Lk = (imaginaryUnit/(2.0*wavenumber))*(1.0-Lk);
        Nk = (iWavenumber/2.0-Nk);

        if(substituteDPhiWithPhi(columnIndex))
        {
            return - (Lk + couplingParameterHalf) / beta(columnIndex);
        }
        else
        {
            return (couplingParameter * Nk - 0.5) / alpha(columnIndex);
        }
    }
    else
    {
        std::complex<double> Lk = 0.0;
        std::complex<double> Mk = 0.0;
        std::complex<double> Mtk = 0.0;
        std::complex<double> Nk = 0.0;

        Eigen::Vector3d nq = boundaryElements.triangles.at(columnIndex).normal;
        Eigen::Vector3d np = boundaryElements.triangles.at(rowIndex).normal;
        double upNq = np.dot(nq);

        for(int i = 0; i < numberOfRegularWeightsAndAbscissa; i++) //gauss quadrature of standard triangle
        {
            double quadratureWeight = regularWeightsAndAbscissa(i,2);
            double quadratureAbscissaNode1 = regularWeightsAndAbscissa(i,0);
            double quadratureAbscissaNode2 = regularWeightsAndAbscissa(i,1);
            double quadratureAbscissaNode3 = (1.0-quadratureAbscissaNode2-quadratureAbscissaNode1);
            Eigen::Vector3d q = (quadratureAbscissaNode1*boundaryElements.triangles.at(columnIndex).node1) + (quadratureAbscissaNode2*boundaryElements.triangles.at(columnIndex).node2) + (quadratureAbscissaNode3*boundaryElements.triangles.at(columnIndex).node3);
            Eigen::Vector3d rVector = boundaryElements.collocationPoints.at(rowIndex) - q; //p-q
            double rr = rVector.squaredNorm();
            double r = std::sqrt(rr);
            double rnq = - (rVector.dot(nq))/r;
            double rnp = (rVector.dot(np))/r;
            std::complex<double> ikr = iWavenumber*r;
            std::complex<double> greensFunction = (std::exp(ikr))/(PI4*r);
            Lk += quadratureWeight*greensFunction;
            std::complex<double> dGreensDR = (greensFunction/r)*(ikr-1.0);
            Mtk += quadratureWeight*dGreensDR*rnp;
            if(!substituteDPhiWithPhi(columnIndex))
            {
                Mk += quadratureWeight*dGreensDR*rnq;
                std::complex<double> d2GreensDR2 = (greensFunction/rr)*(2.0-2.0*ikr-wavenumberSquared*rr);
                std::complex<double> d2Gkupnq = dGreensDR*(-(upNq+rnp*rnq)/r)+d2GreensDR2*rnp*rnq;
                Nk += quadratureWeight*d2Gkupnq;
            }
        }
        if(substituteDPhiWithPhi(columnIndex))
        {
            return - (Lk + couplingParameter * Mtk ) * boundaryElements.trianglesArea(columnIndex) / beta(columnIndex);
        }
        else
        {
            return (Mk + couplingParameter * Nk) * boundaryElements.trianglesArea(columnIndex) / alpha(columnIndex);
        }
    }
}

std::complex<double> BoundaryElementSolver::implicitDPhiMatrix(long rowIndex, long columnIndex)
{
    if(rowIndex == columnIndex)
    {
        static double xAbszissa[16] = {0.0483076656877383162348126,0.1444719615827964934851864,0.2392873622521370745446032,0.3318686022821276497799168,0.4213512761306353453641194,0.5068999089322293900237475,0.5877157572407623290407455,0.6630442669302152009751152,0.7321821187402896803874267,0.7944837959679424069630973,0.8493676137325699701336930,0.8963211557660521239653072,0.9349060759377396891709191,0.9647622555875064307738119,0.9856115115452683354001750,0.9972638618494815635449811};
        static double xWeights[16] = {0.0965400885147278005667648,0.0956387200792748594190820,0.0938443990808045656391802,0.0911738786957638847128686,0.0876520930044038111427715,0.0833119242269467552221991,0.0781938957870703064717409,0.0723457941088485062253994,0.0658222227763618468376501,0.0586840934785355471452836,0.0509980592623761761961632,0.0428358980222266806568786,0.0342738629130214331026877,0.0253920653092620594557526,0.0162743947309056706051706,0.0070186100094700966004071};
        int numberOfXWeightsAndAbszissas = 16;
        VectorTriangle triangle1(boundaryElements.triangles.at(columnIndex).node1,boundaryElements.triangles.at(columnIndex).node2,boundaryElements.collocationPoints.at(columnIndex));
        VectorTriangle triangle2(boundaryElements.triangles.at(columnIndex).node2,boundaryElements.triangles.at(columnIndex).node3,boundaryElements.collocationPoints.at(columnIndex));
        VectorTriangle triangle3(boundaryElements.triangles.at(columnIndex).node3,boundaryElements.triangles.at(columnIndex).node1,boundaryElements.collocationPoints.at(columnIndex));
        triangle1.normal=boundaryElements.triangles.at(columnIndex).normal;
        triangle2.normal=triangle1.normal;
        triangle3.normal=triangle1.normal;

        std::complex<double> Ltest1,Ltest2,Ltest3;
        std::complex<double> Nktest1,Nktest2,Nktest3;

        Eigen::Vector3d triangleCSide=(triangle1.node3-triangle1.node1);
        Eigen::Vector3d triangleBSide=(triangle1.node3-triangle1.node2);
        Eigen::Vector3d triangleASide=(triangle1.node2-triangle1.node1);
        double bSquared=triangleBSide.dot(triangleBSide);
        double cSquared=triangleCSide.dot(triangleCSide);
        double aSquared=triangleASide.dot(triangleASide);
        double b=std::sqrt(bSquared);
        double c=std::sqrt(cSquared);
//        double a=std::sqrt(aSquared);
        if(c<b)
        {
            double tmp=c;
            c=b;
            b=tmp;
        }
        double A=std::acos((bSquared+cSquared-aSquared)/(2.0*b*c));
        double B=std::atan((b*std::sin(A))/(c-b*std::cos(A)));

        double lower = 0.5*(A-0);
        double upper = 0.5*(A+0) ;
        for(int i=0; i<numberOfXWeightsAndAbszissas; i++)
        {
            double xVal1 = xAbszissa[i]*lower+upper;
            double xVal2 = -xAbszissa[i]*lower+upper;
            double r1 = c*sin(B)/(sin(xVal1+B));
            double r2 = c*sin(B)/(sin(xVal2+B));
            std::complex<double> ikr1 = iWavenumber*r1;
            std::complex<double> ikr2 = iWavenumber*r2;
            Ltest1+=xWeights[i]*(std::exp(ikr1)+std::exp(ikr2));
            Nktest1+=xWeights[i]*(std::exp(ikr1)/r1+std::exp(ikr2)/r2);
        }
        Ltest1 = Ltest1*lower;
        Nktest1 = Nktest1*lower;

        triangleCSide = (triangle2.node3-triangle2.node1);
        triangleBSide = (triangle2.node3-triangle2.node2);
        triangleASide = (triangle2.node2-triangle2.node1);
        bSquared = triangleBSide.dot(triangleBSide);
        cSquared = triangleCSide.dot(triangleCSide);
        aSquared = triangleASide.dot(triangleASide);
        b = std::sqrt(bSquared);
        c = std::sqrt(cSquared);
        if(c<b)
        {
            double tmp = c;
            c = b;
            b = tmp;
        }
        double AA = A;
        A = std::acos((bSquared+cSquared-aSquared)/(2.0*b*c));
        B = std::atan((b*std::sin(A))/(c-b*std::cos(A)));

        lower  =  0.5*(A);
        upper  =  0.5*(A+AA+AA);
        for(int i = 0; i<numberOfXWeightsAndAbszissas; i++)
        {
            double xVal1 = xAbszissa[i]*lower+upper;
            double xVal2 = -xAbszissa[i]*lower+upper;
            double r1 = c*sin(B)/(sin(xVal1-AA+B));
            double r2 = c*sin(B)/(sin(xVal2-AA+B));
            std::complex<double> ikr1 = iWavenumber*r1;
            std::complex<double> ikr2 = iWavenumber*r2;
            Ltest2+=xWeights[i]*(std::exp(ikr1)+std::exp(ikr2));
            Nktest2+=xWeights[i]*(std::exp(ikr1)/r1+std::exp(ikr2)/r2);
        }
        Ltest2 = Ltest2*lower;
        Nktest2 = Nktest2*lower;

        triangleCSide = (triangle3.node3-triangle3.node1);
        triangleBSide = (triangle3.node3-triangle3.node2);
        triangleASide = (triangle3.node2-triangle3.node1);
        bSquared = triangleBSide.dot(triangleBSide);
        cSquared = triangleCSide.dot(triangleCSide);
        aSquared = triangleASide.dot(triangleASide);
        b = std::sqrt(bSquared);
        c = std::sqrt(cSquared);
        if(c<b)
        {
            double tmp = c;
            c = b;
            b = tmp;
        }
        double AAA = AA+A;
        A = std::acos((bSquared+cSquared-aSquared)/(2.0*b*c));
        B = std::atan((b*std::sin(A))/(c-b*std::cos(A)));

        lower  =  0.5*(A);
        upper  =  0.5*(A+AAA+AAA);
        for(int i = 0; i<numberOfXWeightsAndAbszissas; i++)
        {
            double xVal1 = xAbszissa[i]*lower+upper;
            double xVal2 = -xAbszissa[i]*lower+upper;
            double r1 = c*sin(B)/(sin(xVal1-AAA+B));
            double r2 = c*sin(B)/(sin(xVal2-AAA+B));
            std::complex<double> ikr1 = iWavenumber*r1;
            std::complex<double> ikr2 = iWavenumber*r2;
            Ltest3+=xWeights[i]*(std::exp(ikr1)+std::exp(ikr2));
            Nktest3+=xWeights[i]*(std::exp(ikr1)/r1+std::exp(ikr2)/r2);
        }

        Ltest3 = Ltest3*lower;
        Nktest3 = Nktest3*lower;
        std::complex<double> Lk = (Ltest1+Ltest2+Ltest3)/PI2;
        std::complex<double> Nk = (Nktest1+Nktest2+Nktest3)/PI4;

        Lk = (imaginaryUnit/(2.0*wavenumber))*(1.0-Lk);
        Nk = (iWavenumber/2.0-Nk);

        if(substituteDPhiWithPhi(columnIndex))
        {
            return -alpha(columnIndex) / beta(columnIndex) * (Lk + couplingParameterHalf) - (couplingParameter * Nk - 0.5);
        }
        else
        {
            return Lk + couplingParameterHalf + beta(columnIndex)/alpha(columnIndex) * (couplingParameter * Nk - 0.5);
        }
    }
    else
    {
        std::complex<double> Lk = 0.0;
        std::complex<double> Mk = 0.0;
        std::complex<double> Mtk = 0.0;
        std::complex<double> Nk = 0.0;

        Eigen::Vector3d np = boundaryElements.triangles.at(rowIndex).normal;
        Eigen::Vector3d nq = boundaryElements.triangles.at(columnIndex).normal;
        double upNq = np.dot(nq);

        for(int i = 0; i < numberOfRegularWeightsAndAbscissa; i++) //gauss quadrature of standard triangle
        {
            double quadratureWeight = regularWeightsAndAbscissa(i,2);
            double quadratureAbscissaNode1 = regularWeightsAndAbscissa(i,0);
            double quadratureAbscissaNode2 = regularWeightsAndAbscissa(i,1);
            double quadratureAbscissaNode3 = (1.0-quadratureAbscissaNode2-quadratureAbscissaNode1);
            Eigen::Vector3d q = (quadratureAbscissaNode1*boundaryElements.triangles.at(columnIndex).node1) + (quadratureAbscissaNode2*boundaryElements.triangles.at(columnIndex).node2) + (quadratureAbscissaNode3*boundaryElements.triangles.at(columnIndex).node3);
            Eigen::Vector3d rVector = boundaryElements.collocationPoints.at(rowIndex) - q; //p-q
            double rr = rVector.squaredNorm();
            double r = std::sqrt(rr);
            double rnq = - (rVector.dot(nq))/r;
            double rnp = (rVector.dot(np))/r;
            std::complex<double> ikr = iWavenumber*r;
            std::complex<double> greensFunction = (std::exp(ikr))/(PI4*r);
            Lk += quadratureWeight*greensFunction;
            std::complex<double> dGreensDR = (greensFunction/r)*(ikr-1.0);
            Mk += quadratureWeight*dGreensDR*rnq;
            Mtk += quadratureWeight*dGreensDR*rnp;
            std::complex<double> d2GreensDR2 = (greensFunction/rr)*(2.0-2.0*ikr-wavenumberSquared*rr);
            std::complex<double> d2Gkupnq = dGreensDR*(-(upNq+rnp*rnq)/r)+d2GreensDR2*rnp*rnq;
            Nk += quadratureWeight*d2Gkupnq;
        }
        if(substituteDPhiWithPhi(columnIndex))
        {
            return (-alpha(columnIndex) / beta(columnIndex) * (Lk + couplingParameter * Mtk ) - (Mk + couplingParameter * Nk)) * boundaryElements.trianglesArea(columnIndex);
        }
        else
        {
            return (Lk + couplingParameter * Mtk + beta(columnIndex)/alpha(columnIndex) * (Mk + couplingParameter * Nk)) * boundaryElements.trianglesArea(columnIndex);
        }
    }
}

std::complex<double> BoundaryElementSolver::implicitPhiMatrix2(long rowIndex, long columnIndex)
{
    std::complex<double> Lk;
    std::complex<double> Mk;
    std::complex<double> Mtk;
    std::complex<double> Nk;

    BemOperatorsConst(rowIndex,columnIndex,Lk,Mk,Mtk,Nk);
    double onPanel = (rowIndex == columnIndex);
    if(substituteDPhiWithPhi(columnIndex))
    {
        return - (Lk + couplingParameter * Mtk + ((double)onPanel) * couplingParameterHalf) / beta(columnIndex);
    }
    else
    {
        return (Mk + couplingParameter * Nk - ((double)onPanel) * 0.5) / alpha(columnIndex);
    }
}

std::complex<double> BoundaryElementSolver::implicitDPhiMatrix2(long rowIndex, long columnIndex)
{
    std::complex<double> Lk;
    std::complex<double> Mk;
    std::complex<double> Mtk;
    std::complex<double> Nk;

    BemOperatorsConst(rowIndex,columnIndex,Lk,Mk,Mtk,Nk);
    double onPanel = (rowIndex == columnIndex);
    if(substituteDPhiWithPhi(columnIndex))
    {
        return - alpha(columnIndex) / beta(columnIndex) * (Lk + couplingParameter * Mtk + ((double)onPanel) * couplingParameterHalf) - (Mk + couplingParameter * Nk - ((double)onPanel) * 0.5);
    }
    else
    {
        return Lk + couplingParameter * Mtk + ((double)onPanel) * couplingParameterHalf + beta(columnIndex)/alpha(columnIndex) * (Mk + couplingParameter * Nk - ((double)onPanel) * 0.5);
    }
}

std::complex<double> BoundaryElementSolver::implicitPhiMatrixPlusReflected(long rowIndex, long columnIndex, const BoundaryElements &boundaryElements, const BoundaryElements &reflectedElements)
{
    return implicitPhiMatrix2(rowIndex, columnIndex) + implicitPhiMatrixReflected(rowIndex, columnIndex, boundaryElements, reflectedElements);
}

std::complex<double> BoundaryElementSolver::implicitDPhiMatrixPlusReflected(long rowIndex, long columnIndex, const BoundaryElements &boundaryElements, const BoundaryElements &reflectedElements)
{
    return implicitDPhiMatrix2(rowIndex, columnIndex) + implicitDPhiMatrixReflected(rowIndex, columnIndex, boundaryElements, reflectedElements);
}

QVector<std::pair<long,long>> BoundaryElementSolver::getNormalFilteredPivotIndices(BlockCluster* block, double filterAngle)
{
    long rowStartIndex = block->rowStartIndex();
    long columnStartIndex = block->colStartIndex();
    long blockRows = block->rows(); //assumes contiguous ascending indexes
    long blockColumns = block->cols();
    long maxContainerIndex = 20;
    double filterValue = std::cos(filterAngle); // filter cone

    QVector<long> rowIndices = global::createContiguousIndexVector(rowStartIndex, rowStartIndex + blockRows -1);

    QVector<QVector<long>> rowIndicesFilterContainer(maxContainerIndex);
    rowIndicesFilterContainer[0] = rowIndices;
    long numberOfUnfilteredIndices = rowIndices.length();
    long rowContainerIndex = 0;
    long ignoreValue = -1;
    while(numberOfUnfilteredIndices > 1 && rowContainerIndex < maxContainerIndex - 1)
    {
        rowIndicesFilterContainer[rowContainerIndex + 1] = QVector<long>(0);
        long ignoreCounter = 0;
        for(long i = 1; i < rowIndicesFilterContainer.at(rowContainerIndex).length(); i++)
        {
            Eigen::Vector3d referenceNormal = boundaryElements.triangles.at(rowIndicesFilterContainer.at(rowContainerIndex).at(0)).normal;
            Eigen::Vector3d normal = boundaryElements.triangles.at(rowIndicesFilterContainer.at(rowContainerIndex).at(i)).normal;
            if(std::abs(referenceNormal.dot(normal)) < filterValue)
            {
                rowIndicesFilterContainer[rowContainerIndex + 1].append(rowIndicesFilterContainer.at(rowContainerIndex).at(i));
                rowIndicesFilterContainer[rowContainerIndex][i] = ignoreValue;
                ignoreCounter++;
            }
        }
        if(ignoreCounter > 0)
        {
            QVector<long> tmp;
            tmp.reserve(ignoreCounter);
            for(long i = 0; i < rowIndicesFilterContainer.at(rowContainerIndex).length(); i++)
            {
                if(rowIndicesFilterContainer.at(rowContainerIndex).at(i) != ignoreValue)
                {
                    tmp.append(rowIndicesFilterContainer.at(rowContainerIndex).at(i));
                }
            }
            rowIndicesFilterContainer[rowContainerIndex] = tmp;
        }
        numberOfUnfilteredIndices = rowIndicesFilterContainer.at(rowContainerIndex + 1).length();
        if(numberOfUnfilteredIndices >= 1)
        {
            rowContainerIndex++;
        }
    }
    rowIndices.clear();

    QVector<long> columnIndices;
    columnIndices = global::createContiguousIndexVector(columnStartIndex, columnStartIndex + blockColumns -1);

    QVector<QVector<long>> columnIndicesFilterContainer(maxContainerIndex);
    columnIndicesFilterContainer[0] = columnIndices;
    numberOfUnfilteredIndices = columnIndices.length();
    long columnContainerIndex = 0;
    ignoreValue = -1;
    while(numberOfUnfilteredIndices > 1 && columnContainerIndex < maxContainerIndex - 1)
    {
        columnIndicesFilterContainer[columnContainerIndex + 1] = QVector<long>(0);
        long ignoreCounter = 0;
        for(long i = 1; i < columnIndicesFilterContainer.at(columnContainerIndex).length(); i++)
        {
            Eigen::Vector3d referenceNormal = boundaryElements.triangles.at(columnIndicesFilterContainer.at(columnContainerIndex).at(0)).normal;
            Eigen::Vector3d normal = boundaryElements.triangles.at(columnIndicesFilterContainer.at(columnContainerIndex).at(i)).normal;
            if(std::abs(referenceNormal.dot(normal)) < filterValue)
            {
                columnIndicesFilterContainer[columnContainerIndex + 1].append(columnIndicesFilterContainer.at(columnContainerIndex).at(i));
                columnIndicesFilterContainer[columnContainerIndex][i] = ignoreValue;
                ignoreCounter++;
            }
        }
        if(ignoreCounter > 0)
        {
            QVector<long> tmp;
            tmp.reserve(ignoreCounter);
            for(long i = 0; i < columnIndicesFilterContainer.at(columnContainerIndex).length(); i++)
            {
                if(columnIndicesFilterContainer.at(columnContainerIndex).at(i) != ignoreValue)
                {
                    tmp.append(columnIndicesFilterContainer.at(columnContainerIndex).at(i));
                }
            }
            columnIndicesFilterContainer[columnContainerIndex] = tmp;
        }
        numberOfUnfilteredIndices = columnIndicesFilterContainer.at(columnContainerIndex + 1).length();
        if(numberOfUnfilteredIndices >= 1)
        {
            columnContainerIndex++;
        }
    }
    columnIndices.clear();

    long numberOfRowContainers = rowContainerIndex + 1;
    long numberOfColumnContainers = columnContainerIndex + 1;

    QVector<std::pair<long,long>> pivotIndices(numberOfRowContainers * numberOfColumnContainers);
    for(long i = 0; i < numberOfRowContainers ; i++)
    {
        for(long j = 0; j < numberOfColumnContainers ; j++)
        {
            long indexForRowIndex = randomGenerator.bounded((qint32) rowIndicesFilterContainer.at(i).size());
//            rowIndices[i * numberOfColumnContainers + j] = rowIndicesFilterContainer.at(i).at(indexForRowIndex);
            long indexForColumnIndex = randomGenerator.bounded((qint32) columnIndicesFilterContainer.at(j).size());
//            columnIndices[i * numberOfColumnContainers + j] = columnIndicesFilterContainer.at(j).at(indexForColumnIndex);
            pivotIndices[i * numberOfColumnContainers + j] = {rowIndicesFilterContainer.at(i).at(indexForRowIndex), columnIndicesFilterContainer.at(j).at(indexForColumnIndex)};
        }
    }
    rowIndicesFilterContainer.clear();
    columnIndicesFilterContainer.clear();

    for(long i = 0; i < pivotIndices.length(); i++) /// make absolute indices relative
    {
        pivotIndices[i].first -= rowStartIndex;
        pivotIndices[i].second -= columnStartIndex;
    }
    return pivotIndices;
}

void BoundaryElementSolver::calcHBlockRowVector(Eigen::RowVectorXcd &rowVector, const long rowIndex, const long columnStartIndex, std::function<std::complex<double> (long, long)> implicitMatrix)
{
    #pragma omp parallel for // parallelizes vector initialization
    for(long i = 0; i < rowVector.size(); i++)   //assemble rowVectors
    {
        rowVector(i) = implicitMatrix(rowIndex, columnStartIndex + i);
    }
}

void BoundaryElementSolver::calcHBlockColumnVector(Eigen::VectorXcd &columnVector, const long rowStartIndex, const long columnIndex, std::function<std::complex<double> (long, long)> implicitMatrix)
{
    #pragma omp parallel for // parallelizes vector initialization
    for(long i = 0; i < columnVector.size(); i++) //assemble columnVectors
    {
        columnVector(i) = implicitMatrix(rowStartIndex + i, columnIndex);
    }
}

void BoundaryElementSolver::calculateFieldSolution()
{
    std::cout << "Calculating field solution at " << frequency << " Hz." << std::endl;
    wavenumber = (2.0 * frequency * global::PI) / waveSpeed;
    wavenumberSquared = wavenumber * wavenumber;
    wavenumberSquaredHalf = wavenumberSquared / 2.0;
    iWavenumber = wavenumber * imaginaryUnit;

    setUpQuadratureRules();
    int numberOfFieldPoints=0;
    for(int i=0;i<observationFields.size();i++)
    {
        observationFields[i].calculateTrianglesMidpoints();
        numberOfFieldPoints += observationFields.at(i).triangles.size();
    }
    std::cout<<numberOfFieldPoints<<" field points."<<std::endl;
    QString message=QString("Number of field points: "+QString::number(numberOfFieldPoints));
    logStrings::logString.append(message+"\r\n");
    emit updateLog();
    Timer timer;
    timer.start();

    if(hMatFieldSolving) // solve the field via H-matrix technique
    {
        calculateFieldSolutionFast(fieldACARelativeError, fieldACAMaxRank);
    }
    else
    {
        calculateFieldSolutionRegular();
    }

    timer.stop();
//    std::cerr << "Convert to sound pressure." << std::endl;
    for(int i=0; i<observationFields.size(); i++)
    {
        if(observationFields.at(i).soundPressure.size() != observationFields.at(i).phiSolution.size())
        {
            observationFields[i].soundPressure = Eigen::VectorXcd(observationFields.at(i).phiSolution.size());
        }
        if(frequency !=0)
        {
            observationFields[i].soundPressure=(imaginaryUnit*airDensity*PI2 * frequency) * observationFields.at(i).phiSolution;
        }
        else
        {
            observationFields[i].soundPressure=(imaginaryUnit*airDensity*PI2) * observationFields.at(i).phiSolution;
        }
    }
    std::cout << "Runtime of field calculation: " << timer.secs() << " seconds" << std::endl;
    message=QString("Runtime of field calculation: "+QString::number(timer.secs())+" seconds.");
    logStrings::logString.append(message+"\r\n");
    emit updateLog();
}

void BoundaryElementSolver::calculateFieldSolutionRegular()
{
    std::cout<<"Calculating field."<<std::endl;
    QString message = QString("Calculating field.");
    logStrings::logString.append(message+"\r\n");
    emit updateLog();
    std::complex<double> Mk = 0.0;
    std::complex<double> Lk = 0.0;
    for(int i = 0; i < observationFields.size(); i++)
    {
        ObservationField currentObsField = observationFields.at(i);
        long numberOfFieldElements = currentObsField.triangles.size();
        long numberOfBoundaryElements = boundaryElements.triangles.length();

        observationFields[i].phiSolution = Eigen::VectorXcd::Zero(numberOfFieldElements);
//        observationFields[i].dPhiSolution = Eigen::VectorXcd (numberOfFieldElements);
        #pragma omp parallel for private(Mk,Lk)// parallelizes field calculation
        for(long fIndex=0; fIndex<numberOfFieldElements;fIndex++) //field point index
        {
            Eigen::Vector3d observationPoint = currentObsField.triangleMidPoints.at(fIndex);
            std::complex<double> phi = 0.0;
            std::complex<double> sourceObservationValue = 0.0;
            for(int sourceIndex = 0; sourceIndex < pointSources.size(); sourceIndex++)
            {
                PointSource currentPointSource = pointSources.at(sourceIndex);
                sourceObservationValue = sourceObservationValue + sourcePhiTerm(currentPointSource,observationPoint);
            }
            phi = sourceObservationValue;
            for(long bIndex = 0; bIndex < numberOfBoundaryElements; bIndex++) // boundary collocation points index
            {
                BemOperatorField(observationPoint,bIndex,Mk,Lk);
                Lk = Lk*boundaryElements.dPhiSolution(bIndex);
                Mk = Mk*boundaryElements.phiSolution(bIndex);
                phi += Mk-Lk;
            }
            observationFields[i].phiSolution[fIndex] = phi;
        }
    }
}

void BoundaryElementSolver::BemOperatorField(const Eigen::Vector3d observationPoint,const int boundaryTriangleIndex,std::complex<double>& Mk, std::complex<double>& Lk)
{
    Mk = 0.0;
    Lk = 0.0;
    Eigen::Vector3d nq;
    double quadratureWeight;
    double quadratureAbscissaNode1;
    double quadratureAbscissaNode2;
    double quadratureAbscissaNode3;
    Eigen::Vector3d q;
    Eigen::Vector3d rVector; //p-q
    double r;
    double rnq;
    std::complex<double> ikr;
    std::complex<double> greensFunction;
    std::complex<double> dGreensDR;
    for(int i=0; i<numberOfLowOrderWeightsAndAbscissa; i++) //gauss quadrature of standard triangle
    {
        nq = boundaryElements.triangles.at(boundaryTriangleIndex).normal;
        quadratureWeight = lowOrderweightsAndAbscissa(i,2);
        quadratureAbscissaNode1 = lowOrderweightsAndAbscissa(i,0);
        quadratureAbscissaNode2 = lowOrderweightsAndAbscissa(i,1);
        quadratureAbscissaNode3 = (1-quadratureAbscissaNode2-quadratureAbscissaNode1);
        q = (quadratureAbscissaNode1*boundaryElements.triangles.at(boundaryTriangleIndex).node1) + (quadratureAbscissaNode2*boundaryElements.triangles.at(boundaryTriangleIndex).node2) + (quadratureAbscissaNode3*boundaryElements.triangles.at(boundaryTriangleIndex).node3);
        rVector = observationPoint-q; //p-q
        r = rVector.norm();
        rnq = -1.0*(rVector.dot(nq))/r;
        ikr = iWavenumber*r;
        greensFunction = (std::exp(ikr))/(PI4*r);
        Lk += quadratureWeight*greensFunction;
        dGreensDR=(greensFunction/r)*(ikr-1.0);
        Mk += quadratureWeight*dGreensDR*rnq;
    }
    Lk *= boundaryElements.trianglesArea(boundaryTriangleIndex);
    Mk *= boundaryElements.trianglesArea(boundaryTriangleIndex);
}

//std::tuple<std::complex<double>,std::complex<double>,std::complex<double>,std::complex<double>> BoundaryElementSolver::bemOperatorField(const long rowIndex, const long columnIndex, const QVector<Eigen::Vector3d> &observationPoints, const BoundaryElements &boundaryElements)
//{
//    std::complex<double> Mk = 0.0;
//    std::complex<double> Lk = 0.0;
//    std::complex<double> Mtk = 0.0;
//    std::complex<double> Nk = 0.0;

//    Eigen::Vector3d nq;
//    double quadratureWeight;
//    double quadratureAbscissaNode1;
//    double quadratureAbscissaNode2;
//    double quadratureAbscissaNode3;
//    Eigen::Vector3d q;
//    Eigen::Vector3d rVector; //p-q
//    double r;
//    double rnq;
//    std::complex<double> ikr;
//    std::complex<double> greensFunction;
//    std::complex<double> dGreensDR;
//    for(int i=0; i<numberOfLowOrderWeightsAndAbscissa; i++) //gauss quadrature of standard triangle
//    {
//        nq=boundaryElements.triangles.at(columnIndex).normal;
//        quadratureWeight=lowOrderweightsAndAbscissa(i,2);
//        quadratureAbscissaNode1=lowOrderweightsAndAbscissa(i,0);
//        quadratureAbscissaNode2=lowOrderweightsAndAbscissa(i,1);
//        quadratureAbscissaNode3=(1-quadratureAbscissaNode2-quadratureAbscissaNode1);
//        q=(quadratureAbscissaNode1*boundaryElements.triangles.at(columnIndex).node1) + (quadratureAbscissaNode2*boundaryElements.triangles.at(columnIndex).node2) + (quadratureAbscissaNode3*boundaryElements.triangles.at(columnIndex).node3);
//        rVector=observationPoints.at(rowIndex)-q; //p-q
//        r=rVector.norm();
//        rnq=-1.0*(rVector.dot(nq))/r;
//        ikr=iWavenumber*r;
//        greensFunction=(std::exp(ikr))/(PI4*r);
//        Lk=Lk+quadratureWeight*greensFunction;
//        dGreensDR=(greensFunction/r)*(ikr-1.0);
//        Mk=Mk+quadratureWeight*dGreensDR*rnq;
//    }
//    Lk=Lk*boundaryElements.trianglesArea.at(columnIndex);
//    Mk=Mk*boundaryElements.trianglesArea.at(columnIndex);

//    return std::tie(Mk, (std::complex<double>) 0.0, Lk, (std::complex<double>) 0.0);
//}

void BoundaryElementSolver::calculateFieldSolutionFast(double relativeError, long maxRank)
{
//    std::cerr<<"Fast field calculation."<<std::endl;
    QString message=QString("Fast field calculation.");
    logStrings::logString.append(message+"\r\n");
    emit updateLog();

    // give each triangle its solution value; otherwise the association would be lost during reordering
    global::setTrianglesPhi(boundaryElements.triangles, boundaryElements.phiSolution);
    global::setTrianglesDPhi(boundaryElements.triangles, boundaryElements.dPhiSolution);
    ClusterTree boundaryClusterTree(&boundaryElements);

    // the creation of the clustertree has reordered the triangles -> extract the boundary solution vectors from the triangle vectors
    boundaryElements.phiSolution = global::getTrianglesPhi(boundaryElements.triangles);
    boundaryElements.dPhiSolution = global::getTrianglesDPhi(boundaryElements.triangles);

//    std::complex<double> storedCouplingParameter = couplingParameter;
    for(int i = 0; i < observationFields.size(); i++)
    {
        ObservationField* currentObsField = &observationFields[i];
        long numberOfFieldElements = currentObsField->triangles.length();
        std::cout<< "currentObsField->triangles.length(): " << numberOfFieldElements << std::endl;
        currentObsField->phiSolution = Eigen::VectorXcd::Zero(numberOfFieldElements);

        ClusterTree fieldClusterTree(&(currentObsField->triangles));
        currentObsField->calculateTrianglesMidpoints();

        HMatrix blockClusterTree;
        blockClusterTree.setClusterTrees(&fieldClusterTree, &boundaryClusterTree);
        blockClusterTree.populateBlockClusterTree();
        QVector<BlockCluster*> partition = blockClusterTree.getMinPartition();

        // calculate source terms
        #pragma omp parallel for
        for(long fIndex=0; fIndex<numberOfFieldElements; fIndex++) //field point index
        {
            Eigen::Vector3d observationPoint = currentObsField->triangleMidPoints.at(fIndex);
            std::complex<double> sourceObservationValue = 0.0;
            for(int sourceIndex = 0; sourceIndex<pointSources.size(); sourceIndex++)
            {
                PointSource currentPointSource = pointSources.at(sourceIndex);
                sourceObservationValue += sourcePhiTerm(currentPointSource,observationPoint);
            }
//            fieldPhiSolution(fIndex) = sourceObservationValue;
            currentObsField->phiSolution(fIndex) = sourceObservationValue;
        }

        // calculate boundary terms
//        #pragma omp parallel for
        #pragma omp parallel master
        for(long i = 0; i < partition.length(); i++)  // we dont assemble the matrices; we discard assembled blocks directly after the block-vector product
        {
            #pragma omp task
            {
                BlockCluster block = *partition.at(i);

                long rowStartIndex = block.rowStartIndex();
                long numberOfRows = block.rows();
                long columnStartIndex = block.colStartIndex();
                long numberOfColumns = block.cols();

                Eigen::VectorXcd tmp;
                if(block.isAdmissible) // use ACA for low-rank block assembly
                {
                    partialPivotACA(&block, maxRank, relativeError, std::bind(&BoundaryElementSolver::implicitPhiMatrixField, this, std::placeholders::_1, std::placeholders::_2, currentObsField->triangleMidPoints, boundaryElements));
                    tmp.noalias() = block.UMat * (block.singularValues.asDiagonal() * (block.VAdjMat * boundaryElements.phiSolution.segment(columnStartIndex, numberOfColumns)));

                    partialPivotACAextra(&block, maxRank, relativeError, std::bind(&BoundaryElementSolver::getNormalFilteredPivotIndicesForField, this, std::placeholders::_1, std::placeholders::_2), std::bind(&BoundaryElementSolver::implicitDPhiMatrixField, this, std::placeholders::_1, std::placeholders::_2, currentObsField->triangleMidPoints, boundaryElements));
                    tmp.noalias() -= block.UMat * (block.singularValues.asDiagonal() * (block.VAdjMat * boundaryElements.dPhiSolution.segment(columnStartIndex, numberOfColumns)));
                }
                else // use full-rank block assembly
                {
                    Eigen::MatrixXcd phiBlock(numberOfRows,numberOfColumns);
                    Eigen::MatrixXcd dPhiBlock(numberOfRows,numberOfColumns);
                    Eigen::VectorXcd tmpColumn(numberOfRows);
                    for(long i = 0; i < numberOfColumns; i++)
                    {
                        calcHBlockColumnVector(tmpColumn, rowStartIndex, columnStartIndex + i, std::bind(&BoundaryElementSolver::implicitPhiMatrixField, this, std::placeholders::_1, std::placeholders::_2, currentObsField->triangleMidPoints, boundaryElements));
                        phiBlock.col(i) = tmpColumn;
                        calcHBlockColumnVector(tmpColumn, rowStartIndex, columnStartIndex + i, std::bind(&BoundaryElementSolver::implicitDPhiMatrixField, this, std::placeholders::_1, std::placeholders::_2, currentObsField->triangleMidPoints, boundaryElements));
                        dPhiBlock.col(i) = tmpColumn;
                    }
                    tmp.noalias() = phiBlock * boundaryElements.phiSolution.segment(columnStartIndex, numberOfColumns) - dPhiBlock * boundaryElements.dPhiSolution.segment(columnStartIndex, numberOfColumns);
                }
                #pragma omp critical
                {
                    currentObsField->phiSolution.segment(rowStartIndex, numberOfRows) += tmp;
                }
            }
        }
        // calculate boundary terms for the reflected geometry
        if(boundaryElements.impedancePlanes.length() == 1) // calculate field for a single impedance plane // we dont assemble the matrices; we discard assembled blocks directly after the block-vector product
        {
            BoundaryElements reflectedBoundaryElements = boundaryElements;
            reflectedBoundaryElements.reflectGeometry(boundaryElements.impedancePlanes.at(0));

            ClusterTree reflectedClusterTree(boundaryClusterTree);  //copy the original clustertree
            reflectedClusterTree.updateMinCuboids(&reflectedBoundaryElements);

            HMatrix halfSpaceMat;
            halfSpaceMat.setClusterTrees(&fieldClusterTree, &reflectedClusterTree);
            halfSpaceMat.populateBlockClusterTree();

            QVector<BlockCluster*> partition = halfSpaceMat.getMinPartition();
            std::cout<< "halfSpaceMat partition.size(): " << partition.size() << std::endl;
//            #pragma omp parallel for
            #pragma omp parallel master
            for(long i = 0; i < partition.length(); i++)
            {
                #pragma omp task
                {
                    BlockCluster block = *partition.at(i);

                    long rowStartIndex = block.rowStartIndex();
                    long numberOfRows = block.rows();
                    long columnStartIndex = block.colStartIndex();
                    long numberOfColumns = block.cols();

                    Eigen::VectorXcd tmp;
                    if(block.isAdmissible) // use ACA for low-rank block assembly
                    {
                        partialPivotACA(&block, maxRank, relativeError, std::bind(&BoundaryElementSolver::implicitPhiMatrixField, this, std::placeholders::_1, std::placeholders::_2, currentObsField->triangleMidPoints, reflectedBoundaryElements));
                        tmp.noalias() = block.UMat * (block.singularValues.asDiagonal() * (block.VAdjMat * boundaryElements.phiSolution.segment(columnStartIndex, numberOfColumns)));

                        partialPivotACAextra(&block, maxRank, relativeError, std::bind(&BoundaryElementSolver::getNormalFilteredPivotIndicesForField, this, std::placeholders::_1, std::placeholders::_2), std::bind(&BoundaryElementSolver::implicitDPhiMatrixField, this, std::placeholders::_1, std::placeholders::_2, currentObsField->triangleMidPoints, reflectedBoundaryElements));
                        tmp.noalias() -= block.UMat * (block.singularValues.asDiagonal() * (block.VAdjMat * boundaryElements.dPhiSolution.segment(columnStartIndex, numberOfColumns)));
                    }
                    else // use full-rank block assembly
                    {
                        Eigen::MatrixXcd phiBlock(numberOfRows,numberOfColumns);
                        Eigen::MatrixXcd dPhiBlock(numberOfRows,numberOfColumns);

                        Eigen::VectorXcd tmpColumn(numberOfRows);

                        for(long i = 0; i < numberOfColumns; i++)
                        {
                            calcHBlockColumnVector(tmpColumn, rowStartIndex, columnStartIndex + i, std::bind(&BoundaryElementSolver::implicitPhiMatrixField, this, std::placeholders::_1, std::placeholders::_2, currentObsField->triangleMidPoints, reflectedBoundaryElements));
                            phiBlock.col(i) = tmpColumn;
                            calcHBlockColumnVector(tmpColumn, rowStartIndex, columnStartIndex + i, std::bind(&BoundaryElementSolver::implicitDPhiMatrixField, this, std::placeholders::_1, std::placeholders::_2, currentObsField->triangleMidPoints, reflectedBoundaryElements));
                            dPhiBlock.col(i) = tmpColumn;
                        }

                        tmp.noalias() = phiBlock * boundaryElements.phiSolution.segment(columnStartIndex, numberOfColumns) - dPhiBlock * boundaryElements.dPhiSolution.segment(columnStartIndex, numberOfColumns);
                    }
                    #pragma omp critical
                    {
                        currentObsField->phiSolution.segment(rowStartIndex, numberOfRows).noalias() += tmp;
                    }
                }
            }
            halfSpaceMat.clear(true);
            reflectedClusterTree.clear();
        }
        else if(boundaryElements.impedancePlanes.length() >= 2) // calculate field for two parallel impedance planes. Here we actually assemble the reflection matrices.
        {
            HMatrix halfSpaceReflMatPhi;
            HMatrix halfSpaceReflMatDPhi;
            calcReflectionMatricesField(halfSpaceReflMatPhi, halfSpaceReflMatDPhi, maxRank, relativeError, currentObsField->triangleMidPoints, boundaryElements, fieldClusterTree, boundaryClusterTree);
            HArithm::MVM(currentObsField->phiSolution, halfSpaceReflMatPhi, boundaryElements.phiSolution);
            HArithm::MVM(currentObsField->phiSolution, halfSpaceReflMatDPhi, -boundaryElements.dPhiSolution);
            halfSpaceReflMatPhi.clear();
            halfSpaceReflMatDPhi.clear(true);
        }
//        currentObsField->phiSolution = fieldPhiSolution;
        fieldClusterTree.clear();
        blockClusterTree.clear();
    }
    boundaryClusterTree.clear();
}

void BoundaryElementSolver::hBlockAssemblyField(BlockCluster* phiBlock, BlockCluster* dPhiBlock, const long maxRank, const double relativeError, const QVector<Eigen::Vector3d> &observationPoints, const BoundaryElements &boundaryElements)
{
    if(phiBlock -> isAdmissible /*&& 0*/) // Adaptive Cross Approximation
    {
        partialPivotACA(phiBlock, maxRank, relativeError, std::bind(&BoundaryElementSolver::implicitPhiMatrixField, this, std::placeholders::_1, std::placeholders::_2, observationPoints, boundaryElements));
        partialPivotACAextra(dPhiBlock, maxRank, relativeError, std::bind(&BoundaryElementSolver::getNormalFilteredPivotIndicesForField, this, std::placeholders::_1, std::placeholders::_2), std::bind(&BoundaryElementSolver::implicitDPhiMatrixField, this, std::placeholders::_1, std::placeholders::_2, observationPoints, boundaryElements));
    }
    else // nearfield block full matrix assembly
    {
        long rowStartIndex = phiBlock->rowStartIndex();
        long columnStartIndex = phiBlock->colStartIndex();
        long blockRows = phiBlock->rows();
        long blockColumns = phiBlock->cols();

        phiBlock->fullMat = Eigen::MatrixXcd(blockRows, blockColumns);
        dPhiBlock->fullMat = Eigen::MatrixXcd(blockRows, blockColumns);

        Eigen::VectorXcd tmpColumn(blockRows);

        #pragma omp parallel for private(tmpColumn)
        for(long columnIndex = 0; columnIndex < blockColumns; columnIndex++)
        {
            calcHBlockColumnVector(tmpColumn, rowStartIndex, columnStartIndex + columnIndex, std::bind(&BoundaryElementSolver::implicitPhiMatrixField, this, std::placeholders::_1, std::placeholders::_2, observationPoints, boundaryElements));
            phiBlock->fullMat.col(columnIndex) = tmpColumn;
            calcHBlockColumnVector(tmpColumn, rowStartIndex, columnStartIndex + columnIndex, std::bind(&BoundaryElementSolver::implicitDPhiMatrixField, this, std::placeholders::_1, std::placeholders::_2, observationPoints, boundaryElements));
            dPhiBlock->fullMat.col(columnIndex) = tmpColumn;
        }
    }
}

void BoundaryElementSolver::calcReflectionMatricesField(HMatrix &reflMatPhi, HMatrix &reflMatDPhi, const long maxRank, const double relativeError, const QVector<Eigen::Vector3d> &observationPoints, const BoundaryElements &elements, ClusterTree &obsClusterTree, ClusterTree &elementsClusterTree, int lastPlaneIndex)
{
    for(int planeIndex=0; planeIndex<elements.impedancePlanes.length() && planeIndex < 2; planeIndex++)
    {
        if(lastPlaneIndex == planeIndex)
        {
            continue;
        }
        BoundaryElements reflectedBoundaryElements = elements;
        reflectedBoundaryElements.reflectGeometry(elements.impedancePlanes.at(planeIndex));
//        reflectedBoundaryElements.calculateTriangleMidPoints();


        ClusterTree reflectedClusterTree(elementsClusterTree);  //copy the original clustertree
        reflectedClusterTree.updateMinCuboids(&reflectedBoundaryElements);

        HMatrix newReflMatPhi;
        newReflMatPhi.setClusterTrees(&obsClusterTree, &reflectedClusterTree);
        newReflMatPhi.populateBlockClusterTree();

        HMatrix newReflMatDPhi;
        newReflMatDPhi.setClusterTrees(&obsClusterTree, &reflectedClusterTree);
        newReflMatDPhi.populateBlockClusterTree();

        QVector<BlockCluster*> newReflMatPhiPartition = newReflMatPhi.getMinPartition();
        QVector<BlockCluster*> newReflMatDPhiPartition = newReflMatDPhi.getMinPartition();

//        std::cout<< "Started the block assembly of the reflection Matrix" << std::endl;
        Timer timer;
        timer.start();
//        std::cout<< "newReflMatDPhiPartition.size(): " << newReflMatDPhiPartition.size() << std::endl;
//        std::cout<< "reflectedDPhiHMatPartition.size(): " << reflectedDPhiHMatPartition.size() << std::endl;

        #pragma omp parallel for
        for(long blockIndex = 0; blockIndex < newReflMatPhiPartition.size(); blockIndex++)
        {
            hBlockAssemblyField(newReflMatPhiPartition.at(blockIndex), newReflMatDPhiPartition.at(blockIndex), maxRank, relativeError, observationPoints, reflectedBoundaryElements);
        }
//        std::cout<< "Finished the block assembly of the reflection Matrix" << std::endl;
        timer.stop();
        std::cout << "runtime of hBlockAssemblyReflectedMatrix: " << timer.secs() << std::endl;
        timer.reset();

        HArithm::compressHMat(newReflMatPhi, maxRank, 0.1 * relativeError);
        HArithm::compressHMat(newReflMatDPhi, maxRank, 0.1 * relativeError);

        double newPhiNorm = newReflMatPhi.norm();
        double newDPhiNorm = newReflMatDPhi.norm();
        std::cout << "newPhiNorm: " << newPhiNorm << std::endl;

        if(lastPlaneIndex == -1 && planeIndex == 0)
        {
            reflMatPhi = newReflMatPhi;
            reflMatDPhi = newReflMatDPhi;

//            calcReflectionMatrices(reflMatPhi, reflMatDPhi, maxRank, relativeError, elements, reflectedBoundaryElements, clusterTree, planeIndex);
            calcReflectionMatricesField(reflMatPhi, reflMatDPhi, maxRank, relativeError, observationPoints, reflectedBoundaryElements, obsClusterTree, elementsClusterTree, planeIndex);
//            reflectedClusterTree.clear(); // can't be deleted here -> still a memory leak
            continue;
        }
        HArithm::recursiveHMatAddition(* reflMatPhi.getRootBlock(), * newReflMatPhi.getRootBlock(), maxRank, 0.1 * relativeError);
        HArithm::recursiveHMatAddition(* reflMatDPhi.getRootBlock(), * newReflMatDPhi.getRootBlock(), maxRank, 0.1 * relativeError);
//        HArithm::recursiveHMatSubstraction(reflMatPhi.getRootBlock(), newReflMatPhi.getRootBlock(), maxRank, 0.1 * relativeError);
//        HArithm::recursiveHMatSubstraction(reflMatDPhi.getRootBlock(), newReflMatDPhi.getRootBlock(), maxRank, 0.1 * relativeError);

        newReflMatDPhi.clear(); // newReflMatDPhi in not needed anymore
        newReflMatPhi.clear(); // newReflMatPhi in not needed anymore
        reflectedClusterTree.clear();

        double reflMatPhiNorm = reflMatPhi.norm();
        double reflMatDPhiNorm = reflMatDPhi.norm();
        std::cout << "reflMatPhiNorm: " << reflMatPhiNorm << std::endl;

        std::cout << "newPhiNorm / reflMatPhiNorm: " << newPhiNorm / reflMatPhiNorm << std::endl;
        std::cout << "newDPhiNorm / reflMatDPhiNorm: " << newDPhiNorm / reflMatDPhiNorm << std::endl;

        if(newPhiNorm / reflMatPhiNorm < relativeError && newDPhiNorm / reflMatDPhiNorm < relativeError)
        {
            continue;
        }
        else
        {
//            calcReflectionMatrices(reflMatPhi, reflMatDPhi, maxRank, relativeError, elements, reflectedBoundaryElements, clusterTree, planeIndex);
            calcReflectionMatricesField(reflMatPhi, reflMatDPhi, maxRank, relativeError, observationPoints, reflectedBoundaryElements, obsClusterTree, elementsClusterTree, planeIndex);
        }
    }
    if(lastPlaneIndex == -1)
    {
        global::trimMemory();
    }
}

std::complex<double> BoundaryElementSolver::implicitPhiMatrixField(const long rowIndex, const long columnIndex, const QVector<Eigen::Vector3d> &observationPoints, const BoundaryElements &boundaryElements)
{
    Eigen::Vector3d observationPoint = observationPoints.at(rowIndex);
    std::complex<double> Mk = 0.0;
//        std::complex<double> Lk =0.0;
    Eigen::Vector3d nq;
    double quadratureWeight;
    double quadratureAbscissaNode1;
    double quadratureAbscissaNode2;
    double quadratureAbscissaNode3;
    Eigen::Vector3d q;
    Eigen::Vector3d rVector; //p-q
    double r;
    double rnq;
    std::complex<double> ikr;
    std::complex<double> greensFunction;
    std::complex<double> dGreensDR;
    for(int i = 0; i < numberOfLowOrderWeightsAndAbscissa; i++) //gauss quadrature of standard triangle
    {
        nq = boundaryElements.triangles.at(columnIndex).normal;
        quadratureWeight = lowOrderweightsAndAbscissa(i,2);
        quadratureAbscissaNode1 = lowOrderweightsAndAbscissa(i,0);
        quadratureAbscissaNode2 = lowOrderweightsAndAbscissa(i,1);
        quadratureAbscissaNode3 = (1-quadratureAbscissaNode2-quadratureAbscissaNode1);
        q = (quadratureAbscissaNode1*boundaryElements.triangles.at(columnIndex).node1) + (quadratureAbscissaNode2*boundaryElements.triangles.at(columnIndex).node2) + (quadratureAbscissaNode3*boundaryElements.triangles.at(columnIndex).node3);
        rVector = observationPoint-q; //p-q
        r = rVector.norm();
        rnq = -(rVector.dot(nq))/r;
        ikr = iWavenumber*r;
        greensFunction = (std::exp(ikr))/(PI4*r);
//            Lk += quadratureWeight*greensFunction;
        dGreensDR = (greensFunction/r)*(ikr-1.0);
        Mk += quadratureWeight*dGreensDR*rnq;
    }
    return  Mk * boundaryElements.trianglesArea(columnIndex);
}

std::complex<double> BoundaryElementSolver::implicitDPhiMatrixField(const long rowIndex, const long columnIndex, const QVector<Eigen::Vector3d> &observationPoints, const BoundaryElements &boundaryElements)
{
    Eigen::Vector3d observationPoint = observationPoints.at(rowIndex);
//        std::complex<double> Mk = 0.0;
    std::complex<double> Lk =0.0;
    Eigen::Vector3d nq;
    double quadratureWeight;
    double quadratureAbscissaNode1;
    double quadratureAbscissaNode2;
    double quadratureAbscissaNode3;
    Eigen::Vector3d q;
    Eigen::Vector3d rVector; //p-q
    double r;
//        double rnq;
    std::complex<double> ikr;
    std::complex<double> greensFunction;
//        std::complex<double> dGreensDR;
    for(int i = 0; i < numberOfLowOrderWeightsAndAbscissa; i++) //gauss quadrature of standard triangle
    {
        nq=boundaryElements.triangles.at(columnIndex).normal;
        quadratureWeight = lowOrderweightsAndAbscissa(i,2);
        quadratureAbscissaNode1 = lowOrderweightsAndAbscissa(i,0);
        quadratureAbscissaNode2 = lowOrderweightsAndAbscissa(i,1);
        quadratureAbscissaNode3 = (1-quadratureAbscissaNode2-quadratureAbscissaNode1);
        q = (quadratureAbscissaNode1*boundaryElements.triangles.at(columnIndex).node1) + (quadratureAbscissaNode2*boundaryElements.triangles.at(columnIndex).node2) + (quadratureAbscissaNode3*boundaryElements.triangles.at(columnIndex).node3);
        rVector = observationPoint-q; //p-q
        r = rVector.norm();
//            rnq = -1.0*(rVector.dot(nq))/r;
        ikr = iWavenumber*r;
        greensFunction = (std::exp(ikr))/(PI4*r);
        Lk += quadratureWeight*greensFunction;
//            dGreensDR=(greensFunction/r)*(ikr-1.0);
//            Mk += quadratureWeight*dGreensDR*rnq;
    }
    return Lk * boundaryElements.trianglesArea(columnIndex);
}

QVector<std::pair<long,long>> BoundaryElementSolver::getNormalFilteredPivotIndicesForField(BlockCluster* block, double filterAngle)
{
//    long rowStartIndex = block->rowStartIndex();
    long columnStartIndex = block->colStartIndex();
    long blockRows = block->rows(); //assumes contiguous ascending indexes
    long blockColumns = block->cols();
    long maxContainerIndex = 20;
    double filterValue = std::cos(filterAngle); // filter cone

    QVector<long> columnIndices;
    columnIndices = global::createContiguousIndexVector(columnStartIndex, columnStartIndex + blockColumns -1);

    QVector<QVector<long>> columnIndicesFilterContainer(maxContainerIndex);
    columnIndicesFilterContainer[0] = columnIndices;
    long numberOfUnfilteredIndices = columnIndices.length();
    long columnContainerIndex = 0;
    long ignoreValue = -1;
    while(numberOfUnfilteredIndices > 1 && columnContainerIndex < maxContainerIndex - 1)
    {
//            std::cerr << "columnIndicesFilterContainer[columnContainerIndex]: ";
//            global::printLongVector(columnIndicesFilterContainer[columnContainerIndex]);
//            std::cerr << std::endl;

        columnIndicesFilterContainer[columnContainerIndex + 1] = QVector<long>(0);
        long ignoreCounter = 0;
        for(long i = 1; i < columnIndicesFilterContainer.at(columnContainerIndex).length(); i++)
        {
            Eigen::Vector3d referenceNormal = boundaryElements.triangles.at(columnIndicesFilterContainer.at(columnContainerIndex).at(0)).normal;
            Eigen::Vector3d normal = boundaryElements.triangles.at(columnIndicesFilterContainer.at(columnContainerIndex).at(i)).normal;
            if(std::abs(referenceNormal.dot(normal)) < filterValue)
            {
//                    std::cerr << "column std::abs(referenceNormal.dot(normal)) < filterValue" << std::endl;
                columnIndicesFilterContainer[columnContainerIndex + 1].append(columnIndicesFilterContainer.at(columnContainerIndex).at(i));
                columnIndicesFilterContainer[columnContainerIndex][i] = ignoreValue;
                ignoreCounter++;
            }
        }
        if(ignoreCounter > 0)
        {
            QVector<long> tmp;
            tmp.reserve(ignoreCounter);
            for(long i = 0; i < columnIndicesFilterContainer.at(columnContainerIndex).length(); i++)
            {
                if(columnIndicesFilterContainer.at(columnContainerIndex).at(i) != ignoreValue)
                {
                    tmp.append(columnIndicesFilterContainer.at(columnContainerIndex).at(i));
                }
            }
            columnIndicesFilterContainer[columnContainerIndex] = tmp;
        }
        numberOfUnfilteredIndices = columnIndicesFilterContainer.at(columnContainerIndex + 1).length();
        if(numberOfUnfilteredIndices >= 1)
        {
//                std::cerr << "numberOfUnfilteredIndices >= 1" << std::endl;
            columnContainerIndex++;
        }
    }
    columnIndices.clear();
    long numberOfColumnContainers = columnContainerIndex + 1;

    QVector<long> rowIndices(numberOfColumnContainers); // = global::createContiguousIndexVector(rowStartIndex, rowStartIndex + blockRows -1);
    for(long i = 0; i < numberOfColumnContainers ; i++)
    {
        rowIndices[i] = /*rowStartIndex +*/ randomGenerator.bounded((qint32) blockRows);
    }

//    long numberOfRowContainers = rowContainerIndex + 1;

    QVector<std::pair<long,long>> pivotIndices(rowIndices.size() * numberOfColumnContainers);
    for(long i = 0; i < rowIndices.size() ; i++)
    {
        for(long j = 0; j < numberOfColumnContainers ; j++)
        {
            long indexForColumnIndex = randomGenerator.bounded((qint32) columnIndicesFilterContainer.at(j).size());
            pivotIndices[i * numberOfColumnContainers + j] = {rowIndices.at(i), columnIndicesFilterContainer.at(j).at(indexForColumnIndex)};
        }
    }
    rowIndices.clear();
    columnIndicesFilterContainer.clear();

    for(long i = 0; i < pivotIndices.length(); i++) /// make absolute indices relative
    {
//        pivotIndices[i].first -= rowStartIndex;
        pivotIndices[i].second -= columnStartIndex;
    }
    return pivotIndices;
}

//void BoundaryElementSolver::BemOperatorsLinear(const long collocationPointIndex, const long domainIndex, Eigen::Vector3d shapeFuncNodeWeights, std::complex<double>& Lk, std::complex<double>& Mk, std::complex<double>& Mtk, std::complex<double>& Nk, const LinearBoundaryElements &elements, const bool singularity)
//{
//    if(!singularity)
//    {
//        Lk = 0.0;
//        Mk = 0.0;
//        Mtk = 0.0;
//        Nk = 0.0;

//        LinearTriangle domain = elements.triangles.at(domainIndex);

//        Eigen::Vector3d node1 = elements.nodes.at(domain.node1Index).coordinates;
//        Eigen::Vector3d node2 = elements.nodes.at(domain.node2Index).coordinates;
//        Eigen::Vector3d node3 = elements.nodes.at(domain.node3Index).coordinates;

//        double triangleArea = elements.trianglesArea(domainIndex);

//        Eigen::Vector3d p = elements.collocationPoints.at(collocationPointIndex);
//        Eigen::Vector3d np = elements.triangles.at(collocationPointIndex).normal; // assume here that the collocation point is on the triangle with collocationPointIndex
//        Eigen::Vector3d nq = domain.normal;
//        double upNq = np.dot(nq);

//        for(int i = 0; i < numberOfRegularWeightsAndAbscissa; i++) //gauss quadrature of standard triangle
//        {
//            double quadratureWeight = regularWeightsAndAbscissa(i,2);
//            double quadratureAbscissaNode1 = regularWeightsAndAbscissa(i,0);
//            double quadratureAbscissaNode2 = regularWeightsAndAbscissa(i,1);
//            double quadratureAbscissaNode3 = (1.0-quadratureAbscissaNode2-quadratureAbscissaNode1);

//            double shapeFunction = shapeFuncNodeWeights(0) * quadratureAbscissaNode1 + shapeFuncNodeWeights(1) * quadratureAbscissaNode2 +shapeFuncNodeWeights(2) * quadratureAbscissaNode3;

//            Eigen::Vector3d q = (quadratureAbscissaNode1*node1) + (quadratureAbscissaNode2*node2) + (quadratureAbscissaNode3*node3);
//            Eigen::Vector3d rVector = p - q;
//            double rr = rVector.squaredNorm();
//            double r = std::sqrt(rr);
//            double rnq = - (rVector.dot(nq))/r;
//            double rnp = (rVector.dot(np))/r;
//            std::complex<double> ikr = iWavenumber*r;
//            std::complex<double> greensFunction = (std::exp(ikr))/(PI4*r);
//            Lk += quadratureWeight * greensFunction * shapeFunction;
//            std::complex<double> dGreensDR = (greensFunction/r)*(ikr-1.0);
//            Mk += quadratureWeight*dGreensDR*rnq * shapeFunction;
//            Mtk += quadratureWeight*dGreensDR*rnp * shapeFunction;
//            std::complex<double> d2GreensDR2 = (greensFunction/rr)*(2.0-2.0*ikr-wavenumberSquared*rr);
//            std::complex<double> d2Gkupnq = dGreensDR*(-(upNq+rnp*rnq)/r)+d2GreensDR2*rnp*rnq;
//            Nk += quadratureWeight*d2Gkupnq * shapeFunction;
//        }
//        Lk *= triangleArea;
//        Mk *= triangleArea;
//        Mtk *= triangleArea;
//        Nk *= triangleArea;
//    }
//    else
//    {
//        Lk = 0.0;
//        Mk = 0.0;
//        Mtk = 0.0;
//        Nk = 0.0;

//        //observation Point p is on the same Panel as source points(q)

//        static double xAbszissa[16] = {0.0483076656877383162348126,0.1444719615827964934851864,0.2392873622521370745446032,0.3318686022821276497799168,0.4213512761306353453641194,0.5068999089322293900237475,0.5877157572407623290407455,0.6630442669302152009751152,0.7321821187402896803874267,0.7944837959679424069630973,0.8493676137325699701336930,0.8963211557660521239653072,0.9349060759377396891709191,0.9647622555875064307738119,0.9856115115452683354001750,0.9972638618494815635449811};
//        static double xWeights[16] = {0.0965400885147278005667648,0.0956387200792748594190820,0.0938443990808045656391802,0.0911738786957638847128686,0.0876520930044038111427715,0.0833119242269467552221991,0.0781938957870703064717409,0.0723457941088485062253994,0.0658222227763618468376501,0.0586840934785355471452836,0.0509980592623761761961632,0.0428358980222266806568786,0.0342738629130214331026877,0.0253920653092620594557526,0.0162743947309056706051706,0.0070186100094700966004071};
//        int numberOfXWeightsAndAbszissas = 16;


//        LinearTriangle domain = elements.triangles.at(domainIndex);

//        Eigen::Vector3d node1 = elements.nodes.at(domain.node1Index).coordinates;
//        Eigen::Vector3d node2 = elements.nodes.at(domain.node2Index).coordinates;
//        Eigen::Vector3d node3 = elements.nodes.at(domain.node3Index).coordinates;

//        Eigen::Vector3d p = elements.collocationPoints.at(collocationPointIndex);
//        double shapeFunctionOnP = 1/3.0; // assumes that the collocation point p lies in the center of the triangle

//        QVector<VectorTriangle> splitTriangles(3);
//        splitTriangles[0] = VectorTriangle(p, node1, node2);
//        splitTriangles[1] = VectorTriangle(p, node2, node3);
//        splitTriangles[2] = VectorTriangle(p, node3, node1);

//        Eigen::Vector3d np = elements.triangles.at(collocationPointIndex).normal; // assume here that the collocation point is on the triangle with collocationPointIndex
//        Eigen::Vector3d rotationAxis = {np(1), - np(0), 0}; // rotationAxis = np cross z-unit vector
//        rotationAxis.normalize();
//        double angleToZAxis = std::acos(np(2));
////        Eigen::Matrix3Xd rotationMatrix(2,3); // rotation matrix to convert from global coordinates to local the coordinates in the triangle plane
////        Eigen::Transform t(Eigen::AngleAxis(angleToZAxis,rotationAxis));
//        Eigen::Matrix3d t;
//        t = Eigen::AngleAxis(angleToZAxis,rotationAxis);
//        std::cout << "Normal: " << np << std::endl;
//        std::cout << "Transformed Normal: " << t *np << std::endl;

////        for(int i = 0; i < splitTriangles.size(); i++) // quadrature over each triangle
////        {
////            std::complex<double> tmpNk = 0.0;
////            std::complex<double> tmpLk = 0.0;
////            VectorTriangle currentTriangle = splitTriangles.at(i);

////            Eigen::Vector3d triangle13Side = (currentTriangle.node3-currentTriangle.node1);
////            Eigen::Vector3d triangle23Side = (currentTriangle.node3-currentTriangle.node2);
////            Eigen::Vector3d triangle12Side = (currentTriangle.node2-currentTriangle.node1);
////            double side13Squared = triangle13Side.squaredNorm();
////            double side23Squared = triangle23Side.squaredNorm();
////            double side12Squared = triangle12Side.squaredNorm();
////            double side13Length = std::sqrt(side13Squared);
////            double side23Length = std::sqrt(side23Squared);
////            double side12Length = std::sqrt(side12Squared);
////            double angleAtProjectionPoint = std::acos((side13Squared+side12Squared-side23Squared)/(2.0*side13Length*side12Length));
////            double angleAtNode2 = std::acos((side23Squared+side12Squared-side13Squared)/(2.0*side23Length*side12Length));

////            ////////////// for the conversion to global coordinates
////            Eigen::Matrix3Xd rotationMatrix(3,2); // rotation matrix to convert from coordinates in triangle plane to the global coordinates
////            rotationMatrix.col(0) = (triangle12Side.normalized());
////            rotationMatrix.col(1) = (triangle13Side - rotationMatrix.col(0).dot(triangle13Side) * rotationMatrix.col(0)).normalized();

////            for(int i = 0; i < numberOfWeightsAndAbscissaLineQuadrature; i++)
////            {
////                double quadratureAbszissaForAngle = weightsAndAbscissaLineQuadrature(i,0);
////                double quadratureWeightForAngle = weightsAndAbscissaLineQuadrature(i,1);

////                double angle1 = (angleAtProjectionPoint + quadratureAbszissaForAngle * angleAtProjectionPoint)/2.0;
////                double maxRForAngle1 = side12Length*std::sin(angleAtNode2)/(std::sin(angle1 + angleAtNode2));

////                for(int j = 0; j < numberOfWeightsAndAbscissaLineQuadrature; j++)
////                {
////                    double quadratureWeightForS = weightsAndAbscissaLineQuadrature(j,1);

////                    double rs = (maxRForAngle1 + weightsAndAbscissaLineQuadrature(j,0) * maxRForAngle1) / 2.0;

////                    double x = rs * std::cos(angle1);
////                    double y = rs * std::sin(angle1);

////                    Eigen::Vector3d q = currentTriangle.node1 + x * rotationMatrix.col(0) + y * rotationMatrix.col(1);

////                    double factorAndWeights = maxRForAngle1/2.0 * rs * quadratureWeightForAngle * quadratureWeightForS;

////                    Eigen::Vector3d rVector = p - q; //p-q
////                    double rr = rVector.squaredNorm();
////                    double r = std::sqrt(rr);
////                    double rnq = - (rVector.dot(nq))/r;
////                    double rnp = (rVector.dot(np))/r;
////                    std::complex<double> ikr = iWavenumber*r;
////                    std::complex<double> greensFunction = (std::exp(ikr))/(PI4*r);
////                    tmpLk += factorAndWeights*greensFunction;
////                    std::complex<double> dGreensDR = (greensFunction/r)*(ikr-1.0);
////                    std::complex<double> d2GreensDR2 = (greensFunction/rr)*(2.0-2.0*ikr-wavenumberSquared*rr);
////                    std::complex<double> d2Gkupnq = dGreensDR*(-(upNq+rnp*rnq)/r)+d2GreensDR2*rnp*rnq;
////                    tmpNk += factorAndWeights*d2Gkupnq;
////                }
////            }
////            Lk += angleAtProjectionPoint / 2.0 * tmpLk;
////            Mk += angleAtProjectionPoint / 2.0 * tmpMk;
////            Mtk += angleAtProjectionPoint / 2.0 * tmpMtk;
////            Nk += angleAtProjectionPoint / 2.0 * tmpNk;
////        }

////        std::complex<double> LkTri1,LkTri2,LkTri3;
////        std::complex<double> LkNode1Tri1,LkNode1Tri2,LkNode1Tri3;
////        std::complex<double> NkTri1,NkTri2,NkTri3;

////        Eigen::Vector3d triangleCSide = (triangle1.node3-triangle1.node1);
////        Eigen::Vector3d triangleBSide = (triangle1.node3-triangle1.node2);
////        Eigen::Vector3d triangleASide = (triangle1.node2-triangle1.node1);
////        double bSquared = triangleBSide.dot(triangleBSide);
////        double cSquared = triangleCSide.dot(triangleCSide);
////        double aSquared = triangleASide.dot(triangleASide);
////        double b = std::sqrt(bSquared);
////        double c = std::sqrt(cSquared);
////    //        double a=std::sqrt(aSquared);
////        if(c<b)
////        {
////            double tmp=c;
////            c=b;
////            b=tmp;
////        }
////        double A=std::acos((bSquared+cSquared-aSquared)/(2.0*b*c));
////        double B=std::atan((b*std::sin(A))/(c-b*std::cos(A)));


////        double lower = 0.5*(A-0);
////        double upper = 0.5*(A+0) ;
////        for(int i=0; i<numberOfXWeightsAndAbszissas; i++)
////        {
////            double xVal1 = xAbszissa[i]*lower+upper;
////            double xVal2 = -xAbszissa[i]*lower+upper;
////            double r1 = c*sin(B)/(sin(xVal1+B));
////            double r2 = c*sin(B)/(sin(xVal2+B));
////            std::complex<double> ikr1 = iWavenumber*r1;
////            std::complex<double> ikr2 = iWavenumber*r2;
////            LkTri1+=xWeights[i]*(std::exp(ikr1)+std::exp(ikr2));
////            NkTri1+=xWeights[i]*(std::exp(ikr1)/r1+std::exp(ikr2)/r2);
////        }
////        LkTri1 = LkTri1*lower;
////        NkTri1 = NkTri1*lower;

////        triangleCSide = (triangle2.node3-triangle2.node1);
////        triangleBSide = (triangle2.node3-triangle2.node2);
////        triangleASide = (triangle2.node2-triangle2.node1);
////        bSquared = triangleBSide.dot(triangleBSide);
////        cSquared = triangleCSide.dot(triangleCSide);
////        aSquared = triangleASide.dot(triangleASide);
////        b = std::sqrt(bSquared);
////        c = std::sqrt(cSquared);
////        if(c<b)
////        {
////            double tmp = c;
////            c = b;
////            b = tmp;
////        }
////        double AA = A;
////        A = std::acos((bSquared+cSquared-aSquared)/(2.0*b*c));
////        B = std::atan((b*std::sin(A))/(c-b*std::cos(A)));

////        lower  =  0.5*(A);
////        upper  =  0.5*(A+AA+AA);
////        for(int i = 0; i<numberOfXWeightsAndAbszissas; i++)
////        {
////            double xVal1 = xAbszissa[i]*lower+upper;
////            double xVal2 = -xAbszissa[i]*lower+upper;
////            double r1 = c*sin(B)/(sin(xVal1-AA+B));
////            double r2 = c*sin(B)/(sin(xVal2-AA+B));
////            std::complex<double> ikr1 = iWavenumber*r1;
////            std::complex<double> ikr2 = iWavenumber*r2;
////            LkTri2+=xWeights[i]*(std::exp(ikr1)+std::exp(ikr2));
////            NkTri2+=xWeights[i]*(std::exp(ikr1)/r1+std::exp(ikr2)/r2);
////        }
////        LkTri2 = LkTri2*lower;
////        NkTri2 = NkTri2*lower;

////        triangleCSide = (triangle3.node3-triangle3.node1);
////        triangleBSide = (triangle3.node3-triangle3.node2);
////        triangleASide = (triangle3.node2-triangle3.node1);
////        bSquared = triangleBSide.dot(triangleBSide);
////        cSquared = triangleCSide.dot(triangleCSide);
////        aSquared = triangleASide.dot(triangleASide);
////        b = std::sqrt(bSquared);
////        c = std::sqrt(cSquared);
////        if(c<b)
////        {
////            double tmp = c;
////            c = b;
////            b = tmp;
////        }
////        double AAA = AA+A;
////        A = std::acos((bSquared+cSquared-aSquared)/(2.0*b*c));
////        B = std::atan((b*std::sin(A))/(c-b*std::cos(A)));

////        lower  =  0.5*(A);
////        upper  =  0.5*(A+AAA+AAA);
////        for(int i = 0; i<numberOfXWeightsAndAbszissas; i++)
////        {
////            double xVal1 = xAbszissa[i]*lower+upper;
////            double xVal2 = -xAbszissa[i]*lower+upper;
////            double r1 = c*sin(B)/(sin(xVal1-AAA+B));
////            double r2 = c*sin(B)/(sin(xVal2-AAA+B));
////            std::complex<double> ikr1 = iWavenumber*r1;
////            std::complex<double> ikr2 = iWavenumber*r2;
////            LkTri3+=xWeights[i]*(std::exp(ikr1)+std::exp(ikr2));
////            NkTri3+=xWeights[i]*(std::exp(ikr1)/r1+std::exp(ikr2)/r2);
////        }

////        LkTri3 = LkTri3*lower;
////        NkTri3 = NkTri3*lower;
////        std::complex<double> LTest = (LkTri1+LkTri2+LkTri3)/PI2;
////        std::complex<double> NkTest = (NkTri1+NkTri2+NkTri3)/PI4;

////        LTest = (imaginaryUnit/(2.0*wavenumber))*(1.0-LTest);
////        NkTest = (iWavenumber/2.0-NkTest);

////        Lk = LTest;
////        Nk = NkTest;
//    }
//}

//void BoundaryElementSolver::BemOperatorsLinearNode1(const long collocationPointIndex, const long domainIndex, std::complex<double>& Lk, std::complex<double>& Mk, std::complex<double>& Mtk, std::complex<double>& Nk, const LinearBoundaryElements &elements, const bool singularity)
//{

//    if(!singularity)
//    {
//        Lk = 0.0;
//        Mk = 0.0;
//        Mtk = 0.0;
//        Nk = 0.0;

//        LinearTriangle domain = elements.triangles.at(domainIndex);

//        Eigen::Vector3d node1 = elements.nodes.at(domain.node1Index).coordinates;
//        Eigen::Vector3d node2 = elements.nodes.at(domain.node2Index).coordinates;
//        Eigen::Vector3d node3 = elements.nodes.at(domain.node3Index).coordinates;

//        double triangleArea = elements.trianglesArea(domainIndex);

//        Eigen::Vector3d p = elements.collocationPoints.at(collocationPointIndex);
//        Eigen::Vector3d np = elements.triangles.at(collocationPointIndex).normal; // assume here that the collocation point is on the triangle with collocationPointIndex
//        Eigen::Vector3d nq = domain.normal;
//        double upNq = np.dot(nq);

//        for(int i = 0; i < numberOfRegularWeightsAndAbscissa; i++) //gauss quadrature of standard triangle
//        {
//            double quadratureWeight = regularWeightsAndAbscissa(i,2);
//            double quadratureAbscissaNode1 = regularWeightsAndAbscissa(i,0);
//            double quadratureAbscissaNode2 = regularWeightsAndAbscissa(i,1);
//            double quadratureAbscissaNode3 = (1.0-quadratureAbscissaNode2-quadratureAbscissaNode1);
//            double node1ShapeFunction = 1.0 - quadratureAbscissaNode2 - quadratureAbscissaNode3;
//            Eigen::Vector3d q = (quadratureAbscissaNode1*node1) + (quadratureAbscissaNode2*node2) + (quadratureAbscissaNode3*node3);
//            Eigen::Vector3d rVector = p - q;
//            double rr = rVector.squaredNorm();
//            double r = std::sqrt(rr);
//            double rnq = - (rVector.dot(nq))/r;
//            double rnp = (rVector.dot(np))/r;
//            std::complex<double> ikr = iWavenumber*r;
//            std::complex<double> greensFunction = (std::exp(ikr))/(PI4*r);
//            Lk += quadratureWeight * greensFunction * node1ShapeFunction;
//            std::complex<double> dGreensDR = (greensFunction/r)*(ikr-1.0);
//            Mk += quadratureWeight*dGreensDR*rnq * node1ShapeFunction;
//            Mtk += quadratureWeight*dGreensDR*rnp * node1ShapeFunction;
//            std::complex<double> d2GreensDR2 = (greensFunction/rr)*(2.0-2.0*ikr-wavenumberSquared*rr);
//            std::complex<double> d2Gkupnq = dGreensDR*(-(upNq+rnp*rnq)/r)+d2GreensDR2*rnp*rnq;
//            Nk += quadratureWeight*d2Gkupnq * node1ShapeFunction;
//        }
//        Lk *= triangleArea;
//        Mk *= triangleArea;
//        Mtk *= triangleArea;
//        Nk *= triangleArea;
//    }
//    else
//    {
//        Lk = 0.0;
//        Mk = 0.0;
//        Mtk = 0.0;
//        Nk = 0.0;

//        //observation Point p is on the same Panel as source points(q)

//        static double xAbszissa[16] = {0.0483076656877383162348126,0.1444719615827964934851864,0.2392873622521370745446032,0.3318686022821276497799168,0.4213512761306353453641194,0.5068999089322293900237475,0.5877157572407623290407455,0.6630442669302152009751152,0.7321821187402896803874267,0.7944837959679424069630973,0.8493676137325699701336930,0.8963211557660521239653072,0.9349060759377396891709191,0.9647622555875064307738119,0.9856115115452683354001750,0.9972638618494815635449811};
//        static double xWeights[16] = {0.0965400885147278005667648,0.0956387200792748594190820,0.0938443990808045656391802,0.0911738786957638847128686,0.0876520930044038111427715,0.0833119242269467552221991,0.0781938957870703064717409,0.0723457941088485062253994,0.0658222227763618468376501,0.0586840934785355471452836,0.0509980592623761761961632,0.0428358980222266806568786,0.0342738629130214331026877,0.0253920653092620594557526,0.0162743947309056706051706,0.0070186100094700966004071};
//        int numberOfXWeightsAndAbszissas = 16;


//        LinearTriangle domain = elements.triangles.at(domainIndex);

//        Eigen::Vector3d node1 = elements.nodes.at(domain.node1Index).coordinates;
//        Eigen::Vector3d node2 = elements.nodes.at(domain.node2Index).coordinates;
//        Eigen::Vector3d node3 = elements.nodes.at(domain.node3Index).coordinates;

//        Eigen::Vector3d p = elements.collocationPoints.at(collocationPointIndex);
//        double shapeFunctionOnP = 1/3.0; // assumes that the collocation point p lies in the center of the triangle

//        VectorTriangle triangle1(node1, node2, p);
//        VectorTriangle triangle2(node2, node3, p);
//        VectorTriangle triangle3(node3, node1, p);


//        triangle1.normal = elements.triangles.at(collocationPointIndex).normal;
//        triangle2.normal = triangle1.normal;
//        triangle3.normal = triangle1.normal;

//        std::complex<double> Ltest1,Ltest2,Ltest3;
//        std::complex<double> Nktest1,Nktest2,Nktest3;

//        Eigen::Vector3d triangleCSide = (triangle1.node3-triangle1.node1);
//        Eigen::Vector3d triangleBSide = (triangle1.node3-triangle1.node2);
//        Eigen::Vector3d triangleASide = (triangle1.node2-triangle1.node1);
//        double bSquared = triangleBSide.dot(triangleBSide);
//        double cSquared = triangleCSide.dot(triangleCSide);
//        double aSquared = triangleASide.dot(triangleASide);
//        double b = std::sqrt(bSquared);
//        double c = std::sqrt(cSquared);
//    //        double a=std::sqrt(aSquared);
//        if(c<b)
//        {
//            double tmp=c;
//            c=b;
//            b=tmp;
//        }
//        double A=std::acos((bSquared+cSquared-aSquared)/(2.0*b*c));
//        double B=std::atan((b*std::sin(A))/(c-b*std::cos(A)));


//        double lower = 0.5*(A-0);
//        double upper = 0.5*(A+0) ;
//        for(int i=0; i<numberOfXWeightsAndAbszissas; i++)
//        {
//            double xVal1 = xAbszissa[i]*lower+upper;
//            double xVal2 = -xAbszissa[i]*lower+upper;
//            double r1 = c*sin(B)/(sin(xVal1+B));
//            double r2 = c*sin(B)/(sin(xVal2+B));
//            std::complex<double> ikr1 = iWavenumber*r1;
//            std::complex<double> ikr2 = iWavenumber*r2;
//            Ltest1+=xWeights[i]*(std::exp(ikr1)+std::exp(ikr2));
//            Nktest1+=xWeights[i]*(std::exp(ikr1)/r1+std::exp(ikr2)/r2);
//        }
//        Ltest1 = Ltest1*lower;
//        Nktest1 = Nktest1*lower;

//        triangleCSide = (triangle2.node3-triangle2.node1);
//        triangleBSide = (triangle2.node3-triangle2.node2);
//        triangleASide = (triangle2.node2-triangle2.node1);
//        bSquared = triangleBSide.dot(triangleBSide);
//        cSquared = triangleCSide.dot(triangleCSide);
//        aSquared = triangleASide.dot(triangleASide);
//        b = std::sqrt(bSquared);
//        c = std::sqrt(cSquared);
//        if(c<b)
//        {
//            double tmp = c;
//            c = b;
//            b = tmp;
//        }
//        double AA = A;
//        A = std::acos((bSquared+cSquared-aSquared)/(2.0*b*c));
//        B = std::atan((b*std::sin(A))/(c-b*std::cos(A)));

//        lower  =  0.5*(A);
//        upper  =  0.5*(A+AA+AA);
//        for(int i = 0; i<numberOfXWeightsAndAbszissas; i++)
//        {
//            double xVal1 = xAbszissa[i]*lower+upper;
//            double xVal2 = -xAbszissa[i]*lower+upper;
//            double r1 = c*sin(B)/(sin(xVal1-AA+B));
//            double r2 = c*sin(B)/(sin(xVal2-AA+B));
//            std::complex<double> ikr1 = iWavenumber*r1;
//            std::complex<double> ikr2 = iWavenumber*r2;
//            Ltest2+=xWeights[i]*(std::exp(ikr1)+std::exp(ikr2));
//            Nktest2+=xWeights[i]*(std::exp(ikr1)/r1+std::exp(ikr2)/r2);
//        }
//        Ltest2 = Ltest2*lower;
//        Nktest2 = Nktest2*lower;

//        triangleCSide = (triangle3.node3-triangle3.node1);
//        triangleBSide = (triangle3.node3-triangle3.node2);
//        triangleASide = (triangle3.node2-triangle3.node1);
//        bSquared = triangleBSide.dot(triangleBSide);
//        cSquared = triangleCSide.dot(triangleCSide);
//        aSquared = triangleASide.dot(triangleASide);
//        b = std::sqrt(bSquared);
//        c = std::sqrt(cSquared);
//        if(c<b)
//        {
//            double tmp = c;
//            c = b;
//            b = tmp;
//        }
//        double AAA = AA+A;
//        A = std::acos((bSquared+cSquared-aSquared)/(2.0*b*c));
//        B = std::atan((b*std::sin(A))/(c-b*std::cos(A)));

//        lower  =  0.5*(A);
//        upper  =  0.5*(A+AAA+AAA);
//        for(int i = 0; i<numberOfXWeightsAndAbszissas; i++)
//        {
//            double xVal1 = xAbszissa[i]*lower+upper;
//            double xVal2 = -xAbszissa[i]*lower+upper;
//            double r1 = c*sin(B)/(sin(xVal1-AAA+B));
//            double r2 = c*sin(B)/(sin(xVal2-AAA+B));
//            std::complex<double> ikr1 = iWavenumber*r1;
//            std::complex<double> ikr2 = iWavenumber*r2;
//            Ltest3+=xWeights[i]*(std::exp(ikr1)+std::exp(ikr2));
//            Nktest3+=xWeights[i]*(std::exp(ikr1)/r1+std::exp(ikr2)/r2);
//        }

//        Ltest3 = Ltest3*lower;
//        Nktest3 = Nktest3*lower;
//        std::complex<double> LTest = (Ltest1+Ltest2+Ltest3)/PI2;
//        std::complex<double> NkTest = (Nktest1+Nktest2+Nktest3)/PI4;

//        LTest = (imaginaryUnit/(2.0*wavenumber))*(1.0-LTest);
//        NkTest = (iWavenumber/2.0-NkTest);

//        Lk = LTest;
//        Nk = NkTest;
//    }
//}

//void BoundaryElementSolver::BemOperatorsLinearNode2(const long collocationPointIndex, const long domainIndex, std::complex<double>& Lk, std::complex<double>& Mk, std::complex<double>& Mtk, std::complex<double>& Nk, const LinearBoundaryElements &elements, const bool singularity)
//{
//    if(!singularity)
//    {
//        Lk = 0.0;
//        Mk = 0.0;
//        Mtk = 0.0;
//        Nk = 0.0;

//        LinearTriangle domain = elements.triangles.at(domainIndex);

//        Eigen::Vector3d node1 = elements.nodes.at(domain.node1Index).coordinates;
//        Eigen::Vector3d node2 = elements.nodes.at(domain.node2Index).coordinates;
//        Eigen::Vector3d node3 = elements.nodes.at(domain.node3Index).coordinates;

//        double triangleArea = elements.trianglesArea(domainIndex);

//        Eigen::Vector3d p = elements.collocationPoints.at(collocationPointIndex);
//        Eigen::Vector3d np = elements.triangles.at(collocationPointIndex).normal; // assume here that the collocation point is on the triangle with collocationPointIndex
//        Eigen::Vector3d nq = domain.normal;
//        double upNq = np.dot(nq);

//        for(int i = 0; i < numberOfRegularWeightsAndAbscissa; i++) //gauss quadrature of standard triangle
//        {
//            double quadratureWeight = regularWeightsAndAbscissa(i,2);
//            double quadratureAbscissaNode1 = regularWeightsAndAbscissa(i,0);
//            double quadratureAbscissaNode2 = regularWeightsAndAbscissa(i,1);
//            double quadratureAbscissaNode3 = (1.0-quadratureAbscissaNode2-quadratureAbscissaNode1);
//            double node2ShapeFunction = quadratureAbscissaNode2;
//            Eigen::Vector3d q = (quadratureAbscissaNode1*node1) + (quadratureAbscissaNode2*node2) + (quadratureAbscissaNode3*node3);
//            Eigen::Vector3d rVector = p - q;
//            double rr = rVector.squaredNorm();
//            double r = std::sqrt(rr);
//            double rnq = - (rVector.dot(nq))/r;
//            double rnp = (rVector.dot(np))/r;
//            std::complex<double> ikr = iWavenumber*r;
//            std::complex<double> greensFunction = (std::exp(ikr))/(PI4*r);
//            Lk += quadratureWeight * greensFunction * node2ShapeFunction;
//            std::complex<double> dGreensDR = (greensFunction/r)*(ikr-1.0);
//            Mk += quadratureWeight*dGreensDR*rnq * node2ShapeFunction;
//            Mtk += quadratureWeight*dGreensDR*rnp * node2ShapeFunction;
//            std::complex<double> d2GreensDR2 = (greensFunction/rr)*(2.0-2.0*ikr-wavenumberSquared*rr);
//            std::complex<double> d2Gkupnq = dGreensDR*(-(upNq+rnp*rnq)/r)+d2GreensDR2*rnp*rnq;
//            Nk += quadratureWeight*d2Gkupnq * node2ShapeFunction;
//        }
//        Lk *= triangleArea;
//        Mk *= triangleArea;
//        Mtk *= triangleArea;
//        Nk *= triangleArea;
//    }
//}

//void BoundaryElementSolver::BemOperatorsLinearNode3(const long collocationPointIndex, const long domainIndex, std::complex<double>& Lk, std::complex<double>& Mk, std::complex<double>& Mtk, std::complex<double>& Nk, const LinearBoundaryElements &elements, const bool singularity)
//{
//    if(!singularity)
//    {
//        Lk = 0.0;
//        Mk = 0.0;
//        Mtk = 0.0;
//        Nk = 0.0;

//        LinearTriangle domain = elements.triangles.at(domainIndex);

//        Eigen::Vector3d node1 = elements.nodes.at(domain.node1Index).coordinates;
//        Eigen::Vector3d node2 = elements.nodes.at(domain.node2Index).coordinates;
//        Eigen::Vector3d node3 = elements.nodes.at(domain.node3Index).coordinates;

//        double triangleArea = elements.trianglesArea(domainIndex);

//        Eigen::Vector3d p = elements.collocationPoints.at(collocationPointIndex);
//        Eigen::Vector3d np = elements.triangles.at(collocationPointIndex).normal; // assume here that the collocation point is on the triangle with collocationPointIndex
//        Eigen::Vector3d nq = domain.normal;
//        double upNq = np.dot(nq);

//        for(int i = 0; i < numberOfRegularWeightsAndAbscissa; i++) //gauss quadrature of standard triangle
//        {
//            double quadratureWeight = regularWeightsAndAbscissa(i,2);
//            double quadratureAbscissaNode1 = regularWeightsAndAbscissa(i,0);
//            double quadratureAbscissaNode2 = regularWeightsAndAbscissa(i,1);
//            double quadratureAbscissaNode3 = (1.0-quadratureAbscissaNode2-quadratureAbscissaNode1);
//            double node3ShapeFunction = quadratureAbscissaNode3;
//            Eigen::Vector3d q = (quadratureAbscissaNode1*node1) + (quadratureAbscissaNode2*node2) + (quadratureAbscissaNode3*node3);
//            Eigen::Vector3d rVector = p - q;
//            double rr = rVector.squaredNorm();
//            double r = std::sqrt(rr);
//            double rnq = - (rVector.dot(nq))/r;
//            double rnp = (rVector.dot(np))/r;
//            std::complex<double> ikr = iWavenumber*r;
//            std::complex<double> greensFunction = (std::exp(ikr))/(PI4*r);
//            Lk += quadratureWeight * greensFunction * node3ShapeFunction;
//            std::complex<double> dGreensDR = (greensFunction/r)*(ikr-1.0);
//            Mk += quadratureWeight*dGreensDR*rnq * node3ShapeFunction;
//            Mtk += quadratureWeight*dGreensDR*rnp * node3ShapeFunction;
//            std::complex<double> d2GreensDR2 = (greensFunction/rr)*(2.0-2.0*ikr-wavenumberSquared*rr);
//            std::complex<double> d2Gkupnq = dGreensDR*(-(upNq+rnp*rnq)/r)+d2GreensDR2*rnp*rnq;
//            Nk += quadratureWeight*d2Gkupnq * node3ShapeFunction;
//        }
//        Lk *= triangleArea;
//        Mk *= triangleArea;
//        Mtk *= triangleArea;
//        Nk *= triangleArea;
//    }
//}

void BoundaryElementSolver::calculateBoundaryConditionAndSourceTermVector()
{
    long numberOfElements = boundaryElements.triangles.length();
    if(numberOfElements == 0)
    {
        std::cout<<"No elements!"<<std::endl;
        QString message=QString("No elements!");
        logStrings::logString.append(message+"\r\n");
        return;
    }

    sourceTermVector = Eigen::VectorXcd::Zero(numberOfElements);
    alpha = Eigen::VectorXcd::Zero(numberOfElements);
    beta = Eigen::VectorXcd::Zero(numberOfElements);
    f = Eigen::VectorXcd::Zero(numberOfElements);

    PointSource currentPointSource;
    Eigen::Vector3d currentCollocPoint;
    Eigen::Vector3d currentNormal;
    std::complex<double> currentSourceObservationValue;
    #pragma omp parallel for private(currentPointSource,currentCollocPoint,currentNormal,currentSourceObservationValue)// parallelizes source term and bc initialization
    for(long row = 0; row < numberOfElements; row++)
    {
        alpha(row) = boundaryElements.triangles.at(row).robinBoundaryCondition.a;
        beta(row) = boundaryElements.triangles.at(row).robinBoundaryCondition.b;
        f(row) = boundaryElements.triangles.at(row).robinBoundaryCondition.g;
//        alpha(row)=boundaryElements.robinBoundaryConditions.at(row).a;
//        beta(row)=boundaryElements.robinBoundaryConditions.at(row).b;
//        f(row)=boundaryElements.robinBoundaryConditions.at(row).g;
        currentCollocPoint = boundaryElements.collocationPoints.at(row);
        currentNormal = boundaryElements.triangles.at(row).normal;
        currentSourceObservationValue = 0.0;
        for(int sourceIndex = 0; sourceIndex < pointSources.size(); sourceIndex++)
        {
            currentPointSource = pointSources.at(sourceIndex);
            currentSourceObservationValue += sourceTerm(currentPointSource,currentCollocPoint,currentNormal);
        }
        sourceTermVector(row) = currentSourceObservationValue;
    }
    substituteDPhiWithPhi = (alpha.array().abs() < beta.array().abs()).cast<int>(); // vector holds the substitution rules for the matrix assembly

    if(((alpha.array().abs() < global::tiny ) && (beta.array().abs() < global::tiny )).any())
    {
        std::cerr << "Alpha and beta are small in boundary condition!" << std::endl;
    }
}

std::complex<double> BoundaryElementSolver::sourceTerm(const PointSource source, const Eigen::Vector3d listeningPosition, const Eigen::Vector3d normal)
{
    Eigen::Vector3d rVector = listeningPosition - source.position; //p-s
    double r = rVector.norm();
//    double rup= -(rVector.dot(normal))/r; // in exterior AEBEM3 comparison is -normal required
    double rup = (rVector.dot(normal))/r; // in exterior AEBEM3 comparison is -normal required
    std::complex<double> ikr = iWavenumber*r;
    std::complex<double> greensFunction = source.weight*(std::exp(ikr))/(PI4*r); //e^(ikr-iwt)
    std::complex<double> dGreensDR = (greensFunction/r)*(ikr-1.0);
    return greensFunction + couplingParameter * dGreensDR * rup;
}

std::complex<double> BoundaryElementSolver::sourcePhiTerm(const PointSource source,const Eigen::Vector3d listeningPosition)
{
    Eigen::Vector3d rVector = source.position-listeningPosition; //p-q
    double r = rVector.norm();
    std::complex<double> ikr = iWavenumber*r;
    std::complex<double> greensFunction = source.weight*(std::exp(ikr))/(PI4*r);
    return greensFunction;
}

void BoundaryElementSolver::setUpQuadratureRules()
{
    regularWeightsAndAbscissa = triangleQuadratureRules::weightsandAbscissa(orderRegular);
    highOrderweightsAndAbscissa = triangleQuadratureRules::weightsandAbscissa(orderHigh);
    numberOfRegularWeightsAndAbscissa = regularWeightsAndAbscissa.rows();
    numberOfHighOrderWeightsAndAbscissa = highOrderweightsAndAbscissa.rows();

    weightsAndAbscissaLineQuadrature = lineQuadratureRules::weightsandAbscissa(5);
    numberOfWeightsAndAbscissaLineQuadrature = weightsAndAbscissaLineQuadrature.rows();
}

void BoundaryElementSolver::prepareBoundaryElements()
{
    boundaryElements.calculateCollocationPoints();
    boundaryElements.calculateTrianglesAreaAndAverageDimension();
    boundaryElements.calculateTriangleMidPoints();
}

void BoundaryElementSolver::setCouplingParameterNegative()
{
    couplingParameter = -couplingParameter;
    couplingParameterHalf = -couplingParameterHalf;
}
