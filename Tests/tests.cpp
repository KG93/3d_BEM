#include "tests.h"

QTEST_MAIN(Test)

QVector<Eigen::Vector3d> Test::generateRandomPoints(unsigned long n)
{
    QVector<Eigen::Vector3d> points(n);
    for(unsigned long i=0; i<n; i++)
    {
        points[i] = Eigen::Vector3d::Random();
    }
    return points;
}

void Test::testHMatrixTransposition()
{
//    HArithm::compressHMat(hMatrix, 0, 0);
    HMatrix hMatCopy = HMatrix(hMatrix);
    HArithm::transpose(hMatCopy);
    Eigen::MatrixXcd hMatAsFull = hMatrix.toFullMat().transpose();
    Eigen::MatrixXcd hMatCopyAsFull = hMatCopy.toFullMat();
    double relAccuracy = (hMatAsFull - hMatCopyAsFull).norm() / hMatAsFull.norm();
    std::cout << "Relative accuracy of HMatrix transposition: " << relAccuracy << "." << std::endl;
    QVERIFY(relAccuracy < tiny);
}

void Test::testACA()
{
    // test hMat to full matrix conversion and the relative accuracy of the HMatrix approximation
    Eigen::MatrixXcd hMatAsFull = hMatrix.toFullMat();

    Eigen::MatrixXcd test(n,n);
    for(long col = 0; col < n; col++)
    {
        for(long row = 0; row < n; row++)
        {
            test(row, col) = implicitMatrix(row, col, pointsForRows, pointsForCols);
        }
    }
    double relAccuracy = (test - hMatAsFull).norm() / test.norm();
    std::cout << "Relative accuracy of HMatrix approximation: " << relAccuracy << ". Desired accuracy was set to " << relError << std::endl;
    QVERIFY(relAccuracy < 2 * relError); // relaxed accuracy requirements; the stopping criterion for the ACA asseembly is a heuristic
}

void Test::testMatrixVectorProducts()
{
    Eigen::MatrixXcd hMatAsFull = hMatrix.toFullMat();

    // test matrix-vector product
    const Eigen::VectorXcd x = Eigen::VectorXcd::Random(n);
    Eigen::VectorXcd resultFull = hMatAsFull * x;
    Eigen::VectorXcd resultHMat = Eigen::VectorXcd::Zero(n);
    HArithm::MVM(resultHMat, hMatrix, x);
    double relAccuracy =  (resultHMat - resultFull).norm() / resultFull.norm();
    std::cout << "Relative accuracy of MVM: " << relAccuracy << std::endl;
    QVERIFY(relAccuracy < tiny);

    // test recursive matrix-vector product
    resultFull = hMatAsFull * x;
    resultHMat = Eigen::VectorXcd::Zero(n);
    HArithm::recursiveMatrixVectorPoduct(resultHMat, hMatrix.getRootBlock(),x, 0, 0);
    relAccuracy =  (resultHMat - resultFull).norm() / resultFull.norm();
    std::cout << "Relative accuracy of recursiveMatrixVectorPoduct: " << relAccuracy << std::endl;
    QVERIFY(relAccuracy < tiny);

    // test parallel subtractive matrix-vector product
    resultFull = - hMatAsFull * x;
    resultHMat = Eigen::VectorXcd::Zero(n);
    HArithm::subtractiveParallelMatrixVectorPoduct(resultHMat, hMatrix.getRootBlock(), x, 0, 0);
    relAccuracy =  (resultHMat - resultFull).norm() / resultFull.norm();
    std::cout << "Relative accuracy of subtractiveParallelMatrixVectorPoduct: " << relAccuracy << std::endl;
    QVERIFY(relAccuracy < tiny);

    // test parallel recursive matrix-vector product
    resultFull = hMatAsFull * x;
    resultHMat = Eigen::VectorXcd::Zero(n);
    HArithm::parallelMatrixVectorPoduct(resultHMat, hMatrix.getRootBlock(),x, 0, 0);
    relAccuracy =  (resultHMat - resultFull).norm() / resultFull.norm();
    std::cout << "Relative accuracy of parallelMatrixVectorPoduct: " << relAccuracy << std::endl;
    QVERIFY(relAccuracy < tiny);

    // test vector-matrix product
    const Eigen::RowVectorXcd xT = Eigen::RowVectorXcd::Random(n);
    Eigen::RowVectorXcd resultFullT = xT * hMatAsFull;
    Eigen::RowVectorXcd resultHMatT = Eigen::RowVectorXcd::Zero(n);
    HArithm::VMM(resultHMatT, xT, hMatrix);
    relAccuracy =  (resultHMatT - resultFullT).norm() / resultFullT.norm();
    std::cout << "Relative accuracy of VMM: " << relAccuracy << std::endl;
    QVERIFY(relAccuracy < tiny);

    // test recursive vector-matrix product
    resultFullT = xT * hMatAsFull;
    resultHMatT = Eigen::RowVectorXcd::Zero(n);
    HArithm::recursiveVectorMatrixPoduct(resultHMatT, xT, hMatrix.getRootBlock(), 0, 0);
    relAccuracy =  (resultHMatT - resultFullT).norm() / resultFullT.norm();
    std::cout << "Relative accuracy of recursiveVectorMatrixPoduct: " << relAccuracy << std::endl;
    QVERIFY(relAccuracy < tiny);

    // test parallel subtractive vector-matrix product
    resultFullT = -xT * hMatAsFull;
    resultHMatT = Eigen::RowVectorXcd::Zero(n);
    HArithm::subtractiveParallelVectorMatrixPoduct(resultHMatT, xT, hMatrix.getRootBlock(), 0, 0);
    relAccuracy =  (resultHMatT - resultFullT).norm() / resultFullT.norm();
    std::cout << "Relative accuracy of subtractiveParallelVectorMatrixPoduct: " << relAccuracy << std::endl;
    QVERIFY(relAccuracy < tiny);
}

void Test::testMatrixMatrixProduct()
{
    // test matrix-matrix multiplication
    double multiplyRelTolerance = 0.001;
    HMatrix hMatFromLeft = HMatrix(&columnClustertree, &rowClustertree, maxRank, relError, std::bind(&Test::implicitMatrix, this, std::placeholders::_1, std::placeholders::_2, std::ref(pointsForRows), std::ref(pointsForCols)));
    Eigen::MatrixXcd hMatAsFull = hMatrix.toFullMat();
    Eigen::MatrixXcd hMatFromLeftAsFull = hMatFromLeft.toFullMat();
    HMatrix pruductHMat = HMultiply::multiplyHMat(hMatFromLeft, hMatrix, 0, multiplyRelTolerance);
    Eigen::VectorXcd x = Eigen::VectorXcd::Random(n);
    Eigen::VectorXcd resultF = hMatFromLeftAsFull * (hMatAsFull * x);
    Eigen::VectorXcd resultH = Eigen::VectorXcd::Zero(n);
    HArithm::MVM(resultH, pruductHMat, x);
    double relAccuracy = (resultF - resultH).norm() / resultF.norm();
    std::cout << "Relative accuracy of multiplyHMat: " << relAccuracy << ". Rounding accuracy was set to " << multiplyRelTolerance << std::endl;
    pruductHMat.clear();
    hMatFromLeft.clear();
    QVERIFY(relAccuracy < 10 * multiplyRelTolerance + 100 * tiny); // relaxed accuracy requirements due to multiple truncations
}

void Test::testMatrixSubtraction()
{
    // test recursiveHMatSubstraction
    double subtrRelError = 0.01;
    Eigen::MatrixXcd hMatAsFull = hMatrix.toFullMat();
    HMatrix hMatCopy1 = HMatrix(hMatrix);
    HMatrix hMatCopy2 = HMatrix(hMatrix);
    HArithm::multiplyHMatByMinusOne(hMatCopy2);
    HArithm::recursiveHMatSubstraction(hMatCopy1.getRootBlock(), hMatCopy2.getRootBlock(), 0, subtrRelError);

    Eigen::MatrixXcd resultF = hMatCopy1.toFullMat();
    double relAccuracy = (resultF - 2*hMatAsFull).norm() / (2*hMatAsFull.norm());

    std::cout << "Relative accuracy of recursiveHMatSubstraction: " << relAccuracy << ". Rounding accuracy was set to " << subtrRelError << std::endl;
    hMatCopy1.clear();
    hMatCopy2.clear();
    QVERIFY(relAccuracy < 10 * subtrRelError + tiny); // relaxed accuracy requirements due to multiple truncations
}

void Test::testMatValForwardSubstitutions()
{
    // test matrix-valued forward substitution
    double relTruncTol = 0; // relative truncation tolerance
    // hMatSymmetricClusterTrees has the same ClusterTrees for the rows and the columns and is therefore amenable for H-LU-factorization
    HMatrix hMatSym = HMatrix(&rowClustertree, &rowClustertree, maxRank, relError, std::bind(&Test::implicitMatrix, this, std::placeholders::_1, std::placeholders::_2, std::ref(pointsForRows), std::ref(pointsForCols)));
    const Eigen::MatrixXcd hMatSymFull = hMatSym.toFullMat(); // for accuracy comparison

    HMatrix Z = HMatrix(hMatSym); // copy the hMatSym matrix for the forwSubsMatVal operation
    HMatrix copyForLu = HMatrix(hMatSym); // copy the hMatSym matrix for the LUDecomposition
    std::pair<HMatrix,HMatrix> LUPair = HArithm::LUDecomposition(copyForLu, 0, relTruncTol); // calculate the LU decomposition of hMatSym

    HMatrix X = HArithm::forwSubsMatVal(LUPair.first, Z, 0, relTruncTol); // Solve L * X = Z for X, L is lower triangular H-matrix; X is the return value
    HMatrix deltaHMatrix = HMultiply::multiplyHMat(LUPair.first, X, 0, 0);
    HMatrix hMatSymCopy = HMatrix(hMatSym);
    HArithm::recursiveHMatSubstraction(hMatSymCopy.getRootBlock(), deltaHMatrix.getRootBlock(), 0, 0); // after the operation hMatSymCopy is the residual hMatSym - L * X
    double relAccuracy = hMatSymCopy.norm() / hMatSym.norm();
    std::cout << "Relative accuracy of forwSubsMatVal: " << relAccuracy << ". Rounding accuracy was set to " << relTruncTol << std::endl;

    Eigen::MatrixXcd matSolve;
    double relAccuracyEigen;
    const Eigen::MatrixXcd fLowerTri = LUPair.first.toFullMat();
    if(inversionInLUDecomp || qrInLUDecomp)
    {
        matSolve = fLowerTri.partialPivLu().solve(hMatSymFull);
        relAccuracyEigen = (fLowerTri * matSolve - hMatSymFull).norm() / hMatSymFull.norm();
    }
    else
    {
        matSolve = fLowerTri.triangularView<Eigen::Lower>().solve(hMatSymFull);
        relAccuracyEigen = (fLowerTri.triangularView<Eigen::Lower>() * matSolve - hMatSymFull).norm() / hMatSymFull.norm();
    }
    std::cout << "Relative accuracy of eigens matrix valued forward substitution: " << relAccuracyEigen << std::endl;

    const Eigen::MatrixXcd fUpperTri = LUPair.second.toFullMat();
    Z = HMatrix(hMatSym); // copy hMatSym again; was destroyed during forwSubsMatVal
    X = HArithm::forwSubsMatValTransposed(LUPair.second, Z, 0, relTruncTol); // Solve X * U = Z for X, U is upper triangular H-matrix, Z is destroyed
    deltaHMatrix = HMultiply::multiplyHMat(X, LUPair.second, 0, 0);
    hMatSymCopy = HMatrix(hMatSym);
    HArithm::recursiveHMatSubstraction(hMatSymCopy.getRootBlock(), deltaHMatrix.getRootBlock(), 0, 0); // after the operation hMatSymCopy is the residual hMatSym - X * U
    relAccuracy = hMatSymCopy.norm() / hMatSym.norm();
    std::cout << "Relative accuracy of forwSubsMatValTransposed: " << relAccuracy << ". Rounding accuracy was set to " << relTruncTol << std::endl;

    //clear H-matrices
    deltaHMatrix.clear();
    hMatSymCopy.clear();
    X.clear();
    Z.clear(); // redundant
    hMatSym.clear();
    LUPair.first.clear();
    LUPair.second.clear();
    QVERIFY(relAccuracy < 100 * relTruncTol + 100 * tiny); // relaxed accuracy requirements due to multiple truncations
}

void Test::normApproximations()
{
    // test matrix-valued forward substitution
    double relTruncTol = 0; // relative truncation tolerance
    // hMatSymmetricClusterTrees has the same ClusterTrees for the rows and the columns and is therefore amenable for H-LU-factorization
    HMatrix hMatSym = HMatrix(&rowClustertree, &rowClustertree, 0, relError, std::bind(&Test::implicitMatrix, this, std::placeholders::_1, std::placeholders::_2, std::ref(pointsForRows), std::ref(pointsForCols)));
    const Eigen::MatrixXcd hMatSymFull = hMatSym.toFullMat(); // for accuracy comparison

    double fullNorm = hMatSymFull.norm();
    double hNorm = hMatSym.norm();

    double relAccuracy = std::abs(hNorm - fullNorm)/ fullNorm;
    std::cout << "Relative accuracy of HMatrix norm() function: " << relAccuracy << std::endl;

    QVERIFY(relAccuracy < 100 * tiny);

    HMatrix copyForLu = HMatrix(hMatSym); // copy the hMatSym matrix for the LUDecomposition
    std::pair<HMatrix,HMatrix> LUPair = HArithm::LUDecomposition(copyForLu, 0, relTruncTol, true); // calculate the LU decomposition of hMatSym

    Eigen::BDCSVD<Eigen::MatrixXcd> svd(hMatSymFull);
    double fullConditionNumber = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size()-1);
    std::cout << "Accurate condition number from full matrix: " << fullConditionNumber << std::endl;

    double hConditionNumber =  HArithm::spectralNorm(hMatSym) * HArithm::spectralNormFromLU(LUPair.first, LUPair.second);
    std::cout << "Fast approximate condition number: " << hConditionNumber << std::endl;

    relAccuracy = std::abs(fullConditionNumber - hConditionNumber)/ fullConditionNumber;
    std::cout << "Relative accuracy: " << relAccuracy << std::endl;

    hMatSym.clear();
    copyForLu.clear();
    LUPair.first.clear();
    LUPair.second.clear();

    QVERIFY(relAccuracy < 0.2);

}

void Test::testLUFactorization()
{
    // test LU-factorization
    double relTruncTol = 0; // relative truncation tolerance
    // hMatSymmetricClusterTrees has the same ClusterTrees for the rows and the columns and is therefore amenable for H-LU-factorization
    HMatrix hMatSym = HMatrix(&rowClustertree, &rowClustertree, maxRank, relError, std::bind(&Test::implicitMatrix, this, std::placeholders::_1, std::placeholders::_2, std::ref(pointsForRows), std::ref(pointsForCols)));

    const Eigen::MatrixXcd fullCopy = hMatSym.toFullMat();
    std::pair<HMatrix,HMatrix> LUPair = HArithm::LUDecomposition(hMatSym, 0, relTruncTol); // calculate the LU decomposition of hMatSymmetricClusterTrees
    HMatrix pruductHMat = HMultiply::multiplyHMat(LUPair.first, LUPair.second, 0, 0);
    const Eigen::MatrixXcd productToFull = pruductHMat.toFullMat();
    pruductHMat.clear();

    Eigen::MatrixXcd LU = global::LUDecompNoPivoting(fullCopy);
    Eigen::MatrixXcd upper = LU.triangularView<Eigen::Upper>();
    Eigen::MatrixXcd matrFromFullLU = LU.triangularView<Eigen::UnitLower>() * upper;
    double relAccuracy = (matrFromFullLU - fullCopy).norm() / fullCopy.norm();
    std::cout << "Relative accuracy of the LU Decomposition without pivoting: " << relAccuracy << std::endl;

    Eigen::MatrixXcd residuum = productToFull - fullCopy;
    relAccuracy = residuum.norm() / fullCopy.norm();
    std::cout << "Relative accuracy of LUDecomposition: " << relAccuracy << std::endl;
    QVERIFY(relAccuracy < 100 * (relTruncTol + tiny));

    const Eigen::MatrixXcd fLowerTri = LUPair.first.toFullMat();
    const Eigen::MatrixXcd fUpperTri = LUPair.second.toFullMat();
    const Eigen::MatrixXcd triProd = fLowerTri * fUpperTri;
    residuum = triProd - fullCopy;

    relAccuracy = residuum.norm() / fullCopy.norm();
    std::cout << "Relative accuracy of LUDecomposition: " << relAccuracy << std::endl;

//    const Eigen::VectorXcd x = Eigen::VectorXcd::Random(n);
//    Eigen::VectorXcd sol = HArithm::forwardSubstitution(LUPair.first, x);

//    Eigen::VectorXcd deltaX = Eigen::VectorXcd::Zero(n);
//    HArithm::parallelMatrixVectorPoduct(deltaX, LUPair.first.getRootBlock(), sol, 0, 0);
//    relAccuracy = (deltaX - x).norm() / x.norm();
//    std::cout << "Relative accuracy of forwardSubstitution: " << relAccuracy << std::endl;
}
