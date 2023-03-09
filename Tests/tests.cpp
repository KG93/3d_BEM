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

void Test::testHMatrixConstruction()
{
//    std::cout << "clustertree.size(): " << rowClustertree.size() << std::endl;
//    const long maxRank = 10;
//    const double relError = 0.05;
//    hMatrix.assembleBlocks(maxRank, relError, std::bind(&Test::implicitMatrix, this, std::placeholders::_1, std::placeholders::_2, std::ref(pointsForRows), std::ref(pointsForCols)));

    // test hMat to full matrix conversion and the relative accuracy of the HMatrix approximation
    Eigen::MatrixXcd hMatAsFull = HArithm::hMatToFullMat(hMatrix);

    Eigen::MatrixXcd test(n,n);
    for(long row = 0; row < n; row++)
    {
        for(long col = 0; col < n; col++)
        {
            test(row, col) = implicitMatrix(row, col, pointsForRows, pointsForCols);
        }
    }
    double relAccuracy = (test - hMatAsFull).norm() / test.norm();
    std::cout << "Relative accuracy of HMatrix approximation: " << relAccuracy << ". Desired accuracy was set to " << relError << std::endl;
    QVERIFY(relAccuracy < 2 * relError); // relaxed accuracy requirements; the stopping criterion for the ACA asseembly is a heuristic

    // test matrix-vector product
    Eigen::VectorXcd x = Eigen::VectorXcd::Random(n);
    Eigen::VectorXcd resultFull = hMatAsFull * x;
    Eigen::VectorXcd resultHMat = Eigen::VectorXcd::Zero(n);
    HArithm::MVM(resultHMat, hMatrix, x);
    relAccuracy =  (resultHMat - resultFull).norm() / resultFull.norm();
    std::cout << "Relative accuracy of MVM: " << relAccuracy << std::endl;
    QVERIFY(relAccuracy < tiny);

    // test recursive matrix-vector product
    x = Eigen::VectorXcd::Random(n);
    resultFull = hMatAsFull * x;
    resultHMat = Eigen::VectorXcd::Zero(n);
    HArithm::recursiveMatrixVectorPoduct(resultHMat, hMatrix.getRootBlock(),x, 0, 0);
    relAccuracy =  (resultHMat - resultFull).norm() / resultFull.norm();
    std::cout << "Relative accuracy of recursiveMatrixVectorPoduct: " << relAccuracy << std::endl;
    QVERIFY(relAccuracy < tiny);


    // test parallel recursive matrix-vector product
    x = Eigen::VectorXcd::Random(n);
    resultFull = hMatAsFull * x;
    resultHMat = Eigen::VectorXcd::Zero(n);
    HArithm::parallelMatrixVectorPoduct(resultHMat, hMatrix.getRootBlock(),x, 0, 0);
    relAccuracy =  (resultHMat - resultFull).norm() / resultFull.norm();
    std::cout << "Relative accuracy of parallelMatrixVectorPoduct: " << relAccuracy << std::endl;
    QVERIFY(relAccuracy < tiny);

    // test vector-matrix product
    Eigen::RowVectorXcd xT = Eigen::RowVectorXcd::Random(n);
    Eigen::RowVectorXcd resultFullT = xT * hMatAsFull;
    Eigen::RowVectorXcd resultHMatT = Eigen::RowVectorXcd::Zero(n);
    HArithm::VMM(resultHMatT, xT, hMatrix);
    relAccuracy =  (resultHMatT - resultFullT).norm() / resultFullT.norm();
    std::cout << "Relative accuracy of VMM: " << relAccuracy << std::endl;
    QVERIFY(relAccuracy < tiny);

    // test recursive vector-matrix product
    xT = Eigen::RowVectorXcd::Random(n);
    resultFullT = xT * hMatAsFull;
    resultHMatT = Eigen::RowVectorXcd::Zero(n);
    HArithm::recursiveVectorMatrixPoduct(resultHMatT, xT, hMatrix.getRootBlock(), 0, 0);
    relAccuracy =  (resultHMatT - resultFullT).norm() / resultFullT.norm();
    std::cout << "Relative accuracy of recursiveMatrixVectorPoduct: " << relAccuracy << std::endl;
    QVERIFY(relAccuracy < tiny);

    // test parallel subtractive vector-matrix product
    xT = Eigen::RowVectorXcd::Random(n);
    resultFullT = -xT * hMatAsFull;
    resultHMatT = Eigen::RowVectorXcd::Zero(n);
    HArithm::subtractiveParallelVectorMatrixPoduct(resultHMatT, xT, hMatrix.getRootBlock(), 0, 0);
    relAccuracy =  (resultHMatT - resultFullT).norm() / resultFullT.norm();
    std::cout << "Relative accuracy of subtractiveParallelVectorMatrixPoduct: " << relAccuracy << std::endl;
    QVERIFY(relAccuracy < tiny);

    // test matrix-matrix multiplication
    double multiplyRelError = 0.0001;
    HMatrix hMatFromLeft = HMatrix(&columnClustertree, &rowClustertree, maxRank, relError, std::bind(&Test::implicitMatrix, this, std::placeholders::_1, std::placeholders::_2, std::ref(pointsForRows), std::ref(pointsForCols)));
    Eigen::MatrixXcd hMatFromLeftAsFull = HArithm::hMatToFullMat(hMatFromLeft);
    HMatrix pruductHMat = HMultiply::multiplyHMat(hMatFromLeft, hMatrix, 0, multiplyRelError);
    x = Eigen::VectorXcd::Random(n);
    resultFull = hMatFromLeftAsFull * (hMatAsFull * x);
    resultHMat = Eigen::VectorXcd::Zero(n);
    HArithm::MVM(resultHMat, pruductHMat, x);
    relAccuracy = (resultFull - resultHMat).norm() / resultFull.norm();
    std::cout << "Relative accuracy of HMultiply::multiplyHMat: " << relAccuracy << ". Rounding accuracy was set to " << multiplyRelError << std::endl;
    QVERIFY(relAccuracy < 10 * multiplyRelError + 100 * tiny); // relaxed accuracy requirements due to multiple truncations
    hMatFromLeft.clear();
    hMatrix.clear();
    pruductHMat.clear();
}
