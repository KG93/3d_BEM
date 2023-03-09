#ifndef TEST_H
#define TEST_H

#include "vectortriangle.h"
#include "HMatrix/clustertree.h"
#include "HMatrix/hmatrix.h"
#include "HMatrix/harithm.h"

#include <QTest>
#include <functional>

class Test : public QObject
{
    Q_OBJECT
public:
    Test(){}
    ~Test(){}

private:
    static constexpr double tiny = 10.0 * Eigen::NumTraits< std::complex<double> >::epsilon(); // tiny number for accuracy tests

    QVector<Eigen::Vector3d> generateRandomPoints(unsigned long n);
    long n = 2000; // numer of points and therefore also the matrix dimension

    QVector<Eigen::Vector3d> pointsForRows = generateRandomPoints(n);
    QVector<Eigen::Vector3d> pointsForCols = generateRandomPoints(n);

    inline static std::complex<double> helholtzGreensf(const Eigen::Vector3d source, const Eigen::Vector3d observer)
    {
        const double wavenumber = 10;
        const std::complex<double> imaginaryUnit = std::complex<double>(0.0, 1.0);

        double r = (source-observer).norm();
        if(r==0) // prevent Nan's due to the singularity
        {
            //std::cerr << "r==0" << std::endl;
            return 0;
        }
        return (std::exp(wavenumber * imaginaryUnit * r))/(4 * std::numbers::pi * r);
    }
    inline std::complex<double> implicitMatrix(const long row, const long col, const QVector<Eigen::Vector3d> &points, const QVector<Eigen::Vector3d> &points2)
    {
        return helholtzGreensf(points.at(row), points2.at(col));
    }

    ClusterTree rowClustertree = ClusterTree(&pointsForRows);
    ClusterTree columnClustertree = ClusterTree(&pointsForCols);
    const long maxRank = 0;
    const double relError = 0.01;
//    HMatrix hMatrix = HMatrix(&rowClustertree, &columnClustertree, true);
    HMatrix hMatrix = HMatrix(&rowClustertree, &columnClustertree, maxRank, relError, std::bind(&Test::implicitMatrix, this, std::placeholders::_1, std::placeholders::_2, std::ref(pointsForRows), std::ref(pointsForCols)));

private slots:
    void testHMatrixConstruction();
};

#endif // TEST_H
