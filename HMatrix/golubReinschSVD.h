#ifndef GOLUBREINSCHSVD_H
#define GOLUBREINSCHSVD_H

#include <eigen3/Eigen/Geometry>
#include <QVector>
#include <iostream>

/**
* \class GoloubReinschSvd
* \brief The class implements the Golub Reinsch singular value decomposition.
*
* The class implements the Golub Reinsch singular value decomposition found in Gene H. Golub and Charles F. van Loan. Matrix Computations. JHU Press, fourth edition. The method keeps pace with Eigen::BDCSVD up to dimension 200x200.
*/
class GoloubReinschSvd
{
public:

    static void realBidiagonalization(Eigen::MatrixXcd &A, Eigen::MatrixXcd &U, Eigen::VectorXd &diagonal, Eigen::VectorXd &superDiag, Eigen::MatrixXcd &VAdj, const double epsilon); /*!< The function factorizes a complex matrix into two orthogonal matrices and a real bidiagonal matrix. A = U * BD * VAdj */
    static void bidiagonalSVD(long diagStartIndex, long diagLength,  Eigen::MatrixXcd &U, Eigen::VectorXd& diagonal, Eigen::VectorXd& superDiagonal, Eigen::MatrixXcd &VAdj); /*!< The function computes the SVD of a real bidiagonal matrix. */

    static void chaseDiagonalZero(Eigen::MatrixXcd &U, Eigen::VectorXcd& diagonal, Eigen::VectorXcd& upperSecDiagonalDiagonal, Eigen::MatrixXcd &VAdj);
    static void givensRotToBidiagFromLeft(long lowerIndex, Eigen::JacobiRotation<double> givensRotation, double &secondLowerDiagElem, Eigen::VectorXd &diag, Eigen::VectorXd &upperSecDiag,  double &tertUpperDiagElem);
    static void givensRotToBidiagFromRight(long lowerIndex, Eigen::JacobiRotation<double> givensRotation, double &secondLowerDiagElem, Eigen::VectorXd &diag, Eigen::VectorXd &upperSecDiag, double &tertUpperDiagElem);

    static void driveOutDiagonalZero(long index, Eigen::MatrixXcd &U, Eigen::VectorXd &diag, Eigen::VectorXd &upperSecDiag);
    static void svdPostprocessing(Eigen::MatrixXcd &U, Eigen::VectorXd &singularValues, Eigen::MatrixXcd &VAdj, const double tolerance);

    /**
    * \brief The function computes the SVD of a complex matrix.
    */
    static void goloubReinschSVD(Eigen::MatrixXcd &A, Eigen::MatrixXcd &U, Eigen::VectorXd &singularValues, Eigen::MatrixXcd &VAdj);
    static Eigen::MatrixXcd bidiagonalToFullMat(const Eigen::VectorXcd& diagonal, const Eigen::VectorXcd& upperSecDiagonalDiagonal);
    static long minRankforError(const Eigen::VectorXd &singularValues, const long maxRank, const double relError);

//    static void test();

};

#endif // GOLUBREINSCHSVD_H

