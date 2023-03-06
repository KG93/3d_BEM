#include "golubReinschSVD.h"


void GoloubReinschSvd::realBidiagonalization(Eigen::MatrixXcd &A, Eigen::MatrixXcd &U, Eigen::VectorXd &diagonal, Eigen::VectorXd &superDiag, Eigen::MatrixXcd &VAdj, const double epsilon)
{
    long rows = A.rows();
    long cols = A.cols();
    bool moreColsThanRows = cols > rows;
    long diagonalLength = std::min (rows, cols);
    long superDiagLength;
    long numberOfHouseholderTransformsFromLeft;
    long numberOfHouseholderTransformsFromRight;

    if(moreColsThanRows) // long matrix
    {
        superDiagLength = diagonalLength;
        numberOfHouseholderTransformsFromLeft = diagonalLength;
        numberOfHouseholderTransformsFromRight = diagonalLength;
    }
    else // tall or square matrix
    {
        superDiagLength = diagonalLength - 1;
        numberOfHouseholderTransformsFromLeft = diagonalLength;
        numberOfHouseholderTransformsFromRight = superDiagLength;
    }

    Eigen::VectorXcd tmpHouseholderVectorLeft;
    Eigen::VectorXcd householderCoefficientsLeft(numberOfHouseholderTransformsFromLeft);
    Eigen::VectorXcd scalingCoefficientsLeft(numberOfHouseholderTransformsFromLeft);

    Eigen::RowVectorXcd tmpHouseholderVectorRight;
    Eigen::VectorXcd householderCoefficientsRight(numberOfHouseholderTransformsFromRight);
    Eigen::VectorXcd scalingCoefficientsRight(numberOfHouseholderTransformsFromRight);

    diagonal.resize(diagonalLength);
    superDiag.resize(superDiagLength);


    for(long index = 0 ; index < std::max(numberOfHouseholderTransformsFromLeft, numberOfHouseholderTransformsFromRight); index++)
    {
        if(index < numberOfHouseholderTransformsFromLeft)
        {
            // calculate left housholder vector and coefficient
            tmpHouseholderVectorLeft = A.block(index, index, rows - index, 1);
            double colNorm2 = tmpHouseholderVectorLeft.squaredNorm();
            double colNorm = sqrt(colNorm2);

            if(colNorm > epsilon) // avoid division by zero
            {
                double absTmpHouseholderVectorLeft0 = std::abs(tmpHouseholderVectorLeft(0));
                std::complex<double> xi;
                if(absTmpHouseholderVectorLeft0 == 0)
                {
                    xi = colNorm;
                }
                else
                {
                    xi = colNorm * (tmpHouseholderVectorLeft(0) / absTmpHouseholderVectorLeft0);
                }
                tmpHouseholderVectorLeft(0) += xi;

                householderCoefficientsLeft(index) = colNorm2 + colNorm * absTmpHouseholderVectorLeft0;

                //apply householder rotation to A from the left
                A.bottomRightCorner(rows - index, cols - index - 1) = A.bottomRightCorner(rows - index, cols - index - 1) - tmpHouseholderVectorLeft * ((tmpHouseholderVectorLeft.adjoint() / householderCoefficientsLeft(index)) * A.bottomRightCorner(rows - index, cols - index - 1));

                // overwrite A  lower column with the essential householder vector
                A.block(index, index, rows - index, 1) = tmpHouseholderVectorLeft;

                diagonal(index) = std::abs(xi);
                scalingCoefficientsLeft(index) = (- xi / diagonal(index));
                A.block(index, index + 1, 1, cols - index - 1) /= scalingCoefficientsLeft(index);
            }
            else
            {
                diagonal(index) = 0.0;
                householderCoefficientsLeft(index) = 1.0;
                scalingCoefficientsLeft(index) = 1.0;
            }
        }

//        std::cout << "a left transform done, now right" << std::endl;

        // calculate right housholder vector and coefficient
        if(index < numberOfHouseholderTransformsFromRight)
        {
            tmpHouseholderVectorRight = A.block(index, index + 1, 1, cols - index - 1);
            double rowNorm2 = tmpHouseholderVectorRight.squaredNorm();
            double rowNorm = sqrt(rowNorm2);

            if(rowNorm > epsilon) // avoid division by zero
            { ////////////////index had to be skipped -> todo
                double absTmpHouseholderVectorRight0 = std::abs(tmpHouseholderVectorRight(0));

                std::complex<double> xiRight;
                if(absTmpHouseholderVectorRight0 == 0)
                {
                    xiRight = rowNorm;
                }
                else
                {
                    xiRight = rowNorm * (tmpHouseholderVectorRight(0) / absTmpHouseholderVectorRight0);
                }
                tmpHouseholderVectorRight(0) += xiRight;

                householderCoefficientsRight(index) = rowNorm2 + rowNorm * absTmpHouseholderVectorRight0;

                //apply modified householder rotation to A from the right
                A.bottomRightCorner(rows - index - 1, cols - index - 1) = A.bottomRightCorner(rows - index - 1, cols - index - 1) - (A.bottomRightCorner(rows - index - 1, cols - index - 1) * (tmpHouseholderVectorRight.adjoint() / householderCoefficientsRight(index))) * tmpHouseholderVectorRight;

                // overwrite A  right part of row with the essential householder vector
                A.block(index, index + 1, 1, cols - index - 1) = tmpHouseholderVectorRight;

                superDiag(index) =  std::abs(xiRight);
                scalingCoefficientsRight(index) = (- xiRight / superDiag(index));
                A.block(index + 1, index + 1, rows - index - 1, 1) /= scalingCoefficientsRight(index);
            }
            else
            {
//                std::cerr << "here2" << std::endl;
                superDiag(index) = 0.0;
                householderCoefficientsRight(index) = 1.0;
                scalingCoefficientsRight(index) = 1.0;
            }
        }
    }
//    std::cout << "householder transforms done" << std::endl;

    // compute U
    U = Eigen::MatrixXcd::Identity(rows, rows);
    if(numberOfHouseholderTransformsFromLeft >= 1)
    {
        U(numberOfHouseholderTransformsFromLeft - 1, numberOfHouseholderTransformsFromLeft - 1) *= scalingCoefficientsLeft(numberOfHouseholderTransformsFromLeft - 1); //scale U-row, but U is diagonal at this point
        tmpHouseholderVectorLeft = A.block(numberOfHouseholderTransformsFromLeft - 1, numberOfHouseholderTransformsFromLeft - 1, rows - (numberOfHouseholderTransformsFromLeft - 1), 1); // essential part of the householder vector as stored in A earlier
        U.bottomRightCorner(rows - (numberOfHouseholderTransformsFromLeft - 1), rows - (numberOfHouseholderTransformsFromLeft - 1)) = U.bottomRightCorner(rows - (numberOfHouseholderTransformsFromLeft - 1), rows - (numberOfHouseholderTransformsFromLeft - 1)) - tmpHouseholderVectorLeft * ((tmpHouseholderVectorLeft.adjoint() / householderCoefficientsLeft(numberOfHouseholderTransformsFromLeft - 1)) * U.bottomRightCorner(rows - (numberOfHouseholderTransformsFromLeft - 1), rows - (numberOfHouseholderTransformsFromLeft - 1)));

        for(long index = numberOfHouseholderTransformsFromLeft - 2 ; index >= 0; index--)
        {
            U.block(index ,index, 1, rows - index) *= scalingCoefficientsLeft(index); // scale U.row(index); but U is diagonal, except for lower right corner full matrix
            tmpHouseholderVectorLeft = A.block(index, index, rows - index, 1); // essential part of the householder vector as stored in A earlier
            U.bottomRightCorner(rows - index, rows - index) = U.bottomRightCorner(rows - index, rows - index) - tmpHouseholderVectorLeft * ((tmpHouseholderVectorLeft.adjoint() / householderCoefficientsLeft(index)) * U.bottomRightCorner(rows - index, rows - index));
        }
    }

    // compute V^h
    VAdj = Eigen::MatrixXcd::Identity(cols, cols);
    if(numberOfHouseholderTransformsFromRight >= 1)
    {
        VAdj(numberOfHouseholderTransformsFromRight, numberOfHouseholderTransformsFromRight) *= scalingCoefficientsRight(numberOfHouseholderTransformsFromRight - 1); //scale V-col, but U is diagonal at this point
        tmpHouseholderVectorRight = A.block(numberOfHouseholderTransformsFromRight - 1, numberOfHouseholderTransformsFromRight, 1, cols - numberOfHouseholderTransformsFromRight); // essential part of the householder vector as stored in A earlier
        VAdj.bottomRightCorner(cols - numberOfHouseholderTransformsFromRight, cols - numberOfHouseholderTransformsFromRight) = VAdj.bottomRightCorner(cols - numberOfHouseholderTransformsFromRight, cols - numberOfHouseholderTransformsFromRight) - (VAdj.bottomRightCorner(cols - numberOfHouseholderTransformsFromRight, cols - numberOfHouseholderTransformsFromRight) * (tmpHouseholderVectorRight.adjoint() / householderCoefficientsRight(numberOfHouseholderTransformsFromRight - 1))) * tmpHouseholderVectorRight;

        for(long index = numberOfHouseholderTransformsFromRight - 2; index >= 0; index--)
        {
            VAdj.block(index + 1, index + 1, cols - index - 1, 1) *= scalingCoefficientsRight(index); // scale VAdj.col(index + 1); but VAdj is diagonal, except for lower right corner full matrix
            tmpHouseholderVectorRight = A.block(index, index + 1, 1, cols - index - 1); // essential part of the householder vector as stored in A earlier
            VAdj.bottomRightCorner(cols - index - 1, cols - index - 1) = VAdj.bottomRightCorner(cols - index - 1, cols - index - 1)  - (VAdj.bottomRightCorner(cols - index - 1, cols - index - 1) * (tmpHouseholderVectorRight.adjoint() / householderCoefficientsRight(index))) * tmpHouseholderVectorRight;
         }
    }
}

void GoloubReinschSvd::bidiagonalSVD(long diagStartIndex, long diagLength,  Eigen::MatrixXcd &U, Eigen::VectorXd &diagonal, Eigen::VectorXd &superDiagonal, Eigen::MatrixXcd &VAdj)
{
    // it is presumed here, that diagonal.size() == upperSecDiagonalDiagonal.size() + 1 !!!

//    std::cerr << "bidiagonalSVD2" << std::endl;

//    long diagLength = diagonal.size();
//    long upperSecDiagonalLength = upperSecDiagonalDiagonal.size();
    long superDiagLength = diagLength - 1;
    if(diagLength == 1)
    {
        return;
    }

    // compute eigenvalues of trailing 2x2 submatrix of T = B^h * B

    double T11 = std::pow(diagonal(diagStartIndex + diagLength - 2), 2);
    if(diagLength > 2) // for diagLength== 2 there is only 1 superdiagnal element
    {
        T11 += std::pow(superDiagonal(diagStartIndex + superDiagLength - 2), 2);
    }
    double T12 = diagonal(diagStartIndex + diagLength - 2) * superDiagonal(diagStartIndex + superDiagLength - 1);
//    double T21 = T12;
    double T22 = std::pow(diagonal(diagStartIndex + diagLength - 1), 2) + std::pow(superDiagonal(diagStartIndex + superDiagLength - 1), 2);

    double T11PlusT22Half = (T11 + T22)/2.0;
    double pqRest = sqrt(std::pow(T11PlusT22Half, 2) - T11 * T22 + std::pow(T12, 2));

//    //kontrollrechnung
//    Eigen::MatrixXcd fullMat = bidiagonalToFullMat(diagonal.segment(diagStartIndex, diagLength), superDiagonal.segment(diagStartIndex,superDiagLength));
//    Eigen::Matrix2cd control = (fullMat.adjoint() * fullMat).bottomRightCorner(2,2);
//    std::cerr << control(0,0) -  T11<< std::endl;
//    std::cerr << control(0,1) -  T12 << std::endl;
//    std::cerr << control(1,0) -  T21 << std::endl;
//    std::cerr << control(1,1) -  T22 << std::endl;
//    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces;
//    ces.compute(control);
//    std::cout << "The eigenvalues of A are:" << std::endl << ces.eigenvalues() << std::endl;
//    //kontrollrechnung ende

//    //eigenvalues
    double lambda1T = T11PlusT22Half + pqRest;
    double lambda2T = T11PlusT22Half - pqRest;
//    std::cerr << "lambda1T: " << lambda1T << std::endl;
//    std::cerr << "lambda2T: " << lambda2T << std::endl;

//    // choose eigenvalue that is closer to T22 as Wilkinson shift
    double wilkinsonShift;
    if(std::abs(lambda1T - T22) <  std::abs(lambda2T - T22))
    {
        wilkinsonShift = lambda1T;
    }
    else
    {
        wilkinsonShift = lambda2T;
    }
//    std::cout << "wilkinsonShift: " << wilkinsonShift << std::endl;
    double tertiaryUpperDiagElement = 0;
    double secondaryLowerDiagElement = 0;

//    // compute first givens rotation from the right
    Eigen::JacobiRotation<double> givensRotationfromRight;
    if( std::isfinite(wilkinsonShift) ) // avoid overflow
    {
        givensRotationfromRight.makeGivens((std::pow(diagonal(diagStartIndex + 0), 2) - wilkinsonShift), diagonal(diagStartIndex + 0) * superDiagonal(diagStartIndex + 0));
    }
    else
    {
        givensRotationfromRight.makeGivens(diagonal(diagStartIndex + 0), superDiagonal(diagStartIndex + 0));
    }
    givensRotToBidiagFromRight(diagStartIndex + 0, givensRotationfromRight.transpose().adjoint(), secondaryLowerDiagElement, diagonal, superDiagonal,  tertiaryUpperDiagElement);

    VAdj.applyOnTheLeft(diagStartIndex + 0, diagStartIndex + 1, givensRotationfromRight.transpose());

    //chase unwanted offdiagonal elements
    Eigen::JacobiRotation<double> givensRotationFromLeft;
    for(long index = diagStartIndex + 0 ; index < diagStartIndex + (diagLength - 1); index++)
    {
        if(index != diagStartIndex) // no givensRotationfromRight on first iteration
        {

            //apply givens rotation from the right
            givensRotationfromRight.makeGivens(superDiagonal(index - 1), tertiaryUpperDiagElement);
            givensRotToBidiagFromRight(index, givensRotationfromRight /*givensRotationfromRight.transpose().adjoint()*/, secondaryLowerDiagElement, diagonal, superDiagonal,  tertiaryUpperDiagElement);
            VAdj.applyOnTheLeft(/*diagStartIndex + */index , /*diagStartIndex + */index + 1, givensRotationfromRight.transpose());

        }
        //apply givens rotation from the left
        givensRotationFromLeft.makeGivens(diagonal(index), secondaryLowerDiagElement);
        givensRotToBidiagFromLeft(index, givensRotationFromLeft.transpose(), secondaryLowerDiagElement, diagonal, superDiagonal,  tertiaryUpperDiagElement);
        U.applyOnTheRight(/*diagStartIndex + */index , /*diagStartIndex + */index + 1, givensRotationFromLeft);
    }
}

void GoloubReinschSvd::givensRotToBidiagFromLeft(long lowerIndex, Eigen::JacobiRotation<double> givensRotation, double &secondLowerDiagElem, Eigen::VectorXd &diag, Eigen::VectorXd &upperSecDiag,  double &tertUpperDiagElem)
{
    if(lowerIndex < diag.size() - 2)
    {
        Eigen::Vector2d col1(diag(lowerIndex), secondLowerDiagElem);
        col1.applyOnTheLeft(0, 1, givensRotation);
        Eigen::Vector2d col2(upperSecDiag(lowerIndex), diag(lowerIndex + 1));
        col2.applyOnTheLeft(0, 1, givensRotation);
        Eigen::Vector2d col3(0.0, upperSecDiag(lowerIndex + 1));
        col3.applyOnTheLeft(0, 1, givensRotation);

        diag(lowerIndex) = col1(0);
        secondLowerDiagElem = 0.0;

        upperSecDiag(lowerIndex) = col2(0);
        diag(lowerIndex + 1) = col2(1);

        tertUpperDiagElem = col3(0);
        upperSecDiag(lowerIndex + 1) = col3(1);
    }
    else
    {
        Eigen::Vector2d col1(diag(lowerIndex), secondLowerDiagElem);
        col1.applyOnTheLeft(0, 1, givensRotation);
        Eigen::Vector2d col2(upperSecDiag(lowerIndex), diag(lowerIndex + 1));
        col2.applyOnTheLeft(0, 1, givensRotation);

        diag(lowerIndex) = col1(0);
        secondLowerDiagElem = 0.0;

        upperSecDiag(lowerIndex) = col2(0);
        diag(lowerIndex + 1) = col2(1);

        tertUpperDiagElem = 0.0;
    }
}

void GoloubReinschSvd::givensRotToBidiagFromRight(long lowerIndex, Eigen::JacobiRotation<double> givensRotation, double &secondLowerDiagElem, Eigen::VectorXd &diag, Eigen::VectorXd &upperSecDiag, double &tertUpperDiagElem)
{
    if(lowerIndex == 0)
    {
        Eigen::RowVector2d row2(diag(lowerIndex), upperSecDiag(lowerIndex));
        row2.applyOnTheRight(0, 1, givensRotation);
        Eigen::RowVector2d row3(0.0, diag(lowerIndex + 1));
        row3.applyOnTheRight(0, 1, givensRotation);

        diag(lowerIndex) = row2(0);
        upperSecDiag(lowerIndex) = row2(1);

        secondLowerDiagElem = row3(0);
        diag(lowerIndex + 1) = row3(1);
    }
    else
    {
        Eigen::RowVector2d row1(upperSecDiag(lowerIndex - 1), tertUpperDiagElem);
        row1.applyOnTheRight(0, 1, givensRotation);
        Eigen::RowVector2d row2(diag(lowerIndex), upperSecDiag(lowerIndex));
        row2.applyOnTheRight(0, 1, givensRotation);
        Eigen::RowVector2d row3(0.0, diag(lowerIndex + 1));
        row3.applyOnTheRight(0, 1, givensRotation);

        upperSecDiag(lowerIndex - 1) = row1(0);
        tertUpperDiagElem = 0.0;

        diag(lowerIndex) = row2(0);
        upperSecDiag(lowerIndex) = row2(1);

        secondLowerDiagElem = row3(0);
        diag(lowerIndex + 1) = row3(1);
    }
}

void GoloubReinschSvd::driveOutDiagonalZero(long diagIndex, Eigen::MatrixXcd &U, Eigen::VectorXd &diag, Eigen::VectorXd &upperSecDiag)
{
//    std::cout << "driveOutDiagonalZero" << std::endl;

    Eigen::JacobiRotation<double> givensRotationFromLeft;
    double drivenElement = upperSecDiag(diagIndex); // matrix element to be driven towards the right out of the matrix

    for(long i = diagIndex + 1; i < diag.size(); i++)
    {
//        (c, s^)  a = x     -> (-d , e^) f = 0  -> (d , e^) -f = 0
//        (-s, c^) b   0        (e, d^)  g  x       (-e, d^)  g  x
        givensRotationFromLeft.makeGivens(diag(i), drivenElement);

        double c = givensRotationFromLeft.c();
        double s = givensRotationFromLeft.s();

        givensRotationFromLeft = Eigen::JacobiRotation<double>(c, - s).transpose();
        if(i < diag.size() - 1)
        {
            Eigen::Vector2d col1(drivenElement, diag(i));
            col1.applyOnTheLeft(0, 1, givensRotationFromLeft);
            Eigen::Vector2d col2(0.0, upperSecDiag(i));
            col2.applyOnTheLeft(0, 1, givensRotationFromLeft);

            diag(i) = col1(1);
            drivenElement = col2(0);
            upperSecDiag(i) = col2(1);
            U.applyOnTheRight(diagIndex, i, givensRotationFromLeft.adjoint());
        }
        else
        {
            Eigen::Vector2d col1(drivenElement, diag(i));
            col1.applyOnTheLeft(0, 1, givensRotationFromLeft);
            diag(i) = col1(1);
            U.applyOnTheRight(diagIndex, i, givensRotationFromLeft.adjoint());
        }
    }
    upperSecDiag(diagIndex) = 0.0; // first driven element was a super diagonal element which has to be set zo zero
}

void GoloubReinschSvd::svdPostprocessing(Eigen::MatrixXcd &U, Eigen::VectorXd &singularValues, Eigen::MatrixXcd &VAdj, const double tolerance)
{
    long rows = U.rows();
    long cols = VAdj.cols();
    long rank = singularValues.size();

    for(long i = 0 ; i < rank; i++) // make singular values positive
    {
        if(/*std::real*/(singularValues(i)) < 0.0)
        {
            singularValues(i) = - singularValues(i);
            U.col(i) = - U.col(i);
        }
    }

    // sort singular values in descending order
    QVector<long> indices(rank); // indices holds permutation need to sort the singular values

    for(long i = 0 ; i < rank; i++)
    {
        indices[i] = i;
    }

    std::sort(indices.begin(), indices.end(),
              [&](const long& a, const long& b) {
                  return (/*std::abs*/(singularValues(a)) > /*std::abs*/(singularValues(b)));
              });

    long newRank = 0;
    for(long i = 0 ; i < rank; i++) // find number of nonzero singuar values
    {
        if(/*std::abs*/(singularValues(indices.at(i))) > tolerance)
        {
            newRank++;
        }
        else
        {
            break;
        }
    }
    rank = std::max(newRank, (long) 1);

    Eigen::VectorXd sortedSingVal(rank);
    Eigen::MatrixXcd newU(rows, rank);
    Eigen::MatrixXcd newVadj(rank, cols);

    for(long i = 0 ; i < rank; i++) // reorder singular values, U-cols and V-rows
    {
        sortedSingVal(i) = singularValues(indices.at(i));
        newU.col(i) = U.col(indices.at(i));
        newVadj.row(i) = VAdj.row(indices.at(i));
    }
    singularValues = sortedSingVal;
    U = newU;
    VAdj = newVadj;
}

void GoloubReinschSvd::goloubReinschSVD(Eigen::MatrixXcd &A, Eigen::MatrixXcd &U, Eigen::VectorXd &singularValues, Eigen::MatrixXcd &VAdj) // The SVD  Algorithm 8.6.2 from Matrix Computations 4. Edition by Goloub and Van Loan
{
    Eigen::NumTraits< std::complex<double> > numT;
    double epsilon = 10.0 * numT.epsilon();

    long rows = A.rows();
    long cols = A.cols();
    long minDimension = std::min(rows, cols);
    long minDimension2 = minDimension * minDimension;
    long minDimension3 = minDimension * minDimension2;

    bool wideMat = rows < cols; // matrix istn't tall or square -> transpose
    if(wideMat)
    {
        A = A.transpose().eval();
        rows = A.rows();
        cols = A.cols();
    }

    Eigen::VectorXd superDiagonal;

    realBidiagonalization(A, U, singularValues, superDiagonal, VAdj, epsilon); // transform A to bidiagonal Form, Information is contained in diagonal and upperSecDiagonalDiagonal

    U.conservativeResize(Eigen::NoChange, singularValues.size());
    VAdj.conservativeResize(singularValues.size(), Eigen::NoChange);


    long qrIterations = 0;
    long q = 0;
    while(q < cols - 1)
    {
        qrIterations++;
        if(qrIterations == minDimension2) // maybe no convergence due to rounding threshold too low
        {
            epsilon *= 100;
        }
        if(qrIterations == minDimension3) // no convergence
        {
            std::cerr << "No convergence in goloubReinschSVD()." << std::endl;
            break;
        }
        for(long i = 0; i < cols - 1; i++) // set b(i,i+1) to zero if abs(b(i,i+1)) <= epsilon * ( abs(b(i,i)) + abs(b(i+1,i+1)))
        {
            if(std::abs(superDiagonal(i)) <= epsilon * (std::abs(singularValues(i)) + std::abs(singularValues(i+1))))
            {
                superDiagonal(i) = 0.0;
            }
        }
        for(long i = superDiagonal.size() - 1 - q; i >= 0; i--) // find size (q) of largest diagonal matrix in the lower right corner
        {
            if(std::abs(superDiagonal(i)) <= epsilon)
            {
                q++;
            }
            else
            {
                break;
            }
        }
        long p = 0;
        for(long i = superDiagonal.size() - 1 - q; i >= 0; i--) // find size (cols -p -q) of largest bidiagonal matrix in the mid
        {
            if(std::abs(superDiagonal(i)) <= epsilon)
            {
                p = i + 1; // i is index, so size = index + 1
                break;
            }
        }
        if(q < cols)
        {
            bool noZeroOnDiagonal = true;
            for(long i = p; i < cols - q - 1; i++)
            {
                if(std::abs(singularValues(i)) <= epsilon)
                {

                    driveOutDiagonalZero(i, U, singularValues, superDiagonal);
                    noZeroOnDiagonal = false;
                }
            }
            if(noZeroOnDiagonal)
            {
                bidiagonalSVD(p, cols - p - q, U, singularValues, superDiagonal, VAdj);
            }
        }
    }
    svdPostprocessing(U, singularValues, VAdj, epsilon);
    if(wideMat)
    {
        Eigen::MatrixXcd tmp = U;
        U = VAdj.transpose();
        VAdj = tmp.transpose();
    }
}

Eigen::MatrixXcd GoloubReinschSvd::bidiagonalToFullMat(const Eigen::VectorXcd& diagonal, const Eigen::VectorXcd& superDiagonal)
{
    Eigen::MatrixXcd bidiagonalMatrix;
    if(!(superDiagonal.size() + 1 == diagonal.size() || superDiagonal.size() == diagonal.size()))
    {
        std::cerr << "incompatible diagonal and superDiagonal size in bidiagonalToFullMat call" << std::endl;

        return bidiagonalMatrix;
    }
    if(superDiagonal.size() + 1 == diagonal.size())
    {
        bidiagonalMatrix = Eigen::MatrixXcd::Zero(diagonal.size(), diagonal.size());
        bidiagonalMatrix.diagonal() = diagonal;
        for(long index = 0; index < superDiagonal.size(); index++)
        {
            bidiagonalMatrix(index, index + 1) = superDiagonal(index);
        }
    }
    else
    {
        bidiagonalMatrix = Eigen::MatrixXcd::Zero(diagonal.size(), diagonal.size() + 1);
        bidiagonalMatrix.diagonal() = diagonal;
        for(long index = 0; index < diagonal.size(); index++)
        {
            bidiagonalMatrix(index, index + 1) = superDiagonal(index);
        }
    }

    return bidiagonalMatrix;
}

long GoloubReinschSvd::minRankforError(const Eigen::VectorXd &singularValues, const long maxRank, const double relError)
{
    const long nonzeroSingularValues = singularValues.size();
    if(relError <= 0)
    {
        if(maxRank <= 0)
        {
            return nonzeroSingularValues;
        }
        else
        {
            return std::min(nonzeroSingularValues, maxRank);
        }
    }
    else
    {
        long maxIndex = nonzeroSingularValues;
        if(maxRank > 0)
        {
            maxIndex = std::min(maxRank, nonzeroSingularValues);
        }
        double frobeniousNorm2FullMatrix = singularValues.squaredNorm();
        double frobeniousNorm2RkResidue = frobeniousNorm2FullMatrix;
        double errorSquared = std::pow(relError, 2);
        for(long optimalRank = 1; optimalRank <= maxIndex; optimalRank++) //for absolute error
        {
            frobeniousNorm2RkResidue -= std::pow(singularValues(optimalRank - 1), 2);
            if(frobeniousNorm2RkResidue / frobeniousNorm2FullMatrix < errorSquared) //for relative error in frobebinous norm
            {
                return optimalRank;
            }
        }
        return maxIndex;
    }
}
