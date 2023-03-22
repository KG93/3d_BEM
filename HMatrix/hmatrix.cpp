#include "hmatrix.h"


long BlockCluster::rows() const
{
    if(rowCluster != nullptr)
    {
        return rowCluster->indices.last() - rowCluster->indices.first() + 1;
    }
    else
    {
        return 0;
    }
}

long BlockCluster::cols() const
{
    if(columnCluster != nullptr)
    {
        return columnCluster->indices.last() - columnCluster->indices.first() + 1;
    }
    else
    {
        return 0;
    }
}

long BlockCluster::size() const
{
    return rows() * cols();
}

long BlockCluster::rowStartIndex() const
{
    if(rowCluster != nullptr)
    {
        return rowCluster->indices.first();
    }
    else
    {
        return 0;
    }
}

long BlockCluster::colStartIndex() const
{
    if(columnCluster != nullptr)
    {
        return columnCluster->indices.first();
    }
    else
    {
        return 0;
    }
}

double BlockCluster::compression() const
{
    return ((double) size()) / (UMat.size() + singularValues.size() + VAdjMat.size() + fullMat.size());
}

void BlockCluster::clearMatrixInfo()
{
    UMat.resize(0,0);
    singularValues.resize(0);
    VAdjMat.resize(0,0);
    fullMat.resize(0,0);
}

void BlockCluster::clear()
{
    clearMatrixInfo();
    trimBelow();
}

void BlockCluster::trimBelow()
{
    isLeaf = true;
    if(son11 != nullptr)
    {
        son11->clear();
        son11 = nullptr;
    }
    if(son12 != nullptr)
    {
        son12->clear();
        son12 = nullptr;
    }
    if(son21 != nullptr)
    {
        son21->clear();
        son21 = nullptr;
    }
    if(son22 != nullptr)
    {
        son22->clear();
        son22 = nullptr;
    }
}

void BlockCluster::compressAdmissibleBlock(const long rank, const double relError)
{
    long localRank = UMat.cols();
    if(!isAdmissible)
    {
        return;
    }
    else if(relError <= 0 && (localRank <= rank || rank <= 0)) //might be wrong with rk mat in svd form
    {
//        std::cout << "Unneccessary RMatRankReduction() call." << std::endl;
        return;
    }
    else    //block is far field
    {
        long localARank = std::min(localRank, UMat.rows());
        long localBRank = std::min(localRank, VAdjMat.cols());

        Eigen::HouseholderQR<Eigen::MatrixXcd> qrA(UMat /** matrixBlock.singularValues.asDiagonal()*/);
        Eigen::HouseholderQR<Eigen::MatrixXcd> qrB_T(VAdjMat.transpose());
        Eigen::MatrixXcd R_A = qrA.matrixQR().triangularView<Eigen::Upper>().toDenseMatrix().topRows(localARank);
        Eigen::MatrixXcd R_B_Transpose = qrB_T.matrixQR().triangularView<Eigen::Upper>().transpose().toDenseMatrix().leftCols(localBRank);

        Eigen::MatrixXcd QA = qrA.householderQ() * Eigen::MatrixXcd::Identity(UMat.rows(), localARank);//.leftCols(localARank);
        Eigen::MatrixXcd QB = qrB_T.householderQ() * Eigen::MatrixXcd::Identity(VAdjMat.cols(), localBRank);//.leftCols(localBRank);

        Eigen::MatrixXcd A = (R_A * singularValues.head(R_A.cols()).asDiagonal()) * R_B_Transpose.triangularView<Eigen::Lower>();
        Eigen::MatrixXcd U;
        Eigen::VectorXd singVals;
        Eigen::MatrixXcd VAdj;

        GoloubReinschSvd::goloubReinschSVD(A, U, singVals , VAdj);

        long newRank = GoloubReinschSvd::minRankforError(singVals, rank, relError);

        UMat = QA * U.leftCols(newRank);
        singularValues = singVals.head(newRank);
        VAdjMat = VAdj.topRows(newRank) * QB.transpose();
        frobeniusNorm = singularValues.norm();
    }
}

void BlockCluster::compressSVDBlock(const long rank, const double error)
{
    long newRank = GoloubReinschSvd::minRankforError(singularValues.cwiseAbs(), rank, error);

    UMat.conservativeResize(Eigen::NoChange, newRank);
    singularValues.conservativeResize(newRank);
    VAdjMat.conservativeResize(newRank, Eigen::NoChange);
}

bool BlockCluster::hasNan() const
{
    return UMat.hasNaN() || singularValues.hasNaN() || VAdjMat.hasNaN() || fullMat.hasNaN();
}

double BlockCluster::norm()
{
    if(isAdmissible)
    {
        std::complex<double> blockNormSquared = 0;
        for(int l = 0; l < singularValues.size(); l++)
        {
            blockNormSquared += std::abs(std::pow(singularValues(l), 2)) * UMat.col(l).squaredNorm() * VAdjMat.row(l).squaredNorm();
            for(long j = 0; j < l; j++)
            {
                blockNormSquared += 2.0 * (UMat.col(l).adjoint() * UMat.col(j))(0,0) * ((singularValues(j) * VAdjMat.row(j)) * (singularValues(l) *  VAdjMat.row(l)).adjoint())(0,0);
            }
        }
        frobeniusNorm = std::sqrt(std::real(blockNormSquared));
        return frobeniusNorm;
    }
    else
    {
        frobeniusNorm = fullMat.norm();
        return frobeniusNorm;
    }
}

bool BlockCluster::hasFourChildren() const
{
    return (this->son11 != nullptr) && (this->son12 != nullptr) && (this->son21 != nullptr) && (this->son22 != nullptr);
}

void BlockCluster::getPartition(const BlockCluster &block, QVector<const BlockCluster*> &partition)
{
    if(block.isLeaf) // block is in partition
    {
        if(block.fullMat.size() != 0 || block.singularValues.size() != 0)
        {
            #pragma omp critical
            partition.push_back(&block);
        }
    }
    else // block has leaves in blockclustertree
    {
        #pragma omp parallel sections
        {
            #pragma omp section
            {
                if(block.son11 != nullptr)
                {
                    BlockCluster::getPartition(* block.son11, partition);
                }
            }
            #pragma omp section
            {
                if(block.son12 != nullptr)
                {
                    BlockCluster::getPartition(* block.son12, partition);
                }
            }
            #pragma omp section
            {
                if(block.son21 != nullptr)
                {
                    BlockCluster::getPartition(* block.son21, partition);
                }
            }
            #pragma omp section
            {
                if(block.son22 != nullptr)
                {
                    BlockCluster::getPartition(* block.son22, partition);
                }
            }
        }
    }
}

std::unique_ptr<BlockCluster> BlockCluster::returnCopy(BlockCluster* father, long rank, double error, bool originalIsInSVDFormat) // returns a pointer to a copy of the entire subtree of this
{
    std::unique_ptr<BlockCluster> returnCluster = std::make_unique<BlockCluster>();

    returnCluster->isAdmissible = this->isAdmissible;
    returnCluster->isRoot = this->isRoot;
    returnCluster->isLeaf = this->isLeaf;

    returnCluster->rowCluster = this->rowCluster;
    returnCluster->columnCluster = this->columnCluster;

    returnCluster->frobeniusNorm = this->frobeniusNorm;

    returnCluster->fullMat = this->fullMat;
    returnCluster->UMat = this->UMat;
    returnCluster->singularValues = this->singularValues;
    returnCluster->VAdjMat = this->VAdjMat;

    if(returnCluster->isAdmissible && returnCluster->singularValues.size() > 0)
    {
        if(originalIsInSVDFormat)
        {
            returnCluster->compressSVDBlock(rank, error);
        }
        else
        {
            returnCluster->compressAdmissibleBlock(rank, error);
        }
    }

    if(father != nullptr)
    {
        returnCluster->father = father;
    }
    else
    {
        returnCluster->isRoot = true;
    }
    if(!this->isLeaf)
    {
        if(this->son11 != nullptr)
        {
            returnCluster->son11 = this->son11->returnCopy(returnCluster.get(), rank, error, originalIsInSVDFormat);
        }
        else
        {
            std::cerr << "!this->isLeaf && this->son11 == nullptr in returnCopy() call!" << std::endl;
        }
        if(this->son12 != nullptr)
        {
            returnCluster->son12 = this->son12->returnCopy(returnCluster.get(), rank, error, originalIsInSVDFormat);
        }
        else
        {
            std::cerr << "!this->isLeaf && this->son12 == nullptr in returnCopy() call!" << std::endl;
        }
        if(this->son21 != nullptr)
        {
            returnCluster->son21 = this->son21->returnCopy(returnCluster.get(), rank, error, originalIsInSVDFormat);
        }
        else
        {
            std::cerr << "!this->isLeaf && this->son21 == nullptr in returnCopy() call!" << std::endl;
        }
        if(this->son22 != nullptr)
        {
            returnCluster->son22 = this->son22->returnCopy(returnCluster.get(), rank, error, originalIsInSVDFormat);
        }
        else
        {
            std::cerr << "!this->isLeaf && this->son22 == nullptr in returnCopy() call!" << std::endl;
        }
    }
    return returnCluster;
}

Eigen::VectorXcd BlockCluster::row(long rowIndex) const
{
    if(isAdmissible)
    {
        return (UMat.row(rowIndex) * singularValues.asDiagonal()) * VAdjMat;
    }
    else
    {
        return fullMat.row(rowIndex);
    }
}

Eigen::VectorXcd BlockCluster::col(long colIndex) const
{
    if(isAdmissible)
    {
        return UMat * (singularValues.asDiagonal() * VAdjMat.col(colIndex));
    }
    else
    {
        return fullMat.col(colIndex);
    }
}

//Eigen::MatrixXcd BlockCluster::block(long rowIndex, long colIndex, long numberOfRows, long numberOfCols) const
//{
//    if(isAdmissible)
//    {
//        long rank = singularValues.size();
//        if(numberOfRows <= numberOfCols)
//        {
//            return (UMat.block(rowIndex, 0 , numberOfRows, rank) * singularValues.asDiagonal()) * VAdjMat.block(0, colIndex, rank, numberOfCols);
//        }
//        else
//        {
//            return UMat.block(rowIndex, 0 , numberOfRows, rank) * (singularValues.asDiagonal() * VAdjMat.block(0, colIndex, rank, numberOfCols));
//        }
//    }
//    else
//    {
//        return fullMat.block(rowIndex, colIndex, numberOfRows, numberOfCols);
//    }
//}

void BlockCluster::fullRankAssembly(std::function<std::complex<double> (long, long)> implicitMatrix)
{
    const long rowStartIndex = this->rowStartIndex();
    const long columnStartIndex = this->colStartIndex();
    const long blockRows = this->rows();
    const long blockColumns = this->cols();

    fullMat = Eigen::MatrixXcd(blockRows, blockColumns);

    #pragma omp parallel for
    for(long column = 0; column < blockColumns; column++)
    {
        long globalColumnIndex = columnStartIndex + column;
        for(long row = 0; row < blockRows; row++)
        {
            long globalRowIndex = rowStartIndex + row;
            fullMat(row, column) = implicitMatrix(globalRowIndex, globalColumnIndex);
//            if(isnan(std::abs(fullMat(row, column))))
//            {
//                std::cout<< row << " " << column << std::endl;
//                std::cout<< globalRowIndex << " " << globalColumnIndex << std::endl;
//            }
        }
    }
}

void BlockCluster::fullPivotACA(const long rank, const double relativeError, std::function<std::complex<double> (long, long)> implicitMatrix)
{
    VAdjMat.resize(0,0);
    singularValues.resize(0);
    UMat.resize(0,0);

    double norm = 0; // full rank matrix norm approximation by the R1 matrix norm multiplied by relative error
    double residuumNorm = 0;
    const long rowStartIndex = this->rowStartIndex();
    const long columnStartIndex = this->colStartIndex();
    const long blockRows = this->rows(); //assumes contiguous ascending indexes
    const long blockColumns = this->cols();

    long maxRank = std::min(blockRows, blockColumns);
    long initialRankReservation = 5;
    if(rank > 0)
    {
        maxRank = std::min(maxRank, rank);
        UMat.resize(blockRows, maxRank);
        singularValues.resize(maxRank);
        VAdjMat.resize(maxRank, blockColumns);
    }
    else
    {
        long reservationRank = std::min(initialRankReservation, maxRank);
        UMat.resize(blockRows, reservationRank);
        singularValues.resize(reservationRank);
        VAdjMat.resize(reservationRank, blockColumns);
    }
    long tmpRank = 0;
    Eigen::RowVectorXcd rowVector(blockColumns);
    Eigen::MatrixXcd Residuum(blockRows, blockColumns);
    for (int i = 0; i < blockRows ; i++)
    {
        calcHBlockRowFromImplicit(rowVector, rowStartIndex + i, columnStartIndex, implicitMatrix);
        Residuum.row(i) = rowVector;
    }
    norm = Residuum.norm();
    long rowIndex;
    long columnIndex;

    Eigen::VectorXcd columnVector(blockRows);
    while(tmpRank < maxRank)
    {
        Residuum.cwiseAbs().maxCoeff(&rowIndex, &columnIndex);

        calcHBlockRowFromImplicit(rowVector, rowStartIndex + rowIndex, columnStartIndex, implicitMatrix);

        if(tmpRank >= 1)
        {
            rowVector -= (UMat.leftCols(tmpRank).row(rowIndex) * singularValues.head(tmpRank).asDiagonal()) * VAdjMat.topRows(tmpRank);
        }

        calcHBlockColumnFromImplicit(columnVector, rowStartIndex, columnStartIndex + columnIndex,  implicitMatrix);

        std::complex<double> alpha = 1.0 / rowVector(columnIndex);

        if(tmpRank >= 1)
        {
            columnVector -= UMat.leftCols(tmpRank) * (singularValues.head(tmpRank).asDiagonal() * VAdjMat.topRows(tmpRank).col(columnIndex));
        }

        if(singularValues.size() <= tmpRank)
        {
            long newSize = std::min(2*tmpRank, maxRank);
            UMat.conservativeResize(Eigen::NoChange, newSize);
            singularValues.conservativeResize(newSize);
            VAdjMat.conservativeResize(newSize, Eigen::NoChange);
         }

        UMat.col(tmpRank) = columnVector;
        singularValues(tmpRank) = alpha;
        VAdjMat.row(tmpRank) = rowVector;
        Residuum -= UMat.col(tmpRank) * singularValues(tmpRank) * VAdjMat.row(tmpRank);
        tmpRank++; // increment rank counter to actual current rank of the low rank matrix

        residuumNorm = Residuum.norm();

        if(relativeError > 0 && residuumNorm / norm < relativeError)
        {
            break;
        }
    }
    UMat.conservativeResize(Eigen::NoChange, tmpRank);
    singularValues.conservativeResize(tmpRank);
    VAdjMat.conservativeResize(tmpRank, Eigen::NoChange);

    if(tmpRank == 0) // if nothing has been sampled
    {
         UMat = Eigen::MatrixXcd::Zero(blockRows, 1);
         singularValues = Eigen::VectorXcd::Zero(1);
         VAdjMat = Eigen::MatrixXcd::Zero(1, blockColumns);
    }
}

void BlockCluster::partialPivotACA(const long rank, const double relativeError, std::function<std::complex<double> (long, long)> implicitMatrix)
{
    VAdjMat.resize(0,0);
    singularValues.resize(0);
    UMat.resize(0,0);

    std::complex<double> normEstimate = 0; // full rank matrix norm approximation by the R1 matrix norm multiplied by relative error
    double relErrorTarget = 0; // full rank matrix norm approximation by the R1 matrix norm multiplied by relative error
    const long rowStartIndex = this->rowStartIndex();
    const long columnStartIndex = this->colStartIndex();
    const long blockRows = this->rows(); //assumes contiguous ascending indexes
    const long blockColumns = this->cols();

    long maxRank = std::min(blockRows, blockColumns);
    const long initialRankReservation = 5;
    if(rank > 0)
    {
        maxRank = std::min(maxRank, rank);
        UMat.resize(blockRows, maxRank);
        singularValues.resize(maxRank);
        VAdjMat.resize(maxRank, blockColumns);
    }
    else
    {
        long reservationRank = std::min(initialRankReservation, maxRank);
        UMat.resize(blockRows, reservationRank);
        singularValues.resize(reservationRank);
        VAdjMat.resize(reservationRank, blockColumns);
    }
    long tmpRank = 0;

    Eigen::RowVectorXcd rowVector(blockColumns);
    Eigen::VectorXcd columnVector(blockRows);

    long rowIndex;
    long columnIndex;

    columnIndex = QRandomGenerator::system()->bounded((qint32) blockColumns);
    calcHBlockColumnFromImplicit(columnVector, rowStartIndex, columnStartIndex + columnIndex,  implicitMatrix);

    while(tmpRank < maxRank)
    {
        if(tmpRank >= 1)
        {
            columnVector(rowIndex) = 0;
        }
        else //generate new random row index
        {
//            rowIndex = QRandomGenerator::global()->bounded((qint32) blockRows); // find random row index
            rowIndex = QRandomGenerator::system()->bounded((qint32) blockRows); // find random row index
        }
        columnVector.cwiseAbs().maxCoeff(&rowIndex);

        if(std::abs(columnVector(rowIndex)) <= global::tiny && tmpRank >= 1) //the max absolute element is too small -> no suitable column index for rowIndex
        {
            std::cerr << "max(columnVector) <= global::Tiny" << std::endl;
            break;    //change later
        }

        // calculate the block row rowIndex
        calcHBlockRowFromImplicit(rowVector, rowStartIndex + rowIndex, columnStartIndex, implicitMatrix);

        if(tmpRank >= 1)
        {
            rowVector -= (UMat.leftCols(tmpRank).row(rowIndex) * singularValues.head(tmpRank).asDiagonal()) * VAdjMat.topRows(tmpRank);
        }

         // find column index for largest absolute element in row with above rowindex
        double maxRowVal = rowVector.cwiseAbs().maxCoeff(&columnIndex); //  columnIndex is now max column index of phiAbsVec

        if(maxRowVal <= global::tiny && tmpRank >= 1) //the max absolute element is too small -> no suitable column index for rowIndex
        {
            std::cerr << "max(rowAbsVec) <= Tiny" << std::endl;
            break;    //change later
        }

        // calculate the block column columnIndex
        calcHBlockColumnFromImplicit(columnVector, rowStartIndex, columnStartIndex + columnIndex, implicitMatrix);

        std::complex<double> alpha = 1.0 / rowVector(columnIndex);

        if(tmpRank >= 1)
        {
            columnVector -= UMat.leftCols(tmpRank) * (singularValues.head(tmpRank).asDiagonal() * VAdjMat.topRows(tmpRank).col(columnIndex));
        }

        if(singularValues.size() <= tmpRank)
        {
            long newSize = std::min(2*tmpRank, maxRank);
            UMat.conservativeResize(Eigen::NoChange, newSize);
            singularValues.conservativeResize(newSize);
            VAdjMat.conservativeResize(newSize, Eigen::NoChange);
         }

        UMat.col(tmpRank) = columnVector;
        VAdjMat.row(tmpRank) = rowVector;
        singularValues(tmpRank) = alpha;
        tmpRank++; // increment rank counter to actual current rank of the low rank matrix

        normEstimate += std::abs(std::pow(alpha, 2)) *  rowVector.squaredNorm() * columnVector.squaredNorm();
        for(long i = 0; i < tmpRank - 1; i++)
        {
            normEstimate += 2.0 * std::conj(alpha) * singularValues(i) * (columnVector.adjoint() * UMat.col(i))(0,0) * (( VAdjMat.row(i)) * ( rowVector).adjoint())(0,0);
        }
        relErrorTarget = std::real(std::pow(relativeError, 2) * normEstimate);
        if(relativeError > 0 && std::pow(std::abs(alpha), 2) * rowVector.squaredNorm() * columnVector.squaredNorm() < relErrorTarget)
        {
            break;
        }
    }
    UMat.conservativeResize(Eigen::NoChange, tmpRank);
    singularValues.conservativeResize(tmpRank);
    VAdjMat.conservativeResize(tmpRank, Eigen::NoChange);

    frobeniusNorm = std::sqrt(std::abs(normEstimate));

    if(tmpRank == 0) // if the first pivoting element was already do small -> nothing has been sampled
    {
         UMat = Eigen::MatrixXcd::Zero(blockRows, 1);
         singularValues = Eigen::VectorXcd::Zero(1);
         VAdjMat = Eigen::MatrixXcd::Zero(1, blockColumns);
    }
}

void BlockCluster::partialPivotACAextra(const long rank, const double relativeError, std::function<QVector<std::pair<long,long>>(BlockCluster*,double)> getPivotIndices, std::function<std::complex<double> (long, long)> implicitMatrix) /*!< Low-rank assembly of an h-block by ACA with heuristic partial pivoting. The (guaranteed) accuracy is improved for dPhi-blocks. */
{
    VAdjMat.resize(0,0);
    singularValues.resize(0);
    UMat.resize(0,0);

    std::complex<double> normEstimate = 0; // full rank matrix norm approximation by the R1 matrix norm multiplied by relative error
    double relErrorTarget = 0; // full rank matrix norm approximation by the R1 matrix norm multiplied by relative error
    const long rowStartIndex = this->rowStartIndex();
    const long columnStartIndex = this->colStartIndex();
    const long blockRows = rows(); //assumes contiguous ascending indexes
    const long blockColumns = cols();
    long maxRank = std::min(blockRows, blockColumns);
    const long initialRankReservation = 5;
    if(rank > 0)
    {
        maxRank = std::min(maxRank, rank);
        UMat.resize(blockRows, maxRank);
        singularValues.resize(maxRank);
        VAdjMat.resize(maxRank, blockColumns);
    }
    else
    {
        long reservationRank = std::min(initialRankReservation, maxRank);
        UMat.resize(blockRows, reservationRank);
        singularValues.resize(reservationRank);
        VAdjMat.resize(reservationRank, blockColumns);
    }

    QVector<std::pair<long,long>> pivotIndices = getPivotIndices(this, 0.44);

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
        if(std::abs((pivotElement - ((UMat.row(rowIndex).head(tmpRank).transpose().array() * singularValues.head(tmpRank).array()).array() * VAdjMat.col(columnIndex).head(tmpRank).array()).sum()) / pivotElement) < relativeError)
        {
            continue; // pivot element is already sufficiently approximated
        }
        calcHBlockColumnFromImplicit(columnVector, rowStartIndex, columnStartIndex + columnIndex, implicitMatrix);
        if(tmpRank >= 1)
        {
            columnVector -= UMat.leftCols(tmpRank) * (singularValues.head(tmpRank).asDiagonal() * VAdjMat.topRows(tmpRank).col(columnIndex));
        }
        bool firstIterationForCurrentPivotElement = true;
        while(tmpRank < maxRank)
        {
            Eigen::VectorXd colAbsVec;
            colAbsVec  = columnVector.cwiseAbs();
            if(!firstIterationForCurrentPivotElement)
            {
                colAbsVec(rowIndex) = 0; //at this point rowIndex is still the number from the last iteration
            }
            colAbsVec.maxCoeff(&rowIndex); //  columnIndex is now max column index of phiAbsVec
            firstIterationForCurrentPivotElement = false;
            if(std::abs(columnVector(rowIndex)) <= global::tiny && tmpRank >= 1) //the max absolute element is too small -> no suitable column index for rowIndex
            {
                 std::cerr << "max(columnVector) <= global::Tiny" << std::endl;
                 break;    //change later
            }

            // calculate the block row
            calcHBlockRowFromImplicit(rowVector, rowStartIndex + rowIndex, columnStartIndex, implicitMatrix);

            if(tmpRank >= 1)
            {
                rowVector -= (UMat.leftCols(tmpRank).row(rowIndex) * singularValues.head(tmpRank).asDiagonal()) * VAdjMat.topRows(tmpRank);
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
            calcHBlockColumnFromImplicit(columnVector, rowStartIndex, columnStartIndex + columnIndex, implicitMatrix);

            std::complex<double> alpha = 1.0 / rowVector(columnIndex);

            if(tmpRank >= 1)
            {
                columnVector -= UMat.leftCols(tmpRank) * (singularValues.head(tmpRank).asDiagonal() * VAdjMat.topRows(tmpRank).col(columnIndex));
            }

            if(singularValues.size() <= tmpRank)
            {
                long newSize = std::min(2*tmpRank, maxRank);
                UMat.conservativeResize(Eigen::NoChange, newSize);
                singularValues.conservativeResize(newSize);
                VAdjMat.conservativeResize(newSize, Eigen::NoChange);
             }

            UMat.col(tmpRank) = columnVector;
            VAdjMat.row(tmpRank) = rowVector;
            singularValues(tmpRank) = alpha;
            tmpRank++; // increment rank counter to actual current rank of the low rank matrix

            std::complex<double> oldNorm = normEstimate;

            normEstimate += std::abs(std::pow(alpha, 2)) * rowVector.squaredNorm() * columnVector.squaredNorm();
            for(long i = 0; i < tmpRank - 1; i++)
            {
                normEstimate += 2.0 * std::conj(alpha) * singularValues(i) * (columnVector.adjoint() * UMat.col(i))(0,0) * (( VAdjMat.row(i)) * ( rowVector).adjoint())(0,0);
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
    UMat.conservativeResize(Eigen::NoChange, tmpRank);
    singularValues.conservativeResize(tmpRank);
    VAdjMat.conservativeResize(tmpRank, Eigen::NoChange);

    frobeniusNorm = std::sqrt(std::abs(normEstimate));

    if(tmpRank == 0)
    {
         UMat = Eigen::MatrixXcd::Zero(blockRows, 1);
         singularValues = Eigen::VectorXcd::Zero(1);
         VAdjMat = Eigen::MatrixXcd::Zero(1, blockColumns);
    }
}

void BlockCluster::calcHBlockRowFromImplicit(Eigen::RowVectorXcd &rowVector, const long rowIndex, const long columnStartIndex, std::function<std::complex<double> (long, long)> implicitMatrix)
{
    #pragma omp parallel for // parallelizes vector initialization
    for(long i = 0; i < rowVector.size(); i++)   //assemble rowVectors
    {
        rowVector(i) = implicitMatrix(rowIndex, columnStartIndex + i);
    }
}

void BlockCluster::calcHBlockColumnFromImplicit(Eigen::VectorXcd &columnVector, const long rowStartIndex, const long columnIndex, std::function<std::complex<double> (long, long)> implicitMatrix)
{
    #pragma omp parallel for // parallelizes vector initialization
    for(long i = 0; i < columnVector.size(); i++) //assemble columnVectors
    {
        columnVector(i) = implicitMatrix(rowStartIndex + i, columnIndex);
    }
}

//BlockClusterTree::BlockClusterTree()
//{

//}

HMatrix::HMatrix(const HMatrix &original, long rank, double error, bool originalIsInSVDFormat)
{
    this->rowClustertree = original.rowClustertree;
    this->columnClustertree = original.columnClustertree;
    this->rootBlock = original.getRootBlock()->returnCopy(nullptr, rank, error, originalIsInSVDFormat);
    this->updatePartition();
}

void HMatrix::setClusterTrees(std::shared_ptr<ClusterTree> rowClustertree, std::shared_ptr<ClusterTree> columnClustertree)
{
    this->rowClustertree = rowClustertree;
    this->columnClustertree = columnClustertree;
}

QVector<BlockCluster*> HMatrix::getMinPartition()
{
    if(rootBlock == nullptr)
    {
        minPartition.clear();
    }
    return minPartition;
}

void HMatrix::assembleBlocks(const long rank, double relativeError, std::function<std::complex<double> (long, long)> implicitMatrix, bool updatePartition)
{
    if(updatePartition)
    {
        this->updatePartition(); // update the partition vector
    }
    #pragma omp parallel master
    for(long i=0; i<minPartition.size(); i++)
    {
        #pragma omp task
        {
            if(minPartition.at(i)->isAdmissible)
            {
                minPartition.at(i)->partialPivotACA(rank, relativeError, implicitMatrix); // low rank assembly via adaptive cross approximation
            }
            else
            {
                minPartition.at(i)->fullRankAssembly(implicitMatrix); // full rank assembly
            }
            if(minPartition.at(i)->hasNan())
            {
                std::cerr << "NaNs in block " << i << std::endl;
                std::cerr << "Is Block admissible? " << minPartition.at(i)->isAdmissible << std::endl;
                std::cerr << "rowStartIndex " << minPartition.at(i)->rowStartIndex() << std::endl;
                std::cerr << "rows " << minPartition.at(i)->rows() << std::endl;
                std::cerr << "colStartIndex " << minPartition.at(i)->colStartIndex() << std::endl;
                std::cerr << "cols " << minPartition.at(i)->cols() << std::endl;
            }
        }
    }
}

void HMatrix::assembleBlocksExtra(const long rank, double relativeError, /*std::function<std::complex<double> (long, long)> implicitMatrixFull, */ std::function<QVector<std::pair<long,long>>(BlockCluster*,double)> getPivotIndices, std::function<std::complex<double> (long, long)> implicitMatrixLowRank, bool updatePartition)
{
    if(updatePartition)
    {
        this->updatePartition(); // update the partition vector
    }
    #pragma omp parallel master
    for(long i=0; i<minPartition.size(); i++)
    {
        #pragma omp task
        {
            if(minPartition.at(i)->isAdmissible)
            {
                minPartition.at(i)->partialPivotACAextra(rank, relativeError, getPivotIndices, implicitMatrixLowRank); // low rank assembly via adaptive cross approximation with extra pivot indices
            }
            else
            {
                minPartition.at(i)->fullRankAssembly(implicitMatrixLowRank); // full rank assembly
            }
        }
    }
}

QVector<BlockCluster*> HMatrix::getDiagonalBlocks() const
{
    QVector<BlockCluster*> diagonalBlocks;
    if(rootBlock != nullptr)
    {
        recursiveAppendDiagonalBlock(*rootBlock, diagonalBlocks);
    }
    return diagonalBlocks;
}

void HMatrix::clear(bool callMalloc_trim)
{
    if(rootBlock != nullptr)
    {
        rootBlock->clear();
        rootBlock = nullptr;
        minPartition.clear();
    }
    if(callMalloc_trim)
    {
        global::trimMemory();
    }
}

void HMatrix::recursiveClear(BlockCluster &blockClusterNode)
{
//    std::cerr << "recursiveClear" << std::endl;
    if(blockClusterNode.son11 != nullptr)
    {
        recursiveClear(*blockClusterNode.son11);
    }
    if(blockClusterNode.son12 != nullptr)
    {
        recursiveClear(*blockClusterNode.son12);
    }
    if(blockClusterNode.son21 != nullptr)
    {
        recursiveClear(*blockClusterNode.son21);
    }
    if(blockClusterNode.son22 != nullptr)
    {
        recursiveClear(*blockClusterNode.son22);
    }
    blockClusterNode.son11 = nullptr;
    blockClusterNode.son12 = nullptr;
    blockClusterNode.son21 = nullptr;
    blockClusterNode.son22 = nullptr;

    blockClusterNode.clearMatrixInfo();

//    delete blockClusterNode;
//    blockClusterNode = nullptr;
}

void HMatrix::setRootBlock(std::unique_ptr<BlockCluster> &root)
{
    rootBlock = std::move(root);
}

BlockCluster* HMatrix::getRootBlock() const
{
    if(rootBlock == nullptr)
    {
        std::cerr << "getRootBlock() called on empty BlockClusterTree" << std::endl;
    }
    return rootBlock.get();
}

//void HMatrix::setRootBlock(BlockCluster* rootBlockCluster)
//{
//    this->rootBlockCluster = rootBlockCluster;
//    this->rootBlockCluster->isRoot = true;
//}

void HMatrix::populateBlockClusterTree(bool randBlocks)
{
//    std::cout<<"populateBlockClusterTree called!"<<std::endl;
//    delete rootBlockCluster;
    rootBlock = nullptr;

    createRootNode();
//    std::cout<<"RootBlockClusterNode created!"<<std::endl;

    createBlockClusterChildren(*rootBlock);
//    std::cout<<"BlockClusterTree populated!"<<std::endl;

    if(randBlocks)
    {
        setUpRandomMatrices(maxRank);
        std::cout << "Memory for Matrix Blocks alloccated." << std::endl;
    }
}

void HMatrix::createRootNode()
{
    minPartition.clear();

    if(rowClustertree->getRootCluster() == nullptr || columnClustertree->getRootCluster() == nullptr)
    {
        std::cerr << "nullptr returned by getRootCluster() in createRootNode call. rootBlockCluster set to nullptr." << std::endl;
        rootBlock = nullptr;
        return;
    }
    Cluster* rootCluster1 = rowClustertree->getRootCluster();
    Cluster* rootCluster2 = columnClustertree->getRootCluster();
//    rootBlockCluster = new BlockCluster(rootCluster1, rootCluster2);
    rootBlock = std::make_unique<BlockCluster>(rootCluster1, rootCluster2);
    rootBlock->isRoot = true;
}

void HMatrix::createBlockClusterChildren(BlockCluster &block)
{
    if(isBlockAdmissible(block) /*&& (!blockClusterNode->rowCluster->isLeaf && !blockClusterNode->columnCluster->isLeaf)*/)
    {
        block.isAdmissible = true;
        block.isLeaf = true;

        block.son11 = nullptr;
        block.son12 = nullptr;
        block.son21 = nullptr;
        block.son22 = nullptr;

        minPartition.append(&block);  //TreeLeaves constitute the Blockpartition

        return;
    }
    else if(block.rowCluster->isLeaf || block.columnCluster->isLeaf)
    {
        block.isAdmissible = false;
        block.isLeaf = true;

        block.son11 = nullptr;
        block.son12 = nullptr;
        block.son21 = nullptr;
        block.son22 = nullptr;

        minPartition.append(&block);  //TreeLeaves constitute the Blockpartition

        return;
    }
    else
    {
        block.isAdmissible = false;
        block.isLeaf = false;

        block.son11 = std::make_unique<BlockCluster>(block.rowCluster->son1.get(), block.columnCluster->son1.get(), &block);
        createBlockClusterChildren(*block.son11);

        block.son12 = std::make_unique<BlockCluster>(block.rowCluster->son1.get(), block.columnCluster->son2.get(), &block);
        createBlockClusterChildren(*block.son12);

        block.son21 = std::make_unique<BlockCluster>(block.rowCluster->son2.get(), block.columnCluster->son1.get(), &block);
        createBlockClusterChildren(*block.son21);

        block.son22 = std::make_unique<BlockCluster>(block.rowCluster->son2.get(), block.columnCluster->son2.get(), &block);
        createBlockClusterChildren(*block.son22);

        return;
    }
}

void HMatrix::setUpRandomMatrices(long maxRank)
{
    for(int i = 0; i < minPartition.length(); i++)
    {
        if(minPartition.at(i)->rowCluster == nullptr || minPartition.at(i)->columnCluster == nullptr)
        {
            std::cerr << "cluster == nullpointr for block in setUpMatrices call." << std::endl;
            return;
        }
        else
        {
            long blockRows = minPartition.at(i)->rowCluster->indices.last() - minPartition.at(i)->rowCluster->indices.first() + 1; //assumes contiguous ascending indexes
            long blockColumns = minPartition.at(i)->columnCluster->indices.last() - minPartition.at(i)->columnCluster->indices.first() + 1;
            if(minPartition.at(i)->isAdmissible)
            {
                Eigen::HouseholderQR<Eigen::MatrixXcd> qrA(Eigen::MatrixXcd::Random(blockRows, maxRank));
                Eigen::HouseholderQR<Eigen::MatrixXcd> qrB_T(Eigen::MatrixXcd::Random(maxRank, blockColumns).transpose());

                long localARank = std::min(maxRank, blockRows);
                long localBRank = std::min(maxRank, blockColumns);
                long rank = std::min(localARank, localBRank);
                Eigen::MatrixXcd R_A = qrA.matrixQR().triangularView<Eigen::Upper>().toDenseMatrix().topRows(localARank);
                Eigen::MatrixXcd R_B_Transpose = qrB_T.matrixQR().triangularView<Eigen::Upper>().transpose().toDenseMatrix().leftCols(localBRank);

                Eigen::MatrixXcd QA = qrA.householderQ() * Eigen::MatrixXcd::Identity(blockRows, localARank);//.leftCols(localARank);
                Eigen::MatrixXcd QB = qrB_T.householderQ() * Eigen::MatrixXcd::Identity(blockColumns, localBRank);//.leftCols(localBRank);
                Eigen::BDCSVD<Eigen::MatrixXcd,Eigen::ComputeThinU|Eigen::ComputeThinV> svd( R_A * R_B_Transpose );

                Eigen::MatrixXcd U = QA * svd.matrixU().leftCols(rank);
                Eigen::MatrixXcd V = svd.matrixV().adjoint().topRows(rank) * QB.transpose();

                minPartition[i]->UMat = U;
                minPartition[i]->singularValues = svd.singularValues().head(rank);
                minPartition[i]->VAdjMat = V;
            }
            else
            {
                minPartition[i]->fullMat = Eigen::MatrixXcd::Random(blockRows, blockColumns);
            }
        }
    }
}

bool HMatrix::isBlockAdmissible(const BlockCluster &block)
{
    double clusterDist = clusterDistance(block.rowCluster, block.columnCluster);
    if(clusterDist == 0)
    {
        return false;
    }
    else if(std::min(clusterSize(block.rowCluster), clusterSize(block.columnCluster)) / clusterDist < admissibilityConstant)
//    else if(std::max(clusterSize(blockCluster->rowCluster), clusterSize(blockCluster->columnCluster)) / clusterDist < admissibilityConstant)
    {
        return true;
    }
    else
    {
        return false;
    }
}

double HMatrix::clusterDistance(Cluster* cluster1, Cluster* cluster2)
{
    Eigen::Vector3d minP1 = cluster1->minCuboid.minPoint;
    Eigen::Vector3d maxP1 = cluster1->minCuboid.maxPoint;
    Eigen::Vector3d minP2 = cluster2->minCuboid.minPoint;
    Eigen::Vector3d maxP2 = cluster2->minCuboid.maxPoint;

    bool xOverlap = (minP1.x() < maxP2.x() && maxP1.x() > minP2.x()) || (minP2.x() < maxP1.x() && maxP2.x() > minP1.x());
    bool yOverlap = (minP1.y() < maxP2.y() && maxP1.y() > minP2.y()) || (minP2.y() < maxP1.y() && maxP2.y() > minP1.y());
    bool zOverlap = (minP1.z() < maxP2.z() && maxP1.z() > minP2.z()) || (minP2.z() < maxP1.z() && maxP2.z() > minP1.z());
    int numberOfOverlaps = xOverlap + yOverlap + zOverlap;

    if(numberOfOverlaps == 3)   // quboids overlap
    {
        return 0;
    }
    else if(numberOfOverlaps == 2) //min distance is between faces
    {
        if(!xOverlap)
        {
            return std::min(std::abs(minP1.x() - maxP2.x()), std::abs(minP2.x() - maxP1.x()));
        }
        else if(!yOverlap)
        {
            return std::min(std::abs(minP1.y() - maxP2.y()), std::abs(minP2.y() - maxP1.y()));
        }
        else
        {
            return std::min(std::abs(minP1.z() - maxP2.z()), std::abs(minP2.z() - maxP1.z()));
        }
    }
    else if(numberOfOverlaps == 1) //min distance between lines
    {
        if(xOverlap)
        {
            double yDis = std::min(std::abs(minP1.y() - maxP2.y()), std::abs(minP2.y() - maxP1.y()));
            double zDis = std::min(std::abs(minP1.z() - maxP2.z()), std::abs(minP2.z() - maxP1.z()));
            return sqrt(yDis*yDis + zDis*zDis);
        }
        else if(yOverlap)
        {
            double xDis = std::min(std::abs(minP1.x() - maxP2.x()), std::abs(minP2.x() - maxP1.x()));
            double zDis = std::min(std::abs(minP1.z() - maxP2.z()), std::abs(minP2.z() - maxP1.z()));
            return sqrt(xDis*xDis + zDis*zDis);
        }
        else
        {
            double xDis = std::min(std::abs(minP1.x() - maxP2.x()), std::abs(minP2.x() - maxP1.x()));
            double yDis = std::min(std::abs(minP1.y() - maxP2.y()), std::abs(minP2.y() - maxP1.y()));
            return sqrt(xDis*xDis + yDis*yDis);
        }
    }
    else //min distance is between corners of cuboids
    {
        double xDis = std::min(std::abs(minP1.x() - maxP2.x()), std::abs(minP2.x() - maxP1.x()));
        double yDis = std::min(std::abs(minP1.y() - maxP2.y()), std::abs(minP2.y() - maxP1.y()));
        double zDis = std::min(std::abs(minP1.z() - maxP2.z()), std::abs(minP2.z() - maxP1.z()));
        return sqrt(xDis*xDis + yDis*yDis + zDis*zDis);
    }
}

double HMatrix::clusterSize(Cluster* cluster) const
{
    Eigen::Vector3d cuboidSideLengths = cluster->minCuboid.maxPoint - cluster->minCuboid.minPoint;
    return cuboidSideLengths.norm();
}

bool HMatrix::isCrossPartitioned() const
{
    bool onlyCrossPartitions = true; // bool value returned true if all non leaf blocks have four sons

    checkBlockPartition(*rootBlock, onlyCrossPartitions);
    return onlyCrossPartitions;
}

void HMatrix::checkBlockPartition(BlockCluster &block, bool &onlyCrossPartitions) const
{
    if(onlyCrossPartitions == false)
    {
        return;
    }
    else
    {
        if(block.isLeaf == true) //block is not partitioned ->Y break
        {
            return;
        }
        else if(block.son11 != nullptr && block.son12 != nullptr && block.son21 != nullptr && block.son22 != nullptr)
        {  //block is cross partitioned -> check partition of all sons
            checkBlockPartition(*block.son11, onlyCrossPartitions);
            checkBlockPartition(*block.son12, onlyCrossPartitions);
            checkBlockPartition(*block.son21, onlyCrossPartitions);
            checkBlockPartition(*block.son22, onlyCrossPartitions);
            return;
        }
        else // block is not cross partitioned
        {
            onlyCrossPartitions = false;
            return;
        }
    }
}

void HMatrix::updatePartition()
{
    minPartition.clear();
    farFieldBlocks = 0;
    nearFieldBlocks = 0;
    if(rootBlock != nullptr)
    {
        recursiveAppendPartitionBlock(*rootBlock);
    }
    else
    {
        std::cerr<<"updatePartition() called on empty BlockClusterTree." <<std::endl;
    }
}

void HMatrix::recursiveAppendPartitionBlock(BlockCluster &block)
{
    if(block.isLeaf)
    {
        minPartition.append(&block);
        if(block.isAdmissible)
        {
            farFieldBlocks++;
        }
        else
        {
            nearFieldBlocks++;
        }
    }
    else if(block.son11 != nullptr)
    {
        recursiveAppendPartitionBlock(*block.son11);
        recursiveAppendPartitionBlock(*block.son12);
        recursiveAppendPartitionBlock(*block.son21);
        recursiveAppendPartitionBlock(*block.son22);
    }
    else
    {
        std::cerr << "Block with non-leaf designation but without children in recursiveAppendPartitionBlock() call. " << std::endl;
        block.isLeaf = true;
    }
}

void HMatrix::updatePartitionWithNullptrCheck()
{
    minPartition.clear();
    farFieldBlocks = 0;
    nearFieldBlocks = 0;
    int depth = 0;
    if(rootBlock != nullptr)
    {
        recAppendPartitionBlockNullptrCheck(*rootBlock, depth);
    }
    else
    {
        std::cerr<<"updatePartition() called on empty BlockClusterTree." <<std::endl;
    }
}

void HMatrix::recAppendPartitionBlockNullptrCheck(BlockCluster &block, int depth)
{
    ++depth;
//    std::cout<<"depth: "<<depth<<std::endl;
    if(block.isLeaf)
    {
        minPartition.append(&block);
        if(block.son11 == nullptr)
        {
            if(block.isAdmissible)
            {
                farFieldBlocks++;
            }
            else
            {
                nearFieldBlocks++;
            }
        }
    }
    if(block.son11 != nullptr)
    {
        recAppendPartitionBlockNullptrCheck(*block.son11, depth);
        recAppendPartitionBlockNullptrCheck(*block.son12, depth);
        recAppendPartitionBlockNullptrCheck(*block.son21, depth);
        recAppendPartitionBlockNullptrCheck(*block.son22, depth);
    }
}

void HMatrix::updatePartitionWithMatrixInfo()
{
    minPartition.clear();
    farFieldBlocks = 0;
    nearFieldBlocks = 0;
    if(rootBlock!=nullptr)
    {
        QVector<BlockCluster*> partition;
        recAppendPartitionBlockMatInfo(*rootBlock, partition);
        minPartition = partition;
    }
    else
    {
        std::cerr<<"updatePartition() called on empty BlockClusterTree." <<std::endl;
    }
}

void HMatrix::recAppendPartitionBlockMatInfo(BlockCluster &block, QVector<BlockCluster*> &subBlocksWithInformation)
{
    if(block.UMat.cols() != 0)
    {
        if(block.fullMat.cols() != 0)
        {
            std::cout << "blockCluster.AMat.cols() != 0 && blockCluster.fullMat.cols() != 0" << std::endl;
        }

        block.isAdmissible = true;
        block.isLeaf = true;
        subBlocksWithInformation.append(&block);
//        ++farFieldBlocks;
    }
    else if(block.fullMat.cols() != 0)
    {
        block.isAdmissible = false;
        block.isLeaf = true;
        subBlocksWithInformation.append(&block);
//        ++nearFieldBlocks;
    }
    else if(block.hasFourChildren())
    {
        block.isLeaf = false;
    }
    else
    {
        block.isLeaf = true;
    }

    if(block.son11 != nullptr)
    {
        recAppendPartitionBlockMatInfo(*block.son11, subBlocksWithInformation);
    }
    if(block.son12 != nullptr)
    {
        recAppendPartitionBlockMatInfo(*block.son12, subBlocksWithInformation);
    }
    if(block.son21 != nullptr)
    {
        recAppendPartitionBlockMatInfo(*block.son21, subBlocksWithInformation);
    }
    if(block.son22 != nullptr)
    {
        recAppendPartitionBlockMatInfo(*block.son22, subBlocksWithInformation);
    }
}

void HMatrix::recursiveAppendDiagonalBlock(BlockCluster &block, QVector<BlockCluster*> &diagonalBlocks) const
{
    if(block.son11 != nullptr && block.son22 != nullptr)
    {
        HMatrix::recursiveAppendDiagonalBlock(*block.son11, diagonalBlocks);
        HMatrix::recursiveAppendDiagonalBlock(*block.son22, diagonalBlocks);
    }
    else
    {
        if(block.isLeaf == true)
        {
            diagonalBlocks.append(&block);
        }
        else
        {
            std::cerr << " blockCluster.son11 == nullptr && blockCluster.son22 == nullptr and  blockCluster.isLeaf == false in recursiveAppendDiagonalBlock call." << std::endl;
        }
    }
}


QVector<BlockCluster*> HMatrix::getSubBlocksWithMatrixInformation(BlockCluster &block)
{
    QVector<BlockCluster*> subBlocksWithInformation;

    if(block.son11 != nullptr)
    {
        recAppendPartitionBlockMatInfo(*block.son11, subBlocksWithInformation);
    }
    if(block.son12 != nullptr)
    {
        recAppendPartitionBlockMatInfo(*block.son12, subBlocksWithInformation);
    }
    if(block.son21 != nullptr)
    {
        recAppendPartitionBlockMatInfo(*block.son21, subBlocksWithInformation);
    }
    if(block.son22 != nullptr)
    {
        recAppendPartitionBlockMatInfo(*block.son22, subBlocksWithInformation);
    }
    return subBlocksWithInformation;
}

long HMatrix::rows() const
{
    if(rootBlock != nullptr)
    {
        return rootBlock->rows();
    }
    else
    {
        return 0;
    }
}

long HMatrix::cols() const
{
    if(rootBlock != nullptr)
    {
        return rootBlock->cols();
    }
    else
    {
        return 0;
    }
}

long HMatrix::size() const
{
    return rows() * cols();
}

long HMatrix::rowStartIndex() const
{
    if(rootBlock != nullptr)
    {
        return rootBlock->rowStartIndex();
    }
    else
    {
        return 0;
    }
}

long HMatrix::colStartIndex() const
{
    if(rootBlock != nullptr)
    {
        return rootBlock->colStartIndex();
    }
    else
    {
        return 0;
    }
}

double HMatrix::getCompressionRatio()
{
    long fullMatMemory = rows() * cols();
    long hMatMemory = 0;
    updatePartition();

    for(long partitionIndex = 0 ; partitionIndex < minPartition.length(); partitionIndex++)
    {
        BlockCluster* blockcluster = minPartition.at(partitionIndex);

        if(blockcluster -> isAdmissible) //farfield block with Rk-Matrix
        {
            hMatMemory += (blockcluster->rows() + blockcluster->cols() + 1) * blockcluster->singularValues.size();
        }
        else //nearfield block with full matrix
        {
            hMatMemory += blockcluster->rows() * blockcluster->cols();
        }
    }
    return ((double) fullMatMemory) / hMatMemory;
}

bool HMatrix::consistencyCheck() const
{
    if(rootBlock != nullptr)
    {
        return recursiveConsistencyCheck(*rootBlock);
    }
    else
    {
        return false;
    }
}

bool HMatrix::recursiveConsistencyCheck(BlockCluster &block) const
{
    bool useCerr = true;
    if(block.isLeaf)
    {
        if(block.isAdmissible)
        {
            if(block.fullMat.size() > 0 || block.UMat.size() <= 0)
            {
                if(useCerr)
                {
                    std::cerr << "blockCluster is admissible;  blockCluster.fullMat.size() > 0 || blockCluster.UMat.size() <= 0" << std::endl;
                }
                return false;
            }
            if(block.UMat.rows() != block.rows() || block.VAdjMat.cols() != block.cols())
            {
                if(useCerr)
                {
                    std::cerr << "blockCluster is admissible;  blockCluster.UMat.rows() != blockCluster.rows() || blockCluster.VAdjMat.cols() != blockCluster.cols()" << std::endl;
                }
                return false;
            }
            if(block.UMat.cols() != block.singularValues.size() || block.VAdjMat.rows() != block.singularValues.size())
            {
                if(useCerr)
                {
                    std::cerr << "blockCluster is admissible;  blockCluster.UMat.cols() != blockCluster.singularValues.size() || blockCluster.VAdjMat.rows() != blockCluster.singularValues.size()" << std::endl;
                }
                return false;
            }
        }
        else    // !blockCluster.isAdmissible
        {
            if(block.fullMat.size() <= 0 || block.UMat.size() > 0)
            {
                if(useCerr)
                {
                    std::cerr << "blockCluster is not admissible;  blockCluster.fullMat.size() <= 0 || blockCluster.UMat.size() > 0" << std::endl;
                }
                return false;
            }
            if(block.fullMat.cols() != block.cols() || block.fullMat.rows() != block.rows())
            {
                if(useCerr)
                {
                    std::cerr << "blockCluster is not admissible;  blockCluster.fullMat.cols() != blockCluster.cols() || blockCluster.fullMat.rows() != blockCluster.rows()" << std::endl;
                }
                return false;
            }
        }
        if(block.son11 != nullptr || block.son12 != nullptr || block.son21 != nullptr || block.son22 != nullptr)
        {
            if(useCerr)
            {
                std::cerr << "blockCluster is leaf;  blockCluster.son11 != nullptr || blockCluster.son12 != nullptr || blockCluster.son21 != nullptr || blockCluster.son22 != nullptr" << std::endl;
            }
            return false;
        }
    }
    else // !blockCluster.isLeaf
    {
        if( !(block.son11 != nullptr && block.son12 != nullptr && block.son21 != nullptr && block.son22 != nullptr) )
        {
            if(useCerr)
            {
                std::cerr << "blockCluster is not leaf;  !(blockCluster.son11 == nullptr && blockCluster.son12 == nullptr && blockCluster.son21 == nullptr && blockCluster.son22 == nullptr)" << std::endl;
            }
            return false;
        }
        if(block.fullMat.size() > 0 || block.UMat.size() > 0)
        {
            if(useCerr)
            {
                std::cerr << "blockCluster is not leaf;  blockCluster.fullMat.size() > 0 || blockCluster.UMat.size() > 0" << std::endl;
            }
            return false;
        }
        if(block.rowStartIndex() != block.son11->rowStartIndex() || block.colStartIndex() != block.son11->colStartIndex())
        {
            if(useCerr)
            {
                std::cerr << "blockCluster.rowStartindex() != blockCluster.son11->rowStartindex() || blockCluster.colStartindex() != blockCluster.son11->colStartindex()" << std::endl;
            }
            return false;
        }
        if((block.rowStartIndex() + block.son11->rows()) != block.son21->rowStartIndex() || (block.colStartIndex() + block.son11->cols()) != block.son12->colStartIndex())
        {
            if(useCerr)
            {
                std::cerr << "(blockCluster.rowStartindex() + blockCluster.son11->rows()) != blockCluster.son21->rowStartindex() || (blockCluster.colStartindex() + blockCluster.son11->cols()) != blockCluster.son12->colStartindex()" << std::endl;
            }
            return false;
        }
        if((block.rowStartIndex() + block.son12->rows()) != block.son22->rowStartIndex() || (block.colStartIndex() + block.son21->cols()) != block.son22->colStartIndex())
        {
            if(useCerr)
            {
                std::cerr << "(blockCluster.rowStartindex() + blockCluster.son12->rows()) != blockCluster.son22->rowStartindex() || (blockCluster.colStartindex() + blockCluster.son21->cols()) != blockCluster.son22->colStartindex()" << std::endl;
            }
            return false;
        }
        if(block.rows() != block.son11->rows() + block.son21->rows() || block.rows() != block.son12->rows() + block.son22->rows())
        {
            if(useCerr)
            {
                std::cerr << "blockCluster.rows() != blockCluster.son11->rows() + blockCluster.son21->rows() || blockCluster.rows() != blockCluster.son12->rows() + blockCluster.son22->rows()" << std::endl;
            }
            return false;
        }
        if(block.cols() != block.son11->cols() + block.son12->cols() || block.cols() != block.son21->cols() + block.son22->cols())
        {
            if(useCerr)
            {
                std::cerr << "blockCluster.cols() != blockCluster.son11->cols() + blockCluster.son12->cols() || blockCluster.cols() != blockCluster.son21->cols() + blockCluster.son22->cols()" << std::endl;
            }
            return false;
        }
    }

    if(block.son11 != nullptr)
    {
        if(recursiveConsistencyCheck(*block.son11) == false)
        {
            return false;
        }
    }

    if(block.son12 != nullptr)
    {
        if(recursiveConsistencyCheck(*block.son12) == false)
        {
            return false;
        }
    }

    if(block.son21 != nullptr)
    {
        if(recursiveConsistencyCheck(*block.son21) == false)
        {
            return false;
        }
    }
    if(block.son22 != nullptr)
    {
        if(recursiveConsistencyCheck(*block.son22) == false)
        {
            return false;
        }
    }
    return true;
}

long HMatrix::calculateMaxBlockRank()
{
    long maxBlockRank = 0;
    for( int i = 0; i < minPartition.size(); i++)
    {
        if(minPartition.at(i)->isAdmissible)
        {
            long localBlockRank = minPartition.at(i)->singularValues.size();
            if(localBlockRank > maxBlockRank)
            {
                maxBlockRank = localBlockRank;
            }
        }
    }
    return maxBlockRank;
}

long HMatrix::calculateSizeOfAdmissibleBlocks()
{
    long size = 0;
    for( int i = 0; i < minPartition.size(); i++)
    {
        if(minPartition.at(i)->isAdmissible)
        {
            size += minPartition.at(i)->size();
        }
    }
    return size;
}

double HMatrix::norm(bool calcMinPartition)
{
    if(calcMinPartition)
    {
        updatePartition();
    }
    double norm = 0;
    #pragma omp parallel for
    for( int i = 0; i < minPartition.size(); i++)
    {
        if(minPartition.at(i)->isAdmissible)
        {
            std::complex<double> blockNorm = 0;
            for(int l = 0; l < minPartition.at(i)->singularValues.size(); l++)
            {
                blockNorm += std::abs(std::pow(minPartition.at(i)->singularValues(l), 2)) * minPartition.at(i)->UMat.col(l).squaredNorm() * minPartition.at(i)->VAdjMat.row(l).squaredNorm();
                for(long j = 0; j < l; j++)
                {
                    blockNorm += 2.0 * (minPartition.at(i)->UMat.col(l).adjoint() * minPartition.at(i)->UMat.col(j))(0,0) * ((minPartition.at(i)->singularValues(j) * minPartition.at(i)->VAdjMat.row(j)) * (minPartition.at(i)->singularValues(l) *  minPartition.at(i)->VAdjMat.row(l)).adjoint())(0,0);
                }
            }
            #pragma omp critical
            norm += std::real(blockNorm);
            minPartition[i]->frobeniusNorm = std::sqrt(std::real(blockNorm));
//            std::cerr << minPartition[i]->frobeniusNorm - (minPartition[i]->UMat * minPartition[i]->singularValues.asDiagonal() * minPartition[i]->VAdjMat).norm() << std::endl;
        }
        else
        {
            double blockNorm = minPartition.at(i)->fullMat.squaredNorm();
            minPartition[i]->frobeniusNorm = std::sqrt(blockNorm);
            #pragma omp critical
            norm += blockNorm;
        }
    }
    return std::sqrt(norm);
}

double HMatrix::normInSVDForm(bool calcMinPartition)
{
    if(calcMinPartition)
    {
        updatePartition();
    }
    double norm = 0;
    #pragma omp parallel for
    for( int i = 0; i < minPartition.size(); i++)
    {
        if(minPartition.at(i)->isAdmissible)
        {
            #pragma omp critical
            norm += minPartition.at(i)->singularValues.squaredNorm();
        }
        else
        {
            #pragma omp critical
            norm += minPartition.at(i)->fullMat.squaredNorm();
        }
    }
    return std::sqrt(norm);
}

bool HMatrix::hasNan(bool calcMinPartition)
{
    if(!calcMinPartition)
    {
        updatePartition();
    }
    bool hasNan = false;
    #pragma omp parallel for
    for(long i=0; i<minPartition.size(); i++)
    {
        if(!hasNan && minPartition[i]->hasNan())
        {
            hasNan = true;
        }
    }
    return hasNan;
}

Eigen::MatrixXcd HMatrix::toFullMat() //test function; extremely expensive for large matrices O(n^2)
{
    const long totalNumberOfRows = this->rows();
    const long totalNumberOfColumns = this->cols();
    const long globalRowStartIndex = this->rowStartIndex(); // The H-matrix may be a subblock of another h-matrix
    const long globalColumnStartIndex = this->colStartIndex(); // -> the corresponding the row- and col-indices therefore may not start at zero

    if(this->minPartition.size() == 0)
    {
//        this->updatePartition();
        this->updatePartitionWithMatrixInfo();
    }
    if(this->minPartition.size() == 0)
    {
        std::cerr << "Empty partition in toFullMat call." << std::endl;
    }

    Eigen::MatrixXcd returnMatrix = Eigen::MatrixXcd(totalNumberOfRows, totalNumberOfColumns);

    #pragma omp parallel for
    for(long i = 0; i < this->minPartition.length(); i++)
    {
        const BlockCluster* tmpBlock = this->minPartition.at(i);
        const long rowStartIndex = tmpBlock->rowStartIndex();
        const long numberOfRows = tmpBlock->rows();
        const long columnStartIndex = tmpBlock->colStartIndex();
        const long numberOfColumns = tmpBlock->cols();

        if(tmpBlock->isAdmissible) // far field block
        {
            if(numberOfRows < numberOfColumns)
            {
                returnMatrix.block(rowStartIndex - globalRowStartIndex, columnStartIndex - globalColumnStartIndex, numberOfRows, numberOfColumns).noalias() = (tmpBlock->UMat * tmpBlock->singularValues.asDiagonal()) * tmpBlock->VAdjMat;
            }
            else
            {
                returnMatrix.block(rowStartIndex - globalRowStartIndex, columnStartIndex - globalColumnStartIndex, numberOfRows, numberOfColumns).noalias() = tmpBlock->UMat * (tmpBlock->singularValues.asDiagonal() * tmpBlock->VAdjMat);
            }
        }
        else
        {
//            std::cerr << rowStartIndex - globalRowStartIndex << " " << columnStartIndex - globalColumnStartIndex << std::endl;
            returnMatrix.block(rowStartIndex - globalRowStartIndex, columnStartIndex - globalColumnStartIndex, numberOfRows, numberOfColumns) = tmpBlock->fullMat;
        }
    }
    return returnMatrix;
}
