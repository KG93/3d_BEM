#include "harithm.h"

void HArithm::MVM(Eigen::VectorXcd &y, HMatrix &hmatrix, const Eigen::VectorXcd &x)
{
    QVector<BlockCluster*> minPartition = hmatrix.getMinPartition();
    #pragma omp parallel for
    for(long i = 0; i < minPartition.length(); i++)
    {
        BlockCluster* tmpBlock = minPartition.at(i);
        long rowStartIndex = tmpBlock->rowStartIndex();
        long numberOfRows = tmpBlock->rows();
        long columnStartIndex = tmpBlock->colStartIndex();
        long numberOfColumns = tmpBlock->cols();
        Eigen::VectorXcd tmp;
        if(tmpBlock->isAdmissible) // far field block
        {
            tmp.noalias() = tmpBlock->UMat * (tmpBlock->singularValues.asDiagonal() * (tmpBlock->VAdjMat * x.segment(columnStartIndex, numberOfColumns)));
            #pragma omp critical
            y.segment(rowStartIndex, numberOfRows) += tmp;
        }
        else // near field block
        {
            tmp.noalias() = tmpBlock->fullMat * x.segment(columnStartIndex, numberOfColumns);
            #pragma omp critical
            y.segment(rowStartIndex, numberOfRows) += tmp;
        }
    }
}

void HArithm::VMM(Eigen::RowVectorXcd &y, const Eigen::RowVectorXcd &x, HMatrix &hmatrix) // vector-H-matrix product. y' += x' * hmatrix
{
    QVector<BlockCluster*> minPartition = hmatrix.getMinPartition();
    #pragma omp parallel for
    for(long i = 0; i < minPartition.length(); i++)
    {
        BlockCluster* tmpBlock = minPartition.at(i);
        long rowStartIndex = tmpBlock->rowStartIndex();
        long numberOfRows = tmpBlock->rows();
        long columnStartIndex = tmpBlock->colStartIndex();
        long numberOfColumns = tmpBlock->cols();
        Eigen::RowVectorXcd tmp;
        if(tmpBlock->isAdmissible) // far field block
        {
            tmp.noalias() = ((x.segment(rowStartIndex, numberOfRows) * tmpBlock->UMat) * tmpBlock->singularValues.asDiagonal()) * tmpBlock->VAdjMat;
            #pragma omp critical
            y.segment(columnStartIndex, numberOfColumns) += tmp;
        }
        else // near field block
        {
            tmp.noalias() = x.segment(rowStartIndex, numberOfRows) * tmpBlock->fullMat;
            #pragma omp critical
            y.segment(columnStartIndex, numberOfColumns) += tmp;
        }
    }
}

void HArithm::compressHMat(HMatrix &hmatrix, const long maxRank, const double error)
{
    hmatrix.updatePartition();
    QVector<BlockCluster*> partition = hmatrix.getMinPartition();

    #pragma omp parallel for
//    #pragma omp parallel master
    for(long i = 0; i < partition.length(); i++)
    {
//        #pragma omp task
        {
            BlockCluster* tmpBlock = partition.at(i);
            if(tmpBlock->isAdmissible)
            {
    //            #pragma omp task
                RkMatRankReduction(*tmpBlock, maxRank, error);
            }
        }
    }
}

long HArithm::minRankforError(const Eigen::BDCSVD<Eigen::MatrixXcd> &svd, const long maxRank, const double relError)
{
    const long nonzeroSingularValues = svd.nonzeroSingularValues();
    if(relError <= 0)
    {
        if(maxRank <= 0)
        {
            return nonzeroSingularValues;
        }
        else
        {
            return std::min(svd.nonzeroSingularValues(), maxRank);
        }
    }
    else
    {
        long maxIndex = nonzeroSingularValues;
        if(maxRank > 0)
        {
            maxIndex = std::min(maxRank, nonzeroSingularValues);
        }
        double frobeniousNorm2FullMatrix = svd.singularValues().squaredNorm();
        double frobeniousNorm2RkMatrix = svd.singularValues().squaredNorm();
        double errorSquared = std::pow(relError, 2);
        for(long optimalRank = 1; optimalRank <= maxIndex; optimalRank++) //for absolute error
//        for(long optimalRank = 1; optimalRank <= maxIndex-1; optimalRank++) //for relative error
        {
//            if(svd.singularValues()(optimalRank - 1) < errorSquared)  //for absolute error
//            if(svd.singularValues()(optimalRank ) < error * svd.singularValues()(0)) //for relative error
            frobeniousNorm2RkMatrix -= std::pow(svd.singularValues()(optimalRank - 1), 2);
            if(frobeniousNorm2RkMatrix / frobeniousNorm2FullMatrix < errorSquared) //for relative error in frobebious norm
            {
                return optimalRank;
            }
        }
        return maxIndex;
    }
}

void HArithm::svdOnBlock(BlockCluster &block)
{
    std::cerr << "svdOnBlock" << std::endl;
    long localARank = block.UMat.rows();
    long localBRank = block.VAdjMat.cols();

    Eigen::HouseholderQR<Eigen::MatrixXcd> qrA(block.UMat /** matrixBlock.singularValues.asDiagonal()*/);
    Eigen::HouseholderQR<Eigen::MatrixXcd> qrB_T(block.VAdjMat.transpose());
    Eigen::MatrixXcd R_A = qrA.matrixQR().triangularView<Eigen::Upper>().toDenseMatrix().topRows(localARank);
    Eigen::MatrixXcd R_B_Transpose = qrB_T.matrixQR().triangularView<Eigen::Upper>().transpose().toDenseMatrix().leftCols(localBRank);

    Eigen::MatrixXcd QA = qrA.householderQ() * Eigen::MatrixXcd::Identity(block.UMat.rows(), localARank);//.leftCols(localARank);
    Eigen::MatrixXcd QB = qrB_T.householderQ() * Eigen::MatrixXcd::Identity(block.VAdjMat.cols(), localBRank);//.leftCols(localBRank);


    if constexpr(useEigenSvd)
    {
//        Eigen::JacobiSVD<Eigen::MatrixXcd> svd( (R_A * matrixBlock.singularValues.head(R_A.cols()).asDiagonal()) * R_B_Transpose.triangularView<Eigen::Lower>());
        Eigen::BDCSVD<Eigen::MatrixXcd> svd( (R_A * block.singularValues.head(R_A.cols()).asDiagonal()) * R_B_Transpose.triangularView<Eigen::Lower>(), Eigen::ComputeThinU | Eigen::ComputeThinV);

        block.frobeniusNorm = svd.singularValues().norm();
        block.UMat.noalias() = QA * svd.matrixU();
        block.singularValues = svd.singularValues();
        block.VAdjMat.noalias() = svd.matrixV().adjoint() * QB.transpose();
    }
    else
    {
        Eigen::MatrixXcd A = (R_A * block.singularValues.head(R_A.cols()).asDiagonal()) * R_B_Transpose.triangularView<Eigen::Lower>();
        Eigen::MatrixXcd U;
        Eigen::VectorXd singVals;
        Eigen::MatrixXcd VAdj;

        GoloubReinschSvd::goloubReinschSVD(A, U, singVals , VAdj);

        block.UMat.noalias() = QA * U;
        block.singularValues = singVals;
        block.VAdjMat.noalias() = VAdj * QB.transpose();
        block.frobeniusNorm = singVals.norm();
    }
}

void HArithm::RkMatRankReduction(BlockCluster &block, long rank, const double relError)
{
    long localRank = block.singularValues.size();

    if((localRank <= rank && relError == 0) || (rank <= 0 && relError <= 0)) //might be wrong with rk mat in svd form
    {
        return;
    }
    else    //block is far field
    {
        long localARank = std::min(localRank, block.UMat.rows());
        long localBRank = std::min(localRank, block.VAdjMat.cols());

        Eigen::HouseholderQR<Eigen::MatrixXcd> qrA(block.UMat);
        Eigen::HouseholderQR<Eigen::MatrixXcd> qrB_T(block.VAdjMat.transpose());

        Eigen::MatrixXcd R_A = qrA.matrixQR().triangularView<Eigen::Upper>().toDenseMatrix().topRows(localARank);
        Eigen::MatrixXcd R_B_Transpose = qrB_T.matrixQR().triangularView<Eigen::Upper>().transpose().toDenseMatrix().leftCols(localBRank);

        if constexpr(useEigenSvd)
        {
//            Eigen::JacobiSVD<Eigen::MatrixXcd,Eigen::ComputeThinU|Eigen::ComputeThinV> svd( (R_A * matrixBlock.singularValues.head(R_A.cols()).asDiagonal()) * R_B_Transpose.triangularView<Eigen::Lower>());
            Eigen::BDCSVD<Eigen::MatrixXcd> svd( (R_A * block.singularValues.head(R_A.cols()).asDiagonal()) * R_B_Transpose.triangularView<Eigen::Lower>(), Eigen::ComputeThinU | Eigen::ComputeThinV);
            rank = minRankforError(svd, rank, relError);
            block.frobeniusNorm = svd.singularValues().head(rank).norm();
            Eigen::MatrixXcd tmpForU(block.rows(), rank);
            tmpForU << svd.matrixU().leftCols(rank), Eigen::MatrixXcd::Zero(block.rows() - localARank, rank);
            block.UMat = qrA.householderQ() * tmpForU;
            block.singularValues = svd.singularValues().head(rank);
            Eigen::MatrixXcd tmpForVAdj(rank, block.cols());
            tmpForVAdj << svd.matrixV().adjoint().topRows(rank), Eigen::MatrixXcd::Zero(rank, block.cols() - localBRank);
            block.VAdjMat = tmpForVAdj * qrB_T.householderQ().transpose();
        }
        else
        {
//            Eigen::MatrixXcd full = matrixBlock.UMat * matrixBlock.singularValues.asDiagonal() * matrixBlock.VAdjMat;

            Eigen::MatrixXcd A = (R_A * block.singularValues.head(R_A.cols()).asDiagonal()) * R_B_Transpose.triangularView<Eigen::Lower>();

            Eigen::MatrixXcd tmpU;
            Eigen::VectorXd tmpSingularValues;
            Eigen::MatrixXcd tmpVAdj;

            GoloubReinschSvd::goloubReinschSVD(A, tmpU, tmpSingularValues , tmpVAdj);

            rank = GoloubReinschSvd::minRankforError(tmpSingularValues, rank, relError);

            tmpSingularValues.conservativeResize(rank);

            Eigen::MatrixXcd tmpForU(block.rows(), rank);
            tmpForU << tmpU.leftCols(rank), Eigen::MatrixXcd::Zero(block.rows() - localARank, rank);
            block.UMat.noalias() = qrA.householderQ() * tmpForU;

            block.singularValues = tmpSingularValues;

            Eigen::MatrixXcd tmpForVAdj(rank, block.cols());
            tmpForVAdj << tmpVAdj.topRows(rank), Eigen::MatrixXcd::Zero(rank, block.cols() - localBRank);

            block.VAdjMat.noalias() = tmpForVAdj * qrB_T.householderQ().transpose();
            block.frobeniusNorm = tmpSingularValues.norm();

//            double relErrorReal = (full- matrixBlock.UMat * matrixBlock.singularValues.asDiagonal() * matrixBlock.VAdjMat).norm() / full.norm();
//            if(relErrorReal > relError)
//            {
//                std::cerr << "svd rel error: " << relErrorReal << std::endl;
//            }
        }
        block.isAdmissible = true;
    }
}

void HArithm::fullMatRankReduction(BlockCluster &block, const long rank, const double relError)
{
    bool rKtoFullFirst = true;

    if(rKtoFullFirst && block.UMat.size() != 0) //convert existing rk matrix and add to full matrix
    {
        addrkMatToFull(block.fullMat, block.UMat, block.singularValues, block.VAdjMat);
        clearRkMatrix(block);
    }
    if(!block.isAdmissible || block.fullMat.size() != 0) //block is far field
    {        
        if constexpr(useEigenSvd)
        {
//            Eigen::JacobiSVD<Eigen::MatrixXcd,Eigen::ComputeThinU|Eigen::ComputeThinV> svd( matrixBlock.fullMat);
            Eigen::BDCSVD<Eigen::MatrixXcd> svd( block.fullMat, Eigen::ComputeThinU | Eigen::ComputeThinV);
            long localRank = minRankforError(svd, rank, relError);

            block.fullMat.resize(0, 0);
            block.isAdmissible = true;
            horizontalJoinMatricesInToFirstOne(block.UMat, svd.matrixU().leftCols(localRank));
            joinVectorsInToFirstOne(block.singularValues, svd.singularValues().head(localRank));
            verticalJoinMatricesInToFirstOne(block.VAdjMat, svd.matrixV().adjoint().topRows(localRank));
            block.frobeniusNorm = svd.singularValues().head(localRank).norm();
        }
        else    // use own implementation of the golub-reinsch svd algorithm
        {
            Eigen::MatrixXcd U;
            Eigen::VectorXd singVals;
            Eigen::MatrixXcd VAdj;
            GoloubReinschSvd::goloubReinschSVD(block.fullMat, U, singVals , VAdj); // goloubReinschSVD overwrites matrixBlock.fullMat
            long localRank = GoloubReinschSvd::minRankforError(singVals, rank, relError);
            block.fullMat.resize(0, 0);
            block.isAdmissible = true;

            horizontalJoinMatricesInToFirstOne(block.UMat, U.leftCols(localRank));
            joinVectorsInToFirstOne(block.singularValues, singVals.head(localRank));
            verticalJoinMatricesInToFirstOne(block.VAdjMat, VAdj.topRows(localRank));
            block.frobeniusNorm = singVals.head(localRank).norm();
        }
    }
    else
    {
        std::cerr << "fullMatRankReduction called on not-near-field block";
    }

    if(!rKtoFullFirst) //alternative method; collect rk matrices first and then reduce rank
    {
        RkMatRankReduction(block, rank, relError);
    }
}

void HArithm::addHMat2ToHMat1SamePartition(HMatrix &matrix1, HMatrix &matrix2, const long rank, const double relError)
{
    QVector<BlockCluster*> partition = matrix1.getMinPartition();
    QVector<BlockCluster*> partition2 = matrix2.getMinPartition();
    if(partition.length() != partition2.length())
    {
        std::cerr << "Unequal partition sizes in addHMat2ToHMat1() call!" << std::endl;
        return;
    }

    #pragma omp parallel for
    for(int i = 0; i < partition.length(); i++)
    {
        if(partition.at(i)->isAdmissible) //low rank addition for far field blocks
        {
            horizontalJoinMatricesInToFirstOne(partition.at(i)->UMat, partition2.at(i)->UMat);
            joinVectorsInToFirstOne(partition.at(i)->singularValues, partition2.at(i)->singularValues);
            verticalJoinMatricesInToFirstOne(partition.at(i)->VAdjMat, partition2.at(i)->VAdjMat);
        }
        else //full rank addition for near field blocks
        {
            long rows1 = partition.at(i)->fullMat.rows();
            long rows2 = partition2.at(i)->fullMat.rows();
            long columns1 = partition.at(i)->fullMat.cols();
            long columns2 = partition2.at(i)->fullMat.cols();

            if(rows1 != rows2 || columns1 != columns2)
            {
                std::cerr << "Incompatible matrix dimensions in full-mat addition in addHMat2ToHMat1SamePartition() call! " <<rows1<<"x"<< columns1<<", " <<rows2<<"x"<< columns2<<std::endl;
            }

            partition.at(i)->fullMat.noalias() += partition2.at(i)->fullMat;
        }
    }

    #pragma omp parallel for
    for(long i = 0; i < partition.length(); i++) //reduce rank of farfield blocks
    {

        if(partition.at(i)->isAdmissible)
        {
            HArithm::RkMatRankReduction(*partition.at(i), rank, relError); //also takes care of norm calculation
        }
    }
}

void HArithm::matByBlock(const Eigen::MatrixXcd &lFactor, const BlockCluster &rFactorBlock, long columnOffsetFactor1, long columnOffsetFactor2, Eigen::MatrixXcd &product)
{
    if(rFactorBlock.isLeaf)
    {
        long rowStartIndex = rFactorBlock.rowStartIndex();
        long numberOfRows = rFactorBlock.rows();
        long columnStartIndex = rFactorBlock.colStartIndex();
        long numberOfColumns = rFactorBlock.cols();
        if(rFactorBlock.isAdmissible) // far field block
        {
            product.block(0, columnStartIndex-columnOffsetFactor2, lFactor.rows(), numberOfColumns).noalias() += (lFactor.block(0, rowStartIndex-columnOffsetFactor1, lFactor.rows(), numberOfRows) * rFactorBlock.UMat )* rFactorBlock.singularValues.asDiagonal() * rFactorBlock.VAdjMat ;
        }
        else // near field block
        {
            product.block(0, columnStartIndex-columnOffsetFactor2, lFactor.rows(), numberOfColumns).noalias() += lFactor.block(0, rowStartIndex-columnOffsetFactor1, lFactor.rows(), numberOfRows) * rFactorBlock.fullMat;
        }
    }
    else
    {
        matByBlock(lFactor, * rFactorBlock.son11, columnOffsetFactor1, columnOffsetFactor2, product);
        matByBlock(lFactor, * rFactorBlock.son12, columnOffsetFactor1, columnOffsetFactor2, product);
        matByBlock(lFactor, * rFactorBlock.son21, columnOffsetFactor1, columnOffsetFactor2, product);
        matByBlock(lFactor, * rFactorBlock.son22, columnOffsetFactor1, columnOffsetFactor2, product);
    }
}

void HArithm::blockByMat(const BlockCluster &lFactorBlock, const Eigen::MatrixXcd &rFactor, long rowOffsetFactor1, long rowOffsetFactor2, Eigen::MatrixXcd &product)
{
    if(lFactorBlock.isLeaf)
    {
        long rowStartIndex = lFactorBlock.rowStartIndex();
        long numberOfRows = lFactorBlock.rows();
        long columnStartIndex = lFactorBlock.colStartIndex();
        long numberOfColumns = lFactorBlock.cols();
        if(lFactorBlock.isAdmissible) // far field block
        {
            product.block(rowStartIndex-rowOffsetFactor1, 0, numberOfRows, rFactor.cols()).noalias() += lFactorBlock.UMat * lFactorBlock.singularValues.asDiagonal() * (lFactorBlock.VAdjMat * rFactor.block(columnStartIndex-rowOffsetFactor2,0,numberOfColumns,rFactor.cols()));
        }
        else // near field block
        {
            product.block(rowStartIndex-rowOffsetFactor1, 0, numberOfRows, rFactor.cols()).noalias() += lFactorBlock.fullMat * rFactor.block(columnStartIndex-rowOffsetFactor2,0,numberOfColumns,rFactor.cols());
        }
    }
    else
    {
        blockByMat(* lFactorBlock.son11, rFactor, rowOffsetFactor1, rowOffsetFactor2, product);
        blockByMat(* lFactorBlock.son12, rFactor, rowOffsetFactor1, rowOffsetFactor2, product);
        blockByMat(* lFactorBlock.son21, rFactor, rowOffsetFactor1, rowOffsetFactor2, product);
        blockByMat(* lFactorBlock.son22, rFactor, rowOffsetFactor1, rowOffsetFactor2, product);
    }
}

void HArithm::roundedAddRMatToBlock(BlockCluster& block, const Eigen::MatrixXcd &UMat, const Eigen::VectorXcd &singVals, const Eigen::MatrixXcd &VAdjMat, const long maxRank, const double relError)
{
    long initialBlockRank = block.singularValues.size();
    addRMatToBlock(block, UMat, singVals, VAdjMat);
    if(initialBlockRank > 0)
    {
        RkMatRankReduction(block, maxRank, relError);
    }
}

void HArithm::addRMatToBlock(BlockCluster& block, const Eigen::MatrixXcd &UMat, const Eigen::VectorXcd &singVals, const Eigen::MatrixXcd &VAdjMat)
{
    horizontalJoinMatricesInToFirstOne(block.UMat, UMat);
    joinVectorsInToFirstOne(block.singularValues, singVals);
    verticalJoinMatricesInToFirstOne(block.VAdjMat, VAdjMat);
}

void HArithm::horizontalJoinMatricesInToFirstOne(Eigen::MatrixXcd &firstMat, const Eigen::MatrixXcd &secondMat)
{
    long rows1 = firstMat.rows();
    long rows2 = secondMat.rows();
    long columns1 = firstMat.cols();
    long columns2 = secondMat.cols();

    if( rows1 != rows2 && ((rows1 != 0) && (rows2 != 0)) )
    {
        std::cerr << "Incompatible matrix dimensions in r-mat addition in horizontalJoinMatricesInToFirstOne() call! " <<rows1<<"x"<< columns1<<", " <<rows2<<"x"<< columns2<<std::endl;
    }
    if(rows1 == 0)
    {
        firstMat = secondMat;
        return;
    }
    if(rows2 == 0)
    {
        return;
    }

    Eigen::MatrixXcd joinMatrix(rows1, columns1 + columns2);
    joinMatrix << firstMat, secondMat;
    firstMat = joinMatrix;
}

void HArithm::verticalJoinMatricesInToFirstOne(Eigen::MatrixXcd &firstMat, const Eigen::MatrixXcd &secondMat)
{

    long rows1 = firstMat.rows();
    long rows2 = secondMat.rows();
    long columns1 = firstMat.cols();
    long columns2 = secondMat.cols();


    if( columns1 != columns2 && ((rows1 !=0) && (rows2 !=0)) )
    {
        std::cerr << "Incompatible matrix dimensions in r-mat addition in verticalJoinMatricesInToFirstOne() call! " <<rows1<<"x"<< columns1<<", " <<rows2<<"x"<< columns2<<std::endl;
    }
    if(rows1 == 0)
    {
        firstMat = secondMat;
        return;
    }
    if(rows2 == 0)
    {
        return;
    }

    Eigen::MatrixXcd joinMatrix = Eigen::MatrixXcd(rows1 + rows2, columns1);
    joinMatrix << firstMat, secondMat;
    firstMat = joinMatrix;
}

void HArithm::joinVectorsInToFirstOne(Eigen::VectorXcd &firstVector, const Eigen::VectorXcd &secondVector)
{
    long size1 = firstVector.size();
    long size2 = secondVector.size();

    if(size1 == 0)
    {
        firstVector = secondVector;
        return;
    }
    if(size2 == 0)
    {
        return;
    }

    Eigen::VectorXcd joinVector = Eigen::VectorXcd(size1 + size2);
    joinVector << firstVector, secondVector;
    firstVector = joinVector;
}

void HArithm::addFullMatrixInToFirstOne(Eigen::MatrixXcd &firstMat, const Eigen::MatrixXcd &secondMat)
{
    if(firstMat.size() != 0)
    {
        firstMat += secondMat;
    }
    else
    {
        firstMat = secondMat;
    }
}

void HArithm::addrkMatToFull(Eigen::MatrixXcd &fullMat, const Eigen::MatrixXcd &UMat, const Eigen::VectorXcd &singVals, const Eigen::MatrixXcd &VAdjMat) /*!< Efficiently add low rank matrix into a full matrix. */
{
    long rows1 = fullMat.rows();
    long rows2 = UMat.rows();
    long columns1 = fullMat.cols();
    long columns2 = VAdjMat.cols();

    if( (columns1 != columns2 || (rows1 != rows2)) && ((rows1 !=0) && (rows2 !=0)) )
    {
        std::cerr << "Incompatible matrix dimensions in addrkMatToFull call! " <<rows1<<"x"<< columns1<<", " <<rows2<<"x"<< columns2<<std::endl;
    }
    if(fullMat.size() != 0)
    {
        if(UMat.rows() < VAdjMat.cols())
        {
            fullMat.noalias() += (UMat * singVals.asDiagonal()) * VAdjMat;
        }
        else
        {
            fullMat.noalias() += UMat * (singVals.asDiagonal() * VAdjMat);
        }
    }
    else
    {
        if(UMat.rows() < VAdjMat.cols())
        {
            fullMat.noalias() = (UMat * singVals.asDiagonal()) * VAdjMat;
        }
        else
        {
            fullMat.noalias() = UMat * (singVals.asDiagonal() * VAdjMat);
        }
    }
}

void HArithm::transportSubblockInToFullAncestorBlock(BlockCluster &ancestorBlock, BlockCluster &subBlock)
{
    long rows = ancestorBlock.rows();
    long cols = ancestorBlock.cols();

    if(subBlock.fullMat.size() != 0 || subBlock.UMat.size() != 0)
    {
        if(subBlock.UMat.size() != 0)
        {
            addrkMatToFull(subBlock.fullMat, subBlock.UMat, subBlock.singularValues, subBlock.VAdjMat);
        }
        if(ancestorBlock.fullMat.size() == 0)
        {
            ancestorBlock.fullMat = Eigen::MatrixXcd::Zero(rows, cols);
        }

        long subRows = subBlock.rows();
        long subCols = subBlock.cols();
        long subBlockRelativeRowStartIndex = subBlock.rowStartIndex() - ancestorBlock.rowStartIndex();
        long subBlockRelativeColStartIndex = subBlock.colStartIndex() - ancestorBlock.colStartIndex();

        #pragma omp critical
        ancestorBlock.fullMat.block(subBlockRelativeRowStartIndex, subBlockRelativeColStartIndex, subRows, subCols) += subBlock.fullMat;

        subBlock.clearMatrixInfo();
        ancestorBlock . isLeaf = true;
        ancestorBlock . isAdmissible = false;
    }
    subBlock . isLeaf = false;
}

void HArithm::setRkMatrixZero(BlockCluster &block)
{
    block.UMat = Eigen::MatrixXcd::Zero(block.rows(), 1);
    block.singularValues = Eigen::VectorXcd::Zero(1);
    block.VAdjMat = Eigen::MatrixXcd::Zero(1, block.cols());
}

void HArithm::clearRkMatrix(BlockCluster &block)
{
    block.UMat.resize(0,0);
    block.singularValues.resize(0);
    block.VAdjMat.resize(0,0);
}

void HArithm::fullRankDirektAgglomerationWithTrimming(BlockCluster &block)
{
    fullRankDirektAgglomeration(block);
    block.trimBelow();
}

void HArithm::fullRankDirektAgglomeration(BlockCluster &block)
{
    if(block.UMat.size() != 0)
    {
        addrkMatToFull(block.fullMat, block.UMat, block.singularValues, block.VAdjMat);
        clearRkMatrix(block);
    }
    QVector<BlockCluster*> subBlocksWithInformation = HMatrix::getSubBlocksWithMatrixInformation(block);
    #pragma omp parallel for
    for(long blockIndex = 0; blockIndex < subBlocksWithInformation.length(); blockIndex++)
    {
        // add subblocks to block.fullmat
        transportSubblockInToFullAncestorBlock(block, * subBlocksWithInformation.at(blockIndex));
    }
    block.isAdmissible = false;
    if(block.fullMat.size() != 0)
    {
        block . isLeaf = true;
    }
}

void HArithm::fullRankAgglomerationWithTrimming(BlockCluster &block)
{
    fullRankAgglomeration(block);
    block.trimBelow();
}


void HArithm::fullRankAgglomeration(BlockCluster &block) //agglomerates all subblocks into the full matrix of the argument block
{
    long rowIndexStart = block.rowStartIndex();
    long numberOfRows = block.rows();
    long columnIndexStart = block.colStartIndex();
    long numberOfColumns = block.cols();


    if(block.UMat.size() != 0)
    {
        addrkMatToFull(block.fullMat, block.UMat, block.singularValues, block.VAdjMat);
        clearRkMatrix(block);
    }
    if(block.fullMat.size() == 0)
    {
        block.fullMat = Eigen::MatrixXcd::Zero(numberOfRows, numberOfColumns);
    }
    block.isLeaf = true;
    block.isAdmissible = false;

    if(block.son11 != nullptr)
    {
        long sonRowIndexStart = block.son11->rowCluster->indices.first();
        long sonNumberOfRows = block.son11->rows();
        long sonColumnIndexStart = block.son11->columnCluster->indices.first();
        long sonNumberOfColumns = block.son11->cols();
        block.fullMat.block(sonRowIndexStart - rowIndexStart, sonColumnIndexStart - columnIndexStart, sonNumberOfRows, sonNumberOfColumns) += fullRankAgglomerationRecursion(* block.son11);
    }
    if(block.son12 != nullptr)
    {
        long sonRowIndexStart = block.son12->rowCluster->indices.first();
        long sonNumberOfRows = block.son12->rows();
        long sonColumnIndexStart = block.son12->columnCluster->indices.first();
        long sonNumberOfColumns = block.son12->cols();
        block.fullMat.block(sonRowIndexStart - rowIndexStart, sonColumnIndexStart - columnIndexStart, sonNumberOfRows, sonNumberOfColumns) += fullRankAgglomerationRecursion(* block.son12);
    }
    if(block.son21 != nullptr)
    {
        long sonRowIndexStart = block.son21->rowCluster->indices.first();
        long sonNumberOfRows = block.son21->rows();
        long sonColumnIndexStart = block.son21->columnCluster->indices.first();
        long sonNumberOfColumns = block.son21->cols();
        block.fullMat.block(sonRowIndexStart - rowIndexStart, sonColumnIndexStart - columnIndexStart, sonNumberOfRows, sonNumberOfColumns) += fullRankAgglomerationRecursion(* block.son21);
    }
    if(block.son22 != nullptr)
    {
        long sonRowIndexStart = block.son22->rowCluster->indices.first();
        long sonNumberOfRows = block.son22->rows();
        long sonColumnIndexStart = block.son22->columnCluster->indices.first();
        long sonNumberOfColumns = block.son22->cols();
        block.fullMat.block(sonRowIndexStart - rowIndexStart, sonColumnIndexStart - columnIndexStart, sonNumberOfRows, sonNumberOfColumns) += fullRankAgglomerationRecursion(* block.son22);
    }
}

Eigen::MatrixXcd HArithm::fullRankAgglomerationRecursion(BlockCluster &block) //agglomerates all subblocks into the full return matrix
{
    Eigen::MatrixXcd son11Matrix, son12Matrix, son21Matrix, son22Matrix;
    if(block.son11 != nullptr)
    {
        son11Matrix = fullRankAgglomerationRecursion(* block.son11);
    }
    if(block.son12 != nullptr)
    {
        son12Matrix = fullRankAgglomerationRecursion(* block.son12);
    }
    if(block.son21 != nullptr)
    {
        son21Matrix = fullRankAgglomerationRecursion(* block.son21);
    }
    if(block.son22 != nullptr)
    {
       son22Matrix = fullRankAgglomerationRecursion(* block.son22);
    }

    long rowIndexStart = block.rowStartIndex();
    long numberOfRows = block.rows();
    long columnIndexStart = block.colStartIndex();
    long numberOfColumns = block.cols();
    Eigen::MatrixXcd returnMatrix = Eigen::MatrixXcd::Zero(numberOfRows, numberOfColumns);

    if(block.son11 != nullptr)
    {
        long sonRowIndexStart = block.son11->rowStartIndex();
        long sonNumberOfRows = block.son11->rows();
        long sonColumnIndexStart = block.son11->colStartIndex();
        long sonNumberOfColumns = block.son11->cols();
        returnMatrix.block(sonRowIndexStart - rowIndexStart, sonColumnIndexStart - columnIndexStart, sonNumberOfRows, sonNumberOfColumns) += son11Matrix;
    }
    if(block.son12 != nullptr)
    {
        long sonRowIndexStart = block.son12->rowStartIndex();
        long sonNumberOfRows = block.son12->rows();
        long sonColumnIndexStart = block.son12->colStartIndex();
        long sonNumberOfColumns = block.son12->cols();
        returnMatrix.block(sonRowIndexStart - rowIndexStart, sonColumnIndexStart - columnIndexStart, sonNumberOfRows, sonNumberOfColumns) += son12Matrix;
    }
    if(block.son21 != nullptr)
    {
        long sonRowIndexStart = block.son21->rowStartIndex();
        long sonNumberOfRows = block.son21->rows();
        long sonColumnIndexStart = block.son21->colStartIndex();
        long sonNumberOfColumns = block.son21->cols();
        returnMatrix.block(sonRowIndexStart - rowIndexStart, sonColumnIndexStart - columnIndexStart, sonNumberOfRows, sonNumberOfColumns) += son21Matrix;
    }
    if(block.son22 != nullptr)
    {
        long sonRowIndexStart = block.son22->rowStartIndex();
        long sonNumberOfRows = block.son22->rows();
        long sonColumnIndexStart = block.son22->colStartIndex();
        long sonNumberOfColumns = block.son22->cols();
        returnMatrix.block(sonRowIndexStart - rowIndexStart, sonColumnIndexStart - columnIndexStart, sonNumberOfRows, sonNumberOfColumns) += son22Matrix;
    }

    if(block.UMat.size() != 0)
    {
        addrkMatToFull(returnMatrix, block.UMat, block.singularValues, block.VAdjMat);
        clearRkMatrix(block);
    }

    if(block.fullMat.size() != 0)
    {
        returnMatrix += block.fullMat;
        block.fullMat.resize(0,0);
    }
    block.isLeaf = false;

    return returnMatrix;
}

void HArithm::reducedRankAgglomerationWithTrimming(BlockCluster &block, const long rank, const double relError)
{
    reducedRankAgglomeration(block, rank, relError);
    block.trimBelow();
}

void HArithm::reducedRankAgglomeration(BlockCluster &block, const long rank, const double relError) //agglomerates all subblocks into the RK matrix of the argument block
{
//    long rowIndexStart = block.rowStartIndex();
    long numberOfRows = block.rows();
//    long columnIndexStart = block.colStartIndex();
    long numberOfColumns = block.cols();

    block.isLeaf = true;
    block.isAdmissible = true;

    if(block.fullMat.size() != 0)
    {
        fullMatRankReduction(block, rank, relError); // function moves full matrix information into rk matrix
    }

    std::tuple<Eigen::MatrixXcd, Eigen::VectorXcd, Eigen::MatrixXcd> tmp;
    if(block.son11 != nullptr)
    {
//        long sonRowIndexStart = block.son11->rowStartIndex();
        long sonNumberOfRows = block.son11->rows();
//        long sonColumnIndexStart = block.son11->colStartIndex();
        long sonNumberOfColumns = block.son11->cols();
        tmp = reducedRankAgglomerationRecursion(* block.son11, rank, relError);
        if( std::abs(std::get<1>(tmp)(0)) != 0) // if first singular value is zero -> whole matrix is zero
        {
            long sonRank = std::get<0>(tmp).cols();

            Eigen::MatrixXcd  tmpUMat = Eigen::MatrixXcd(numberOfRows, sonRank);
            tmpUMat << std::get<0>(tmp), Eigen::MatrixXcd::Zero(numberOfRows - sonNumberOfRows, sonRank);
            horizontalJoinMatricesInToFirstOne(block.UMat, tmpUMat);

            joinVectorsInToFirstOne(block.singularValues, std::get<1>(tmp));

            Eigen::MatrixXcd  tmpVAdjMat = Eigen::MatrixXcd(sonRank, numberOfColumns);
            tmpVAdjMat << std::get<2>(tmp), Eigen::MatrixXcd::Zero(sonRank, numberOfColumns - sonNumberOfColumns);
            verticalJoinMatricesInToFirstOne(block.VAdjMat, tmpVAdjMat);
            RkMatRankReduction(block, rank, relError);
        }
    }
    if(block.son12 != nullptr)
    {
//        long sonRowIndexStart = block.son12->rowStartIndex();
        long sonNumberOfRows = block.son12->rows();
//        long sonColumnIndexStart = block.son12->colStartIndex();
        long sonNumberOfColumns = block.son12->cols();
        tmp = reducedRankAgglomerationRecursion(* block.son12, rank, relError);
        if( std::abs(std::get<1>(tmp)(0)) != 0)
        {
            long sonRank = std::get<0>(tmp).cols();

            Eigen::MatrixXcd  tmpUMat = Eigen::MatrixXcd(numberOfRows, sonRank);
            tmpUMat << std::get<0>(tmp), Eigen::MatrixXcd::Zero(numberOfRows - sonNumberOfRows, sonRank);
            horizontalJoinMatricesInToFirstOne(block.UMat, tmpUMat);

            joinVectorsInToFirstOne(block.singularValues, std::get<1>(tmp));

            Eigen::MatrixXcd  tmpVAdjMat = Eigen::MatrixXcd(sonRank, numberOfColumns);
            tmpVAdjMat << Eigen::MatrixXcd::Zero(sonRank, numberOfColumns - sonNumberOfColumns), std::get<2>(tmp);
            verticalJoinMatricesInToFirstOne(block.VAdjMat, tmpVAdjMat);
            RkMatRankReduction(block, rank, relError);
        }
    }
    if(block.son21 != nullptr)
    {
//        long sonRowIndexStart = block.son21->rowStartIndex();
        long sonNumberOfRows = block.son21->rows();
//        long sonColumnIndexStart = block.son21->colStartIndex();
        long sonNumberOfColumns = block.son21->cols();
        tmp = reducedRankAgglomerationRecursion(* block.son21, rank, relError);
        if( std::abs(std::get<1>(tmp)(0)) != 0)
        {
            long sonRank = std::get<0>(tmp).cols();

            Eigen::MatrixXcd  tmpUMat = Eigen::MatrixXcd(numberOfRows, sonRank);
            tmpUMat << Eigen::MatrixXcd::Zero(numberOfRows - sonNumberOfRows, sonRank), std::get<0>(tmp);
            horizontalJoinMatricesInToFirstOne(block.UMat, tmpUMat);

            joinVectorsInToFirstOne(block.singularValues, std::get<1>(tmp));

            Eigen::MatrixXcd  tmpVAdjMat = Eigen::MatrixXcd(sonRank, numberOfColumns);
            tmpVAdjMat << std::get<2>(tmp), Eigen::MatrixXcd::Zero(sonRank, numberOfColumns - sonNumberOfColumns);
            verticalJoinMatricesInToFirstOne(block.VAdjMat, tmpVAdjMat);
            RkMatRankReduction(block, rank, relError);
        }
    }
    if(block.son22 != nullptr)
    {
//        long sonRowIndexStart = block.son22->rowStartIndex();
        long sonNumberOfRows = block.son22->rows();
//        long sonColumnIndexStart = block.son22->colStartIndex();
        long sonNumberOfColumns = block.son22->cols();
        tmp = reducedRankAgglomerationRecursion(* block.son22, rank, relError);
        if( std::abs(std::get<1>(tmp)(0)) != 0)
        {
            long sonRank = std::get<0>(tmp).cols();

            Eigen::MatrixXcd  tmpUMat = Eigen::MatrixXcd(numberOfRows, sonRank);
            tmpUMat << Eigen::MatrixXcd::Zero(numberOfRows - sonNumberOfRows, sonRank), std::get<0>(tmp);
            horizontalJoinMatricesInToFirstOne(block.UMat, tmpUMat);

            joinVectorsInToFirstOne(block.singularValues, std::get<1>(tmp));

            Eigen::MatrixXcd  tmpVAdjMat = Eigen::MatrixXcd(sonRank, numberOfColumns);
            tmpVAdjMat << Eigen::MatrixXcd::Zero(sonRank, numberOfColumns - sonNumberOfColumns), std::get<2>(tmp);
            verticalJoinMatricesInToFirstOne(block.VAdjMat, tmpVAdjMat);
            RkMatRankReduction(block, rank, relError);
        }
    }
    if(block.UMat.cols() != 0 || block.fullMat.size() != 0)
    {
        block.isLeaf = true;
        if(block.fullMat.size() != 0)
        {
            block.isAdmissible = false;
        }
        else
        {
            block.isAdmissible = true;
        }
    }
    else
    {
        block.isLeaf = false;
    }
}

std::tuple<Eigen::MatrixXcd, Eigen::VectorXcd, Eigen::MatrixXcd> HArithm::reducedRankAgglomerationRecursion(BlockCluster &block, const long rank, const double relError) //agglomerates all subblocks into the full return matrix
{
    std::tuple<Eigen::MatrixXcd, Eigen::VectorXcd, Eigen::MatrixXcd> son11RKMatrix, son12RKMatrix, son21RKMatrix, son22RKMatrix; //reduced rank matrices for the sons of the block
    #pragma omp parallel
    {
        #pragma omp sections
        {
            #pragma omp section
            {
                if(block.son11 != nullptr)
                {
                    son11RKMatrix = reducedRankAgglomerationRecursion(* block.son11, rank, relError);
                }
            }
            #pragma omp section
            {
                if(block.son12 != nullptr)
                {
                    son12RKMatrix = reducedRankAgglomerationRecursion(* block.son12, rank, relError);
                }
            }
            #pragma omp section
            {
                if(block.son21 != nullptr)
                {
                    son21RKMatrix = reducedRankAgglomerationRecursion(* block.son21, rank, relError);
                }
            }
            #pragma omp section
            {
                if(block.son22 != nullptr)
                {
                   son22RKMatrix = reducedRankAgglomerationRecursion(* block.son22, rank, relError);
                }
            }
        }
    }
//    long rowIndexStart = block.rowStartIndex();
    long numberOfRows = block.rows();
//    long columnIndexStart = block.colStartIndex();
    long numberOfColumns = block.cols();

    if(block.son11 != nullptr && std::abs(std::get<1>(son11RKMatrix)(0)) != 0)
    {
//        long sonRowIndexStart = block.son11->rowCluster->triangleIndexes.first();
        long sonNumberOfRows = block.son11->rows();
//        long sonColumnIndexStart = block.son11->columnCluster->triangleIndexes.first();
        long sonNumberOfColumns = block.son11->cols();
        long sonRank = std::get<0>(son11RKMatrix).cols();
        Eigen::MatrixXcd  tmpUMat = Eigen::MatrixXcd(numberOfRows, sonRank);
        tmpUMat << std::get<0>(son11RKMatrix), Eigen::MatrixXcd::Zero(numberOfRows - sonNumberOfRows, sonRank);
        horizontalJoinMatricesInToFirstOne(block.UMat, tmpUMat);

        joinVectorsInToFirstOne(block.singularValues, std::get<1>(son11RKMatrix));

        Eigen::MatrixXcd  tmpVAdjMat = Eigen::MatrixXcd(sonRank, numberOfColumns);
        tmpVAdjMat << std::get<2>(son11RKMatrix), Eigen::MatrixXcd::Zero(sonRank, numberOfColumns - sonNumberOfColumns);
        verticalJoinMatricesInToFirstOne(block.VAdjMat, tmpVAdjMat);

        RkMatRankReduction(block, rank, relError);
    }
    if(block.son12 != nullptr && std::abs(std::get<1>(son12RKMatrix)(0)) != 0)
    {
//        long sonRowIndexStart = block.son12->rowCluster->triangleIndexes.first();
        long sonNumberOfRows = block.son12->rows();
//        long sonColumnIndexStart = block.son12->columnCluster->triangleIndexes.first();
        long sonNumberOfColumns = block.son12->cols();
        long sonRank = std::get<0>(son12RKMatrix).cols();

        Eigen::MatrixXcd  tmpUMat = Eigen::MatrixXcd(numberOfRows, sonRank);
        tmpUMat << std::get<0>(son12RKMatrix), Eigen::MatrixXcd::Zero(numberOfRows - sonNumberOfRows, sonRank);
        horizontalJoinMatricesInToFirstOne(block.UMat, tmpUMat);

        joinVectorsInToFirstOne(block.singularValues, std::get<1>(son12RKMatrix));

        Eigen::MatrixXcd  tmpVAdjMat = Eigen::MatrixXcd(sonRank, numberOfColumns);
        tmpVAdjMat << Eigen::MatrixXcd::Zero(sonRank, numberOfColumns - sonNumberOfColumns), std::get<2>(son12RKMatrix);
        verticalJoinMatricesInToFirstOne(block.VAdjMat, tmpVAdjMat);

        RkMatRankReduction(block, rank, relError);
    }
    if(block.son21 != nullptr && std::abs(std::get<1>(son21RKMatrix)(0)) != 0)
    {
//        long sonRowIndexStart = block.son21->rowCluster->triangleIndexes.first();
        long sonNumberOfRows = block.son21->rows();
//        long sonColumnIndexStart = block.son21->columnCluster->triangleIndexes.first();
        long sonNumberOfColumns = block.son21->cols();
        long sonRank = std::get<0>(son21RKMatrix).cols();

        Eigen::MatrixXcd  tmpUMat = Eigen::MatrixXcd(numberOfRows, sonRank);
        tmpUMat << Eigen::MatrixXcd::Zero(numberOfRows - sonNumberOfRows, sonRank), std::get<0>(son21RKMatrix);
        horizontalJoinMatricesInToFirstOne(block.UMat, tmpUMat);

        joinVectorsInToFirstOne(block.singularValues, std::get<1>(son21RKMatrix));

        Eigen::MatrixXcd  tmpVAdjMat = Eigen::MatrixXcd(sonRank, numberOfColumns);
        tmpVAdjMat << std::get<2>(son21RKMatrix), Eigen::MatrixXcd::Zero(sonRank, numberOfColumns - sonNumberOfColumns);
        verticalJoinMatricesInToFirstOne(block.VAdjMat, tmpVAdjMat);

        RkMatRankReduction(block, rank, relError);
    }
    if(block.son22 != nullptr && std::abs(std::get<1>(son22RKMatrix)(0)) != 0)
    {
//        long sonRowIndexStart = block.son22->rowCluster->triangleIndexes.first();
        long sonNumberOfRows = block.son22->rows();
//        long sonColumnIndexStart = block.son22->columnCluster->triangleIndexes.first();
        long sonNumberOfColumns = block.son22->cols();
        long sonRank = std::get<0>(son22RKMatrix).cols();

        Eigen::MatrixXcd  tmpUMat = Eigen::MatrixXcd(numberOfRows, sonRank);
        tmpUMat << Eigen::MatrixXcd::Zero(numberOfRows - sonNumberOfRows, sonRank), std::get<0>(son22RKMatrix);
        horizontalJoinMatricesInToFirstOne(block.UMat, tmpUMat);

        joinVectorsInToFirstOne(block.singularValues, std::get<1>(son22RKMatrix));

        Eigen::MatrixXcd  tmpVAdjMat = Eigen::MatrixXcd(sonRank, numberOfColumns);
        tmpVAdjMat << Eigen::MatrixXcd::Zero(sonRank, numberOfColumns - sonNumberOfColumns), std::get<2>(son22RKMatrix);
        verticalJoinMatricesInToFirstOne(block.VAdjMat, tmpVAdjMat);

        RkMatRankReduction(block, rank, relError);
    }

    if(block.fullMat.size() != 0)
    {
        fullMatRankReduction(block, rank, relError);
    }

    if(block.UMat.size() == 0)
    {
        setRkMatrixZero(block);
    }
    std::tuple<Eigen::MatrixXcd, Eigen::VectorXcd, Eigen::MatrixXcd> returnRKMatrix = {block.UMat, block.singularValues, block.VAdjMat};
    clearRkMatrix(block);

    block.isLeaf = false;
    return returnRKMatrix;
}

void HArithm::blockSplitting(BlockCluster &block, const long rank, const double error)
{
    if(block.fullMat.size() == 0 && block.UMat.size() == 0) // no work to do on block
    {
        return;
    }
    else if(block.son11 == nullptr) // block is leaf node
    {
         if(block.fullMat.size() != 0 && block.UMat.size()!= 0)
         {
            block.fullMat.noalias() += block.UMat * block.singularValues.asDiagonal() * block.VAdjMat;
            clearRkMatrix(block);
            block.isAdmissible = false;
         }
         else if( block.UMat.size()!= 0) //block is rk-matrix
         {
            block.isAdmissible = true;
         }
         else //block.fullMat.size()!= 0
         {
            block.isAdmissible = false;
         }
         block.isLeaf = true;
         return;
    }
    else // block isn't tree leaf
    {
        block.isLeaf = false;
        if(block.fullMat.size() != 0)
        {
//            #pragma omp parallel
//            {
//                #pragma omp sections
//                {
                    #pragma omp task shared(block)
                    {
//                        long sonRowIndexStart = block.son11->rowStartIndex();
                        long sonNumberOfRows = block.son11->rows();
//                        long sonColumnIndexStart = block.son11->colStartIndex();
                        long sonNumberOfColumns = block.son11->cols();
                        addFullMatrixInToFirstOne(block.son11->fullMat, block.fullMat.topLeftCorner(sonNumberOfRows, sonNumberOfColumns));
                    }
                    #pragma omp task shared(block)
                    {
//                        long sonRowIndexStart = block.son12->rowStartIndex();
                        long sonNumberOfRows = block.son12->rows();
//                        long sonColumnIndexStart = block.son12->colStartIndex();
                        long sonNumberOfColumns = block.son12->cols();
                        addFullMatrixInToFirstOne(block.son12->fullMat, block.fullMat.topRightCorner(sonNumberOfRows, sonNumberOfColumns));
                    }
                    #pragma omp task shared(block)
                    {
//                        long sonRowIndexStart = block.son21->rowStartIndex();
                        long sonNumberOfRows = block.son21->rows();
//                        long sonColumnIndexStart = block.son21->colStartIndex();
                        long sonNumberOfColumns = block.son21->cols();
                        addFullMatrixInToFirstOne(block.son21->fullMat, block.fullMat.bottomLeftCorner(sonNumberOfRows, sonNumberOfColumns));
                     }
//                    #pragma omp task
                    {
//                        long sonRowIndexStart = block.son22->rowStartIndex();
                        long sonNumberOfRows = block.son22->rows();
//                        long sonColumnIndexStart = block.son22->colStartIndex();
                        long sonNumberOfColumns = block.son22->cols();
                        addFullMatrixInToFirstOne(block.son22->fullMat, block.fullMat.bottomRightCorner(sonNumberOfRows, sonNumberOfColumns));
                    }
//                }
//            }
            #pragma omp taskwait
            block.fullMat.resize(0,0);
        }
        if(block.UMat.size() != 0)
        {
            if(block.fullMat.size() > 0)
            {
                std::cerr << "fullMat.size() > 0 " << std::endl;
            }
//            #pragma omp parallel
//            {
//                #pragma omp single
//                {
                    #pragma omp task shared(block)
                    {
//                        long sonRowIndexStart = block.son11->rowStartIndex();
                        long sonNumberOfRows = block.son11->rows();
//                        long sonColumnIndexStart = block.son11->colStartIndex();
                        long sonNumberOfColumns = block.son11->cols();
                        roundedAddRMatToBlock(*block.son11, block.UMat.topRows(sonNumberOfRows), block.singularValues, block.VAdjMat.leftCols(sonNumberOfColumns), rank, error);
                    }
                    #pragma omp task shared(block)
                    {
//                         long sonRowIndexStart = block.son12->rowStartIndex();
                         long sonNumberOfRows = block.son12->rows();
//                         long sonColumnIndexStart = block.son12->colStartIndex();
                         long sonNumberOfColumns = block.son12->cols();
                         roundedAddRMatToBlock(*block.son12, block.UMat.topRows(sonNumberOfRows), block.singularValues, block.VAdjMat.rightCols(sonNumberOfColumns), rank, error);
                    }
                    #pragma omp task shared(block)
                    {
//                        long sonRowIndexStart = block.son21->rowStartIndex();
                        long sonNumberOfRows = block.son21->rows();
//                        long sonColumnIndexStart = block.son21->colStartIndex();
                        long sonNumberOfColumns = block.son21->cols();
                        roundedAddRMatToBlock(*block.son21, block.UMat.bottomRows(sonNumberOfRows), block.singularValues, block.VAdjMat.leftCols(sonNumberOfColumns), rank, error);
                    }
//                    #pragma omp task
                    {
//                        long sonRowIndexStart = block.son22->rowStartIndex();
                        long sonNumberOfRows = block.son22->rows();
//                        long sonColumnIndexStart = block.son22->colStartIndex();
                        long sonNumberOfColumns = block.son22->cols();
                        roundedAddRMatToBlock(*block.son22, block.UMat.bottomRows(sonNumberOfRows), block.singularValues, block.VAdjMat.rightCols(sonNumberOfColumns), rank, error);
                    }
//                }
//            }
            #pragma omp taskwait
            clearRkMatrix(block);
            block.isAdmissible = false;
        }
//        #pragma omp parallel
//        {
//            #pragma omp single nowait
//            {
                #pragma omp task shared(block)
                {
                    blockSplitting(*block.son11, rank, error);
                }
                #pragma omp task shared(block)
                {
                    blockSplitting(*block.son12, rank, error);
                }
                #pragma omp task shared(block)
                {
                    blockSplitting(*block.son21, rank, error);
                }
//                #pragma omp task
                {
                    blockSplitting(*block.son22, rank, error);
                }
                #pragma omp taskwait

//            }
//        }
    }
}

Eigen::VectorXcd HArithm::LUSolve(HMatrix &matrix, Eigen::VectorXcd &rightHandSide, const long rank, const double relError)
{
    const long vectorLength = matrix.cols();
    const long startIndex = matrix.colStartIndex();

    Eigen::VectorXcd solutionVector = Eigen::VectorXcd::Zero(vectorLength);
    Eigen::VectorXcd tmpVector = Eigen::VectorXcd::Zero(vectorLength);
    if(matrix.rows() != matrix.cols() || matrix.rowStartIndex() != matrix.colStartIndex())
    {
        std::cerr << "Matrix in LUSolve isn't square." << std::endl;
        return solutionVector;
    }
    if(vectorLength != rightHandSide.size())
    {
        std::cerr << "Matrix in LUSolve and rightHandSide vector dimensions aren't compatible." << std::endl;
        return solutionVector;
    }

    std::pair<HMatrix,HMatrix> LUPair = LUDecomposition(matrix, rank, relError, true);
    forwardSubstitution(*LUPair.first.getRootBlock(), tmpVector, rightHandSide, startIndex);
    backwardSubstitution(*LUPair.second.getRootBlock(), solutionVector, tmpVector, startIndex);

    LUPair.first.clear();
    LUPair.second.clear(true);

    return solutionVector;
}

Eigen::VectorXcd HArithm::LUSubstitutionSolve(HMatrix &L, HMatrix &U, const Eigen::VectorXcd &rightHandSide)
{
    long vectorLength = L.rows();

    Eigen::VectorXcd solutionVector = Eigen::VectorXcd::Zero(vectorLength);
    Eigen::VectorXcd tmpVector = Eigen::VectorXcd::Zero(vectorLength);
    Eigen::VectorXcd rightHandSideCopy = rightHandSide;
    if(L.rows() != L.cols() || U.rows() != U.cols() || L.rows() != U.rows())
    {
        std::cerr << "Wrong matrix dimensions in LUSubstitutionSolve() call." << std::endl;
        return solutionVector;
    }
    if(vectorLength != rightHandSideCopy.size())
    {
        std::cerr << "matrix in LUSolve and rightHandSide vector dimensions aren't compatible." << std::endl;
        return solutionVector;
    }

    forwardSubstitution(*L.getRootBlock(), tmpVector, rightHandSideCopy);
    backwardSubstitution(*U.getRootBlock(), solutionVector, tmpVector);

    return solutionVector;
}

std::pair<HMatrix, HMatrix> HArithm::LUDecomposition(HMatrix &matrix, const long rank, const double relError, const bool zeroBlocksUninitialized)
{
    HMatrix L(matrix.getRowClustertree(), matrix.getColumnClustertree());
    L.createRootNode();
    HMatrix U(matrix.getRowClustertree(), matrix.getColumnClustertree());
    U.createRootNode();

    recursiveLUDecomposition(*L.getRootBlock(), *U.getRootBlock(), *matrix.getRootBlock(), rank, relError, zeroBlocksUninitialized);

    return std::pair<HMatrix,HMatrix> (L,U);
}

void HArithm::recursiveLUDecomposition(BlockCluster &LBlock, BlockCluster &UBlock, BlockCluster &ABlock, const long rank, const double relError, const bool zeroBlocksUninitialized) // L*U=A
{
    if(ABlock.isLeaf)
    {
        if(ABlock.isAdmissible)
        {
            fullRankAgglomerationWithTrimming(ABlock);
            std::cerr<<"ABlock->isAdmissible in recursiveLUDecomposition"<<std::endl;
        }
        long rows =  ABlock.fullMat.rows();
        long cols =  ABlock.fullMat.cols();
        if constexpr(inversionInLUDecomp)
        {
            LBlock.fullMat = ABlock.fullMat.inverse();
            UBlock.fullMat = Eigen::MatrixXcd::Identity(rows, cols);

//            if(ABlock->fullMat.fullPivLu().isInvertible())
//            {
//                std::cerr << "Noninvertable block." << std::endl;
//            }
//            if(LBlock->fullMat.hasNaN())
//            {
//                std::cerr << "NaNs encoundered during the inversion step in the recursiveLUDecomposition() routine!" << std::endl;
//            }
//            if(!LBlock->fullMat.allFinite())
//            {
//                std::cerr << "Infs encoundered during the inversion step in the recursiveLUDecomposition() routine!" << std::endl;
//            }
//            if((LBlock->fullMat * ABlock->fullMat - Eigen::MatrixXcd::Identity(rows, cols)).norm() > global::Tiny)
//            {
//                std::cerr << "Large error in recursiveLUDecomposition." << std::endl;
//            }
        }
        else if constexpr(qrInLUDecomp)
        {
            Eigen::HouseholderQR<Eigen::MatrixXcd> Qr(ABlock.fullMat);
            LBlock.fullMat = Qr.householderQ();
            UBlock.fullMat = Qr.matrixQR().triangularView<Eigen::Upper>().toDenseMatrix();
        }
        else
        {
            Eigen::MatrixXcd LU = global::LUDecompNoPivoting(ABlock.fullMat);
            LBlock.fullMat = LU.leftCols(std::min(rows,cols)).triangularView<Eigen::UnitLower>();
            UBlock.fullMat = LU.topRows(std::min(rows,cols)).triangularView<Eigen::Upper>();
        }
        LBlock.isLeaf = true;
        LBlock.isAdmissible = false;
        UBlock.isLeaf = true;
        UBlock.isAdmissible = false;
    }
    else
    {
        LBlock.son11 = std::make_unique<BlockCluster>(ABlock.rowCluster->son1.get(), ABlock.columnCluster->son1.get(), &LBlock);
        LBlock.son12 = std::make_unique<BlockCluster>(ABlock.rowCluster->son1.get(), ABlock.columnCluster->son2.get(), &LBlock);

        if(!zeroBlocksUninitialized)
        {
            convertToAdmissibleZeroBlock(*LBlock.son12);
        }
        else
        {
            LBlock.son12->isLeaf = true;
        }
        LBlock.son21 = std::make_unique<BlockCluster>(ABlock.rowCluster->son2.get(), ABlock.columnCluster->son1.get(), &LBlock);
        LBlock.son22 = std::make_unique<BlockCluster>(ABlock.rowCluster->son2.get(), ABlock.columnCluster->son2.get(), &LBlock);

        UBlock.son11 = std::make_unique<BlockCluster>(ABlock.rowCluster->son1.get(), ABlock.columnCluster->son1.get(), &UBlock);
        UBlock.son12 = std::make_unique<BlockCluster>(ABlock.rowCluster->son1.get(), ABlock.columnCluster->son2.get(), &UBlock);
        UBlock.son21 = std::make_unique<BlockCluster>(ABlock.rowCluster->son2.get(), ABlock.columnCluster->son1.get(), &UBlock);

        if(!zeroBlocksUninitialized)
        {
            convertToAdmissibleZeroBlock(*UBlock.son21);
        }
        else
        {
            UBlock.son21->isLeaf = true;
        }
        UBlock.son22 = std::make_unique<BlockCluster>(ABlock.rowCluster->son2.get(), ABlock.columnCluster->son2.get(), &UBlock);

        recursiveLUDecomposition(*LBlock.son11, *UBlock.son11, *ABlock.son11, rank, relError, zeroBlocksUninitialized);

        forwSubsMatVal(*LBlock.son11, *UBlock.son12, *ABlock.son12, rank, relError);
        forwSubsMatValTransposed(*LBlock.son21, *UBlock.son11, *ABlock.son21, rank, relError);

        HMatrix product22 = HMultiply::multiplyHMat(*LBlock.son21, *UBlock.son12, rank, relError);
        recursiveHMatSubstraction(*ABlock.son22, *product22.getRootBlock(), rank, relError);
        product22.clear();

        recursiveLUDecomposition(*LBlock.son22, *UBlock.son22, *ABlock.son22, rank, relError, zeroBlocksUninitialized);
    }
    ABlock.clear();
}

Eigen::VectorXcd HArithm::forwardSubstitution(const HMatrix &L, Eigen::VectorXcd b) // solve L * y = b for the unknown y, L is lower triangular H-matrix
{
    const long rowStartIndex = L.rowStartIndex();
    const long colStartIndex = L.colStartIndex();
    const long rows = L.rows();
    const long cols = L.cols();
    if(rows != cols || rowStartIndex != colStartIndex || rows != b.size())
    {
        std::cerr << "Incompatible matrix or right-hand side in forwardSubstitution() call." << std::endl;
        return Eigen::VectorXcd::Zero(cols);
    }
    Eigen::VectorXcd y(rows);
    HArithm::forwardSubstitution(*L.getRootBlock(), y, b, rowStartIndex);
    return y;
}

void HArithm::forwardSubstitution(const BlockCluster &LBlock, Eigen::VectorXcd &solution, Eigen::VectorXcd &b, const long vectorStartIndex) // solve L * y = b for the unknown y, L is lower triangular H-matrix
{
    if(LBlock.isLeaf)
    {
        const long startIndex = LBlock.colStartIndex();
        const long vectorLength = LBlock.rows();
        if constexpr(inversionInLUDecomp)
        {
            solution.segment(startIndex - vectorStartIndex, vectorLength) = LBlock.fullMat * b.segment(startIndex - vectorStartIndex, vectorLength);
        }
        else if constexpr(qrInLUDecomp)
        {
            solution.segment(startIndex - vectorStartIndex, vectorLength) = LBlock.fullMat.adjoint() * b.segment(startIndex - vectorStartIndex, vectorLength);
        }
        else
        {
            solution.segment(startIndex - vectorStartIndex, vectorLength) = LBlock.fullMat.triangularView<Eigen::Lower>().solve(b.segment(startIndex - vectorStartIndex, vectorLength));
//            solution.segment(startIndex - vectorStartIndex, vectorLength) = global::forwardSubstitution(LBlock.fullMat, b.segment(startIndex - vectorStartIndex, vectorLength));
        }
    }
    else
    {
        forwardSubstitution(*LBlock.son11, solution, b, vectorStartIndex);
        subtractiveParallelMatrixVectorPoduct(b, *LBlock.son21, solution, vectorStartIndex, vectorStartIndex);
        forwardSubstitution(*LBlock.son22, solution, b, vectorStartIndex);
    }
}

void HArithm::forwardSubstitutionTransposed(Eigen::RowVectorXcd &solution, const BlockCluster &UBlock, Eigen::RowVectorXcd &b, const long vectorStartIndex) // solve y^T * U = b^T for y, U is upper triangular H-matrix
{
    if(UBlock.isLeaf)
    {
        const long startIndex = UBlock.rowStartIndex();
        const long vectorLength = UBlock.cols();
        if constexpr(inversionInLUDecomp)
        {
            solution.segment(startIndex - vectorStartIndex, vectorLength) = b.segment(startIndex - vectorStartIndex, vectorLength);
        }
        else if constexpr(qrInLUDecomp)
        {
            solution.segment(startIndex - vectorStartIndex, vectorLength) = UBlock.fullMat.transpose().triangularView<Eigen::Lower>().solve(b.segment(startIndex - vectorStartIndex, vectorLength).transpose());
        }
        else
        {
            solution.segment(startIndex - vectorStartIndex, vectorLength) = UBlock.fullMat.transpose().triangularView<Eigen::Lower>().solve(b.segment(startIndex - vectorStartIndex, vectorLength).transpose());
        }
    }
    else
    {
        forwardSubstitutionTransposed(solution, *UBlock.son11, b, vectorStartIndex);
//        subtractiveRecursiveVectorMatrixPoduct(b, solution, UBlock.son12, vectorStartIndex, vectorStartIndex); // faster than recursiveVectorMatrixPoduct(b, -solution, UBlock.son12, vectorStartIndex, vectorStartIndex);
        subtractiveParallelVectorMatrixPoduct(b, solution, *UBlock.son12, vectorStartIndex, vectorStartIndex);
        forwardSubstitutionTransposed(solution, *UBlock.son22, b, vectorStartIndex);
    }
}

Eigen::VectorXcd HArithm::backwardSubstitution(const HMatrix &U, Eigen::VectorXcd y) // \brief Solve U * x = y for x; U is upper triangular H-matrix.
{
    const long rowStartIndex = U.rowStartIndex();
    const long colStartIndex = U.colStartIndex();
    const long rows = U.rows();
    const long cols = U.cols();
    if(rows != cols || rowStartIndex != colStartIndex || rows != y.size())
    {
        std::cerr << "Incompatible matrix or right-hand side in backwardSubstitution() call." << std::endl;
        return Eigen::VectorXcd::Zero(cols);
    }
    Eigen::VectorXcd x(rows);
    HArithm::backwardSubstitution(*U.getRootBlock(), x, y, rowStartIndex);
    return x;
}

void HArithm::backwardSubstitution(const BlockCluster &UBlock, Eigen::VectorXcd &solution, Eigen::VectorXcd &y, const long vectorStartIndex) // solve U * x = y, U is upper triangular H-matrix
{
    if(UBlock.isLeaf)
    {
        const long startIndex = UBlock.colStartIndex();
        const long vectorLength = UBlock.rows();
        if constexpr(inversionInLUDecomp)
        {
            solution.segment(startIndex - vectorStartIndex, vectorLength) = y.segment(startIndex - vectorStartIndex, vectorLength);
        }
        else if constexpr(qrInLUDecomp)
        {
            solution.segment(startIndex - vectorStartIndex, vectorLength) = UBlock.fullMat.triangularView<Eigen::Upper>().solve(y.segment(startIndex - vectorStartIndex, vectorLength));
        }
        else
        {
            solution.segment(startIndex - vectorStartIndex, vectorLength) = UBlock.fullMat.triangularView<Eigen::Upper>().solve(y.segment(startIndex - vectorStartIndex, vectorLength));
//            solution.segment(startIndex - vectorStartIndex, vectorLength) = global::backwardSubstitution(UBlock.fullMat, y.segment(startIndex - vectorStartIndex, vectorLength));
        }
    }
    else
    {
        backwardSubstitution(*UBlock.son22, solution, y, vectorStartIndex);
        subtractiveParallelMatrixVectorPoduct(y, *UBlock.son12, solution, vectorStartIndex, vectorStartIndex);
        backwardSubstitution(*UBlock.son11, solution, y, vectorStartIndex);
    }
}

void HArithm::backwardSubstitutionTransposed(Eigen::RowVectorXcd &solution, const BlockCluster &LBlock,  Eigen::RowVectorXcd &y, const long vectorStartIndex) // solve x^T * L = y^T for y, L is lower triangular H-matrix
{
    if(LBlock.isLeaf)
    {
        const long startIndex = LBlock.rowStartIndex();
        const long vectorLength = LBlock.rows();
        if constexpr(inversionInLUDecomp)
        {
            solution.segment(startIndex - vectorStartIndex, vectorLength) = y.segment(startIndex - vectorStartIndex, vectorLength) * LBlock.fullMat;
        }
        else if constexpr(qrInLUDecomp)
        {
            solution.segment(startIndex - vectorStartIndex, vectorLength) = y.segment(startIndex - vectorStartIndex, vectorLength) * LBlock.fullMat.adjoint();
        }
        else
        {
            solution.segment(startIndex - vectorStartIndex, vectorLength) = LBlock.fullMat.transpose().triangularView<Eigen::Upper>().solve(y.segment(startIndex - vectorStartIndex, vectorLength).transpose());
        }
    }
    else
    {
        backwardSubstitutionTransposed(solution, *LBlock.son22, y, vectorStartIndex);
        subtractiveParallelVectorMatrixPoduct(y, solution, *LBlock.son21, vectorStartIndex, vectorStartIndex);
        backwardSubstitutionTransposed(solution, *LBlock.son11, y, vectorStartIndex);
    }
}

HMatrix HArithm::forwSubsMatVal(const HMatrix &L, HMatrix &Z, const long rank, const double relError)
{
    HMatrix X(Z.getRowClustertree(), Z.getColumnClustertree());
    X.createRootNode();
    forwSubsMatVal(*L.getRootBlock(), *X.getRootBlock(), *Z.getRootBlock(), rank, relError);
    Z.clear();
    return X;
}

void HArithm::forwSubsMatVal(const BlockCluster &LBlock, BlockCluster &XBlock, BlockCluster &ZBlock, const long rank, const double relError) // solve L * X = Z for X, L is lower triangular h-matrix
{
    if(!ZBlock.isLeaf  && LBlock.isLeaf)
    {
        reducedRankAgglomerationWithTrimming(ZBlock, rank, relError);
    }
    if(ZBlock.isLeaf )
    {       
        if(ZBlock.isAdmissible)
        {
            XBlock.UMat = Eigen::MatrixXcd(ZBlock.UMat.rows(), ZBlock.UMat.cols());
            #pragma omp parallel for
            for(long j = 0; j < ZBlock.UMat.cols(); j++)
            {
                Eigen::VectorXcd xVector = Eigen::VectorXcd::Zero(ZBlock.UMat.rows());
                Eigen::VectorXcd zVector = ZBlock.UMat.col(j);
                long vectorStartIndex = ZBlock.rowStartIndex();

                forwardSubstitution(LBlock, xVector, zVector, vectorStartIndex);
                XBlock.UMat.col(j) = xVector;
                ZBlock.UMat.col(j) = zVector;
            }
            XBlock.singularValues = ZBlock.singularValues;
            XBlock.VAdjMat = ZBlock.VAdjMat;
            XBlock.isLeaf = true;
            XBlock.isAdmissible = true;
        }
        else
        {
            XBlock.fullMat = Eigen::MatrixXcd(ZBlock.fullMat.rows(), ZBlock.fullMat.cols());
            #pragma omp parallel for
            for(long j = 0; j < ZBlock.fullMat.cols(); j++)
            {
                Eigen::VectorXcd xVector = Eigen::VectorXcd::Zero(ZBlock.fullMat.rows());
                Eigen::VectorXcd zVector = ZBlock.fullMat.col(j);
                long vectorStartIndex = ZBlock.rowStartIndex();

                forwardSubstitution(LBlock, xVector, zVector, vectorStartIndex);
                XBlock.fullMat.col(j) = xVector;
                ZBlock.fullMat.col(j) = zVector;
            }
            XBlock.isLeaf = true;
            XBlock.isAdmissible = false;
        }
    }
    else
    {
        XBlock.son11 = std::make_unique<BlockCluster>(ZBlock.rowCluster->son1.get(), ZBlock.columnCluster->son1.get(), &XBlock);
        XBlock.son12 = std::make_unique<BlockCluster>(ZBlock.rowCluster->son1.get(), ZBlock.columnCluster->son2.get(), &XBlock);
        XBlock.son21 = std::make_unique<BlockCluster>(ZBlock.rowCluster->son2.get(), ZBlock.columnCluster->son1.get(), &XBlock);
        XBlock.son22 = std::make_unique<BlockCluster>(ZBlock.rowCluster->son2.get(), ZBlock.columnCluster->son2.get(), &XBlock);

        forwSubsMatVal(*LBlock.son11, *XBlock.son11, *ZBlock.son11, rank, relError);

        HMatrix product21 = HMultiply::multiplyHMat(*LBlock.son21, *XBlock.son11, rank, relError);
        recursiveHMatSubstraction(*ZBlock.son21, *product21.getRootBlock(), rank, relError);
        product21.clear();

        forwSubsMatVal(*LBlock.son11, *XBlock.son12, *ZBlock.son12, rank, relError);

        HMatrix product22 = HMultiply::multiplyHMat(*LBlock.son21, *XBlock.son12, rank, relError);
        recursiveHMatSubstraction(*ZBlock.son22, *product22.getRootBlock(), rank, relError);
        product22.clear();

        forwSubsMatVal(*LBlock.son22, *XBlock.son21, *ZBlock.son21, rank, relError);
        forwSubsMatVal(*LBlock.son22, *XBlock.son22, *ZBlock.son22, rank, relError);
    }
}

HMatrix HArithm::forwSubsMatValTransposed(const HMatrix &U, HMatrix &Z, const long rank, const double relError)
{
    HMatrix X(Z.getRowClustertree(), Z.getColumnClustertree());
    X.createRootNode();
    forwSubsMatValTransposed(*X.getRootBlock(), *U.getRootBlock(), *Z.getRootBlock(), rank, relError);
    Z.clear();
    return X;
}

void HArithm::forwSubsMatValTransposed(BlockCluster &XBlock, const BlockCluster &UBlock, BlockCluster &ZBlock, const long rank, const double relError)
{
    if(!ZBlock.isLeaf  && UBlock.isLeaf)    /// if the h-matrix multiplication uses agglomeration a lot, the Z mmatrix can be less deep than the U Matrix
    {
        reducedRankAgglomerationWithTrimming(ZBlock, rank, relError);
    }
    if(ZBlock.isLeaf)
    {
        if(ZBlock.isAdmissible)
        {
            XBlock.VAdjMat = Eigen::MatrixXcd(ZBlock.VAdjMat.rows(), ZBlock.VAdjMat.cols());
            #pragma omp parallel for
            for(long i = 0; i < ZBlock.VAdjMat.rows(); i++)
            {
                Eigen::RowVectorXcd xVector = Eigen::VectorXcd::Zero(ZBlock.VAdjMat.cols());
                Eigen::RowVectorXcd zVector = ZBlock.VAdjMat.row(i);
                long vectorStartIndex = ZBlock.colStartIndex();

                forwardSubstitutionTransposed(xVector, UBlock, zVector, vectorStartIndex);
                XBlock.VAdjMat.row(i) = xVector;
                ZBlock.VAdjMat.row(i) = zVector;
            }
            XBlock.UMat = ZBlock.UMat;
            XBlock.singularValues = ZBlock.singularValues;
            XBlock.isLeaf = true;
            XBlock.isAdmissible = true;
        }
        else
        {
            XBlock.fullMat = Eigen::MatrixXcd(ZBlock.fullMat.rows(), ZBlock.fullMat.cols());
            #pragma omp parallel for
            for(long i = 0; i < ZBlock.fullMat.rows(); i++)
            {
                Eigen::RowVectorXcd xVector = Eigen::VectorXcd::Zero(ZBlock.fullMat.cols());
                Eigen::RowVectorXcd zVector = ZBlock.fullMat.row(i);
                long vectorStartIndex = ZBlock.colStartIndex();

                forwardSubstitutionTransposed(xVector, UBlock, zVector, vectorStartIndex);
                XBlock.fullMat.row(i) = xVector;
                ZBlock.fullMat.row(i) = zVector;
            }
            XBlock.isLeaf = true;
            XBlock.isAdmissible = false;
        }
    }
    else
    {
        XBlock.son11 = std::make_unique<BlockCluster>(ZBlock.rowCluster->son1.get(), ZBlock.columnCluster->son1.get(), &XBlock);
        XBlock.son12 = std::make_unique<BlockCluster>(ZBlock.rowCluster->son1.get(), ZBlock.columnCluster->son2.get(), &XBlock);
        XBlock.son21 = std::make_unique<BlockCluster>(ZBlock.rowCluster->son2.get(), ZBlock.columnCluster->son1.get(), &XBlock);
        XBlock.son22 = std::make_unique<BlockCluster>(ZBlock.rowCluster->son2.get(), ZBlock.columnCluster->son2.get(), &XBlock);

        forwSubsMatValTransposed(*XBlock.son11, *UBlock.son11, *ZBlock.son11, rank, relError);

        HMatrix product12 = HMultiply::multiplyHMat(*XBlock.son11, *UBlock.son12, rank, relError);
        recursiveHMatSubstraction(*ZBlock.son12, *product12.getRootBlock(), rank, relError);
        product12.clear();

        forwSubsMatValTransposed(*XBlock.son21, *UBlock.son11, *ZBlock.son21, rank, relError);

        HMatrix product22 = HMultiply::multiplyHMat(*XBlock.son21, *UBlock.son12, rank, relError);
        recursiveHMatSubstraction(*ZBlock.son22, *product22.getRootBlock(), rank, relError);
        product22.clear();

        forwSubsMatValTransposed(*XBlock.son12, *UBlock.son22, *ZBlock.son12, rank, relError);
        forwSubsMatValTransposed(*XBlock.son22, *UBlock.son22, *ZBlock.son22, rank, relError);
    }
}

void HArithm::recursiveMatrixVectorPoduct(Eigen::VectorXcd &product, const BlockCluster &block, const Eigen::VectorXcd &x, const long prodStartIndex, const long xStartIndex) // y += hmatrix * x
{
    if(block.isLeaf)
    {
        long rowStartIndex = block.rowStartIndex();
        long numberOfRows = block.rows();
        long columnStartIndex = block.colStartIndex();
        long numberOfColumns = block.cols();

        if(block.isAdmissible) // far field block
        {
            product.segment(rowStartIndex - prodStartIndex, numberOfRows).noalias() += block.UMat * (block.singularValues.asDiagonal() * (block.VAdjMat * x.segment(columnStartIndex - xStartIndex, numberOfColumns)));
        }
        else // near field block
        {
            product.segment(rowStartIndex - prodStartIndex, numberOfRows).noalias() += block.fullMat * x.segment(columnStartIndex - xStartIndex, numberOfColumns);
        }
    }
    else
    {
        #pragma omp parallel sections
        {
            #pragma omp section //prevent aliasing by putting the same y rows in one section
            {                   //otherwise use omp critical above
                recursiveMatrixVectorPoduct(product, *block.son11, x, prodStartIndex, xStartIndex);
                recursiveMatrixVectorPoduct(product, *block.son12, x, prodStartIndex, xStartIndex);
            }
            #pragma omp section
            {
                recursiveMatrixVectorPoduct(product, *block.son21, x, prodStartIndex, xStartIndex);
                recursiveMatrixVectorPoduct(product, *block.son22, x, prodStartIndex, xStartIndex);
            }
        }
    }
}

void HArithm::parallelMatrixVectorPoduct(Eigen::VectorXcd &product, const BlockCluster &block, const Eigen::VectorXcd &x, const long prodStartIndex, const long xStartIndex) // y += hmatrix * x
{
    QVector<const BlockCluster*> minPartition;
    BlockCluster::getPartition(block, minPartition);

    #pragma omp parallel for
    for(long i = 0; i < minPartition.length(); i++)
    {
        const BlockCluster* tmpBlock = minPartition.at(i);
        long rowStartIndex = tmpBlock->rowStartIndex();
        long numberOfRows = tmpBlock->rows();
        long columnStartIndex = tmpBlock->colStartIndex();
        long numberOfColumns = tmpBlock->cols();
    //        std::cout << "rowStartIndex: " << rowStartIndex << " columnStartIndex: " << columnStartIndex << std::endl;
    //        std::cout << "numberOfRows: " << numberOfRows << " numberOfColumns: " << numberOfColumns << std::endl;
        Eigen::VectorXcd tmp;
        if(tmpBlock->isAdmissible) // far field block
        {
            tmp.noalias() = tmpBlock->UMat * (tmpBlock->singularValues.asDiagonal() * (tmpBlock->VAdjMat * x.segment(columnStartIndex - xStartIndex, numberOfColumns)));
            #pragma omp critical
            product.segment(rowStartIndex - prodStartIndex, numberOfRows) += tmp;
        }
        else // near field block
        {
            tmp.noalias() = tmpBlock->fullMat * x.segment(columnStartIndex - xStartIndex, numberOfColumns);
            #pragma omp critical
            product.segment(rowStartIndex - prodStartIndex, numberOfRows) += tmp;
        }
    }
}

void HArithm::subtractiveRecursiveMatrixVectorPoduct(Eigen::VectorXcd &product, const BlockCluster &block, const Eigen::VectorXcd &x, const long prodStartIndex, const long xStartIndex) // y -= hmatrix * x
{
    if(block.isLeaf)
    {
        long rowStartIndex = block.rowStartIndex();
        long numberOfRows = block.rows();
        long columnStartIndex = block.colStartIndex();
        long numberOfColumns = block.cols();

        if(block.isAdmissible) // far field block
        {
            product.segment(rowStartIndex - prodStartIndex, numberOfRows).noalias() -= block.UMat * (block.singularValues.asDiagonal() * (block.VAdjMat * x.segment(columnStartIndex - xStartIndex, numberOfColumns)));
        }
        else // near field block
        {
            product.segment(rowStartIndex - prodStartIndex, numberOfRows).noalias() -= block.fullMat * x.segment(columnStartIndex - xStartIndex, numberOfColumns);
        }
    }
    else
    {
        #pragma omp parallel sections
        {
            #pragma omp section //prevent aliasing by putting the same y rows in one section
            {                   //otherwise use omp critical above
                subtractiveRecursiveMatrixVectorPoduct(product, *block.son11, x, prodStartIndex, xStartIndex);
                subtractiveRecursiveMatrixVectorPoduct(product, *block.son12, x, prodStartIndex, xStartIndex);
            }
            #pragma omp section
            {
                subtractiveRecursiveMatrixVectorPoduct(product, *block.son21, x, prodStartIndex, xStartIndex);
                subtractiveRecursiveMatrixVectorPoduct(product, *block.son22, x, prodStartIndex, xStartIndex);
            }
        }
    }
}

void HArithm::subtractiveParallelMatrixVectorPoduct(Eigen::VectorXcd &product, const BlockCluster &block, const Eigen::VectorXcd &x, const long prodStartIndex, const long xStartIndex) // y -= hmatrix * x
{
    QVector<const BlockCluster*> minPartition;
    BlockCluster::getPartition(block, minPartition);

    #pragma omp parallel for
    for(long i = 0; i < minPartition.length(); i++)
    {
        const BlockCluster* tmpBlock = minPartition.at(i);
        long rowStartIndex = tmpBlock->rowStartIndex();
        long numberOfRows = tmpBlock->rows();
        long columnStartIndex = tmpBlock->colStartIndex();
        long numberOfColumns = tmpBlock->cols();

        Eigen::VectorXcd tmp;
        if(tmpBlock->isAdmissible) // far field block
        {
            tmp.noalias() = tmpBlock->UMat * (tmpBlock->singularValues.asDiagonal() * (tmpBlock->VAdjMat * x.segment(columnStartIndex - xStartIndex, numberOfColumns)));
            #pragma omp critical
            product.segment(rowStartIndex - prodStartIndex, numberOfRows) -= tmp;
        }
        else // near field block
        {
            tmp.noalias() = tmpBlock->fullMat * x.segment(columnStartIndex - xStartIndex, numberOfColumns);
            #pragma omp critical
            product.segment(rowStartIndex - prodStartIndex, numberOfRows) -= tmp;
        }
    }
}

void HArithm::recursiveVectorMatrixPoduct(Eigen::RowVectorXcd &product, const Eigen::RowVectorXcd &x, const BlockCluster &block, const long prodStartIndex, const long xStartIndex) // y += x * hmatrix
{
    if(block.isLeaf)
    {
        long rowStartIndex = block.rowStartIndex();
        long numberOfRows = block.rows();
        long columnStartIndex = block.colStartIndex();
        long numberOfColumns = block.cols();

        if(block.isAdmissible) // far field block
        {
            product.segment(columnStartIndex - prodStartIndex, numberOfColumns).noalias() += ((x.segment(rowStartIndex - xStartIndex, numberOfRows) * block.UMat) * block.singularValues.asDiagonal()) * block.VAdjMat;
        }
        else // near field block
        {
            product.segment(columnStartIndex - prodStartIndex, numberOfColumns).noalias() += x.segment(rowStartIndex - xStartIndex, numberOfRows) * block.fullMat;
        }
    }
    else
    {
        #pragma omp parallel sections
        {
            #pragma omp section //prevent aliasing by putting the same product columns in one section
            {                   //otherwise use omp critical above
                recursiveVectorMatrixPoduct(product, x, *block.son11, prodStartIndex, xStartIndex);
                recursiveVectorMatrixPoduct(product, x, *block.son21, prodStartIndex, xStartIndex);
            }
            #pragma omp section
            {
                recursiveVectorMatrixPoduct(product, x, *block.son12, prodStartIndex, xStartIndex);
                recursiveVectorMatrixPoduct(product, x, *block.son22, prodStartIndex, xStartIndex);
            }
        }
    }
}

void HArithm::subtractiveParallelVectorMatrixPoduct(Eigen::RowVectorXcd &product, const Eigen::RowVectorXcd &x, const BlockCluster &block, const long prodStartIndex, const long xStartIndex) // y -= x * hmatrix
{
    QVector<const BlockCluster*> minPartition;
    BlockCluster::getPartition(block, minPartition);

    #pragma omp parallel for
    for(long i = 0; i < minPartition.length(); i++)
    {
        const BlockCluster* tmpBlock = minPartition.at(i);
        long rowStartIndex = tmpBlock->rowStartIndex();
        long numberOfRows = tmpBlock->rows();
        long columnStartIndex = tmpBlock->colStartIndex();
        long numberOfColumns = tmpBlock->cols();

        Eigen::RowVectorXcd tmp;
        if(tmpBlock->isAdmissible) // far field block
        {
            tmp.noalias() = ((x.segment(rowStartIndex - xStartIndex, numberOfRows) * tmpBlock->UMat) * tmpBlock->singularValues.asDiagonal()) * tmpBlock->VAdjMat;
            #pragma omp critical
            product.segment(columnStartIndex - prodStartIndex, numberOfColumns).noalias() -= tmp;
        }
        else // near field block
        {
            tmp.noalias() = x.segment(rowStartIndex - xStartIndex, numberOfRows) * tmpBlock->fullMat;
            #pragma omp critical
            product.segment(columnStartIndex - prodStartIndex, numberOfColumns).noalias() -= tmp;
        }
    }
}

void HArithm::subtractiveRecursiveVectorMatrixPoduct(Eigen::RowVectorXcd &product, const Eigen::RowVectorXcd &x, const BlockCluster &block, const long prodStartIndex, const long xStartIndex) // y -= x * hmatrix
{
    if(block.isLeaf)
    {
        long rowStartIndex = block.rowStartIndex();
        long numberOfRows = block.rows();
        long columnStartIndex = block.colStartIndex();
        long numberOfColumns = block.cols();

        if(block.isAdmissible) // far field block
        {
            product.segment(columnStartIndex - prodStartIndex, numberOfColumns).noalias() -= ((x.segment(rowStartIndex - xStartIndex, numberOfRows) * block.UMat) * block.singularValues.asDiagonal()) * block.VAdjMat;
        }
        else // near field block
        {
            product.segment(columnStartIndex - prodStartIndex, numberOfColumns).noalias() -= x.segment(rowStartIndex - xStartIndex, numberOfRows) * block.fullMat;
        }
    }
    else
    {
        #pragma omp parallel sections
        {
            #pragma omp section //prevent aliasing by putting the same product columns in one section
            {                   //otherwise use omp critical above
                subtractiveRecursiveVectorMatrixPoduct(product, x, *block.son11, prodStartIndex, xStartIndex);
                subtractiveRecursiveVectorMatrixPoduct(product, x, *block.son21, prodStartIndex, xStartIndex);
            }
            #pragma omp section
            {
                subtractiveRecursiveVectorMatrixPoduct(product, x, * block.son12, prodStartIndex, xStartIndex);
                subtractiveRecursiveVectorMatrixPoduct(product, x, * block.son22, prodStartIndex, xStartIndex);
            }
        }
    }
}

void HArithm::recursiveHMatAddition(BlockCluster &mat1Block, const BlockCluster &mat2Block, const long rank, const double relError)
{
    if(mat1Block.isLeaf && mat2Block.isLeaf)
    {
        if(mat1Block.isAdmissible && mat2Block.isAdmissible)
        {
            roundedAddRMatToBlock(mat1Block, mat2Block.UMat, mat2Block.singularValues, mat2Block.VAdjMat, rank, relError);
        }
        else
        {
            if(mat1Block.isAdmissible)
            {
                addrkMatToFull(mat1Block.fullMat, mat1Block.UMat, mat1Block.singularValues, mat1Block.VAdjMat);
                clearRkMatrix(mat1Block);
                mat1Block.isAdmissible = false;
                addFullMatrixInToFirstOne(mat1Block.fullMat, mat2Block.fullMat);
            }
            else if(mat2Block.isAdmissible)
            {
//                addFullMatrixInToFirstOne(mat1Block.fullMat, mat2Block.UMat * mat2Block . singularValues.asDiagonal() * mat2Block.VAdjMat);
                addrkMatToFull(mat1Block.fullMat, mat2Block.UMat, mat2Block.singularValues, mat2Block.VAdjMat);
            }
            else
            {
                addFullMatrixInToFirstOne(mat1Block.fullMat, mat2Block.fullMat);
            }
        }
    }
    else if(mat1Block.isLeaf)
    {
        Cluster *rowCluster =  mat1Block.rowCluster;
        Cluster *colCluster =  mat1Block.columnCluster;

        #pragma omp task shared(mat1Block, mat2Block)
        mat1Block.son11 = copyBlock(*mat2Block.son11, rowCluster->son1.get(), colCluster->son1.get()); // copy blocks from matrix 2 into matrix 1,
        #pragma omp task shared(mat1Block, mat2Block)
        mat1Block.son12 = copyBlock(*mat2Block.son12, rowCluster->son1.get(), colCluster->son2.get()); // direct linking would break independence
        #pragma omp task shared(mat1Block, mat2Block)
        mat1Block.son21 = copyBlock(*mat2Block.son21, rowCluster->son2.get(), colCluster->son1.get());
//        #pragma omp task shared(mat1Block, mat2Block)
        mat1Block.son22 = copyBlock(*mat2Block.son22, rowCluster->son2.get(), colCluster->son2.get());
        #pragma omp taskwait

//        if(mat1Block.isAdmissible)
//        {
//            reducedRankAgglomerationWithTrimming(mat1Block, rank, relError);
//        }
//        else
//        {
////            fullRankAgglomerationWithTrimming(mat1Block);
//            fullRankDirektAgglomerationWithTrimming(mat1Block);
//        }
        blockSplitting(mat1Block, rank, relError);
    }
    else if (mat2Block.isLeaf)
    {
        if(mat2Block.isAdmissible)
        {
//            addRMatToBlock(*mat1Block, mat2Block.UMat, mat2Block.singularValues, mat2Block.VAdjMat);
//            reducedRankAgglomerationWithTrimming(mat1Block, rank, relError);
            roundedAddRMatToBlock(mat1Block, mat2Block.UMat, mat2Block.singularValues, mat2Block.VAdjMat, rank, relError);
        }
        else
        {
            addFullMatrixInToFirstOne(mat1Block.fullMat, mat2Block.fullMat);
//            fullRankDirektAgglomerationWithTrimming(mat1Block);
        }
        blockSplitting(mat1Block, rank, relError);
    }
    else
    {
        #pragma omp task shared(mat1Block, mat2Block)
        recursiveHMatAddition(*mat1Block.son11, *mat2Block.son11, rank, relError);
        #pragma omp task shared(mat1Block, mat2Block)
        recursiveHMatAddition(*mat1Block.son12, *mat2Block.son12, rank, relError);
        #pragma omp task shared(mat1Block, mat2Block)
        recursiveHMatAddition(*mat1Block.son21, *mat2Block.son21, rank, relError);
//        #pragma omp task
        recursiveHMatAddition(*mat1Block.son22, *mat2Block.son22, rank, relError);
        #pragma omp taskwait
    }
}

void HArithm::recursiveHMatSubstraction(BlockCluster &mat1Block, const BlockCluster &mat2Block, const long rank, const double relError) // mat1Block - mat2Block
{
    if(mat1Block.isLeaf && mat2Block.isLeaf)
    {
        if(mat1Block.isAdmissible && mat2Block.isAdmissible)
        {
            roundedAddRMatToBlock(mat1Block, mat2Block.UMat, -mat2Block.singularValues, mat2Block.VAdjMat, rank, relError);
        }
        else
        {
            if(mat1Block.isAdmissible)
            {
               addrkMatToFull(mat1Block.fullMat, mat1Block.UMat, mat1Block.singularValues, mat1Block.VAdjMat);
               clearRkMatrix(mat1Block);
               mat1Block.isAdmissible = false;
               addFullMatrixInToFirstOne(mat1Block.fullMat, -mat2Block.fullMat);
            }
            else if(mat2Block.isAdmissible)
            {
//                addFullMatrixInToFirstOne(mat1Block.fullMat, mat2Block.UMat * ( -mat2Block.singularValues).asDiagonal() * mat2Block.VAdjMat);
                addrkMatToFull(mat1Block.fullMat, mat2Block.UMat, -mat2Block.singularValues, mat2Block.VAdjMat);
            }
            else
            {
                addFullMatrixInToFirstOne(mat1Block.fullMat, -mat2Block.fullMat);
            }
        }
    }
    else if(mat1Block.isLeaf)
    {
        Cluster* rowCluster =  mat1Block.rowCluster;
        Cluster* colCluster =  mat1Block.columnCluster;

//        #pragma omp parallel sections
//        {
            #pragma omp task shared(mat1Block, mat2Block)
            {
                mat1Block.son11 = copyBlock(*mat2Block.son11, rowCluster->son1.get(), colCluster->son1.get()); // copy blocks from matrix 2 into matrix 1,
                mat1Block.son11->father = &mat1Block;
                recursiveMultiplyHMatByMinusOne(*mat1Block.son11);
            }
            #pragma omp task shared(mat1Block, mat2Block)
            {
                mat1Block.son12 = copyBlock(*mat2Block.son12, rowCluster->son1.get(), colCluster->son2.get()); // direct linking would break independence
                mat1Block.son12->father = &mat1Block;
                recursiveMultiplyHMatByMinusOne(*mat1Block.son12);
            }
            #pragma omp task shared(mat1Block, mat2Block)
            {
                mat1Block.son21 = copyBlock(*mat2Block.son21, rowCluster->son2.get(), colCluster->son1.get());
                mat1Block.son21->father = &mat1Block;
                recursiveMultiplyHMatByMinusOne(*mat1Block.son21);
            }
//            #pragma omp section
            {
                mat1Block.son22 = copyBlock(*mat2Block.son22, rowCluster->son2.get(), colCluster->son2.get());
                mat1Block.son22->father = &mat1Block;
                recursiveMultiplyHMatByMinusOne(*mat1Block.son22);
            }
//        }
        #pragma omp taskwait
        blockSplitting(mat1Block, rank, relError);
    }
    else if (mat2Block.isLeaf)
    {
        if(mat2Block.isAdmissible)
        {
            roundedAddRMatToBlock(mat1Block, mat2Block.UMat, -mat2Block.singularValues, mat2Block.VAdjMat, rank, relError);
        }
        else
        {
            addFullMatrixInToFirstOne(mat1Block.fullMat, -mat2Block.fullMat);
        }
        blockSplitting(mat1Block, rank, relError);
    }
    else
    {
//        #pragma omp parallel sections
//        {
            #pragma omp task shared(mat1Block, mat2Block)
            recursiveHMatSubstraction(*mat1Block.son11, *mat2Block.son11, rank, relError);
            #pragma omp task shared(mat1Block, mat2Block)
            recursiveHMatSubstraction(*mat1Block.son12, *mat2Block.son12, rank, relError);
            #pragma omp task shared(mat1Block, mat2Block)
            recursiveHMatSubstraction(*mat1Block.son21, *mat2Block.son21, rank, relError);

            recursiveHMatSubstraction(*mat1Block.son22, *mat2Block.son22, rank, relError);
            #pragma omp taskwait
//        }
    }
}

void HArithm::multiplyHMatByMinusOne(HMatrix &hMat)
{
    if(hMat.getRootBlock() != nullptr)
    {
        recursiveMultiplyHMatByMinusOne(*hMat.getRootBlock());
    }
}

void HArithm::recursiveMultiplyHMatByMinusOne(BlockCluster &matBlock)
{
    if(matBlock.isLeaf)
    {
        if(matBlock.isAdmissible)
        {
            matBlock.singularValues *= -1;
        }
        else
        {
            matBlock.fullMat *= -1;
        }
    }
    else
    {
        #pragma omp parallel sections
        {
            #pragma omp section
            recursiveMultiplyHMatByMinusOne(*matBlock.son11);
            #pragma omp section
            recursiveMultiplyHMatByMinusOne(*matBlock.son12);
            #pragma omp section
             recursiveMultiplyHMatByMinusOne(*matBlock.son21);
            #pragma omp section
            recursiveMultiplyHMatByMinusOne(*matBlock.son22);
        }
    }
}

void HArithm::convertToAdmissibleZeroBlock(BlockCluster &matBlock)
{
    setRkMatrixZero(matBlock);
    matBlock.isLeaf = true;
    matBlock.isAdmissible = true;
}

std::unique_ptr<BlockCluster> HArithm::copyBlock(const BlockCluster &block, Cluster* rowCluster, Cluster* colCluster)
{
    std::unique_ptr<BlockCluster> returnBlock = std::make_unique<BlockCluster>();;

    returnBlock->isAdmissible = block.isAdmissible;
    returnBlock->isRoot = block.isRoot;
    returnBlock->isLeaf = block.isLeaf;

    returnBlock->frobeniusNorm = block.frobeniusNorm;

    returnBlock->fullMat = block.fullMat;
    returnBlock->UMat = block.UMat;
    returnBlock->singularValues = block.singularValues;
    returnBlock->VAdjMat = block.VAdjMat;

    returnBlock->rowCluster = rowCluster;
    returnBlock->columnCluster = colCluster;
    if(!returnBlock->isLeaf)
    {
        if(rowCluster->isLeaf || rowCluster->isLeaf)
        {
            std::cerr << "if(rowCluster.isLeaf || rowCluster.isLeaf) in copyBlock() call" << std::endl;
        }
        #pragma omp parallel sections
        {
            #pragma omp section
            {
                returnBlock->son11 = copyBlock(*block.son11, rowCluster->son1.get(), colCluster->son1.get());
                returnBlock->son11->father = returnBlock.get();
            }
            #pragma omp section
            {
                returnBlock->son12 = copyBlock(*block.son12, rowCluster->son1.get(), colCluster->son2.get());
                returnBlock->son12->father = returnBlock.get();
            }
            #pragma omp section
            {
                returnBlock->son21 = copyBlock(*block.son21, rowCluster->son2.get(), colCluster->son1.get());
                returnBlock->son21->father = returnBlock.get();
            }
            #pragma omp section
            {
                returnBlock->son22 = copyBlock(*block.son22, rowCluster->son2.get(), colCluster->son2.get());
                returnBlock->son22->father = returnBlock.get();
            }
        }
    }
    return returnBlock;
}

Eigen::VectorXcd HArithm::getRankKBlockDiagonal(const BlockCluster &matBlock, long rank)
{
    long diagonalLength = std::min(matBlock.rows(), matBlock.cols());

    if(matBlock.isAdmissible)
    {
        if(matBlock.singularValues.size() == 0)
        {
            return Eigen::VectorXcd::Zero(diagonalLength);
        }

        if(rank <= 0)
        {
            rank = matBlock.singularValues.size();
        }

        Eigen::VectorXcd diagonal(diagonalLength);
        for(long i = 0; i < diagonalLength; i++)
        {
            diagonal(i) = ((matBlock.UMat.row(i).head(rank).transpose().array() * matBlock.singularValues.head(rank).array()).array() * matBlock.VAdjMat.col(i).head(rank).array()).sum();
        }
        return diagonal;
    }
    else
    {
        if(matBlock.fullMat.size() != 0)
        {
            return matBlock.fullMat.diagonal();
        }
        else
        {
            return Eigen::VectorXcd::Zero(diagonalLength);
        }
    }
}

Eigen::VectorXcd HArithm::getRankKBlockIndexedElements(BlockCluster &matBlock, QVector<long> &rowIndices, QVector<long> &colIndices, long rank)
{
    long returnVectorLenght = std::min(rowIndices.size(), colIndices.size());
    if(matBlock.isAdmissible)
    {
        if(matBlock.singularValues.size() == 0)
        {
            return Eigen::VectorXcd::Zero(returnVectorLenght);
        }

        if(rank <= 0)
        {
            rank = matBlock.singularValues.size();
        }

        Eigen::VectorXcd returnVector(returnVectorLenght);

        for(long i = 0; i < returnVectorLenght; i++)
        {
            returnVector(i) = ((matBlock.UMat.row(rowIndices.at(i)).head(rank).transpose().array() * matBlock.singularValues.head(rank).array()).array() * matBlock.VAdjMat.col(colIndices.at(i)).head(rank).array()).sum();
        }
        return returnVector;
    }
    else
    {
        Eigen::VectorXcd returnVector(returnVectorLenght);

        if(matBlock.fullMat.size() != 0)
        {
            for(long i = 0; i < returnVectorLenght; i++)
            {
                returnVector(i) = matBlock.fullMat(rowIndices.at(i), colIndices.at(i));
            }
            return returnVector;
        }
        else
        {
            return Eigen::VectorXcd::Zero(returnVectorLenght);
        }
    }
}

double HArithm::spectralNorm(HMatrix &A, const unsigned long maxIterations) // Approximates the spectral norm a matrix via power iteration.
{
// for reference: J. Kuczyński and H. Woźniakowski, "Estimating the largest eigenvalue by the power and Lanczos algorithms with a random start"
    double normEstimate = 0;
    if(A.cols() != A.rows())
    {
        std::cerr << "Non-square matrix in spectralNorm() routine." << std::endl;
        return 0;
    }
    long dim = A.cols();
    Eigen::RowVectorXcd u = Eigen::RowVectorXcd::Zero(dim);
    Eigen::VectorXcd v = Eigen::VectorXcd::Random(dim);
    v.normalize();

    unsigned long iterCounter = 0;
    while(iterCounter < maxIterations)
    {
        Eigen::RowVectorXcd uTmp = Eigen::VectorXcd::Zero(dim);
        VMM(uTmp, v.conjugate(), A);
        u = uTmp;
        v = Eigen::VectorXcd::Zero(dim)/dim;
        MVM(v, A, u.adjoint());
        normEstimate = std::sqrt(v.norm());
        v.normalize();
        iterCounter++;
    }
    return normEstimate;
}

double HArithm::spectralNormFromLU(HMatrix &L, HMatrix &U, const unsigned long maxIterations) // Approximates the spectral norm of the inverse of a matrix via the LU-decomposition of the matrix. The method is a power iteration with forward- and backward substitution.
{
    if(L.rows() != L.cols() || L.rows() != U.rows() || L.rows() != U.cols())
    {
        std::cerr << "Incompatible matrix dimensions in spectralNormFromLU() call! " << std::endl;
        return 0;
    }
    double normEstimate = 0;
    long dim = L.cols();
    Eigen::RowVectorXcd u = Eigen::RowVectorXcd::Zero(dim);
    Eigen::VectorXcd v = Eigen::VectorXcd::Random(dim);
    v.normalize();

    unsigned long iterCounter = 0;
    while(iterCounter < maxIterations)
    {
        Eigen::RowVectorXcd uTmp = Eigen::RowVectorXcd::Zero(dim);
        Eigen::RowVectorXcd vAdj = v.adjoint();
        HArithm::forwardSubstitutionTransposed(uTmp, *U.getRootBlock(), vAdj); // Solve uTmp * U = v^* for uTmp1; U is upper triangular H-matrix.
        u = Eigen::RowVectorXcd::Zero(dim);
        HArithm::backwardSubstitutionTransposed(u, *L.getRootBlock(), uTmp); // Solve u * L = uTmp for u, L is lower triangular H-matrix

        Eigen::VectorXcd vTmp = Eigen::VectorXcd::Zero(dim);
        Eigen::VectorXcd uAdj = u.adjoint();
        HArithm::forwardSubstitution(*L.getRootBlock(), vTmp, uAdj); // Solve L * vTmp = u^* for vTmp, L is lower triangular H-matrix
        v = Eigen::VectorXcd::Zero(dim);
        HArithm::backwardSubstitution(*U.getRootBlock(), v, vTmp); // Solve U * v = vTmp for v, U is upper triangular H-matrix

        normEstimate = std::sqrt(v.norm());
        v.normalize();
        iterCounter++;

//        std::cout << "normEstimate: " << normEstimate << std::endl;
//        std::cout << "iterCounter: " << iterCounter << std::endl;
    }
//    Eigen::MatrixXcd mL = HArithm::hMatToFullMat(L);
//    Eigen::MatrixXcd mU = HArithm::hMatToFullMat(U);
//    std::cout << "real spectral lu norm: " << (mL*mU).inverse().operatorNorm() << std::endl;
//    std::cout << "spectral lu norm estimate: " << normEstimate << std::endl;
    return normEstimate;
}

double HArithm::spectralNormLeastSignificant(HMatrix &A, const unsigned long maxIterations) // Approximates the spectral norm of a matrix, where only the least significant rank-1 matrix of each block is taken into account. Can be used as an approximation of the spectral norm of the residual matrix after ACA assembly
{
    double normEstimate = 0;
    if(A.cols() != A.rows())
    {
        std::cerr << "Non-square matrix in spectralNorm() routine." << std::endl;
        return 0;
    }
    long dim = A.cols();
    Eigen::RowVectorXcd u = Eigen::RowVectorXcd::Zero(dim);
    Eigen::VectorXcd v = Eigen::VectorXcd::Random(dim);
    v.normalize();

    unsigned long iterCounter = 0;
    while(iterCounter < maxIterations)
    {
        Eigen::RowVectorXcd uTmp = Eigen::VectorXcd::Zero(dim);
        VMMLeastSignificant(uTmp, v.conjugate(), A);
        u = uTmp;
        v = Eigen::VectorXcd::Zero(dim)/dim;
        MVMLeastSignificant(v, A, u.adjoint());
        normEstimate = std::sqrt(v.norm());
        v.normalize();
        iterCounter++;
    }
    return normEstimate;
}

void HArithm::MVMLeastSignificant(Eigen::VectorXcd &y, HMatrix &hmatrix, const Eigen::VectorXcd &x)// vector-H-matrix product. y += hmatrix * x; only least significant rank-1 matrix of each admissible block is accounted for
{
    QVector<BlockCluster*> minPartition = hmatrix.getMinPartition();
    #pragma omp parallel for
    for(long i = 0; i < minPartition.length(); i++)
    {
        BlockCluster* tmpBlock = minPartition.at(i);
        long rowStartIndex = tmpBlock->rowStartIndex();
        long numberOfRows = tmpBlock->rows();
        long columnStartIndex = tmpBlock->colStartIndex();
        long numberOfColumns = tmpBlock->cols();
        long rank = tmpBlock->singularValues.size();
        Eigen::VectorXcd tmp;
        if(tmpBlock->isAdmissible) // far field block
        {
            tmp.noalias() = tmpBlock->UMat.col(rank-1) * (tmpBlock->singularValues(rank-1) * (tmpBlock->VAdjMat.row(rank-1) * x.segment(columnStartIndex, numberOfColumns)));
            #pragma omp critical
            y.segment(rowStartIndex, numberOfRows) += tmp;
        }
        else // near field block
        {
            continue;
        }
    }
}

void HArithm::VMMLeastSignificant(Eigen::RowVectorXcd &y, const Eigen::RowVectorXcd &x, HMatrix &hmatrix) // vector-H-matrix product. y' += x' * hmatrix; only least significant rank-1 matrix of each admissible block is accounted for
{
    QVector<BlockCluster*> minPartition = hmatrix.getMinPartition();
    #pragma omp parallel for
    for(long i = 0; i < minPartition.length(); i++)
    {
        BlockCluster* tmpBlock = minPartition.at(i);
        long rowStartIndex = tmpBlock->rowStartIndex();
        long numberOfRows = tmpBlock->rows();
        long columnStartIndex = tmpBlock->colStartIndex();
        long numberOfColumns = tmpBlock->cols();
        long rank = tmpBlock->singularValues.size();
        Eigen::RowVectorXcd tmp;
        if(tmpBlock->isAdmissible) // far field block
        {
            tmp.noalias() = ((x.segment(rowStartIndex, numberOfRows) * tmpBlock->UMat.col(rank-1)) * tmpBlock->singularValues(rank-1)) * tmpBlock->VAdjMat.row(rank-1);
            #pragma omp critical
            y.segment(columnStartIndex, numberOfColumns) += tmp;
        }
        else // near field block
        {
            continue;
        }
    }
}

double HArithm::frobeniusNormFromLU(HMatrix &L, HMatrix &U, const unsigned long maxIterations) // Approximates the frobenius norm of the inverse of a matrix via the LU-decomposition of the matrix. Low accuracy.
{
    if(L.rows() != L.cols() || L.rows() != U.rows() || L.rows() != U.cols())
    {
        std::cerr << "Incompatible matrix dimensions in frobeniusNormFromLU() call! " << std::endl;
        return 0;
    }
    long dim = L.cols();

    // We approximate (LU)^-1;
    // (LU)^-1 ~= VAdjMat * singularValues * UMat
    // We assemble VAdjMat and UMat via cross approximation; the row- and columnvectors are calculated via forward- and backwardsubstitution
    Eigen::MatrixXcd VAdjMat;
    Eigen::MatrixXcd UMat;
    Eigen::VectorXcd singularValues;

    std::complex<double> lowRankNormSquared = 0; // full rank matrix norm approximation by the R1 matrix norm multiplied by relative error
    long maxRank;

    if(maxIterations > 0)
    {
        maxRank = maxIterations;
        UMat.resize(dim, maxRank);
        singularValues.resize(maxRank);
        VAdjMat.resize(maxRank, dim);
    }
    else
    {
        std::cout << "frobeniusNormFromLU called without explicit iteration limit." << std::endl;
        maxRank = dim;
        long reservationRank = std::min((long)5, maxRank);
        UMat.resize(dim, reservationRank);
        singularValues.resize(reservationRank);
        VAdjMat.resize(reservationRank, dim);
    }
    long tmpRank = 0;

    Eigen::RowVectorXcd rowVector = Eigen::RowVectorXcd::Zero(dim);
    Eigen::VectorXcd columnVector = Eigen::VectorXcd::Zero(dim);

    long rowIndex;
    long columnIndex;

    Eigen::RowVectorXcd unitRowVector(dim);
    Eigen::VectorXcd unitColumnVector = Eigen::VectorXcd::Zero(dim);
    columnIndex = QRandomGenerator::system()->bounded((qint32) dim);
    unitColumnVector(columnIndex) = 1;


    Eigen::VectorXcd tmpCol = Eigen::VectorXcd::Zero(dim);
    HArithm::forwardSubstitution(* L.getRootBlock(), tmpCol, unitColumnVector); // Solve L * tmpCol = unitColumnVector for tmpCol, L is lower triangular H-matrix
    columnVector = Eigen::VectorXcd::Zero(dim);
    HArithm::backwardSubstitution(* U.getRootBlock(), columnVector, tmpCol); // Solve U * columnVector = tmpCol for columnVector, U is upper triangular H-matrix

    std::complex<double> lastGain = 0;
    while(tmpRank < maxRank)
    {
        if(tmpRank >= 1)
        {
            columnVector(rowIndex) = 0;
        }
        else //generate new random row index
        {
//            rowIndex = QRandomGenerator::global()->bounded((qint32) blockRows); // find random row index
            rowIndex = QRandomGenerator::system()->bounded((qint32) dim); // find random row index
        }
        columnVector.cwiseAbs().maxCoeff(&rowIndex);

        if(std::abs(columnVector(rowIndex)) <= global::tiny && tmpRank >= 1) //the max absolute element is too small -> no suitable column index for rowIndex
        {
            std::cerr << "max(columnVector) <= global::Tiny" << std::endl;
            break;    //change later
        }

        // calculate the row rowIndex of the inverse of L*U
        unitRowVector = Eigen::RowVectorXcd::Zero(dim);
        unitRowVector(rowIndex) = 1;
        Eigen::RowVectorXcd tmpRow = Eigen::RowVectorXcd::Zero(dim);
        HArithm::forwardSubstitutionTransposed(tmpRow, * U.getRootBlock(), unitRowVector); // Solve tmpRow * U = unitRowVector for tmpRow; U is upper triangular H-matrix.
        rowVector = Eigen::RowVectorXcd::Zero(dim);
        HArithm::backwardSubstitutionTransposed(rowVector, * L.getRootBlock(), tmpRow); // Solve u * L = uTmp for u, L is lower triangular H-matrix

        if(tmpRank >= 1)
        {
            rowVector -= (UMat.leftCols(tmpRank).row(rowIndex) * singularValues.head(tmpRank).asDiagonal()) * VAdjMat.topRows(tmpRank);
        }

         // find column index for largest absolute element in row with above rowindex
        double maxRowVal = rowVector.cwiseAbs().maxCoeff(&columnIndex); //  columnIndex is now max column index of phiAbsVec

        if(maxRowVal <= global::tiny && tmpRank >= 1) //the max absolute element is too small -> no suitable column index for rowIndex
        {
            std::cerr << "max(rowAbsVec) <= Tiny" << std::endl;
            break;
        }

        // calculate the block column columnIndex
        unitColumnVector = Eigen::VectorXcd::Zero(dim);
        unitColumnVector(columnIndex) = 1;
        tmpCol = Eigen::VectorXcd::Zero(dim);
        HArithm::forwardSubstitution(* L.getRootBlock(), tmpCol, unitColumnVector); // Solve L * tmpCol = unitColumnVector for tmpCol, L is lower triangular H-matrix
        columnVector = Eigen::VectorXcd::Zero(dim);
        HArithm::backwardSubstitution(* U.getRootBlock(), columnVector, tmpCol); // Solve U * columnVector = tmpCol for columnVector, U is upper triangular H-matrix

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

        std::complex<double> oldNormSquared2 = lowRankNormSquared;
        lowRankNormSquared += std::pow(std::abs(alpha), 2) *  rowVector.squaredNorm() * columnVector.squaredNorm();
        for(long i = 0; i < tmpRank - 1; i++)
        {
            lowRankNormSquared += 2.0 * std::conj(alpha) * singularValues(i) * (columnVector.adjoint() * UMat.col(i))(0,0) * (( VAdjMat.row(i)) * (rowVector).adjoint())(0,0);
        }
        lastGain = lowRankNormSquared - oldNormSquared2;
    }

    double finalNormEstimate = std::sqrt(std::abs(std::real(lowRankNormSquared + (double)(dim - tmpRank) * lastGain)));

    std::cout << "frobenius norm estimate of the inverse: " << finalNormEstimate << std::endl;

    // test
//    Eigen::MatrixXcd fullMat = HArithm::hMatToFullMat(L) * HArithm::hMatToFullMat(U);
//    std::cout << "frobenius norm estimate of the inverse (accurate): " << fullMat.inverse().norm() << std::endl;

    return finalNormEstimate; // return frobenius norm of low rank approximation of (L*U)^-1
}

void HArithm::transpose(HMatrix &matrix)
{
    std::shared_ptr<ClusterTree> rowClustertree = matrix.getRowClustertree();
    std::shared_ptr<ClusterTree> columnClustertree = matrix.getColumnClustertree();

    matrix.setClusterTrees(columnClustertree, rowClustertree);
    if(matrix.getRootBlock() != nullptr)
    {
        transpose(* matrix.getRootBlock());
    }
}

void HArithm::transpose(BlockCluster &matBlock)
{
    std::swap(matBlock.rowCluster, matBlock.columnCluster);

    if(matBlock.isLeaf)
    {
        if(matBlock.isAdmissible)
        {
            std::swap(matBlock.UMat,matBlock.VAdjMat);
            matBlock.UMat.transposeInPlace();
            matBlock.VAdjMat.transposeInPlace();
        }
        else
        {
            matBlock.fullMat.transposeInPlace();
        }
    }
    else
    {
        std::swap(matBlock.son21, matBlock.son12);
        transpose(* matBlock.son11);
        transpose(* matBlock.son12);
        transpose(* matBlock.son21);
        transpose(* matBlock.son22);
    }
}

