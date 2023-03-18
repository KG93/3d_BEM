#include "hmultiply.h"

//long intermBlCl::rows() const
//{
//    if(rowCluster != nullptr)
//    {
//        return rowCluster->triangleIndexes.last() - rowCluster->triangleIndexes.first() + 1;
//    }
//    else
//    {
//        return 0;
//    }
//}

//long intermBlCl::cols() const
//{
//    if(columnCluster != nullptr)
//    {
//        return columnCluster->triangleIndexes.last() - columnCluster->triangleIndexes.first() + 1;
//    }
//    else
//    {
//        return 0;
//    }
//}

//long intermBlCl::size() const
//{
//    return rows() * cols();
//}

//long intermBlCl::rowStartIndex() const
//{
//    if(rowCluster != nullptr)
//    {
//        return rowCluster->triangleIndexes.first();
//    }
//    else
//    {
//        return 0;
//    }
//}

//long intermBlCl::colStartIndex() const
//{
//    if(columnCluster != nullptr)
//    {
//        return columnCluster->triangleIndexes.first();
//    }
//    else
//    {
//        return 0;
//    }
//}

//void intermBlCl::clearMatrixInfo()
//{
//    UMat.resize(0,0);
//    singularValues.resize(0);
//    VAdjMat.resize(0,0);
//    fullMat.resize(0,0);
//}

void intermBlCl::clear()
{
    clearMatrixInfo();

    fullFactorsList.clear();
    rKFactorsList.clear();
    fullFlushMergeTargets.clear();
    rKFlushTargets.clear();
    if(intSon11 != nullptr)
    {
        intSon11->clear();
        delete intSon11;
        intSon11 = nullptr;
    }
    if(intSon12 != nullptr)
    {
        intSon12->clear();
        delete intSon12;
        intSon12 = nullptr;
    }
    if(intSon21 != nullptr)
    {
        intSon21->clear();
        delete intSon21;
        intSon21 = nullptr;
    }
    if(intSon22 != nullptr)
    {
        intSon22->clear();
        delete intSon22;
        intSon22 = nullptr;
    }
}

void intermBlCl::trimBelow()
{
    isLeaf = true;
    isEventualLeaf = true;
    if(intSon11 != nullptr)
    {
        intSon11->clear();
        delete intSon11;
        intSon11 = nullptr;
    }
    if(intSon12 != nullptr)
    {
        intSon12->clear();
        delete intSon12;
        intSon12 = nullptr;
    }
    if(intSon21 != nullptr)
    {
        intSon21->clear();
        delete intSon21;
        intSon21 = nullptr;
    }
    if(intSon22 != nullptr)
    {
        intSon22->clear();
        delete intSon22;
        intSon22 = nullptr;
    }
}

void intermBlCl::setRKFactors(const QVector<std::pair<BlockCluster *, BlockCluster *> > &newRKFactors)
{
    rKFactorsList = newRKFactors;
}

HMultiply::HMultiply()
{

}

HMatrix HMultiply::multiplyHMat(BlockCluster &factor1, BlockCluster &factor2, const long rank, const double relError)
{
    intermBlCl* productRootBlockCluster = new intermBlCl(factor1.rowCluster, factor2.columnCluster, true);
    MM(factor1, factor2, productRootBlockCluster); // create product tree and save factor blocks into the product nodes
    recursProdPartition(productRootBlockCluster, &factor1, &factor2);  // determine the eventual leave nodes in the product tree

    QVector<intermBlCl*> fullTransitBlocks;
    findFullFlushTargets(productRootBlockCluster, fullTransitBlocks); // connect all product tree nodes with information with their respective eventual leaf nodes

    findRkFlushTargets(productRootBlockCluster);

    QVector<intermBlCl*> independentWorkBlocks;
    findIndependentWorkBlocks(productRootBlockCluster, independentWorkBlocks);

    startWorkOnIndependentBlocks(independentWorkBlocks, rank, relError);

    // convert the intermittent blocks into regular blocks
    BlockCluster* rootBlockCluster = intermBlClToBlockCluster(productRootBlockCluster);
    rootBlockCluster->isRoot = true;
    HMatrix product(rootBlockCluster);

    if(product.consistencyCheck() == false)
    {
        std::cerr << "product.consistencyCheck() == false" << std::endl;
    }
    product.updatePartition();
    return product;
}

HMatrix HMultiply::multiplyHMat(HMatrix &factor1, HMatrix &factor2, const long rank, const double relError)
{
    if( !factor1.isCrossPartitioned() || !factor2.isCrossPartitioned())
    {
        std::cerr<<"At least one factor in multiplyHMat() is not cross partitioned!" << std::endl;
    }
    if(factor1.getRootBlock() == nullptr || factor2.getRootBlock() == nullptr)
    {
        std::cerr<<"At least one factor in multiplyHMat() is an empty h-matrix!" << std::endl;
    }
    if(factor1.getRootBlock()->columnCluster->indices.first() != factor2.getRootBlock()->rowCluster->indices.first() || factor1.getRootBlock()->columnCluster->indices.last() != factor2.getRootBlock()->rowCluster->indices.last())
    {
        std::cerr<<"Factor matrices of incompatible sizes in multiplyHMat() call!" << std::endl;
    }
    if(factor1.consistencyCheck() == false)
    {
        std::cerr << "factor1.consistencyCheck() == false" << std::endl;
    }
    if(factor2.consistencyCheck() == false)
    {
        std::cerr << "factor2.consistencyCheck() == false" << std::endl;
    }

    HMatrix product = multiplyHMat(* factor1.getRootBlock(), * factor2.getRootBlock(), rank, relError);
    product.setClusterTrees(factor1.getRowClustertree(), factor2.getColumnClustertree());

    return product;
}

void HMultiply::MM(/*QVector<intermBlCl*> &blocksWithMatrixLoad,*/ BlockCluster &factorBlock1, BlockCluster &factorBlock2, intermBlCl* productBlock)
{
    if( !factorBlock1.isLeaf && !factorBlock2.isLeaf) // no factor block is neither full matrix nor reduced-rank matrix
    {
        if(factorBlock1.colStartIndex() != factorBlock2.rowStartIndex() || factorBlock1.cols() != factorBlock2.rows())
        {
            std::cerr<<"Incompatible blocks in MM() call!" << std::endl;
            return;
        }
        if(productBlock->intSon11 == nullptr) //son pointers shall not be overriden in following MM() call
        {
            productBlock->intSon11 = new intermBlCl(factorBlock1.rowCluster->son1, factorBlock2.columnCluster->son1, productBlock);
//        }
//        if(productBlock->son12 == nullptr) //son pointers shall not be overriden in following MM() call
//        {
            productBlock->intSon12 = new intermBlCl(factorBlock1.rowCluster->son1, factorBlock2.columnCluster->son2, productBlock);
//        }
//        if(productBlock->son21 == nullptr) //son pointers shall not be overriden in following MM() call
//        {
            productBlock->intSon21 = new intermBlCl(factorBlock1.rowCluster->son2, factorBlock2.columnCluster->son1, productBlock);
//        }
//        if(productBlock->son22 == nullptr) //son pointers shall not be overriden in following MM() call
//        {
            productBlock->intSon22 = new intermBlCl(factorBlock1.rowCluster->son2, factorBlock2.columnCluster->son2, productBlock);
        }

//        #pragma omp parallel
        {
            //omp_set_nested(0);

//            #pragma set OMP_NESTED = FALSE;

//            #pragma omp sections
            {
//                #pragma omp section
                {
                    MM(* factorBlock1.son11, * factorBlock2.son11, productBlock->intSon11);
                    MM(* factorBlock1.son12, * factorBlock2.son21, productBlock->intSon11);
                }
//                #pragma omp section
                {
                    MM(* factorBlock1.son11, * factorBlock2.son12, productBlock->intSon12);
                    MM(* factorBlock1.son12, * factorBlock2.son22, productBlock->intSon12);
                }
//                #pragma omp section
                {
                    MM(* factorBlock1.son21, * factorBlock2.son11, productBlock->intSon21);
                    MM(* factorBlock1.son22, * factorBlock2.son21, productBlock->intSon21);
                }
//                #pragma omp section
                {
                    MM(* factorBlock1.son21, * factorBlock2.son12, productBlock->intSon22);
                    MM(* factorBlock1.son22, * factorBlock2.son22, productBlock->intSon22);
                }
            }
        }
    }
    else
    {
        if((factorBlock1.isLeaf && factorBlock1.isAdmissible) || (factorBlock2.isLeaf && factorBlock2.isAdmissible)) // maybe change to &&
        { // the factor pair results in a reduced rank matrix
            productBlock->rKFactorsList.append( std::pair<BlockCluster*, BlockCluster*> {&factorBlock1, &factorBlock2});
            productBlock->holdsRKInformation = true;
        }
        else // the factor pair results in a full matrix
        {
            productBlock->fullFactorsList.append( std::pair<BlockCluster*, BlockCluster*> {&factorBlock1, &factorBlock2});
            if(!productBlock->isAdmissible)
            {
                productBlock->holdsFullInformation = true;
            }
            else
            {
                productBlock->holdsRKInformation = true;
            }
        }
    }
}

void HMultiply::findFullFlushTargets(intermBlCl* block, QVector<intermBlCl*> fullTransitBlocks)
{
    if(block->isEventualLeaf)
    {
        for(int i = 0; i < fullTransitBlocks.length(); i++)
        {
            fullTransitBlocks[i]->fullFlushMergeTargets.append(block);
        }
        if(block->intSon11 != nullptr) //block has information below, that has to be transported up
        {
            setFullMergeTarget(block, block->intSon11);
            setFullMergeTarget(block, block->intSon12);
            setFullMergeTarget(block, block->intSon21);
            setFullMergeTarget(block, block->intSon22);
        }
    }
    else
    {
        if(block->holdsFullInformation) //block holds information
        {
            fullTransitBlocks.append(block); //the block also needs a target block
        }
        if(block->intSon11 != nullptr)
        {
            findFullFlushTargets(block->intSon11, fullTransitBlocks);
            findFullFlushTargets(block->intSon12, fullTransitBlocks);
            findFullFlushTargets(block->intSon21, fullTransitBlocks);
            findFullFlushTargets(block->intSon22, fullTransitBlocks);
        }
        else
        {
            std::cerr << "block is leaf in intermediate tree but not in eventual partition in findTargetBlocksBelow call" << std::endl;
        }
    }
}

void HMultiply::findRkFlushTargets(intermBlCl* block, intermBlCl* transitBlock)
{
    if(transitBlock != nullptr) // there is already information travelling downwards the tree
    {
        if(block->isEventualLeaf) // the downwards transitioning block information has found it's final target
        {
            transitBlock->rKFlushTargets.append(block);
        }
        else if(block->holdsRKInformation)
        {
            transitBlock->rKFlushTargets.append(block); // the downwards transitioning block information has found it's intermediate target
            findRkFlushTargets(block->intSon11, block);    // -> the information of transitBlock will be added to block
            findRkFlushTargets(block->intSon12, block);    // -> the the result of above operation has to be flushed as well
            findRkFlushTargets(block->intSon21, block);    // -> block is the new transitBlock
            findRkFlushTargets(block->intSon22, block);
        }
        else // the downwards-travelling information has to go down further
        {
            findRkFlushTargets(block->intSon11, transitBlock);
            findRkFlushTargets(block->intSon12, transitBlock);
            findRkFlushTargets(block->intSon21, transitBlock);
            findRkFlushTargets(block->intSon22, transitBlock);
        }
    }
    else
    {
        if(block->isEventualLeaf)
        {
            return; // one could start organizing the merge operation for a possible subtree of the block here
        }
        else
        {
            if(block->holdsRKInformation) // block is not leaf and holds an rk-martix -> block needs to be flushed downwards
            {
                findRkFlushTargets(block->intSon11, block);
                findRkFlushTargets(block->intSon12, block);
                findRkFlushTargets(block->intSon21, block);
                findRkFlushTargets(block->intSon22, block);
            }
            else // *transitBlock is nullpointer -> had no information; block also has no rK information -> call findRkFlushTargets on subtree without transit block
            {
                findRkFlushTargets(block->intSon11, nullptr);
                findRkFlushTargets(block->intSon12, nullptr);
                findRkFlushTargets(block->intSon21, nullptr);
                findRkFlushTargets(block->intSon22, nullptr);
            }
        }
    }
}

void HMultiply::setFullMergeTarget(intermBlCl* targetBlock, intermBlCl* transitBlock)
{
    if(transitBlock->holdsFullInformation) // transit block holds full matrix information
    {
        transitBlock->fullFlushMergeTargets.resize(1);
        transitBlock->fullFlushMergeTargets[0] = targetBlock;
    }
    if(transitBlock->intSon11 != nullptr)
    {
        setFullMergeTarget(targetBlock, transitBlock->intSon11);
        setFullMergeTarget(targetBlock, transitBlock->intSon12);
        setFullMergeTarget(targetBlock, transitBlock->intSon21);
        setFullMergeTarget(targetBlock, transitBlock->intSon22);
    }
}

void HMultiply::findIndependentWorkBlocks(intermBlCl* block, QVector<intermBlCl*> &indepWorkBlocks)
{
    if(block->holdsFullInformation || block->holdsRKInformation || block->isEventualLeaf)
    {
        indepWorkBlocks.append(block);
    }
    else
    {
        if(block->intSon11 != nullptr)
        {
            findIndependentWorkBlocks(block->intSon11, indepWorkBlocks);
            findIndependentWorkBlocks(block->intSon12, indepWorkBlocks);
            findIndependentWorkBlocks(block->intSon21, indepWorkBlocks);
            findIndependentWorkBlocks(block->intSon22, indepWorkBlocks);
        }
        else
        {
            std::cerr << "block->son11 == nullptr && !block->isEventualLeaf in findIndependentWorkBlocks() call." << std::endl;
        }
    }
}

void HMultiply::recursProdPartition(intermBlCl* productBlock, BlockCluster* factorBlock1, BlockCluster* factorBlock2)
{
    if(factorBlock1 -> isLeaf || factorBlock2 -> isLeaf)
    {
        productBlock -> isEventualLeaf = true;
        return;
    }
    else if(productBlock -> intSon11 == nullptr) // if a block is partitioned, it is cross partitioned -> only one son check neccessary
    {
        productBlock -> isEventualLeaf = true;
        return;
    }
    else
    {
        recursProdPartition(productBlock->intSon11, factorBlock1->son11, factorBlock2->son11);
        recursProdPartition(productBlock->intSon12, factorBlock1->son12, factorBlock2->son12);
        recursProdPartition(productBlock->intSon21, factorBlock1->son21, factorBlock2->son21);
        recursProdPartition(productBlock->intSon22, factorBlock1->son22, factorBlock2->son22);
    }
}

void HMultiply::startWorkOnIndependentBlocks(QVector<intermBlCl*> &blocksWithMatrixLoad, const long rank, const double relError)
{
    // can parallelize over blocks, because the blocks (and their submatrices) are disjunct
    #pragma omp parallel for
    for(long i = 0; i < blocksWithMatrixLoad.length(); i++)
    {
        processBlocksDownWardRecursion(blocksWithMatrixLoad[i], rank, relError);
    }
}

void HMultiply::processBlocksDownWardRecursion(intermBlCl* productBlock, const long rank, const double relError)
{
    if(productBlock->holdsFullInformation || productBlock->holdsRKInformation)
    {
        processFactorPairs(productBlock, rank, relError);
    }
    if(productBlock->isEventualLeaf) // block is in the partition; block is a leaf node in the final product tree
    {
        // start merge operation
        if(productBlock->intSon11 != nullptr) // the (eventual) leaf node/block still has a nontrivial subtree
        {
            mergeSubtree(productBlock->intSon11, rank, relError); // start an agglomeration for the subtree of the sons  of productBlock
            mergeSubtree(productBlock->intSon12, rank, relError);
            mergeSubtree(productBlock->intSon21, rank, relError);
            mergeSubtree(productBlock->intSon22, rank, relError);

            int numberOfSonsWithRkMatrices = productBlock->intSon11->holdsRKInformation + productBlock->intSon12->holdsRKInformation + productBlock->intSon21->holdsRKInformation + productBlock->intSon22->holdsRKInformation;
            if(numberOfSonsWithRkMatrices >= 1)
            {
                roundedAgglomerateSonsIntoBlock(productBlock, rank, relError); // agglomerates above resuts into the eventual leaf (partition) block
            }
            productBlock->trimBelow();
            if(productBlock->intSon11 != nullptr)
            {
                std::cerr << "trimming hasn't worked." << std::endl;
            }
        }

        if(!productBlock->isAdmissible && productBlock->singularValues.size() > 0)
        {
            HArithm::addrkMatToFull(productBlock->fullMat, productBlock->UMat, productBlock->singularValues, productBlock->VAdjMat);
            clearRkMat(productBlock);
        }

        if(productBlock->isAdmissible && productBlock->fullMat.size() > 0)
        {
            long rows = productBlock->rows();
            long cols = productBlock->cols();
            if(rows < cols)
            {
                HMultiply::roundedAddRMatToBlock(productBlock, Eigen::MatrixXcd::Identity(rows, rows), Eigen::VectorXd::Ones(rows), productBlock->fullMat, rank, relError);
            }
            else
            {
                HMultiply::roundedAddRMatToBlock(productBlock, productBlock->fullMat, Eigen::VectorXd::Ones(cols), Eigen::MatrixXcd::Identity(cols, cols), rank, relError);
            }

            productBlock->fullMat.resize(0,0);
            productBlock->holdsFullInformation = false;
        }
    }
    else  //flush information downwards
    {
        if(productBlock->holdsRKInformation)
        {
            for(int i = 0; i < productBlock->rKFlushTargets.length(); i++)
            {
                long rows = productBlock->rKFlushTargets[i]->rows();
                long cols = productBlock->rKFlushTargets[i]->cols();
                long rowStartIndex = productBlock->rKFlushTargets[i]->rowStartIndex() - productBlock->rowStartIndex();
                long colStartIndex = productBlock->rKFlushTargets[i]->colStartIndex() - productBlock->colStartIndex();

                HMultiply::roundedAddRMatToBlock(productBlock->rKFlushTargets[i], productBlock->UMat.middleRows(rowStartIndex, rows), productBlock->singularValues, productBlock->VAdjMat.middleCols(colStartIndex, cols), rank, relError);
                productBlock->rKFlushTargets[i]->holdsRKInformation = true;
            }
        }

        if(productBlock->holdsFullInformation)
        {
            if(productBlock->fullMat.size() != 0)
            {
                for(int i = 0; i < productBlock->fullFlushMergeTargets.length(); i++)
                {
                    long rows = productBlock->fullFlushMergeTargets[i]->rows();
                    long cols = productBlock->fullFlushMergeTargets[i]->cols();
                    long rowStartIndex = productBlock->fullFlushMergeTargets[i]->rowStartIndex() - productBlock->rowStartIndex();
                    long colStartIndex = productBlock->fullFlushMergeTargets[i]->colStartIndex() - productBlock->colStartIndex();

                    HArithm::addFullMatrixInToFirstOne(productBlock->fullFlushMergeTargets[i]->fullMat, productBlock->fullMat.block(rowStartIndex, colStartIndex, rows, cols));
                    productBlock->fullFlushMergeTargets[i]->holdsFullInformation = true;
                }
            }
        }
        clearBlock(productBlock);

        processBlocksDownWardRecursion(productBlock->intSon11, rank, relError);
        processBlocksDownWardRecursion(productBlock->intSon12, rank, relError);
        processBlocksDownWardRecursion(productBlock->intSon21, rank, relError);
        processBlocksDownWardRecursion(productBlock->intSon22, rank, relError);
    }
}

void HMultiply::mergeSubtree(intermBlCl* productBlock/*, intermBlCl *nextHighestRkBlock*/, const long rank, const double relError)
{
    if(productBlock->intSon11 != nullptr) // block has non-trivial subtree -> process sons first, so that their results can be agglomerated
    {
        mergeSubtree(productBlock->intSon11, rank, relError);
        mergeSubtree(productBlock->intSon12, rank, relError);
        mergeSubtree(productBlock->intSon21, rank, relError);
        mergeSubtree(productBlock->intSon22, rank, relError);
    }
    if(productBlock->holdsFullInformation || productBlock->holdsRKInformation)
    {
        processFactorPairs(productBlock, rank, relError); // calculate the product full matrix and rk-matrix
        if(productBlock->holdsFullInformation)
        {
            intermBlCl *targetForFull = productBlock->fullFlushMergeTargets[0];
            if(targetForFull->isAdmissible) // full matrix will need to be rk-matrix eventually -> convert it already here; will be aglomerated later
            {
                Eigen::MatrixXcd U;
                Eigen::VectorXd singVals;
                Eigen::MatrixXcd VAdj;

                if(productBlock->rows() < productBlock->cols())
                {
                    U = Eigen::MatrixXcd::Identity(productBlock->rows(), productBlock->rows());
                    singVals = Eigen::VectorXd::Ones(productBlock->rows());
                    VAdj = productBlock->fullMat;
                }
                else
                {
                    U = productBlock->fullMat;
                    singVals = Eigen::VectorXd::Ones(productBlock->cols());
                    VAdj = Eigen::MatrixXcd::Identity(productBlock->cols(), productBlock->cols());
                }

                HMultiply::roundedAddRMatToBlock(productBlock, U, singVals, VAdj, rank, relError);
                productBlock->fullMat.resize(0,0);
                productBlock->holdsFullInformation = false;
                productBlock->holdsRKInformation = true; // the block has an rk-matrix (now)
            }
            else // target block is inadmissible -> accepts the full matrix directly
            {
                if(targetForFull->fullMat.size() == 0) // initialize target full matrix, if it hasn't been already
                {
                    targetForFull->fullMat = Eigen::MatrixXcd::Zero(targetForFull->rows(), targetForFull->cols());
                }

                long rows = productBlock->rows();
                long cols = productBlock->cols();
                long rowStartIndex = targetForFull->rowStartIndex() - productBlock->rowStartIndex();
                long colStartIndex = targetForFull->colStartIndex() - productBlock->colStartIndex();

                targetForFull->fullMat.block(rowStartIndex, colStartIndex, rows, cols) += productBlock->fullMat; // add full matrix to target fullmatrix
                productBlock->fullMat.resize(0,0);
                productBlock->holdsFullInformation = false;
            }
        }
    }
    if(productBlock->intSon11 != nullptr) // block has non-trivial subtree -> agglomerated sons
    {
        int numberOfSonsWithRkMatrices = productBlock->intSon11->holdsRKInformation + productBlock->intSon12->holdsRKInformation + productBlock->intSon21->holdsRKInformation + productBlock->intSon22->holdsRKInformation;
        if(numberOfSonsWithRkMatrices >= 1)
        {
            roundedAgglomerateSonsIntoBlock(productBlock, rank, relError);
        }
    }
}

void HMultiply::processFactorPairs(intermBlCl* productBlock, const long rank, const double relError)
{
    if(productBlock->holdsFullInformation && productBlock->fullMat.size() == 0 && !productBlock->isAdmissible)
    {
        productBlock->fullMat = Eigen::MatrixXcd::Zero(productBlock->rows(), productBlock->cols());
    }

    for(int i = 0; i < productBlock->fullFactorsList.length(); i++) // process factor pairs that result in a full matrix
    {
        factorsToFull(productBlock, productBlock->fullFactorsList[i], rank, relError);
    }    
    for(int i = 0; i < productBlock->rKFactorsList.length(); i++) // process factor pairs that result in a reduced rank matrix
    {
        factorsToRK(productBlock, productBlock->rKFactorsList[i], rank, relError);
    }
    if(productBlock->holdsRKInformation && productBlock->singularValues.size() == 0)
    {
        std::cerr << "productBlock->holdsRKInformation && productBlock->singularValues.size() == 0" << std::endl;
    }
    else if(!productBlock->holdsRKInformation && productBlock->rKFactorsList.length() > 0)
    {
        std::cerr << "!productBlock->holdsRKInformation && productBlock->rKFactors.length() > 0 " << productBlock->rKFactorsList.length() << std::endl;
    }
}

void HMultiply::factorsToFull(intermBlCl* productBlock, const std::pair<BlockCluster*, BlockCluster*> fullFactors, const long rank, const double relError)
{
    if(!productBlock->isAdmissible)
    {
        if(fullFactors.first->isLeaf)  // first block is in the partition -> block holds martix information and is not further subdivided
        {
            HArithm::matByBlock(fullFactors.first->fullMat, *fullFactors.second, fullFactors.first->colStartIndex(), fullFactors.second->colStartIndex(), productBlock->fullMat);
        }
        else // second factor is leaf
        {
            HArithm::blockByMat(*fullFactors.first, fullFactors.second->fullMat, fullFactors.first->rowStartIndex(), fullFactors.second->rowStartIndex(), productBlock->fullMat);
        }
    }
    else
    {
        productBlock->holdsRKInformation = true;
        if(fullFactors.first->isLeaf)  // first block is in the partition -> block holds martix information and is not further subdivided
        {
            Eigen::MatrixXcd productMat = Eigen::MatrixXcd::Zero( fullFactors.first->cols(), fullFactors.second->cols() );
            Eigen::MatrixXcd UMat = fullFactors.first->fullMat;
            Eigen::MatrixXcd singularValues =  Eigen::VectorXcd::Ones(fullFactors.first->cols());
            Eigen::MatrixXcd VAdjMat =  Eigen::MatrixXcd::Identity(fullFactors.first->cols(), fullFactors.first->cols());

            HArithm::matByBlock(VAdjMat, *fullFactors.second, fullFactors.first->colStartIndex(), fullFactors.second->colStartIndex(), productMat);
            HMultiply::roundedAddRMatToBlock(productBlock, UMat, singularValues, productMat, rank, relError);
        }
        else // second factor is leaf
        {
            Eigen::MatrixXcd productMat = Eigen::MatrixXcd::Zero(fullFactors.first->rows(), fullFactors.second->rows());
            Eigen::MatrixXcd UMat = Eigen::MatrixXcd::Identity(fullFactors.second->rows(), fullFactors.second->rows());
            Eigen::MatrixXcd singularValues =  Eigen::VectorXcd::Ones(fullFactors.second->rows());
            Eigen::MatrixXcd VAdjMat =  fullFactors.second->fullMat;

            HArithm::blockByMat(*fullFactors.first, UMat, fullFactors.first->rowStartIndex(), fullFactors.second->rowStartIndex(), productMat);
            HMultiply::roundedAddRMatToBlock(productBlock, productMat, singularValues, VAdjMat, rank, relError);
        }
    }
}

void HMultiply::factorsToRK(intermBlCl* productBlock, const std::pair<BlockCluster*, BlockCluster*> rKFactors, const long rank, const double relError)
{
    // the factor pair results in a reduced rank matrix
    if(rKFactors.first->isLeaf && rKFactors.first->isAdmissible)
    {
        Eigen::MatrixXcd productMat = Eigen::MatrixXcd::Zero( rKFactors.first->VAdjMat.rows(), rKFactors.second->cols() );
        HArithm::matByBlock(rKFactors.first->VAdjMat, *rKFactors.second, rKFactors.first->colStartIndex(), rKFactors.second->colStartIndex(), productMat);
        HMultiply::roundedAddRMatToBlock(productBlock, rKFactors.first->UMat, rKFactors.first->singularValues, productMat, rank, relError);
    }
    else // (rKFactors.second->isLeaf && rKFactors.second->isAdmissible)
    {
        Eigen::MatrixXcd productMat = Eigen::MatrixXcd::Zero(rKFactors.first->rows(), rKFactors.second->UMat.cols());
        HArithm::blockByMat(*rKFactors.first, rKFactors.second ->UMat, rKFactors.first->rowStartIndex(), rKFactors.second->rowStartIndex(), productMat);
        HMultiply::roundedAddRMatToBlock(productBlock, productMat, rKFactors.second->singularValues, rKFactors.second->VAdjMat, rank, relError);
    }
}

void HMultiply::addRMatToBlock(intermBlCl* block, const Eigen::MatrixXcd &UMat, const Eigen::VectorXcd &singVals, const Eigen::MatrixXcd &VAdjMat)
{
    HArithm::horizontalJoinMatricesInToFirstOne(block->UMat, UMat);
    HArithm::joinVectorsInToFirstOne(block->singularValues, singVals);
    HArithm::verticalJoinMatricesInToFirstOne(block->VAdjMat, VAdjMat);
}

void HMultiply::roundedAddRMatToBlock(intermBlCl* block, const Eigen::MatrixXcd &UMat, const Eigen::VectorXcd &singVals, const Eigen::MatrixXcd &VAdjMat, const long maxRank, const double relError)
{
    addRMatToBlock(block, UMat, singVals, VAdjMat);
    RkMatRankReduction(block->UMat, block->singularValues, block->VAdjMat, maxRank, relError);
}

void HMultiply::clearRkMat(intermBlCl* block)
{
    block->UMat.resize(0,0);
    block->singularValues.resize(0);
    block->VAdjMat.resize(0,0);
    block->holdsRKInformation = false;
}

void HMultiply::clearBlock(intermBlCl* block)
{
    block->fullMat.resize(0,0);
    block->holdsFullInformation = false;

    clearRkMat(block);
}

void HMultiply::RkMatRankReduction(Eigen::MatrixXcd &UMat, Eigen::VectorXcd &singularValues, Eigen::MatrixXcd &VAdjMat, long rank, const double relError)
{
    long localRank = singularValues.size();
    long minDim = std::min(UMat.rows(), VAdjMat.cols());
    if(((localRank <= rank && relError == 0) || (rank <= 0 && relError <= 0)) && localRank <= minDim) //might be wrong with rk mat in svd form
    {
//        std::cout << "Unneccessary RMatRankReduction() call." << std::endl;
        return;
    }
    else
    {
//        Eigen::MatrixXcd full = UMat * singularValues.asDiagonal() * VAdjMat;

        long localARank = std::min(localRank, UMat.rows());
        long localBRank = std::min(localRank, VAdjMat.cols());

//        long localARank = localRank;
//        long localBRank = localRank;

        Eigen::HouseholderQR<Eigen::MatrixXcd> qrA(UMat);
        Eigen::HouseholderQR<Eigen::MatrixXcd> qrB_T(VAdjMat.transpose());

        Eigen::MatrixXcd R_A = qrA.matrixQR().triangularView<Eigen::Upper>().toDenseMatrix().topRows(localARank);
        Eigen::MatrixXcd R_B_Transpose = qrB_T.matrixQR().triangularView<Eigen::Upper>().transpose().toDenseMatrix().leftCols(localBRank);

        Eigen::MatrixXcd A = (R_A * singularValues.head(R_A.cols()).asDiagonal()) * R_B_Transpose.triangularView<Eigen::Lower>();

        if constexpr(useEigenSvd)
        {
//            Eigen::JacobiSVD<Eigen::MatrixXcd,Eigen::ComputeThinU|Eigen::ComputeThinV> svd(A);
            Eigen::BDCSVD<Eigen::MatrixXcd,Eigen::ComputeThinU|Eigen::ComputeThinV> svd(A);
            rank = HArithm::minRankforError(svd, rank, relError);
            Eigen::MatrixXcd tmpForU(UMat.rows(), rank);
            tmpForU << svd.matrixU().leftCols(rank), Eigen::MatrixXcd::Zero(UMat.rows() - localARank, rank);
            UMat.noalias() = qrA.householderQ() * tmpForU;
            singularValues = svd.singularValues().head(rank);
            Eigen::MatrixXcd tmpForVAdj(rank, VAdjMat.cols());
            tmpForVAdj << svd.matrixV().adjoint().topRows(rank), Eigen::MatrixXcd::Zero(rank, VAdjMat.cols() - localBRank);
            VAdjMat.noalias() = tmpForVAdj * qrB_T.householderQ().transpose();
        }
        else
        {
            Eigen::MatrixXcd tmpU;
            Eigen::VectorXd tmpSingularValues;
            Eigen::MatrixXcd tmpVAdj;

            GoloubReinschSvd::goloubReinschSVD(A, tmpU, tmpSingularValues , tmpVAdj);

            rank = GoloubReinschSvd::minRankforError(tmpSingularValues, rank, relError);

            Eigen::MatrixXcd tmpForU(UMat.rows(), rank);
            tmpForU << tmpU.leftCols(rank), Eigen::MatrixXcd::Zero(UMat.rows() - localARank, rank);
            UMat.noalias() = qrA.householderQ() * tmpForU;

            singularValues = tmpSingularValues.head(rank);

            Eigen::MatrixXcd tmpForVAdj(rank, VAdjMat.cols());
            tmpForVAdj << tmpVAdj.topRows(rank), Eigen::MatrixXcd::Zero(rank, VAdjMat.cols() - localBRank);
            VAdjMat.noalias() = tmpForVAdj * qrB_T.householderQ().transpose();
        }
//        double relErrorReal = (full- UMat * singularValues.asDiagonal() * VAdjMat).norm() / full.norm();
//        if(relErrorReal > relError)
//        {
//            std::cerr << "svd rel error: " << (full- UMat * singularValues.asDiagonal() * VAdjMat).norm() / full.norm() << std::endl;
//        }
    }
}

void HMultiply::roundedAgglomerateSonsIntoBlock(intermBlCl* block, const long rank, const double relError)
{
    std::tuple<Eigen::MatrixXcd, Eigen::VectorXcd, Eigen::MatrixXcd> agglomeratedRkMatrix = agglomerateSons(block, rank, relError);
    //    HArithm::roundedAddRMatToBlock(*block, std::get<0>(agglomeratedRkMatrix), std::get<1>(agglomeratedRkMatrix), std::get<2>(agglomeratedRkMatrix), rank, relError);
    HMultiply::roundedAddRMatToBlock(block, std::get<0>(agglomeratedRkMatrix), std::get<1>(agglomeratedRkMatrix), std::get<2>(agglomeratedRkMatrix), rank, relError);

    clearRkMat(block->intSon11);
    clearRkMat(block->intSon12);
    clearRkMat(block->intSon21);
    clearRkMat(block->intSon22);

    block->holdsRKInformation = true;
}

std::tuple<Eigen::MatrixXcd, Eigen::VectorXcd, Eigen::MatrixXcd> HMultiply::agglomerateSons(intermBlCl* block, const long rank, const double relError) // only call function on block, whose sons hold actual rk information
{
    if(!block->intSon11->holdsRKInformation && block->intSon11->singularValues.size() >= 1)
    {
        std::cerr << "!block->son11->holdsRKInformation && block->son11->singularValues.size() = " << block->intSon11->singularValues.size() << std::endl;
    }
    if(!block->intSon12->holdsRKInformation && block->intSon12->singularValues.size() >= 1)
    {
        std::cerr << "!block->son12->holdsRKInformation && block->son12->singularValues.size() = " << block->intSon12->singularValues.size() << std::endl;
    }
    if(!block->intSon21->holdsRKInformation && block->intSon21->singularValues.size() >= 1)
    {
        std::cerr << "!block->son21->holdsRKInformation && block->son21->singularValues.size() = " << block->intSon21->singularValues.size() << std::endl;
    }
    if(!block->intSon22->holdsRKInformation && block->intSon22->singularValues.size() >= 1)
    {
        std::cerr << "!block->son22->holdsRKInformation && block->son22->singularValues.size() = " << block->intSon22->singularValues.size() << std::endl;
    }

    Eigen::MatrixXcd U;
    Eigen::VectorXcd singularValues;
    Eigen::MatrixXcd VAdj;

    if(block->intSon11 != nullptr)
    {
        long son11Rows = block->intSon11->rows();
        long son11Cols = block->intSon11->cols();

        long son22Rows = block->intSon22->rows();
        long son22Cols = block->intSon22->cols();

        Eigen::MatrixXcd URow1;
        Eigen::VectorXcd singularValuesRow1;
        Eigen::MatrixXcd VAdjRow1;

        if(block->intSon11->holdsRKInformation || block->intSon12->holdsRKInformation) // at least one of the upper sons holds information
        {
            long son12Cols = block->intSon12->cols();

            long son11Rank = block->intSon11->singularValues.size();
            long son12Rank = block->intSon12->singularValues.size();

            URow1 = Eigen::MatrixXcd/*::Zero*/(son11Rows, son11Rank + son12Rank);
            singularValuesRow1 = Eigen::VectorXcd/*::Zero*/(son11Rank + son12Rank);
            VAdjRow1 = Eigen::MatrixXcd/*::Zero*/(son11Rank + son12Rank, son11Cols + son12Cols);

            if(block->intSon11->holdsRKInformation)
            {
                URow1.leftCols(son11Rank) = block->intSon11->UMat;
                singularValuesRow1.head(son11Rank) = block->intSon11->singularValues;
                VAdjRow1.topLeftCorner(son11Rank, son11Cols) = block->intSon11->VAdjMat;
                VAdjRow1.topRightCorner(son11Rank, son12Cols) = Eigen::MatrixXcd::Zero(son11Rank, son12Cols);
            }

            if(block->intSon12->holdsRKInformation)
            {
                URow1.rightCols(son12Rank) = block->intSon12->UMat;
                singularValuesRow1.tail(son12Rank) = block->intSon12->singularValues;
                VAdjRow1.bottomLeftCorner(son12Rank, son11Cols) = Eigen::MatrixXcd::Zero(son12Rank, son11Cols);
                VAdjRow1.bottomRightCorner(son12Rank, son12Cols) = block->intSon12->VAdjMat;
            }

            if(block->intSon11->holdsRKInformation && block->intSon12->holdsRKInformation)
            {
                RkMatRankReduction(URow1, singularValuesRow1, VAdjRow1, rank, relError);
            }
        }

        Eigen::MatrixXcd URow2;
        Eigen::VectorXcd singularValuesRow2;
        Eigen::MatrixXcd VAdjRow2;

        if(block->intSon21->holdsRKInformation || block->intSon22->holdsRKInformation) // at least one of the lowir sons holds information
        {
            long son21Cols = block->intSon21->cols();

            long son21Rank = block->intSon21->singularValues.size();
            long son22Rank = block->intSon22->singularValues.size();

            URow2 = Eigen::MatrixXcd/*::Zero*/(son22Rows, son21Rank + son22Rank);
            singularValuesRow2 = Eigen::VectorXcd/*::Zero*/(son21Rank + son22Rank);
            VAdjRow2 = Eigen::MatrixXcd/*::Zero*/(son21Rank + son22Rank, son21Cols + son22Cols);

            if(block->intSon21->holdsRKInformation)
            {
                URow2.leftCols(son21Rank) = block->intSon21->UMat;
                singularValuesRow2.head(son21Rank) = block->intSon21->singularValues;
                VAdjRow2.topLeftCorner(son21Rank, son21Cols) = block->intSon21->VAdjMat;
                VAdjRow2.topRightCorner(son21Rank, son22Cols) = Eigen::MatrixXcd::Zero(son21Rank, son22Cols);
            }
            if(block->intSon22->holdsRKInformation)
            {
                URow2.rightCols(son22Rank) = block->intSon22->UMat;
                singularValuesRow2.tail(son22Rank) = block->intSon22->singularValues;
                VAdjRow2.bottomLeftCorner(son22Rank, son21Cols) = Eigen::MatrixXcd::Zero(son22Rank, son21Cols);
                VAdjRow2.bottomRightCorner(son22Rank, son22Cols) = block->intSon22->VAdjMat;
            }

            if(block->intSon21->holdsRKInformation && block->intSon22->holdsRKInformation)
            {
                RkMatRankReduction(URow2, singularValuesRow2, VAdjRow2, rank, relError);
            }
        }

        if( (block->intSon11->holdsRKInformation || block->intSon12->holdsRKInformation) || (block->intSon21->holdsRKInformation || block->intSon22->holdsRKInformation) ) // merge the upper and the lower sons
        {
            U = Eigen::MatrixXcd/*::Zero*/(block->rows(), singularValuesRow1.size() + singularValuesRow2.size());
            singularValues = Eigen::VectorXcd/*::Zero*/(singularValuesRow1.size() + singularValuesRow2.size());
            VAdj = Eigen::MatrixXcd/*::Zero*/(singularValuesRow1.size() + singularValuesRow2.size(), block->cols());

            if((block->intSon11->holdsRKInformation || block->intSon12->holdsRKInformation))
            {
                U.topLeftCorner(son11Rows, singularValuesRow1.size()) = URow1;
                U.bottomLeftCorner(son22Rows, singularValuesRow1.size()) = Eigen::MatrixXcd::Zero(son22Rows, singularValuesRow1.size());
                singularValues.head(singularValuesRow1.size()) = singularValuesRow1;
                VAdj.topRows(singularValuesRow1.size()) = VAdjRow1;
            }
            if(block->intSon21->holdsRKInformation || block->intSon22->holdsRKInformation)
            {
                U.topRightCorner(son11Rows, singularValuesRow2.size()) = Eigen::MatrixXcd::Zero(son11Rows, singularValuesRow2.size());
                U.bottomRightCorner(son22Rows, singularValuesRow2.size()) = URow2;
                singularValues.tail(singularValuesRow2.size()) = singularValuesRow2;
                VAdj.bottomRows(singularValuesRow2.size()) = VAdjRow2;
            }

            if( (block->intSon11->holdsRKInformation || block->intSon12->holdsRKInformation) && (block->intSon21->holdsRKInformation || block->intSon22->holdsRKInformation) ) // upper and lower sons hold information  -> recompress the sum
            {
                RkMatRankReduction(U, singularValues, VAdj, rank, relError);
            }
            return {U, singularValues, VAdj};
        }
    }
    return {Eigen::MatrixXcd::Zero(block->rows(), 1), Eigen::VectorXcd::Zero(1), Eigen::MatrixXcd::Zero(1, block->cols())};
}

//void HMultiply::intermBlClToBlockCluster(BlockCluster* blockClusterNode, intermBlCl* block)
//{
//    if(block->isEventualLeaf)
//    {
//        if(block->isAdmissible)
//        {
//            blockClusterNode->UMat = block->UMat;
//            blockClusterNode->singularValues = block->singularValues;
//            blockClusterNode->VAdjMat = block->VAdjMat;
//        }
//        else
//        {
//            blockClusterNode->fullMat = block->fullMat;
//        }
//        blockClusterNode->isLeaf = true;
//    }
//    else
//    {
//        blockClusterNode->isLeaf = false;
//    }

//    blockClusterNode->isAdmissible = block->isAdmissible;

//    blockClusterNode->rowCluster = block->rowCluster;
//    blockClusterNode->columnCluster = block->columnCluster;

//    if(/*block->son11 != nullptr && */!block->isEventualLeaf)
//    {
//        blockClusterNode->son11 = new BlockCluster(blockClusterNode->rowCluster->son1, blockClusterNode->columnCluster->son1, blockClusterNode);
//        intermBlClToBlockCluster(blockClusterNode->son11, block->intSon11);

//        blockClusterNode->son12 = new BlockCluster(blockClusterNode->rowCluster->son1, blockClusterNode->columnCluster->son2, blockClusterNode);
//        intermBlClToBlockCluster(blockClusterNode->son12, block->intSon12);

//        blockClusterNode->son21 = new BlockCluster(blockClusterNode->rowCluster->son2, blockClusterNode->columnCluster->son1, blockClusterNode);
//        intermBlClToBlockCluster(blockClusterNode->son21, block->intSon21);

//        blockClusterNode->son22 = new BlockCluster(blockClusterNode->rowCluster->son2, blockClusterNode->columnCluster->son2, blockClusterNode);
//        intermBlClToBlockCluster(blockClusterNode->son22, block->intSon22);
//    }
//}

BlockCluster* HMultiply::intermBlClToBlockCluster(intermBlCl* block) /*!< Convert the intermediary product blockcluster into a ordinary block cluster. */
{
    if(block->intSon11 != nullptr && !block->isEventualLeaf)
    {
        #pragma omp parallel sections
        {
            #pragma omp section
            block->son11 = intermBlClToBlockCluster(block->intSon11);
            #pragma omp section
            block->son12 = intermBlClToBlockCluster(block->intSon12);
            #pragma omp section
            block->son21 = intermBlClToBlockCluster(block->intSon21);
            #pragma omp section
            block->son22 = intermBlClToBlockCluster(block->intSon22);
        }
    }
    if(!block->isEventualLeaf)
    {
        block->isAdmissible = false;
        block->isLeaf = false;
    }
    else
    {
        block->isLeaf = true;
    }
    block->fullFactorsList.clear();
    block->rKFactorsList.clear();
    block->fullFlushMergeTargets.clear();
    block->rKFlushTargets.clear();

    block->intSon11 = nullptr;
    block->intSon12 = nullptr;
    block->intSon21 = nullptr;
    block->intSon22 = nullptr;

    BlockCluster* blockPointer = static_cast<BlockCluster*> (block);
    if(block->intSon11 != nullptr && !block->isEventualLeaf)
    {
        block->son11->father = blockPointer;
        block->son12->father = blockPointer;
        block->son21->father = blockPointer;
        block->son22->father = blockPointer;
    }
    return blockPointer;
}
