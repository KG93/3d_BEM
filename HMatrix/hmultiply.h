#ifndef HMULTIPLY_H
#define HMULTIPLY_H
#include "harithm.h"
//#include <tuple>

/**
 * \class intermBlCl
 *
 * \brief Intermediate blockcluster class for the H-matrix-matrix multiplication.
 */
class intermBlCl : public BlockCluster
{
public:
    intermBlCl(Cluster* rowCluster, Cluster* columnCluster) /*!< \brief Construct blockcluster. */
    {
        this->rowCluster = rowCluster;
        this->columnCluster = columnCluster;
        if(std::min(rows(), cols()) < ClusterTree::triangleThreshold)
        {
            this->isAdmissible = false;
        }
        else
        {
            this->isAdmissible = true;
        }
    }
    intermBlCl(Cluster* rowCluster, Cluster* columnCluster, bool isRoot) /*!< \brief Construct root blockcluster. */
    {
        this->rowCluster = rowCluster;
        this->columnCluster = columnCluster;
        this->isRoot = isRoot;
        if(std::min(rows(), cols()) < ClusterTree::triangleThreshold)
        {
            this->isAdmissible = false;
        }
        else
        {
            this->isAdmissible = true;
        }
    }
    intermBlCl(Cluster* rowCluster, Cluster* columnCluster, intermBlCl &father) /*!< \brief Construct blockcluster. */
    {
        this->rowCluster = rowCluster;
        this->columnCluster = columnCluster;
        this->intFather = &father;
        if(std::min(rows(), cols()) < ClusterTree::triangleThreshold)
        {
            this->isAdmissible = false;
        }
        else
        {
            this->isAdmissible = true;
        }
    }
    virtual ~intermBlCl(){}

//    long rows() const; /*!< Get number of rows of the block. */
//    long cols() const; /*!< Get number of colums of the block. */
//    long size() const; /*!< Get the size of the block. */
//    long rowStartIndex() const; /*!< Get the global index of the first block row. */
//    long colStartIndex() const; /*!< Get the global index of the first block column. */
////    double compression() const;
//    void clearMatrixInfo(); /*!< Clear the actual matrices in the block. */
    void clear(); /*!< \brief Clear the entire subtree (the block and all subblocks). */
    void trimBelow(); /*!< \brief Trim the subtree of the block. */

//    Cluster* rowCluster = nullptr;
//    Cluster* columnCluster = nullptr;
    intermBlCl* intFather = nullptr;
//    intermBlCl* intSon11 = nullptr;
//    intermBlCl* intSon12 = nullptr;
//    intermBlCl* intSon21 = nullptr;
//    intermBlCl* intSon22 = nullptr;
    std::unique_ptr<intermBlCl> intSon11 = nullptr;
    std::unique_ptr<intermBlCl> intSon12 = nullptr;
    std::unique_ptr<intermBlCl> intSon21 = nullptr;
    std::unique_ptr<intermBlCl> intSon22 = nullptr;

//    bool isAdmissible = true;
//    bool isRoot = false;
    bool holdsRKInformation = false;
    bool holdsFullInformation = false;
    bool isEventualLeaf = false;

//    //    bool alreadyinMMr = false;
//    Eigen::MatrixXcd fullMat;
//    Eigen::MatrixXcd UMat;
//    Eigen::VectorXcd singularValues;
//    Eigen::MatrixXcd VAdjMat;
//    double frobeniusNorm = -1;
    QVector<std::pair<BlockCluster*, BlockCluster*>> fullFactorsList;
    QVector<std::pair<BlockCluster*, BlockCluster*>> rKFactorsList;
    QVector<intermBlCl*> fullFlushMergeTargets; //if the intermittent block is not in product partition, then the block's information will be moved to actual partition blocks (the targets)
    QVector<intermBlCl*> rKFlushTargets; //if the intermittent block is not in product partition, then the block's information will be moved to actual partition blocks (the targets)

    //EIGEN_MAKE_ALIGNED_OPERATOR_NEW //only needed for fixed size eigen objects
    //    int localRank = 0;
    void setRKFactors(const QVector<std::pair<BlockCluster *, BlockCluster *> > &newRKFactors); /*!< \brief Set the list (sum) of low rank multiplications. */
};

/**
 * \class HMultiply
 *
 * \brief Class for the H-matrix-matrix multiplication.
 */
class HMultiply
{
public:
    HMultiply();
    static HMatrix multiplyHMat(BlockCluster &factor1, BlockCluster &factor2, const long rank, const double relError); /*!< \brief The method returns the product of two H-matrices. */
    static HMatrix multiplyHMat(HMatrix &factor1, HMatrix &factor2, const long rank, const double relError); /*!< \brief The method returns the product of two H-matrices. */

private:
    static void MM(BlockCluster &factorBlock1, BlockCluster &factorBlock2, intermBlCl &productBlock); /*!< \brief The method sets up the product tree and sets up the partial sums in the tree nodes. The actual work is handled by later tasks. */
    static void findFullFlushTargets(intermBlCl &startBlock, QVector<intermBlCl*> transitBlocks); /*!< \brief The method flushes full matrices further down the product tree. */
    static void findRkFlushTargets(intermBlCl &potentialTargetBlock, intermBlCl *transitBlock = nullptr); /*!< \brief The method flushes low rank matrices further down the product tree. */
    static void setFullMergeTarget(intermBlCl &targetBlock, intermBlCl &transitBlock); /*!< \brief The method merges sub blocks (the sub-tree) up into a full matrix in the target block. */
    static void findIndependentWorkBlocks(intermBlCl &block, QVector<intermBlCl*> &indepWorkBlocks); /*!< \brief The method traverses the product tree downward and assambles a list of all independent nodes (blocks) with work. */
    static void recursProdPartition(intermBlCl &productBlock, BlockCluster &factorBlock1, BlockCluster &factorBlock2); /*!< \brief The method determines the eventual leaf nodes in the product tree (and thereby determines the product matrix partition). */
    static void startWorkOnIndependentBlocks(QVector<intermBlCl*> &blocksWithMatrixLoad, const long rank, const double relError); /*!< \brief Start the work on the the independent work nodes. */
    static void processBlocksDownWardRecursion(intermBlCl &productBlock, const long rank, const double relError); /*!< \brief Process the block work and flush the result down or merge up (depending on the leaf nodes.) */
    static void mergeSubtree(intermBlCl &productBlock, const long rank, const double relError); /*!< \brief The method merges the sub tree into the productBlock. */
    static void processFactorPairs(intermBlCl &productBlock, const long rank, const double relError); /*!< \brief Do the multiplication work in the product block. */
    static void factorsToFull(intermBlCl &productBlock, const std::pair<BlockCluster*, BlockCluster*> fullFactors, const long rank, const double relError); /*!< \brief Calculate the sum of factors that result in a full matrix in the block. */
    static void factorsToRK(intermBlCl &productBlock, const std::pair<BlockCluster*, BlockCluster*> rKFactors, const long rank, const double relError); /*!< \brief Calculate the sum of factors that result in a low rank matrix in the block. */
    static void addRMatToBlock(intermBlCl &block, const Eigen::MatrixXcd &UMat, const Eigen::VectorXcd &singVals, const Eigen::MatrixXcd &VAdjMat); /*!< \brief Add a low rank matrix to the block information. */
    static void roundedAddRMatToBlock(intermBlCl &block, const Eigen::MatrixXcd &UMat, const Eigen::VectorXcd &singVals, const Eigen::MatrixXcd &VAdjMat, const long maxRank, const double relError); /*!< \brief Truncated addition of a low rank matrix to the block. */
    static void clearRkMat(intermBlCl &block); /*!< \brief Clear the low rank matrix of the block. */
    static void clearBlock(intermBlCl &block); /*!< \brief Clear the low rank  and full rank matrix of the block. */

    static void RkMatRankReduction(Eigen::MatrixXcd &UMat, Eigen::VectorXcd &singVals, Eigen::MatrixXcd &VAdjMat, long rank, const double relError); /*!< \brief Truncate the low rank matrix (U D VAdj). If rank and relError are unequal to zero, the method truncates to the relative tolerance relError within the rank limitation. */
    static void roundedAgglomerateSonsIntoBlock(intermBlCl &block, const long rank, const double relError); /*!< \brief Truncated agglomeration of the sub tree into the block. */
    static std::tuple<Eigen::MatrixXcd, Eigen::VectorXcd, Eigen::MatrixXcd> agglomerateSons(intermBlCl &block, const long rank, const double relError); /*!< \brief Return the agglomeration low rank matrix of the sub tree. */

//    static void intermBlClToBlockCluster(BlockCluster* factorBlock1, intermBlCl* block); /*!< \brief Convert the intermediary product blockcluster into a ordinary block cluster. */
    static std::unique_ptr<BlockCluster> intermBlClToBlockCluster(std::unique_ptr<intermBlCl> &block); /*!< \brief Convert the intermediary product blockcluster into an ordinary block cluster. */
};

#endif // HMULTIPLY_H
