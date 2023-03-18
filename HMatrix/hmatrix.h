#ifndef HMATRIX_H
#define HMATRIX_H

#include "clustertree.h"
#include <math.h>
#include "golubReinschSVD.h"
#include <functional>

/**
* \class BlockCluster
* \brief A recursive container for a matrix block.
*
* An H-matrix block can itself be cross partitioned into four subblocks. An H-matrix is a quaternary tree made up of BlockClusters as nodes.
*/
class BlockCluster
{
    public:
    BlockCluster(){}
    BlockCluster(Cluster* rowCluster, Cluster* columnCluster)
    {
        this->rowCluster = rowCluster;
        this->columnCluster = columnCluster;
    }
    BlockCluster(Cluster* rowCluster, Cluster* columnCluster, bool isRoot)
    {
        this->rowCluster = rowCluster;
        this->columnCluster = columnCluster;
        this->isRoot = isRoot;
    }
    BlockCluster(Cluster* rowCluster, Cluster* columnCluster, BlockCluster* father)
    {
        this->rowCluster = rowCluster;
        this->columnCluster = columnCluster;
        this->father = father;
    }

    long rows() const; /*!< \brief Get number of rows of the block. */
    long cols() const; /*!< \brief Get number of colums of the block. */
    long size() const; /*!< \brief Get the size of the block. */
    long rowStartIndex() const; /*!< \brief Get the global index of the first block row.*/
    long colStartIndex() const; /*!< \brief Get the global index of the first block column.*/
    double compression() const; /*!< \brief Get the ratio of actual memory requirements of the block compared to the full rank requirements. */
    void clearMatrixInfo(); /*!< \brief Clear the full and low rank matrices in the block. */
    void clear(); /*!< \brief Clear the entire subtree (the block and all subblocks). */
    void trimBelow(); /*!< \brief Trim the subtree of the block. */
    void compressAdmissibleBlock(const long rank, const double error); /*!< \brief Truncate the low rank matrix. */
    void compressSVDBlock(const long rank, const double error); /*!< \brief Truncate the low rank matrix. The matrix is assumed to be already in SVD form. */
    bool hasNan() const; /*!< \brief Check whether the block contains any NaNs. */
    double norm(); /*!< \brief Get the Frobenius norm of the block. */
    bool hasFourChildren() const;  /*!< \brief Check whether the block is cross partitioned. */
    static void getPartition(const BlockCluster &block, QVector<const BlockCluster *> &partition); /*!< \brief Calculate the partition vector that corresponds to the H-matrix/blockcluster tree.  */
    BlockCluster *returnCopy(BlockCluster* father = nullptr, long rank = 0, double error = 0, bool originalIsInSVDFormat = false); /*!< \brief Returns a pointer to a copy of the entire subtree of this. */

    Eigen::VectorXcd row(long colIndex) const; /*!< \brief Get block row. */
    Eigen::VectorXcd col(long colIndex) const; /*!< \brief Get block column. */
//    Eigen::MatrixXcd block(long rowIndex, long colIndex, long numberOfRows, long numberOfCols) const;  /*!< \brief Get sub block of the block as a full matrix. */

    /**
    * \brief Full rank assembly of an H-block by ACA with full pivoting.
    * \param implicitMatrix The implicit matrix used to assemble the block.
    */
    void fullRankAssembly(std::function<std::complex<double> (long, long)> implicitMatrix);

    /**
    * \brief Low rank assembly of an H-block by ACA with full pivoting.
    * This function performs low rank assembly of an H-block by ACA with full pivoting.
    * \param rank The rank of the H-matrix block.
    * \param relativeError The relative approximation error of the adaptive cross approximation.
    * \param implicitMatrix The implicit matrix used in the ACA.
    */
    void fullPivotACA(const long rank, const double relativeError, std::function<std::complex<double> (long, long)> implicitMatrix);

    /**
    * \brief Low rank assembly of an H-block by ACA with partial pivoting.
    * This function performs low rank assembly of an H-block by ACA with partial pivoting.
    * \param rank The rank of the H-matrix block.
    * \param relativeError The relative approximation error of the adaptive cross approximation.
    * \param implicitMatrix The implicit matrix used in the ACA.
    */
    void partialPivotACA(const long rank, const double relativeError, std::function<std::complex<double> (long, long)> implicitMatrix);

    /**
    * \brief Low-rank assembly of an H-block by ACA with heuristic partial pivoting and additional error control heuristic.
    * This function performs low-rank assembly of an H-block by ACA with heuristic partial pivoting.
    * The (guaranteed) accuracy is improved for dPhi-blocks.
    * \param rank The rank of the H-matrix.
    * \param relativeError The relative approximation error of the adaptive cross approximation.
    * \param getPivotIndices The function that returns the pivot indices to use in the ACA.
    * \param implicitMatrix The implicit matrix used in the ACA.
    */
    void partialPivotACAextra(const long rank, const double relativeError, std::function<QVector<std::pair<long,long>>(BlockCluster*,double)> getPivotIndices, std::function<std::complex<double> (long, long)> implicitMatrix);

    /**
    * This function calculates a row vector of an H-matrix block by using the given implicit matrix function. It is used in the ACA algorithms for the low-rank assembly of the H-matrix blocks.
    * \param rowVector The row vector to be calculated.
    * \param rowIndex The index of the row of the H-matrix block.
    * \param columnStartIndex The start index of the columns of the H-matrix block.
    * \param implicitMatrix The implicit matrix function that is used to evaluate the elements of the row vector.
    */
    void calcHBlockRowFromImplicit(Eigen::RowVectorXcd &rowVector, const long rowIndex, const long columnStartIndex, std::function<std::complex<double> (long, long)> implicitMatrix);

    /**
    * This function calculates a column vector of an H-matrix block by using the given implicit matrix function. It is used in the ACA algorithms for the low-rank assembly of the H-matrix blocks.
    * \param columnVector The column vector to be calculated.
    * \param rowStartIndex The start index of the rows of the H-matrix block.
    * \param columnIndex The index of the column of the H-matrix block.
    * \param implicitMatrix The implicit matrix function that is used to evaluate the elements of the column vector.
    */
    void calcHBlockColumnFromImplicit(Eigen::VectorXcd &columnVector, const long rowStartIndex, const long columnIndex, std::function<std::complex<double> (long, long)> implicitMatrix);

    Cluster* rowCluster = nullptr;
    Cluster* columnCluster = nullptr;
    BlockCluster* father = nullptr;
    BlockCluster* son11 = nullptr;
    BlockCluster* son12 = nullptr;
    BlockCluster* son21 = nullptr;
    BlockCluster* son22 = nullptr;
    bool isAdmissible = false;
    bool isRoot = false;
    bool isLeaf = false;

    Eigen::MatrixXcd fullMat; /*!< \brief Full matrix. */
    Eigen::MatrixXcd UMat; /*!< \brief Left factor of low rank matrix. */
    Eigen::VectorXcd singularValues; /*!< \brief Diagonal matrix factor of low rank matrix. Usually in SVD form. */
    Eigen::MatrixXcd VAdjMat; /*!< \brief Right factor of low rank matrix. */
    double frobeniusNorm = -1;  /*!< \brief Frobenius norm of the matrix. */

    //EIGEN_MAKE_ALIGNED_OPERATOR_NEW //only needed for fixed size eigen objects
//    int localRank = 0;
};

/**
* \class HMatrix
* \brief An H-matrix class.
*
*  An H-matrix is a quaternary tree made up of BlockClusters.
*/
class HMatrix
{
public:
    HMatrix(){};

    HMatrix(BlockCluster *rootBlockCluster) //
    {
        this->rootBlockCluster = rootBlockCluster;
        this->frobeniusNorm = rootBlockCluster->frobeniusNorm;
        this->rootBlockCluster->isRoot = true;
    }

    HMatrix(ClusterTree* rowClustertree, ClusterTree* columnClustertree, bool populate = false)
    {
        this->rowClustertree = rowClustertree;
        this->columnClustertree = columnClustertree;
        if(populate)
        {
            this->populateBlockClusterTree();
        }
    }

    HMatrix(ClusterTree* rowClustertree, ClusterTree* columnClustertree, const long rank, double relativeError, std::function<std::complex<double> (long, long)> implicitMatrix)
    {
        this->rowClustertree = rowClustertree;
        this->columnClustertree = columnClustertree;
        this->populateBlockClusterTree();
        assembleBlocks(rank, relativeError, implicitMatrix);
    }

    HMatrix(const HMatrix &original, long rank = 0, double error = 0, bool originalIsInSVDFormat = false); /*!< \brief Construct a truncated copy of the original HMatrix. If originalIsInSVDFormat the copy truncation is a lot faster. */

    ~HMatrix(){
        minPartition.clear();
    }

//    HMatrix& operator=(const HMatrix &hMat) = default; // use default copy assignment operator
    HMatrix& operator=(const HMatrix &hMat) /*!< \brief The copy assignment operator */
    {
        if(this->rootBlockCluster != nullptr)
        {
            this->rootBlockCluster->clear();
//            std::cout << "clear in copy assign constructor called!" << std::endl;
        }
        this->rootBlockCluster = hMat.rootBlockCluster;

        this->minPartition = hMat.minPartition;
        this->nearFieldBlocks = hMat.nearFieldBlocks;
        this->farFieldBlocks = hMat.farFieldBlocks;
        this->rowClustertree = hMat.rowClustertree;
        this->columnClustertree = hMat.columnClustertree;
        this->frobeniusNorm = hMat.frobeniusNorm;

        return *this;
    }

    void createRootNode(); /*!< \brief Create the root block for the blockclustertree. */

    void populateBlockClusterTree(bool randBlocks = false); /*!< \brief Set up the blockcluster tree structure. */
    void setClusterTrees(ClusterTree* rowClustertree, ClusterTree* columnClustertree); /*!< \brief Set the clustertrees for the blockcluster tree. */
    QVector<BlockCluster*> getMinPartition(); /*!< \brief Assemble a vector of the partition blocks. */

    void assembleBlocks(const long rank, double relativeError, std::function<std::complex<double> (long, long)> implicitMatrix, bool updatePartition = true);
    void assembleBlocksExtra(const long rank, double relativeError,  std::function<QVector<std::pair<long,long>>(BlockCluster*,double)> getPivotIndices, std::function<std::complex<double> (long, long)> implicitMatrixLowRank, bool updatePartition = true);

    QVector<BlockCluster*> getDiagonalBlocks() const; /*!< \brief Assemble a vector of the diagonal partition blocks. */
    BlockCluster* getRootBlock() const; /*!< \brief Get the root block of the H-matrix. */
    void setRootBlock(BlockCluster* rootBlockCluster); /*!< \brief Set the root block of the H-matrix. */
    ClusterTree* getRowClustertree(){return rowClustertree;} /*!< \brief Get the clustertree for the matrix row indices. */
    ClusterTree* getColumnClustertree(){return columnClustertree;} /*!< \brief Get the clustertree for the matrix column indices. */
    void clear(bool callMalloc_trim = false); /*!< \brief Clear the entire blockcluster tree. */
    void recursiveClear(BlockCluster* blockClusterNode);  /*!< \brief Recursively clear the blockcluster tree. */
    bool isCrossPartitioned() const; /*!< \brief Check whether the blocks in the blockclustertree are only cross partitioned. */
    bool isBlockAdmissible(const BlockCluster *blockCluster); /*!< \brief Check whether the block is admissible. (Has suitable low rank representation.) */
    void updatePartition(); /*!< \brief Update the partition vector. */
    void updatePartitionWithNullptrCheck(); /*!< \brief Update the partition vector. Check for nullpointers. */
    void updatePartitionWithMatrixInfo(); /*!< \brief Update the partition vector with all blocks that contain actual matrix information. */
    static QVector<BlockCluster*> getSubBlocksWithMatrixInformation(BlockCluster* blockCluster); /*!< \brief Get all tree blocks that contain actual matrix information. */

    long rows() const; /*!< \brief Get number of rows of the matrix. */
    long cols() const; /*!< \brief Get number of colums of the matrix. */
    long size() const; /*!< \brief Get the size of the matrix. */
    long rowStartIndex() const; /*!< \brief Get the global index of the first matrix row.*/
    long colStartIndex() const; /*!< \brief Get the global index of the first matrix column.*/
    double getCompressionRatio(); /*!< \brief Get the ratio of actual (low rank) memory requirements compared to the full rank requirements. */
    bool consistencyCheck() const; /*!< \brief Check the consistency of the blockcluster tree. */
    long calculateMaxBlockRank(); /*!< \brief Get the maximum local block rank (of all low rank blocks). */
    long calculateSizeOfAdmissibleBlocks(); /*!< \brief Calculate the sum of the size of all low rank blocks. */
    double norm(bool calcMinPartition = true); /*!< \brief Get the Frobenius norm of the H-matrix. */
    double normInSVDForm(bool calcMinPartition = true); /*!< \brief Get the Frobenius norm of the H-matrix. Assumes that all low rank blocks are already in SVD representaion. (faster) */
    bool hasNan(bool calcMinPartition = true); /*!< \brief Check whether the H-matrix contains any NaNs. */
    Eigen::MatrixXcd toFullMat(); /*!< \brief Convert the H-matrix to a full matrix. Costs O(n*m). */
//    QVector<BlockCluster*> minPartition;
    int nearFieldBlocks = 0;
    int farFieldBlocks = 0;

private:
//    void createRootNode();
    void createBlockClusterChildren(BlockCluster* blockClusterNode); /*!< \brief Recursively set up the blockcluster tree structure. */
    void setUpRandomMatrices(long maxRank); /*!< \brief Initialize random matrices in the partition blocks. */

    static constexpr double admissibilityConstant = 3; /*!< \brief Constant used in the admissibility condition for a block. */
    double clusterDistance(Cluster* cluster1, Cluster* cluster2); /*!< \brief Calculate the distance between two clusters. */
    double clusterSize(Cluster* cluster) const; /*!< \brief Calculate the size of a cluster. */
    void checkBlockPartition(BlockCluster* blockCluster, bool &onlyCrossPartitions) const; /*!< \brief Check whether all blocks are only cross partitioned. */
    void recursiveAppendPartitionBlock(BlockCluster* blockCluster); /*!< \brief Recursively assemble the partition vector. */
    void recAppendPartitionBlockNullptrCheck(BlockCluster* blockCluster, int depth);
    static void recAppendPartitionBlockMatInfo(BlockCluster* blockCluster, QVector<BlockCluster*> &subBlocksWithInformation);
    void recursiveAppendDiagonalBlock(BlockCluster* blockCluster, QVector<BlockCluster*> &diagonalBlocks) const; /*!< \brief Collect pointers to the diagonal partition blocks in vector. */
    bool recursiveConsistencyCheck(BlockCluster* blockClusterNode) const; /*!< \brief Recursively check the consistency of the blockcluster tree. */

    bool useWeakAdmissibility = false;

    int maxRank = 5; // unused in bem solver
    double frobeniusNorm = -1;
    ClusterTree* rowClustertree = nullptr;
    ClusterTree* columnClustertree = nullptr;
    BlockCluster* rootBlockCluster = nullptr;
    QVector<BlockCluster*> minPartition;
};

#endif // HMATRIX_H
