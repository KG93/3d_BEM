#ifndef HARITHM_H
#define HARITHM_H

#include "hmatrix.h"
#include "Timer.h" //Timer
#include "hmatrixvisuals.h"
#include <eigen3/Eigen/SVD>
#include <eigen3/Eigen/Eigenvalues> // for operator norm
#include "hmultiply.h"

static constexpr bool useEigenSvd = false; /*!< Use Eigens SVD implementation if true or the own Golub Reinsch SVD if false. */
static constexpr bool inversionInLUDecomp = false; /*!< Use the inverse matrix instead of the LU factorization for the diagonal blocks. */
static constexpr bool qrInLUDecomp = true; /*!< Use the QR decomposition instead of the LU factorization for the diagonal blocks. Recommended. */

/**
* \class HArithm
*
* \brief Class containing methods for performing operations on H-matrices.
*
* This class contains methods for performing various operations on H-matrices, such as H-matrix-vector multiplication, vector-H-matrix multiplication, and converting H-matrices to full matrices.
* The class also contains methods for compressing H-matrices, as well as methods for performing rank reduction on H-matrices and their blocks. The class uses the Eigen library for linear algebra operations.
*/
class HArithm
{
public:
    /**
    * \brief H-matrix-vector product. y += hmatrix * x
    * \param y Reference to the vector to which the result will be added.
    * \param hmatrix The h-matrix that will be multiplied.
    * \param x The vector that will be multiplied.
    */
    static void MVM(Eigen::VectorXcd &y, HMatrix &hmatrix, const Eigen::VectorXcd &x);

    /**
    * \brief Vector-H-matrix product.  y' += x' * hmatrix
    * \param y Reference to the row vector to which the result will be added.
    * \param x The row vector that will be multiplied.
    * \param hmatrix The h-matrix that will be multiplied.
    */
    static void VMM(Eigen::RowVectorXcd &y, const Eigen::RowVectorXcd &x, HMatrix &hmatrix);

    /**
    * \brief Compress (truncate) all admissible H-matrix blocks. If either control parameter is zero, the parameter will be ignored.
    * \param hmatrix The h-matrix to be compressed.
    * \param rank The rank to which the blocks will be truncated.
    * \param error The maximum relative error allowed.
    */
    static void compressHMat(HMatrix &hmatrix, const long rank, const double error);

    /**
    * \brief Get the minimum rank that satisfies the relative error.
    * \param svd The SVD object that will be used for the calculation.
    * \param maxRank The maximum rank allowed.
    * \param relError The maximum relative error allowed.
    * \return The minimum rank that satisfies the
    */
    static long minRankforError(const Eigen::BDCSVD<Eigen::MatrixXcd,Eigen::ComputeThinU|Eigen::ComputeThinV> &svd, const long maxRank, const double relError);

    /**
    * \brief Bring low rank matrix in SVD form.
    *
    * \param matrixBlock matrix block to be decomposed
    */
    static void svdOnBlock(BlockCluster &block);

    /**
    * \brief Truncate the low rank matrix of the block.
    *
    * \param matrixBlock matrix block
    * \param rank rank to compress to
    * \param relError maximum relative error
    */
    static void RkMatRankReduction(BlockCluster &block, const long rank, const double relError);

    /**
    * \brief Convert the full matrix of the block to a truncated low rank matrix.
    *
    * This method converts the full matrix of the given block to a truncated low rank matrix, using the specified rank and relative error.
    * \param matrixBlock The matrix block to be converted.
    * \param rank The rank to which the matrix will be truncated.
    * \param error The maximum relative error allowed.
    */
    static void fullMatRankReduction(BlockCluster &block, const long rank, const double error);

    /**
    * \brief Add two H-matrices with the same partition.
    *
    * This method adds two H-matrices with the same partition, and stores the result in the first H-matrixrix.
    * \param matrix1 The first H-matrix, to which the result will be stored.
    * \param matrix2 The second H-matrix, which will be added to the first.
    */
    static void addHMat2ToHMat1SamePartition(HMatrix &matrix1, HMatrix &matrix2, const long rank, const double error);

    static void matByBlock(const Eigen::MatrixXcd &lFactor, const BlockCluster &rFactorBlock, long columnOffsetFactor1, long columnOffsetFactor2, Eigen::MatrixXcd &product);
    static void blockByMat(const BlockCluster &lFactorBlock, const Eigen::MatrixXcd &rFactor, long rowOffsetFactor1, long rowOffsetFactor2, Eigen::MatrixXcd &product);

    /**
    * \brief Add and truncate a low-rank matrix to a block cluster.
    *
    * This method adds a low-rank matrix to a block cluster, truncating the result to a maximum rank and maximum relative error after performing the addition.
    * \param block Reference to the block cluster to which the low-rank matrix will be added.
    * \param UMat The U matrix of the low-rank matrix in SVD form.
    * \param singVals The singular values of the low-rank matrix in SVD form.
    * \param VAdjMat The adjoint V matrix of the low-rank matrix in SVD form.
    * \param maxRank The maximum rank of the truncated low-rank matrix addition.
    * \param relError The maximum relative error allowed for the truncated low-rank matrix addition.
    */
    static void roundedAddRMatToBlock(BlockCluster &block, const Eigen::MatrixXcd &UMat, const Eigen::VectorXcd &singVals, const Eigen::MatrixXcd &VAdjMat, const long maxRank, const double relError);

    /**
    * \brief Add a low-rank matrix to a block cluster.
    *
    * This method adds a low-rank matrix to a block cluster, using the U, singular values, and adjoint V matrices in the low-rank matrix's SVD form.
    * \param block Reference to the block cluster to which the low-rank matrix will be added.
    * \param UMat The U matrix of the low-rank matrix in SVD form.
    * \param singVals The singular values of the low-rank matrix in SVD form.
    * \param VAdjMat The adjoint V matrix of the low-rank matrix in SVD form.
    */
    static void addRMatToBlock(BlockCluster &block, const Eigen::MatrixXcd &UMat, const Eigen::VectorXcd &singVals, const Eigen::MatrixXcd &VAdjMat);

    /**
    * \brief Join two matrices horizontally into the first matrix.
    *
    * This method horizontally concatenates the second matrix to the first matrix, modifying the first matrix in place.
    * \param firstMat Reference to the matrix to which the second matrix will be joined.
    * \param secondMat The matrix that will be joined to the first matrix.
    */
    static void horizontalJoinMatricesInToFirstOne(Eigen::MatrixXcd &firstMat, const Eigen::MatrixXcd &secondMat);

    /**
    * \brief Join two matrices vertically into the first matrix.
    *
    * This method vertically concatenates the second matrix to the first matrix, modifying the first matrix in place.
    * \param firstMat Reference to the matrix to which the second matrix will be joined.
    * \param secondMat The matrix that will be joined to the first matrix.
    */
    static void verticalJoinMatricesInToFirstOne(Eigen::MatrixXcd &firstMat, const Eigen::MatrixXcd &secondMat);

    /**
    * \brief Join two vectors into the first one.
    * This function appends the second vector to the end of the first vector.
    * \param firstVector Reference to the vector to which the second vector will be appended.
    * \param secondVector The vector that will be appended to the first vector.
    */
    static void joinVectorsInToFirstOne(Eigen::VectorXcd &firstVector, const Eigen::VectorXcd &secondVector);

    /**
    * \brief Add full matrix to another full matrix.
    * \param firstMat Reference to the first matrix to which the second matrix will be added.
    * \param secondMat The matrix to be added.
    */
    static void addFullMatrixInToFirstOne(Eigen::MatrixXcd &firstMat, const Eigen::MatrixXcd &secondMat);

    /**
    * \brief Efficiently add low rank matrix into a full matrix.
    *
    * This method efficiently adds a low rank matrix, given in the form of its SVD, into a full matrix.
    * \param fullMat Reference to the full matrix into which the low rank matrix will be added.
    * \param UMat The U matrix in the SVD of the low rank matrix.
    * \param singVals The singular values vector in the SVD of the low rank matrix.
    * \param VAdjMat The adjoint of the V matrix in the SVD of the low rank matrix.
    */
    static void addrkMatToFull(Eigen::MatrixXcd &fullMat, const Eigen::MatrixXcd &UMat, const Eigen::VectorXcd &singVals, const Eigen::MatrixXcd &VAdjMat);

    /**
    * \brief Merge matrices of subtree into the ancestor block.
    *
    * This method merges the matrices of the subblocks of an H-matrix subtree into the full matrix of the ancestor block of the subtree.
    * \param ancestorBlock The ancestor block whose full matrix will be updated.
    * \param subBlock The subblock whose matrix will be added to the ancestor block's full matrix.
    */
    static void transportSubblockInToFullAncestorBlock(BlockCluster &ancestorBlock, BlockCluster &subBlock); /*!< Merge full matrices of subtree into the ancestor block. */

    /**
    * \brief Set low rank matrix to zero.
    *
    * This method sets the U, S, and V matrices of the given block to zero matrices (consisting of zero vectors).
    * \param block The block whose low rank matrix will be set to zero.
    */
    static void setRkMatrixZero(BlockCluster &block);

    /**
    * \brief Clear the low rank matrix of a block.
    * \param block The block whose low rank matrix will be cleared.
    */
    static void clearRkMatrix(BlockCluster &block);

    /**
    * \brief Agglomerate subtree into the full matrix of the block. Then trim subtree below the block.
    *
    * This function performs a full rank direkt agglomeration on a subtree, where the matrices/blocks of the subtree are added to the full matrix of the block. The subtree below the block is then trimmed off/deleted.
    * \param block Pointer to the root of the subtree.
    */
    static void fullRankDirektAgglomerationWithTrimming(BlockCluster &block);
    static void fullRankDirektAgglomeration(BlockCluster &block);
    static void fullRankAgglomerationWithTrimming(BlockCluster &block);
    static void fullRankAgglomeration(BlockCluster &block);
    static Eigen::MatrixXcd fullRankAgglomerationRecursion(BlockCluster &block);

    /**
    * \brief Agglomerate subtree into the low rank matrix of the block. Then trim subtree below the block.
    *
    * This method merges the full matrix of the given block and the matrices/blocks of its descendants into the low rank matrix of the block, using the given rank and error to truncate the resulting matrix. The method then trims the subtree below the block.
    * \param block The block whose subtree will be merged and trimmed.
    * \param rank The rank to which the resulting matrix will be truncated.
    * \param error The maximum error allowed in the truncated matrix.
    */
    static void reducedRankAgglomerationWithTrimming(BlockCluster &block, const long rank, const double relError); /** \brief Agglomerate subtree into the low rank matrix of the block. Then trim subtree below the block. */
    static void reducedRankAgglomeration(BlockCluster &block, const long rank, const double relError); /** \brief Agglomerate subtree into the low rank matrix of the block. */
    static std::tuple<Eigen::MatrixXcd, Eigen::VectorXcd, Eigen::MatrixXcd> reducedRankAgglomerationRecursion(BlockCluster &block, const long rank, const double relError);

    static void blockSplitting(BlockCluster &block, const long rank, const double error);  /*!< \brief Flush the matrices of the block down to the leaf nodes. */

    static Eigen::VectorXcd LUSolve(HMatrix &A, Eigen::VectorXcd &b, const long rank, const double relError);  /*!< \brief Solve  Ax=b  for x with H-LU-factorization. */
    static Eigen::VectorXcd LUSubstitutionSolve(HMatrix &L, HMatrix &U, const Eigen::VectorXcd &rightHandSide);   /*!< \brief Solve  LUx=b for x via forward and backward substitution. */
    static std::pair<HMatrix,HMatrix> LUDecomposition(HMatrix &matrix, const long rank, const double relError, const bool zeroBlocksUninitialized = false);   /*!< \brief Calculate the LU-factorization of an H-matrix. The operation destroys the matrix. */
    static void recursiveLUDecomposition(BlockCluster &LBlock, BlockCluster &UBlock, BlockCluster &ABlock, const long rank, const double relError,  const bool zeroBlocksUninitialized = false);  /*!< \brief Recursive LU factorization. (L*U=A) */
    static Eigen::VectorXcd forwardSubstitution(const HMatrix &L, Eigen::VectorXcd b); /*!< \brief Solve L * y = b for y; L is lower triangular H-matrix. */
    static void forwardSubstitution(const BlockCluster &LBlock, Eigen::VectorXcd &solution, Eigen::VectorXcd &b, const long vectorStartIndex = 0); /*!< \brief Solve L * y = b for y; L is lower triangular H-matrix. */
    static void forwardSubstitutionTransposed(Eigen::RowVectorXcd &solution, const BlockCluster &UBlock, Eigen::RowVectorXcd &b, const long vectorStartIndex = 0); /*!< \brief Solve y^T * U = b^T for y; U is upper triangular H-matrix. */
    static Eigen::VectorXcd backwardSubstitution(const HMatrix &U, Eigen::VectorXcd y); /*!< \brief Solve U * x = y for x; U is upper triangular H-matrix. */
    static void backwardSubstitution(const BlockCluster &UBlock, Eigen::VectorXcd &solution, Eigen::VectorXcd &y, const long vectorStartIndex = 0); /*!< \brief Solve U * x = y for x; U is upper triangular H-matrix. */
    static void backwardSubstitutionTransposed(Eigen::RowVectorXcd &solution, const BlockCluster &LBlock, Eigen::RowVectorXcd &y, const long vectorStartIndex = 0); /*!< \brief Solve x^T * L = y^T for x; L is lower triangular H-matrix*/
    static HMatrix forwSubsMatVal(const HMatrix &L, HMatrix &Z, const long rank, const double relError); /*!< \brief Solve L * X = Z for X, L is lower triangular H-matrix; X is the return value; Z is destroyed */
    static void forwSubsMatVal(const BlockCluster &LBlock, BlockCluster &XBlock, BlockCluster &ZBlock, const long rank, const double relError); /*!< \brief Solve L * X = Z for X, L is lower triangular H-matrix */
    static HMatrix forwSubsMatValTransposed(const HMatrix &U, HMatrix &Z, const long rank, const double relError); /*!< \brief Solve X * U = Z for X, U is upper triangular H-matrix; X is the return value */
    static void forwSubsMatValTransposed(BlockCluster &XBlock, const BlockCluster &UBlock, BlockCluster &ZBlock, const long rank, const double relError); /*!< \brief Solve X * U = Z for X, U is upper triangular H-matrix; Z is destroyed */

    static void recursiveMatrixVectorPoduct(Eigen::VectorXcd &product, const BlockCluster &block, const Eigen::VectorXcd &x, const long prodStartIndex, const long xStartIndex); /*!< \brief Recursively calculate y += hmatrix * x */
    static void parallelMatrixVectorPoduct(Eigen::VectorXcd &product, const BlockCluster &block, const Eigen::VectorXcd &x, const long prodStartIndex, const long xStartIndex); /*!< \brief Calculate y += hmatrix * x; better parallelized. */
    static void subtractiveRecursiveMatrixVectorPoduct(Eigen::VectorXcd &product, const BlockCluster &block, const Eigen::VectorXcd &x, const long prodStartIndex, const long xStartIndex); /*!< \brief Recursively calculate y -= hmatrix * x */
    static void subtractiveParallelMatrixVectorPoduct(Eigen::VectorXcd &product, const BlockCluster &block, const Eigen::VectorXcd &x, const long prodStartIndex, const long xStartIndex); /*!< \brief Calculate y -= hmatrix * x; better parallelized. */
    static void recursiveVectorMatrixPoduct(Eigen::RowVectorXcd &product, const Eigen::RowVectorXcd &x, const BlockCluster &block, const long prodStartIndex, const long xStartIndex); /*!< \brief Recursively calculate y += x * hmatrix; y and x are row-vectors */
    static void subtractiveRecursiveVectorMatrixPoduct(Eigen::RowVectorXcd &product, const Eigen::RowVectorXcd &x, const BlockCluster &block, const long prodStartIndex, const long xStartIndex); /*!< \brief Recursively calculate y -= x * hmatrix; y and x are row-vectors */
    static void subtractiveParallelVectorMatrixPoduct(Eigen::RowVectorXcd &product, const Eigen::RowVectorXcd &x, const BlockCluster &block, const long prodStartIndex, const long xStartIndex); /*!< \brief Recursively calculate y -= x * hmatrix; y and x are row-vectors; better parallelized */
    static void recursiveHMatAddition(BlockCluster &mat1Block, const BlockCluster &mat2Block, const long rank, const double relError); /*!< \brief Recursive addition for H-matrices based on same clustertree; the partitions can be different */
    static void recursiveHMatSubstraction(BlockCluster &mat1Block, const BlockCluster &mat2Block, const long rank, const double relError); /*!< \brief Recursive subtraction for H-matrices based on same clustertree; the partitions can be different */

    static void multiplyHMatByMinusOne(HMatrix &hMat); /*!< \brief Multiply matrix by -1. */
    static void recursiveMultiplyHMatByMinusOne(BlockCluster &matBlock); /*!< \brief Multiply subtree by -1. */
    static void convertToAdmissibleZeroBlock(BlockCluster &matBlock); /*!< \brief Make block zero block in low rank representation. */

    static std::unique_ptr<BlockCluster> copyBlock(const BlockCluster &block, Cluster* rowCluster, Cluster* colCluster); /*!< \brief Make independent copy of the subtree of BlockCluster. */

    static Eigen::VectorXcd getRankKBlockDiagonal(const BlockCluster &matBlock, long rank = 0);
    static Eigen::VectorXcd getRankKBlockIndexedElements(BlockCluster &matBlock, QVector<long> &rowIndices, QVector<long> &colIndices, long rank = 0);

    static double spectralNorm(HMatrix &A, const unsigned long maxIterations = 20); /*!< \brief Approximates the spectral norm a matrix via power iteration. */
    static double spectralNormFromLU(HMatrix &L, HMatrix &U, const unsigned long maxIterations = 20); /*!< \brief Approximates the spectral norm of the inverse of a matrix via the LU-decomposition of the matrix. The method is a power iteration with forward- and backward substitution. */

    static double spectralNormLeastSignificant(HMatrix &A, const unsigned long maxIterations = 20);  /*!< \brief Approximates the spectral norm of an H-matrix, where only the least significant rank-1 matrix of each partition block is taken into account. Can be used as an approximate upper bound of the spectral norm of the residual matrix after ACA assembly. */
    static void MVMLeastSignificant(Eigen::VectorXcd &y, HMatrix &hmatrix, const Eigen::VectorXcd &x);
    static void VMMLeastSignificant(Eigen::RowVectorXcd &y, const Eigen::RowVectorXcd &x, HMatrix &hmatrix);

    static double frobeniusNormFromLU(HMatrix &L, HMatrix &U, const unsigned long maxIterations = 20); /*!< \brief Approximates the frobenius norm of the inverse of a matrix via the LU-decomposition of the matrix. Low accuracy.*/

    static void transpose(HMatrix &matrix); /*!< \brief Transpose the HMatrix in place. */

    static void transpose(BlockCluster &matBlock); /*!< \brief Recursive transposition the subtree of the BlockCluster. */
};
#endif // HARITHM_H
