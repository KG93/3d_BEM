#ifndef GMRES_H
#define GMRES_H

#include "harithm.h"
#include <eigen3/unsupported/Eigen/IterativeSolvers>


class HMatrixWrapper;
using Eigen::SparseMatrix;

namespace Eigen {
    namespace internal {
      template<>
      struct traits<HMatrixWrapper> :  public Eigen::internal::traits<Eigen::SparseMatrix<std::complex<double>> >
      {};
    }
}


/**
* \brief Wrapper class for an HMatrix into an Eigen compatible matrix type.
*/
class HMatrixWrapper : public Eigen::EigenBase<HMatrixWrapper>
{

public:
    typedef std::complex<double> Scalar;
    typedef double RealScalar;
    typedef int StorageIndex;
    enum
    {
        ColsAtCompileTime = Eigen::Dynamic,
        MaxColsAtCompileTime = Eigen::Dynamic,
        IsRowMajor = false
    };

    Eigen::Index rows() const { return hMatrix->rows(); }
    Eigen::Index cols() const { return hMatrix->cols(); }

    template<typename Rhs>
    Eigen::Product<HMatrixWrapper,Rhs,Eigen::AliasFreeProduct> operator*(const Eigen::MatrixBase<Rhs>& x) const {
      return Eigen::Product<HMatrixWrapper,Rhs,Eigen::AliasFreeProduct>(*this, x.derived());
    } /*!< Implements the matrix vector product for the matrix free GMRES api. */

    HMatrixWrapper(HMatrix* hMatrix)
    {
        this->hMatrix = hMatrix;
    }

    HMatrix* hMatrix;
};

namespace Eigen {
    namespace internal {

        template<typename Rhs>
        struct generic_product_impl<HMatrixWrapper, Rhs, SparseShape, DenseShape, GemvProduct> // GEMV stands for matrix-vector
        : generic_product_impl_base<HMatrixWrapper,Rhs,generic_product_impl<HMatrixWrapper,Rhs> >
        {
            typedef typename Product<HMatrixWrapper,Rhs>::Scalar Scalar;

            template<typename Dest>
            static void scaleAndAddTo(Dest& dst, const HMatrixWrapper& lhs, const Rhs& rhs, const Scalar& alpha)
            {
                assert(alpha==Scalar(1) && "scaling is not implemented");
                EIGEN_ONLY_USED_FOR_DEBUG(alpha);

                HArithm::MVM(dst, *lhs.hMatrix, rhs); // y += hmatrix * x
            }
        };
    }
}

class LUPrecondidionedHMatrixWrapper;
using Eigen::SparseMatrix;

namespace Eigen {
namespace internal {
  template<>
  struct traits<LUPrecondidionedHMatrixWrapper> :  public Eigen::internal::traits<Eigen::SparseMatrix<std::complex<double>> >
  {};
}
}

/**
* \brief Wrapper class for an HMatrix and an HLU Matrix for preconditioning into an Eigen compatible matrix type.
*/
class LUPrecondidionedHMatrixWrapper : public Eigen::EigenBase<LUPrecondidionedHMatrixWrapper>
{

public:

    typedef std::complex<double> Scalar;
    typedef double RealScalar;
    typedef int StorageIndex;
    enum
    {
        ColsAtCompileTime = Eigen::Dynamic,
        MaxColsAtCompileTime = Eigen::Dynamic,
        IsRowMajor = false
    };

    Eigen::Index rows() const { return hMatrix->rows(); }
    Eigen::Index cols() const { return hMatrix->cols(); }

    template<typename Rhs>
    Eigen::Product<LUPrecondidionedHMatrixWrapper,Rhs,Eigen::AliasFreeProduct> operator*(const Eigen::MatrixBase<Rhs>& x) const {
      return Eigen::Product<LUPrecondidionedHMatrixWrapper,Rhs,Eigen::AliasFreeProduct>(*this, x.derived());
    } /*!< Implements the preconditioned matrix vector product for the matrix free GMRES api. */

    LUPrecondidionedHMatrixWrapper(HMatrix* hMatrix, HMatrix* hLMatrix, HMatrix* hUMatrix)
    {
        this->hMatrix = hMatrix;
        this->hLMatrix = hLMatrix;
        this->hUMatrix = hUMatrix;
    }

    HMatrix* hMatrix;
    HMatrix* hLMatrix;
    HMatrix* hUMatrix;
};

namespace Eigen
{
    namespace internal
    {
        template<typename Rhs>
        struct generic_product_impl<LUPrecondidionedHMatrixWrapper, Rhs, SparseShape, DenseShape, GemvProduct> // GEMV stands for matrix-vector
        : generic_product_impl_base<LUPrecondidionedHMatrixWrapper,Rhs,generic_product_impl<LUPrecondidionedHMatrixWrapper,Rhs> >
        {
            typedef typename Product<LUPrecondidionedHMatrixWrapper,Rhs>::Scalar Scalar;

            template<typename Dest>
            static void scaleAndAddTo(Dest& dst, const LUPrecondidionedHMatrixWrapper& lhs, const Rhs& rhs, const Scalar& alpha)
            {
                assert(alpha==Scalar(1) && "scaling is not implemented");
                EIGEN_ONLY_USED_FOR_DEBUG(alpha);

                Eigen::VectorXcd copy = dst;
                long rows = lhs.hMatrix->rows();
                HArithm::MVM(copy, *lhs.hMatrix, rhs);
                Eigen::VectorXcd forwSubsSolution = Eigen::VectorXcd::Zero(rows);
                HArithm::forwardSubstitution(*(lhs.hLMatrix->getRootBlock()), forwSubsSolution, copy, 0); // solve L * y = b for y, L is lower triangular h-matrix
                Eigen::VectorXcd backwSubsSolution = Eigen::VectorXcd::Zero(rows);
                HArithm::backwardSubstitution(*(lhs.hUMatrix->getRootBlock()), backwSubsSolution, forwSubsSolution, 0); // U * x = y, U is upper triangular h-matrix
                dst += backwSubsSolution;
            }
        };
    }
}

/**
* \class GMRES
* \brief The GMRES class contains multiple methods for solving linear systems using the GMRES algorithm.
*/
class GMRES
{
public:
    GMRES();

    /**
    * \brief Solves the linear system Ax = b using the GMRES algorithm with an H-matrix
    * \param A is the system H-matrix
    * \param xGuess The starting vector for the iterative solution
    * \param b The right-hand side vector of the linear system
    * \param tolerance The tolerance for the solution of the linear system. An upper bound for the relative residual error.
    * \param iterations The number of iterations to perform before restarting the GMRES algorithm
    * \return The solution vector to the linear system Ax = b
    */
    static Eigen::VectorXcd gmresSolve(const HMatrixWrapper &A, Eigen::VectorXcd xGuess, const Eigen::VectorXcd b, double tolerance = 0.01, long iterations = 300)
    {
        Eigen::GMRES<HMatrixWrapper, Eigen::IdentityPreconditioner> gmres(A);
        gmres.setTolerance(tolerance);
        gmres.set_restart(iterations);

        Eigen::VectorXcd solution = gmres.solveWithGuess(b, xGuess);
        std::cout << "gmres.iterations(): " << gmres.iterations() << std::endl;
        std::cout << "gmres relative residual error estimate: " << gmres.error() << std::endl;
        return solution;
    }

    /**
    * \brief Solves the linear system Ax = b using the GMRES algorithm with an H-matrix and HLU-preconditioning
    * \param AwithLU is the system H-matrix with an HLU preconditioner
    * \param xGuess The starting vector for the iterative solution
    * \param b The right-hand side vector of the linear system
    * \param tolerance The tolerance for the solution of the linear system. An upper bound for the relative residual error.
    * \param iterations The number of iterations to perform before restarting the GMRES algorithm
    * \return The solution vector to the linear system Ax = b
    */
    static Eigen::VectorXcd gmresLUPreconditionedSolve(const LUPrecondidionedHMatrixWrapper &AwithLU, Eigen::VectorXcd xGuess, const Eigen::VectorXcd b, double tolerance = 0.01, long iterations = 300)
    {
        Eigen::GMRES<LUPrecondidionedHMatrixWrapper, Eigen::IdentityPreconditioner> gmres(AwithLU);
        gmres.setTolerance(tolerance); // sets the relative tolerance for the residual
        gmres.set_restart(iterations); // sets the amount of iterations after which the GMRES is restarted

        Eigen::VectorXcd rightHandSideCopy = b;
        Eigen::VectorXcd forwSubsSolution = Eigen::VectorXcd::Zero(AwithLU.rows());
        HArithm::forwardSubstitution(*(AwithLU.hLMatrix->getRootBlock()), forwSubsSolution, rightHandSideCopy); // solve L * y = b for y, L is lower triangular h-matrix
        Eigen::VectorXcd backwSubsSolution = Eigen::VectorXcd::Zero(AwithLU.rows());
        HArithm::backwardSubstitution(*(AwithLU.hUMatrix->getRootBlock()), backwSubsSolution, forwSubsSolution); // U * x = y, U is upper triangular h-matrix

        Eigen::VectorXcd solution = gmres.solveWithGuess(backwSubsSolution, xGuess);
        std::cout << "gmres.iterations(): " << gmres.iterations() << std::endl;
        std::cout << "gmres relative residual error estimate: " << gmres.error() << std::endl;
        return solution;
    }

    /**
    * \brief Solves the linear system Ax = b using the GMRES algorithm
    * \param A is the system matrix
    * \param xGuess The starting vector for the iterative solution
    * \param b The right-hand side vector of the linear system
    * \param tolerance The tolerance for the solution of the linear system. An upper bound for the relative residual error.
    * \param iterations The number of iterations to perform before restarting the GMRES algorithm
    * \return The solution vector to the linear system Ax = b
    */
    static Eigen::VectorXcd gmresSolve(Eigen::MatrixXcd A, Eigen::VectorXcd xGuess, const Eigen::VectorXcd b, double tolerance = 0.01, long iterations = 300)
    {
        Eigen::GMRES<Eigen::MatrixXcd, Eigen::IdentityPreconditioner> gmres(A);
        gmres.setTolerance(tolerance); // sets the relative tolerance for the residual
        gmres.set_restart(iterations); // sets the amount of iterations after which the GMRES is restarted

        Eigen::VectorXcd solution = gmres.solveWithGuess(b, xGuess);
        std::cout << "gmres.iterations(): " << gmres.iterations() << std::endl;
        std::cout << "gmres relative residual error estimate: " << gmres.error() << std::endl;
        return solution;
    }
};

#endif // GMRES_H
