#include "gls.h"

// Ax = By + c with alpha_i*x_i + beta_i * y_i = f_i
void GLS::calculateGLS(Eigen::MatrixXcd &A, Eigen::VectorXcd &x, Eigen::MatrixXcd &B, Eigen::VectorXcd &y, Eigen::VectorXcd &c, Eigen::VectorXcd &alpha, Eigen::VectorXcd &beta, Eigen::VectorXcd &f)
{
    std::cout<<"GLS solver started!"<<std::endl;
    int n = A.rows();
    if(A.rows()!=A.cols())
    {
        std::cout<<"Error! A.rows() != A.cols()"<<std::endl;
        return;
    }
    if(B.rows()!=B.cols())
    {
        std::cout<<"Error! B.rows() != B.cols()"<<std::endl;
        return;
    }
    if(A.rows()!=B.rows())
    {
        std::cout<<"Error! A and B have different dimensions!"<<std::endl;
        return;
    }
    if(A.rows()!=c.size())
    {
        std::cout<<"Error! c has different dimension than the  matrices!"<<std::endl;
        return;
    }
    if(A.rows()!=alpha.size())
    {
        std::cout<<"Error! alpha has different dimension than the  matrices!"<<std::endl;
        return;
    }
    if(A.rows()!=beta.size())
    {
        std::cout<<"Error! beta has different dimension than the  matrices!"<<std::endl;
        return;
    }
    if(A.rows()!=f.size())
    {
        std::cout<<"Error! f has different dimension than the  matrices!"<<std::endl;
        return;
    }
    if(A.rows()!=x.size())
    {
        std::cout<<"Error! x has different dimension than the  matrices!"<<std::endl;
        return;
    }
    if(A.rows()!=y.size())
    {
        std::cout<<"Error! y has different dimension than the  matrices!"<<std::endl;
        return;
    }

    // Compute the 1-norm of both matrices and the colwise 1-norm of both matrices
    double Anorm=A.lpNorm<1>();
    double Bnorm=B.lpNorm<1>();
    Eigen::VectorXd Acolnorms(n);
    Eigen::VectorXd Bcolnorms(n);

    #pragma omp parallel for
    for(int j=0; j<n; j++)
    {
        Acolnorms(j) = A.col(j).lpNorm<1>();
        Bcolnorms(j) = B.col(j).lpNorm<1>();
    }

    // alpha(i) and beta(i) can not both be zero
    for(int i = 0; i < n; i++)
    {
        if(abs(alpha(i)) < global::tiny * Anorm && abs(beta(i)) < global::tiny * Bnorm)
         {
            std::cout<<"alpha(i) and beta(i) are both zero. Abandoning the GLS routine!"<<std::endl;
            return;
        }
    }

    // The following line determines whether the system is initially solved for x_i or y_i
    Eigen::VectorXi solveForX = ((pow(beta.array().abs(),2) * Acolnorms.array()) > (pow(alpha.array().abs(),2) * Bcolnorms.array())).cast<int>();


    // Exchange the matrix columns according to the substitutions
    Eigen::VectorXcd tmp(n);
    #pragma omp parallel for private(tmp)
    for (int j =0; j<n; j++)
    {
        if (!solveForX(j))
        {
           tmp = A.col(j);
           A.col(j) = -B.col(j);
           B.col(j) = -tmp;
        }
    }

    // scale the columns of the B matrix
    std::complex<double> alphaOrBeta;
    #pragma omp parallel for private(alphaOrBeta)
    for (int j =0; j<n; j++)
    {
        if (solveForX(j))
        {
           alphaOrBeta = beta(j);
        }
        else
        {
           alphaOrBeta = alpha(j);
        }

        B.col(j) *= 1.0 / alphaOrBeta;
    }

    // Subtract the terms from the A matrix
    #pragma omp parallel for
    for (int j=0; j<n; j++)
    {
        if (solveForX(j))
        {
            A.col(j) += alpha(j)*B.col(j);
        }
        else
        {
            A.col(j) += beta(j)*B.col(j);
        }
    }

    // Calculate the righthandside vector
    Eigen::VectorXcd rightHandSide;
    rightHandSide = c + B * f;
    B.resize(0,0);

    Eigen::VectorXcd x_tmp;
    Eigen::VectorXcd y_tmp;

    // Solve A*x_tmp = rightHandSide via LU decomposition
    bool fullPivLu = false;
    if(fullPivLu == true)
    {
        std::cout << "LU factorization started!" << std::endl;
        Eigen::FullPivLU<Eigen::MatrixXcd> lu(A);
        x_tmp=lu.solve(rightHandSide);
    }
    else
    {
        std::cout<<"LU factorization started with "<<Eigen::nbThreads( )<<" Threads!"<<std::endl;
        Eigen::PartialPivLU<Eigen::MatrixXcd> luPP(A);
        x_tmp=luPP.solve(rightHandSide);
    }
    A.resize(0,0);
    y_tmp=x_tmp;

    // Resubstitute the missing entries of the solution
    for(int i=0; i<n; i++)
    {
        if (solveForX(i))
        {
            y_tmp(i) = (f(i) - alpha(i) * x_tmp(i)) / beta(i);
        }
        else
        {
            x_tmp(i) = (f(i) - beta(i) * y_tmp(i)) / alpha(i);
        }
    }
    x = x_tmp;
    y = y_tmp;
    std::cout<<"GLS solver finished!"<<std::endl;
}
