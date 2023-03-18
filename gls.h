#ifndef GLS_H
#define GLS_H

#include <QtGlobal>
#include <eigen3/Eigen/Geometry>
#include <iostream>
#include <QVector>
#include "global.h"
// Derived from "Solving the Linear Systems of Equations in the Generalized Direct Boundary Element Method" by Stephen Kirkup


/**
* \brief The GLS class solves a system of the form Ax=By+c with alpha(i)*x(i)+beta(i)y(i)=f(i) for i=0,...,n with A,B n-times-n matrices and x,y and c n-length vectors.
*/
class GLS
{
public:
    GLS();
    //solves Ax=By+c with alpha(i)*x(i)+beta(i)*y(i)=f(i) for i=0,...,n with A,B n-times-n matrices and x,y and c n-length vectors
    /**
    * \brief Solves Ax=By+c with alpha(i)*x(i)+beta(i)*y(i)=f(i) for i=0,...,n with A,B n-times-n matrices and x,y and c n-length vectors.
    * \param A n-times-n matrix
    * \param x n-length vector
    * \param B n-times-n matrix
    * \param y n-length vector
    * \param c n-length vector
    * \param alpha n-length vector
    * \param beta n-length vector
    * \param f n-length vector
    */
    static void calculateGLS(Eigen::MatrixXcd& A, Eigen::VectorXcd& x, Eigen::MatrixXcd& B, Eigen::VectorXcd& y, Eigen::VectorXcd& c, Eigen::VectorXcd& alpha, Eigen::VectorXcd& beta, Eigen::VectorXcd& f);
};

#endif // GLS_H
