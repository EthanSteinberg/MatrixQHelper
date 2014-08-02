#include <iostream>
#include <Eigen/Dense>
#include "R.hpp"
#include "MatrixQ.hpp"
#include <unsupported/Eigen/MatrixFunctions>
#include <Eigen/Eigenvalues> 
using Eigen::MatrixXd;
using Eigen::Matrix4d;
int main()
{
    Matrix4d blah;
    blah << -4.5, 0.5, 0.5000000000000001, 0.5, 
100, -6.0, 1.6666666666666667, 1.5000000000000002 ,
1.5, 2.5, -5.666666666666667, 2.625 ,
2.0, 3.0, 3.5, -4.625 ;

    //std::cout<<blah<<std::endl;


    MatrixQ ok(4,blah);
    //std::cout<<ok.transposedQ<<std::endl;
    //
    MatrixXd regResult = ok.transposedQ.exp();


    Eigen::EigenSolver<MatrixXd> solver(ok.transposedQ,true);

    Eigen::MatrixXd eigenVectors = solver.eigenvectors().real();
    Eigen::MatrixXd eigenInverseVectors = solver.eigenvectors().inverse().real();
    
    Eigen::VectorXd eigenValues = solver.eigenvalues().real();
    //std::cout<<eigenValues.imag().array().abs().sum()<<std::endl;
    //std::cout<<eigenVectors.imag().array().abs().sum()<<std::endl;


    //Eigen::MatrixXd diagonalAfterThing = eigenValues.array().exp().matrix().asDiagonal();
    //std::cout<<eigenValues<<std::endl;
    //std::cout<<diagonalAfterThing<<std::endl;
    //


    //Eigen::MatrixXd newResult = (eigenVectors * diagonalAfterThing * eigenInverseVectors);

    //std::cout<<(newResult - regResult).norm()<<std::endl;

    

    Eigen::VectorXd columnVector = Eigen::VectorXd::Random(ok.transposedQ.cols()).array().abs().matrix();
    
    double sum = 0;
    double sum1 = 0;
    for (int i = 0 ; i < 10000; i++)
    {
        double dist = i/100.0;
    Eigen::MatrixXd diagonalAfterThing = (dist * eigenValues.array()).exp().matrix().asDiagonal();
    //std::cout<<eigenValues<<std::endl;
    //std::cout<<diagonalAfterThing<<std::endl;

    //Eigen::VectorXd newResult = (eigenVectors * (diagonalAfterThing * (eigenInverseVectors * columnVector)));
    //sum += newResult.sum();
    
    Eigen::VectorXd result = (dist * ok.transposedQ).exp() * columnVector;
    sum1 += result.sum();
    }

    std::cout<<sum<<" "<<sum1<<std::endl;

}
