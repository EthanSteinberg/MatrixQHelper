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


    MatrixQ ok(blah,4);
    Eigen::VectorXd columnVector = Eigen::VectorXd::Random(ok.transposedQ.cols()).array().abs().matrix();
    
    for (int i = 0 ; i < 10000; i++)
    {
        Eigen::VectorXd result = ok.multiplyByVector(i/100.0,columnVector);
    }
}
