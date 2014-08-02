#ifndef MATRIX_Q_INCLUDED
#define MATRIX_Q_INCLUDED


#include <Eigen/Dense>
#include "R.hpp"
#include <iostream>
#include <complex>

using Eigen::MatrixXd;
using Eigen::MatrixXcd;
using Eigen::VectorXd;
using Eigen::VectorXcd;

class MatrixQ{
public:
    const int M;
    const int matrixSize;

    MatrixXd transposedQ;
    MatrixXcd realEigenVectors;
    MatrixXcd realEigenInverseVectors;
    VectorXcd realEigenValues;


    typedef Eigen::Matrix4d RateModel;

    MatrixQ(const RateModel& foo, int M, double theta): M(M), matrixSize(R::getMatrixSize(M))
    {
        transposedQ = initializeMatrixQ(M,matrixSize,foo,theta);
        initializeEigenVectorsAndValues();
    }

    MatrixQ(const RateModel& rModel, int M): MatrixQ(rModel,M,1)
    {
    }

    //template <typename Derived>
    //VectorXd multiplyByVector(const Eigen::MatrixBase<Derived>& columnVector,double timeScale)
    VectorXd multiplyByVector(double timeScale,const VectorXd& columnVector)
    {
        int size = columnVector.rows();
        auto diagonalAfterThing = (timeScale * realEigenValues.topRows(size).array()).exp().matrix().asDiagonal();
        return (realEigenVectors.topLeftCorner(size,size) * (diagonalAfterThing * (realEigenInverseVectors.topLeftCorner(size,size) * columnVector.cast< std::complex<double> >()))).real();
    }

    
    void initializeEigenVectorsAndValues()
    {
        Eigen::EigenSolver<MatrixXd> solver(transposedQ,true);

        Eigen::MatrixXcd eigenVectors = solver.eigenvectors();
    
        Eigen::VectorXcd eigenValues = solver.eigenvalues();

        /*double sumImaginaryValues = eigenValues.imag().array().abs().sum();
        double sumImaginaryVectors = eigenVectors.imag().array().abs().sum();

        if (sumImaginaryValues != 0 || sumImaginaryVectors != 0)
        {
            std::cout<<"My assumption about imaginary values was wrong. Sorry"<<std::endl;
            abort();
        }

        assert(sumImaginaryValues == 0 && sumImaginaryVectors == 0);
        */
        realEigenVectors = eigenVectors;
        realEigenValues = eigenValues;
        realEigenInverseVectors = realEigenVectors.inverse();
    }


    static MatrixXd initializeMatrixQ(int M,int matrixSize, const RateModel& rateMatrix,  double theta){
        MatrixXd transposedQ(matrixSize, matrixSize);
        transposedQ.setZero();

        for(int n=1; n<=M; n++){
            for(R & r : loopOver(n)){

                for (int fromType = 0; fromType < R::getNumberOfTypes(); fromType++) {
                    if (r.getNum(fromType) != 0)
                    {
                        for (int toType = 0; toType < R::getNumberOfTypes(); toType++) {
                            if (toType == fromType)
                                continue;

                            transposedQ(r.getIndex(),r.transition(fromType, toType).getIndex()) =  r.getNum(fromType) * rateMatrix(toType, fromType);
                        }
                    }
                }

                if (n != M)
                {
                    for (int splitType = 0; splitType < R::getNumberOfTypes(); splitType++) {
                        transposedQ(r.getIndex(),r.split(splitType).getIndex()) = r.getNum(splitType)*(n+1)/(2.0*theta);
                    }
                }

                double sum = -n*(n-1)/(2.0*theta);

                for (int stayType = 0; stayType < R::getNumberOfTypes(); stayType++) {
                    sum += r.getNum(stayType) * rateMatrix(stayType,stayType);
                }

                transposedQ(r.getIndex(), r.getIndex()) = sum;
            }
        }

        return transposedQ;


    }
};

#endif
