#ifndef MATRIX_Q_INCLUDED
#define MATRIX_Q_INCLUDED


#include <Eigen/Dense>
#include "R.hpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class MatrixQ{
public:
    const int M;

    MatrixXd transposedQ;
    MatrixXd realEigenVectors;
    MatrixXd realEigenInverseVectors;
    VectorXd realEigenValues;


    typedef Eigen::Matrix4d RateModel;

    MatrixQ(int M, const RateModel& foo, double theta): M(M)
    {
        transposedQ = initializeMatrixQ(M,foo,theta);
        initializeEigenVectorsAndValues();
    }

    MatrixQ(int M, const RateModel& rModel): MatrixQ(M,rModel,1)
    {
    }

    template <typename Derived>
    VectorXd multiplyByVector(const MatrixBase<Derived>& columnVector,double timeScale)
    {
         Eigen::MatrixXd diagonalAfterThing = (dist * eigenValues.array()).exp().matrix().asDiagonal();
         return (eigenVectors * (diagonalAfterThing * (eigenInverseVectors * columnVector)));
    }

    
    void initializeEigenVectorsAndValues()
    {
        Eigen::EigenSolver<MatrixXd> solver(transposedQ,true);

        Eigen::MatrixXcd eigenVectors = solver.eigenvectors();
    
        Eigen::VectorXcd eigenValues = solver.eigenvalues();

        double sumImaginaryValues = eigenValues.imag().array().abs().sum();
        double sumImaginaryVectors = eigenVectors.imag().array().abs().sum();

        assert(false);
        assert(sumImaginaryValues == 0 && sumImaginaryVectors == 0);

        realEigenVectors = eigenVectors.real();
        realEigenValues = eigenValues.real();
        realEigenInverseVectors = realEigenVectors.inverse();
    }


    VectorXd 

    static MatrixXd initializeMatrixQ(int M,const RateModel& rateMatrix,  double theta){
        MatrixXd transposedQ(R::getMatrixSize(M), R::getMatrixSize(M));
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
