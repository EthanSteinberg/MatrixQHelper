#include "matrixqhelperlib_JNIMatrixQHelper.h"
#include "MatrixQ.hpp"
#include <iostream>


template<typename T>
struct JArrayHelper{
};

template<>
struct JArrayHelper<double>{
  typedef jdoubleArray type;
};


template<>
struct JArrayHelper<int>{
  typedef jintArray type;
};

template <typename T>
struct JConverter {
public:
  JNIEnv* env;
  typename JArrayHelper<T>::type arr;
  T* actualArr;
  JConverter(JNIEnv *a_env, typename JArrayHelper<T>::type a_arr) :
    env(a_env),
    arr(a_arr),
    actualArr((T*)env->GetPrimitiveArrayCritical(arr, NULL))
  {
  }

  operator T*()
  {
    return actualArr;
  }

  ~JConverter()
  {
    env->ReleasePrimitiveArrayCritical(arr, actualArr, 0);
  }
};

typedef JConverter<double>  JDoubleConverter;
typedef JConverter<int>  JIntConverter;

static_assert(sizeof(jlong) >= sizeof(MatrixQ*),"The handle needs to be the right size to hold the pointers");

JNIEXPORT jlong JNICALL Java_matrixqhelperlib_JNIMatrixQHelper_allocateMatrixQ
  (JNIEnv *env, jclass myself, jdoubleArray rateMatrix, jint M, jdouble theta)
{
    
    JDoubleConverter converter(env,rateMatrix);
    Eigen::Map<Eigen::Matrix4d> actualRateMatrix(converter);

    MatrixQ* result = new MatrixQ(actualRateMatrix.transpose(),M,theta);
    return reinterpret_cast<jlong>(result);
}

JNIEXPORT void JNICALL Java_matrixqhelperlib_JNIMatrixQHelper_releaseMatrixQ
  (JNIEnv *env, jclass myself, jlong handle)
{
    MatrixQ* value = reinterpret_cast<MatrixQ*>(handle);
    delete value;
}

JNIEXPORT void JNICALL Java_matrixqhelperlib_JNIMatrixQHelper_multiplyByVectorInC
  (JNIEnv *env, jclass myself, jlong handle, jdouble timeScale, jdoubleArray inputVector, jdoubleArray outputVector)
{
    MatrixQ* value = reinterpret_cast<MatrixQ*>(handle);

    int inputLength = env->GetArrayLength(inputVector);
    int outputLength = env->GetArrayLength(outputVector);

    if (inputLength != outputLength)
    {
        std::cout<<"The sizes of the passed in vectors were not equal."<<std::endl;
        abort();
    }

    JDoubleConverter convertedInputVector(env,inputVector);
    Eigen::Map<Eigen::VectorXd> mappedInputVector(convertedInputVector,inputLength);

    JDoubleConverter convertedOutputVector(env,outputVector);
    Eigen::Map<Eigen::VectorXd> mappedOutputVector(convertedOutputVector,inputLength);


    mappedOutputVector = value->multiplyByVector(timeScale,mappedInputVector);
}


