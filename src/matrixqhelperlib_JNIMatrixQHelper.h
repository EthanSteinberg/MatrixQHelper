/* DO NOT EDIT THIS FILE - it is machine generated */
#include <jni.h>
/* Header for class matrixqhelperlib_JNIMatrixQHelper */

#ifndef _Included_matrixqhelperlib_JNIMatrixQHelper
#define _Included_matrixqhelperlib_JNIMatrixQHelper
#ifdef __cplusplus
extern "C" {
#endif
/*
 * Class:     matrixqhelperlib_JNIMatrixQHelper
 * Method:    allocateMatrixQ
 * Signature: ([DID)J
 */
JNIEXPORT jlong JNICALL Java_matrixqhelperlib_JNIMatrixQHelper_allocateMatrixQ
  (JNIEnv *, jclass, jdoubleArray, jint, jdouble);

/*
 * Class:     matrixqhelperlib_JNIMatrixQHelper
 * Method:    releaseMatrixQ
 * Signature: (J)V
 */
JNIEXPORT void JNICALL Java_matrixqhelperlib_JNIMatrixQHelper_releaseMatrixQ
  (JNIEnv *, jclass, jlong);

/*
 * Class:     matrixqhelperlib_JNIMatrixQHelper
 * Method:    multiplyByVectorInC
 * Signature: (JD[D[D)V
 */
JNIEXPORT void JNICALL Java_matrixqhelperlib_JNIMatrixQHelper_multiplyByVectorInC
  (JNIEnv *, jclass, jlong, jdouble, jdoubleArray, jdoubleArray);

#ifdef __cplusplus
}
#endif
#endif