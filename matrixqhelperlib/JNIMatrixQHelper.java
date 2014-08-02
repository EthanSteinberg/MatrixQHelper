package matrixqhelperlib;

import java.util.Arrays;
import java.util.Map;
import java.util.HashMap;

public class JNIMatrixQHelper implements AutoCloseable
{

    final long handle;

    public JNIMatrixQHelper(double[] rawRateMatrix, int M)
    {
        this(rawRateMatrix,M,1.0);
    }

    public JNIMatrixQHelper(double[] rawRateMatrix,int M, double theta)
    {
        handle = allocateMatrixQ(rawRateMatrix,M ,theta);
    }

    public double[] getProbabilityForColumn(double timeScale, double[] columnVector)
    {
        double[] result = new double[columnVector.length];
        multiplyByVectorInC(handle,timeScale,columnVector,result);
        return result;
    }


    @Override
    public void close()
    {
        releaseMatrixQ(handle);
    }

    public static void main(String[] args) {

        double[] gtrVector = {-4.5, 0.5, 0.5000000000000001, 0.5, 
1, -6.0, 1.6666666666666667, 1.5000000000000002 ,
1.5, 2.5, -5.666666666666667, 2.625 ,
2.0, 3.0, 3.5, -4.625};

        JNIMatrixQHelper test = new JNIMatrixQHelper(gtrVector,4);

        double[] columnVector = new double[4];
        columnVector[1] = .1;
        System.out.println(Arrays.toString(test.getProbabilityForColumn(4.0,columnVector)));

    }

    private static native long allocateMatrixQ(double[] rateMatrix, int M, double theta);
    private static native void releaseMatrixQ(long handle);

    private static native void multiplyByVectorInC(long handle, double timeScale, double[] columnVector, double[] result);


    static {
        try
        {
            System.loadLibrary("JNIMatrixQHelper");
        }
        catch (Exception e)
        {
                e.printStackTrace();
        }
    }
}
