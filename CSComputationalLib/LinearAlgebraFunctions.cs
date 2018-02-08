using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MatrixLibrary;

namespace ComputationLib
{
    public static class LinearAlgebraFunctions
    {
        public enum DistanceMeasure
        {
            SquaredDifference = 1,
        }
        public enum enumVectorNorm
        {
            L_1 = 0,
            L_2 = 1,
            L_inf = 2,
        }
        public enum enumMatrixNorm
        {
            L_1 = 0,
            L_inf = 1,
            Frobenius = 3, 
        }

        // normalizations
        public static double[] NormalizeToSum1(double[] data)
        {
            double tot = data.Sum();
            double[] result = new double[data.Length];
            for (int i = 0; i < data.Length; i++)
                result[i] = data[i] / tot;

            return result;
        }
        
        public static double[] NormalizeToAnInterval(double[] data, double lBound, double uBound)
        {
            int count = data.Length;
            double[] normalizedData;

            // normalize to a verctor of size 1
            normalizedData = NormalizeToVectorSize1(data);

            double a = (uBound - lBound) / 2;
            double b = (uBound + lBound) / 2;

            for (int i = 0; i < count; ++i)
                normalizedData[i] = a * normalizedData[i] + b;

            return normalizedData;
        }
        public static double[] NormalizeToVectorSize1(double[] data)
        {
            int count = data.Length;
            double[] normalizedData = new double[count];
            double sum2 = 0;
            for (int i = 0; i < count; ++i)
            {
                sum2 += Math.Pow(data[i], 2);
            }
            if (sum2 == 0) return normalizedData;

            double denom = Math.Pow(sum2, 0.5);
            for (int i = 0; i < count; ++i)
            {
                normalizedData[i] = data[i] / denom;
            }

            return normalizedData;
        }
        public static double[,] NormalizeEachColumnToAnInterval(double[,] matrix, double lBound, double uBound)
        {            
            int numOfRows = matrix.GetLength(0);
            int numOfCols = matrix.GetLength(1);
            double[,] normalizedMatrix = new double[numOfRows, numOfCols];

            return normalizedMatrix;
        }
        // dot product of two vectors
        public static double DotProduct(double[] array1, double[] array2)
        {
            double sum = 0;
            for (int i = 0; i < array1.Length; i++)
                sum += array1[i] * array2[i];

            return sum;
        }

        // vector norm
        public static double Norm(double[] array, enumVectorNorm norm)
        {
            int arraySize = array.Length;
            double result = 0;

            switch (norm)
            {
                case enumVectorNorm.L_1:
                    {
                        for (int i = 0; i < arraySize; i++)
                            result += Math.Abs(array[i]);
                    }
                    break;
                case enumVectorNorm.L_2:
                    {
                        for (int i = 0; i < arraySize; i++)
                            result += Math.Pow(array[i], 2);
                        result = Math.Sqrt(result);
                    }
                    break;
                case enumVectorNorm.L_inf:
                    {
                        result = double.MinValue;
                        for (int i = 0; i < arraySize; i++)
                            if (Math.Abs(array[i]) > result)
                                result = Math.Abs(array[i]);
                    }
                    break;
            }
            return result;
        }
        // matrix norm
        public static double Norm(double[,] matrix, enumMatrixNorm norm)
        {
            int numOfRows = matrix.GetLength(0);
            int numOfCols = matrix.GetLength(1);
            double[,] weights = new double[numOfRows, numOfCols];
            SupportFunctions.MakeMatrixEqualTo(ref weights, 1);

            return Norm(matrix, weights, norm);
        }
        // weighted matrix norm
        public static double Norm(double[,] matrix, double[,] weights, enumMatrixNorm norm)
        {
            double result = 0;
            int numOfRows = matrix.GetLength(0);
            int numOfCols = matrix.GetLength(1);

            switch (norm)
            {
                case enumMatrixNorm.L_1:
                    {                        
                        // go over columns
                        for (int j = 0; j < numOfCols; j++)
                        {
                            // find the sum of column j
                            double sum = 0;
                            for (int i = 0; i < numOfRows; i++)
                                sum += weights[i, j] * Math.Abs(matrix[i, j]);
                            // find the max
                            if (sum > result)
                                result = sum;
                        }
                    }
                    break;
                case enumMatrixNorm.L_inf:
                    {
                        // go over rows
                        for (int i = 0; i < numOfRows; i++)
                        {
                            // find the sum of row i
                            double sum = 0;
                            for (int j = 0; i < numOfCols; j++)
                                sum += weights[i, j] * Math.Abs(matrix[i, j]);
                            // find the max
                            if (sum > result)
                                result = sum;
                        }                        
                    }
                    break;                
                case enumMatrixNorm.Frobenius:
                    {
                        double totWeight = 0;
                        for (int i = 0; i < numOfRows; i++)
                            for (int j = 0; j < numOfCols; j++)
                            {
                                if (weights[i, j] > 0 && !double.IsNaN(matrix[i, j]))
                                {
                                    result += weights[i, j] * Math.Pow(matrix[i, j], 2);
                                    totWeight += weights[i, j];
                                }
                            }
                        result = Math.Sqrt(result / totWeight);
                    }
                    break;
            }

            return result;
        }        

        // cosine of the angel between two vectors
        public static double CosineOfAngel(double[] vector1, double[] vector2)
        {
            double nominator = DotProduct(vector1, vector2);
            double denominator = Norm(vector1, enumVectorNorm.L_2) * Norm(vector2, enumVectorNorm.L_2);

            return nominator / denominator;
        }

        // matrix is invertable?
        public static bool Matrix_IfInvertable(double[,] matrix)
        {
            if (Matrix.Det(matrix) == 0)
                return false;

            return true;
        }
        // matrix inverse
        public static double[,] Matrix_Inverse(double[,] matrix)
        {
            return Matrix.Inverse(matrix);            
        }
        // matrix multiply by vector
        public static double[] Matrix_Multiply(double[,] matrix, double[] vector)
        {
            int vectorSize = vector.Length;
            double[] result = new double[vectorSize];
            double[,] V = new double[vectorSize, 1];

            // create a 2 dimential array from vector
            for (int i = 0; i < vectorSize; i++)
                V[i, 0] = vector[i];

            // calculate the multiplication
            Matrix _matX = new Matrix(matrix);
            Matrix _matV = new Matrix(V);
            Matrix multiply = _matX * _matV;

            // the result
            for (int i = 0; i < vectorSize; ++i)                
                result[i] = multiply[i, 0];

            return result;
        }
    }
}
