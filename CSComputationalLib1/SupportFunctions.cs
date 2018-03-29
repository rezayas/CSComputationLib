using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using RandomVariateLib;

namespace ComputationLib
{
    public static class SupportFunctions
    {
        public static void AddToEndOfArray(ref long[] arrInput, long value)
        {
            int sizeOfNewArray;
            long[] arrResult;

            if (arrInput == null)
            {
                sizeOfNewArray = 1;
                arrResult = new long[sizeOfNewArray];
                arrResult[0] = value;
            }
            else
            {
                long[] arrValue = new long[1];
                arrValue[0] = value;

                arrResult = arrInput.Concat(arrValue).ToArray();

                //sizeOfNewArray = arrInput.Length + 1;
                //arrResult = new long[sizeOfNewArray];
                //long[] arrValue = new long[1];
                //arrValue[0] = value;

                //arrInput.CopyTo(arrResult, 0);
                //arrValue.CopyTo(arrResult, arrInput.Length);
            }

            arrInput = (long[])arrResult.Clone();                     
        }
        public static void AddToEndOfArray(ref int[] arrInput, int value)
        {
            int sizeOfNewArray;
            int[] arrResult;

            if (arrInput == null)
            {
                sizeOfNewArray = 1;
                arrResult = new int[sizeOfNewArray];
                arrResult[0] = value;
            }
            else
            {
                int[] arrValue = new int[1];
                arrValue[0] = value;
                arrResult = arrInput.Concat(arrValue).ToArray();

                //sizeOfNewArray = arrInput.Length + 1;
                //arrResult = new int[sizeOfNewArray];
                //int[] arrValue = new int[1];
                //arrValue[0] = value;

                //arrInput.CopyTo(arrResult, 0);
                //arrValue.CopyTo(arrResult, arrInput.Length);
            }

            arrInput = (int[])arrResult.Clone();            
        }
        public static void AddToEndOfArray(ref double[] arrInput, double value)
        {
            int sizeOfNewArray;
            double[] arrResult;

            if (arrInput == null)
            {
                sizeOfNewArray = 1;
                arrResult = new double[sizeOfNewArray];
                arrResult[0] = value;
            }
            else
            {
                double[] arrValue = new double[1];
                arrValue[0] = value;
                arrResult = arrInput.Concat(arrValue).ToArray();

                //sizeOfNewArray = arrInput.Length + 1;
                //arrResult = new double[sizeOfNewArray];
                //double[] arrValue = new double[1];
                //arrValue[0] = value;

                //arrInput.CopyTo(arrResult, 0);
                //arrValue.CopyTo(arrResult, arrInput.Length);
            }

            arrInput = (double[])arrResult.Clone();  

            //int sizeOfNewArray;

            //if (arrInput == null)
            //    sizeOfNewArray = 1;
            //else
            //    sizeOfNewArray = arrInput.Length + 1;


            //double[] newArray = new double[sizeOfNewArray];

            //for (int i = 0; i < sizeOfNewArray - 1; ++i)
            //{
            //    newArray[i] = arrInput[i];
            //}
            //newArray[sizeOfNewArray - 1] = value;

            //arrInput = newArray;
            
        }
        public static void AddToEndOfArray(ref string[] arrInput, string value)
        {
            int sizeOfNewArray;
            string[] arrResult;

            if (arrInput == null)
            {
                sizeOfNewArray = 1;
                arrResult = new string[sizeOfNewArray];
                arrResult[0] = value;
            }
            else
            {
                string[] arrValue = new string[1];
                arrValue[0] = value;
                arrResult = arrInput.Concat(arrValue).ToArray();

                //sizeOfNewArray = arrInput.Length + 1;
                //arrResult = new string[sizeOfNewArray];
                //string[] arrValue = new string[1];
                //arrValue[0] = value;

                //arrInput.CopyTo(arrResult, 0);
                //arrValue.CopyTo(arrResult, arrInput.Length);
            }

            arrInput = (string[])arrResult.Clone();
        }
        public static void AddToEndOfArrayFixedSize(ref long[] arrInput, long value)
        {
            int size = arrInput.Length;
            for (int i = 0; i < size - 1; ++i)
            {
                arrInput[i] = arrInput[i + 1];
            }
            arrInput[size - 1] = value;
        }
        public static void AddToEndOfArrayFixedSize(ref int[] arrInput, int value)
        {
            int size = arrInput.Length;
            for (int i = 0; i < size - 1; ++i)
            {
                arrInput[i] = arrInput[i + 1];
            }
            arrInput[size - 1] = value;
        }
        public static void AddToEndOfArrayFixedSize(ref double[] arrInput, double value)
        {
            int size = arrInput.Length;
            for (int i = 0; i < size - 1; ++i)
            {
                arrInput[i] = arrInput[i + 1];
            }
            arrInput[size - 1] = value;
        }
        public static double[][] ConcatJaggedArray(double[][] array1, double[][] array2)
        {
            return array1.Concat(array2).ToArray();

            //double[][] result = new double[array1.Length + array2.Length][];

            //array1.CopyTo(result, 0);
            //array2.CopyTo(result, array1.Length);

        }
        public static int[][] ConcatJaggedArray(int[][] array1, int[][] array2)
        {
            return array1.Concat(array2).ToArray();

            //int[][] result = new int[array1.Length + array2.Length][];

            //array1.CopyTo(result, 0);
            //array2.CopyTo(result, array1.Length);

            //return result;
        }
        public static string[][] ConcatJaggedArray(string[][] array1, string[][] array2)
        {
            return array1.Concat(array2).ToArray();

            //string[][] result = new string[array1.Length + array2.Length][];

            //array1.CopyTo(result, 0);
            //array2.CopyTo(result, array1.Length);

            //return result;
        }
        public static double[][] ConcatJaggedArray(double[][] array1, double[] array2)
        {
            double[][] tempArray2 = new double[1][];
            tempArray2[0] = new double[array2.Length];

            for (int i = 0; i < array2.Length; i++)
                tempArray2[0][i] = array2[i];

            return ConcatJaggedArray(array1, tempArray2);
        }
        public static string[][] ConcatJaggedArray(string[][] array1, string[] array2)
        {
            string[][] tempArray2 = new string[1][];
            tempArray2[0] = new string[array2.Length];

            for (int i = 0; i < array2.Length; i++)
                tempArray2[0][i] = array2[i];

            return ConcatJaggedArray(array1, tempArray2);
        }
        public static int[][] ConcatJaggedArray(int[][] array1, int[] array2)
        {
            int[][] tempArray2 = new int[1][];
            tempArray2[0] = new int[array2.Length];

            for (int i = 0; i < array2.Length; i++)
                tempArray2[0][i] = array2[i];

            return ConcatJaggedArray(array1, tempArray2);
        }
        public static double[][] RemoveARowFromJaggedArray(double[][] array, int rowIndex)
        {
            double[][] result;
            result = array.Where((el, i) => i != rowIndex).ToArray();
            return result;
        }        

        public static void MakeArrayEqualTo(ref double[] array, double value)
        {            
            for (int i = 0; i < array.Length; ++i)
                array[i] = value;
        }
        public static void MakeArrayEqualTo(ref long[] array, long value)
        {
            for (int i = 0; i < array.Length; ++i)
                array[i] = value;
        }
        public static void MakeArrayEqualTo(ref int[] array, int value)
        {
            for (int i = 0; i < array.Length; ++i)
                array[i] = value;
        }
        public static void MakeMatrixEqualTo(ref double[,] matrix, double value)
        {
            for (int i = 0; i < matrix.GetLength(0); i++)
                for (int j = 0; j < matrix.GetLength(1); j++)
                    matrix[i, j] = value;
        }
        public static bool ConvertYesNoToBool(string input)
        {
            bool result = false;
            switch (input)
            {
                case "No":
                    result = false;
                    break;
                case "Yes":
                    result = true;
                    break;
            }
            return result;
        }

        public static string Left(string text, int length)
        {
            if (length < 0)
                throw new ArgumentOutOfRangeException("length", length, "length must be > 0");
            else if (length == 0 || text.Length == 0)
                return "";
            else if (text.Length <= length)
                return text;
            else
                return text.Substring(0, length);
        }

        public static long[] ConvertArrayToLong(double[] array)
        {
            return Array.ConvertAll<double, long>(array, Convert.ToInt64);
        }
        public static long[] ConvertArrayToLong(int[] array)
        {
            return Array.ConvertAll<int, long>(array, Convert.ToInt64);
        }
        public static int[] ConvertArrayToInt(double[] array)
        {
            return Array.ConvertAll<double, int>(array, Convert.ToInt32);
        }
        public static double[] ConvertArrayToDouble(long[] array)
        {
            return Array.ConvertAll<long, double>(array, Convert.ToDouble);
        }
        public static double[] ConvertArrayToDouble(int[] array)
        {
            return Array.ConvertAll<int, double>(array, Convert.ToDouble);
        }
        public static string[] ConvertArrayToString(double[] array)
        {
            return Array.ConvertAll<double, string>(array, Convert.ToString);
        }
        public static int[] ConvertArrayToInt(string[] array)
        {
            return Array.ConvertAll<string, int>(array, Convert.ToInt32);
        }
        public static int[,] ConvertMatrixToInt(double[,] matrix)
        {
            int[,] newMatrix = new int[matrix.GetLength(0), matrix.GetLength(1)];

            for (int i = 0; i < matrix.GetLength(0); ++i)
                for (int j = 0; j < matrix.GetLength(1); ++j)
                    newMatrix[i,j] = (int)matrix[i, j];

            return (int[,])newMatrix.Clone();
        }

        public static double[,] ConvertFromJaggedArrayToRegularArray(double[][] jaggedArray, int columnSize)
        {
            double[,] result = new double[jaggedArray.GetLength(0), columnSize];
            for (int i = 0; i < jaggedArray.GetLength(0); i++)
                for (int j = 0; j < columnSize; j++)
                    result[i, j] = jaggedArray[i][j];
            return result;
        }
        public static int[,] ConvertFromJaggedArrayToRegularArray(int[][] jaggedArray, int columnSize)
        {
            int[,] result = new int[jaggedArray.GetLength(0), columnSize];
            for (int i = 0; i < jaggedArray.GetLength(0); i++)
                for (int j = 0; j < columnSize; j++)
                    result[i, j] = jaggedArray[i][j];
            return result;
        }
        public static string[,] ConvertFromJaggedArrayToRegularArray(string[][] jaggedArray, int columnSize)
        {
            string[,] result = new string[jaggedArray.GetLength(0), columnSize];
            for (int i = 0; i < jaggedArray.GetLength(0); i++)
                for (int j = 0; j < columnSize; j++)
                    result[i, j] = jaggedArray[i][j];
            return result;
        }

        // convert to base 10 from base 2
        public static int ConvertToBase10FromBase2(int[] arrOfZerosAndOnes)
        {
            int result = 0;
            int numOfBytes = arrOfZerosAndOnes.Length;

            for (int i = 0; i < numOfBytes; ++i)
                result += arrOfZerosAndOnes[i] * (int)Math.Pow(2, numOfBytes - 1 - i);

            return result;
        }
        // convert to base 2 from base 10
        public static int[] ConvertToBase2FromBase10(int number, int numOfBase2Digits)
        {
            int[] tempResult = new int[0];
            int[] result = new int[numOfBase2Digits];
            int remainder = 0;
            int quotient = 0;

            quotient = number;            
            while (quotient > 1)
            {
                quotient = Math.DivRem(quotient, 2, out remainder);
                SupportFunctions.AddToEndOfArray(ref tempResult, remainder);            
            }   
            SupportFunctions.AddToEndOfArray(ref tempResult, quotient);

            int resultIndex = numOfBase2Digits - 1;
            for (int i = 0; i < tempResult.Length; ++i)
                result[resultIndex--] = tempResult[i];                

            return (int[])result.Clone();
        }

        // covert an array to string
        public static string ConvertArrayToString(int[] arr, string splitCharacter)
        {
            return String.Join(splitCharacter, arr.Select(p => p.ToString()).ToArray());
        }
        public static string ConvertArrayToString(long[] arr, string splitCharacter)
        {
            return String.Join(splitCharacter, arr.Select(p => p.ToString()).ToArray());
        }
        public static string ConvertArrayToString(double[] arr, string splitCharacter)
        {
            return String.Join(splitCharacter, arr.Select(p => p.ToString()).ToArray());
        }

        public static int[] ConvertStringToIntArray(string str, char splitCharacter)
        {
            string[] strArray = str.Split(splitCharacter);
            return Array.ConvertAll<string, int>(strArray, Convert.ToInt32);
        }
        public static long[] ConvertStringToLongArray(string str, char splitCharacter)
        {
            string[] strArray = str.Split(splitCharacter);
            return Array.ConvertAll<string, long>(strArray, Convert.ToInt64);
        }

        // suffle an array 
        public static void Shuffle<T>(RNG rng, T[] array)
        {
            int n = array.Length;
            while (n > 1)
            {
                int k = rng.Next(n--);
                T temp = array[n];
                array[n] = array[k];
                array[k] = temp;
            }
        }

        //public static void Shuffle(ThreadSpecificRNG rng, ref double[] array)
        //{
        //    int n = array.Length;
        //    while (n > 1)
        //    {
        //        int k = 1;// rng.Next(n--);
        //        double temp = array[n];
        //        array[n] = array[k];
        //        array[k] = temp;
        //    }
        //}
    }
}
