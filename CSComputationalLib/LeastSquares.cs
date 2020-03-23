using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MatrixLibrary;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearRegression;

namespace ComputationLib
{
    public class LeastSquares
    {

        public int NumOfColumns { get; private set; }
        public Vector<double> Coeff { get; private set; }

        private double _l2Penalty;
        private Vector<double> _y;
        private Matrix<double> _X;
        private Matrix<double> _XTX, _H, _B;

        private int _itr, _minObsToInitializeTraining;


        public LeastSquares(double l2Penalty = 0)
        {
            _l2Penalty = l2Penalty;
        }

        // reset
        public void Reset()
        {
            _itr = 0;
            _X = Matrix<double>.Build.Dense(_minObsToInitializeTraining, NumOfColumns);
            _XTX = Matrix<double>.Build.Dense(NumOfColumns, NumOfColumns);
            _y = Vector<double>.Build.Dense(_minObsToInitializeTraining);
            _H = Matrix<double>.Build.Dense(NumOfColumns, NumOfColumns);
            _B = Matrix<double>.Build.Dense(NumOfColumns, NumOfColumns);
            Coeff = Vector<double>.Build.Dense(NumOfColumns);

        }

        // update coefficient
        public void UpdateCoefficients(double[] coefficients)
        {
            Coeff = Vector<double>.Build.Dense(coefficients);
        }

        // add an L2-Regularization
        public void AddL2Regularization(double penalty)
        {
            _l2Penalty = penalty;
        }

        public void RunRegression(double[,] X, double[] y)
        {
            // coeff = (XT.X)-1.XT.Y
            int numOfObs = y.Length;
            NumOfColumns = X.GetLength(1);

            _X = Matrix<double>.Build.DenseOfArray(X);
            _y = Vector<double>.Build.DenseOfArray(y);
            _XTX = _X.Transpose() * _X;

            // add L2 regularization penalty
            if (_l2Penalty > 0)
                AddL2ToXTX();

            Coeff = _XTX.Cholesky().Solve(_X.TransposeThisAndMultiply(_y));
        }

        private void AddL2ToXTX()
        {
            Matrix<double> I = Matrix<double>.Build.DiagonalIdentity(NumOfColumns);
            I[0, 0] = 0; // no penalty for the intercept
            _XTX += I.Multiply(_l2Penalty);
        }

        public double yValue(double[] x)
        {
            double sum = 0;
            for (int i = 0; i < x.Length; ++i)
                sum += Coeff[i] * x[i];
            return sum;
        }

        public void SetupTraining(int numOfColumns)
        {
            NumOfColumns = numOfColumns;
            _minObsToInitializeTraining = NumOfColumns;
            Reset();
        }

        //Page 249 from Powel (2007) - Page 350 from Powel (2011)
        public void Update(double[] x, double y, double discountRate)
        {
            double gamma;

            _itr += 1;

            if (_itr < _minObsToInitializeTraining)
            {
                // add x to the last line of _arrX0
                for (int i = 0; i < NumOfColumns; i++)
                    _X[_itr - 1, i] = x[i];
                _y[(int)_itr - 1] = y;

            }
            else if (_itr == _minObsToInitializeTraining)
            {
                // add x to the last line of _arrX0
                for (int i = 0; i < NumOfColumns; i++)
                    _X[_itr - 1, i] = x[i];
                _y[(int)_itr - 1] = y;

                // XTX
                _XTX = _X.Transpose() * _X;

                // add L2 regularization penalty
                if (_l2Penalty > 0)
                    AddL2ToXTX();

                // calculate B                
                _B = _XTX.Inverse();
                // calculate coefficient
                Coeff = _B * _X.Transpose() * _y;
            }
            else
            {
                // create a matrix from array x
                Matrix mat_x = new Matrix(_numOfColumns, 1);
                for (int i = 0; i < _numOfColumns; i++)
                    mat_x[i, 0] = x[i];

                // update gamma
                tempMatrix = null;
                tempMatrix = Matrix.Transpose(mat_x) * _matB * mat_x;
                if (double.IsNaN(tempMatrix[0, 0]) || double.IsInfinity(tempMatrix[0, 0]))
                {
                    Console.WriteLine("Error in OLS updating: gamma is either not a number or infinity.");
                    // make the matrix X identity matrix
                    for (int i = 0; i < _numOfColumns; ++i)
                        for (int j = 0; j < _numOfColumns; ++j)
                            if (i == j)
                                _matB[i, j] = 1;
                            else
                                _matB[i, j] = 0;
                    tempMatrix = Matrix.Transpose(mat_x) * _matB * mat_x;
                }
                gamma = discountRate + tempMatrix[0, 0];

                // adjust gamma for stability if needed
                if (Math.Abs(tempMatrix[0, 0]) <= 0.00001)
                    gamma += 0.00001;

                // update H
                _matH = Matrix.ScalarDivide(-gamma, _matB);

                // update epsilon
                tempMatrix = null;
                tempMatrix = Matrix.Transpose(_matCoeff) * mat_x;
                double epsilon = y - tempMatrix[0, 0];

                // update coefficients
                tempMatrix = null;
                tempMatrix = _matH * mat_x;
                _matCoeff = _matCoeff - Matrix.ScalarMultiply(epsilon, tempMatrix);

                // update B
                tempMatrix = null;
                tempMatrix = _matB * mat_x * Matrix.Transpose(mat_x) * _matB;

                tempMatrix2 = null;
                //tempMatrix2 = Matrix.ScalarDivide(_itrNumber,_matB) - Matrix.ScalarDivide(gamma, tempMatrix);
                tempMatrix2 = _matB - Matrix.ScalarDivide(gamma, tempMatrix);

                _matB = Matrix.ScalarDivide(discountRate, tempMatrix2);
            }

            // update the coefficient array
            for (int i = 0; i < _numOfColumns; ++i)
                _arrCoefficients[i] = _matCoeff[i, 0];
        }

            ///// <summary>
            ///// convert a categorial variable to its equivalent binary code 
            ///// </summary>
            ///// <param name="valueOfTheCategorialVariable"> must be between 0 and N - 1, inclusive, where N is the number of possible categories </param>
            ///// <param name="numOfCategories"> number of possible categories</param>
            ///// <returns></returns>
            //static public int[] ConvertACategoricalVariableToItsBinaryCode(int valueOfTheCategorialVariable, int numOfCategories)
            //{
            //    int[] result = new int[numOfCategories - 1];

            //    if (valueOfTheCategorialVariable > 0)
            //        result[valueOfTheCategorialVariable - 1] = 1;

            //    return result;
            //}
        }
}
