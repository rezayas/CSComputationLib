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

                // update gamma 
                Matrix<double> mat_x = Matrix<double>.Build.DenseOfColumnArrays(x);
                Matrix<double> a = mat_x.Transpose() * _B * mat_x;
                gamma = discountRate + a[0, 0];

                // update H
                _H = _B.Divide(-gamma);

                // update epsilon
                
                Matrix<double> b = Coeff.ToRowMatrix() * mat_x;
                double epsilon = y - b[0, 0];

                // update coefficients
                Matrix<double> c = _H * mat_x;
                Matrix<double> mat_coeff = Coeff.ToColumnMatrix();
                mat_coeff = mat_coeff - c.Multiply(epsilon);

                // update B
                Matrix<double> d = _B * mat_x * mat_x.Transpose() * _B;
                Matrix<double> e = _B - d.Divide(gamma);

                _B = e.Divide(discountRate);

                // update the coefficient array
                for (int i = 0; i < NumOfColumns; ++i)
                    Coeff[i] = mat_coeff[i, 0];
            }            
        }

        /// <summary>
        /// convert a categorial variable to its equivalent binary code 
        /// </summary>
        /// <param name="valueOfTheCategorialVariable"> must be between 0 and N - 1, inclusive, where N is the number of possible categories </param>
        /// <param name="numOfCategories"> number of possible categories</param>
        /// <returns></returns>
        static public int[] ConvertACategoricalVariableToItsBinaryCode(int valueOfTheCategorialVariable, int numOfCategories)
        {
            int[] result = new int[numOfCategories - 1];

            if (valueOfTheCategorialVariable > 0)
                result[valueOfTheCategorialVariable - 1] = 1;

            return result;
        }
    }
}
