using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Windows.Forms;
using MatrixLibrary;

namespace ComputationLib
{
    public class LeastSquares
    {
        private string _name;


        double _L2PenaltyParameter;
        // general        
        Matrix _matX, _matY, _matCoeff;
        double[] _arrCoefficients;
        
        // training
        private int _numOfColumns;
        private int _itrNumber, _minObsToInitializeTraining;
        private Matrix _matXTX, _matH, _matB;
        //private Random _RNG;

        public LeastSquares(string name)
        {
            _name = name;
        }
        // return number of columns
        public int numOfColumns
        {
            get {return _numOfColumns;}
        }
        // return coefficients
        public double[] Coefficients
        {
            get { return _arrCoefficients; }
        }
        // reset
        public void Reset()
        {
            _itrNumber = 0;

            _matX = new Matrix(_minObsToInitializeTraining, _numOfColumns);
            _matXTX = new Matrix(_numOfColumns, _numOfColumns);
            _matY = new Matrix(_minObsToInitializeTraining, 1);
            _matH = new Matrix(_numOfColumns, _numOfColumns);
            _matB = new Matrix(_numOfColumns, _numOfColumns);
            _matCoeff = new Matrix(_numOfColumns, 1);
        }
        
        // update coefficient
        public void UpdateCoefficients(double[] arrCoefficients)
        {
            _arrCoefficients = (double[])arrCoefficients.Clone();
            //for (int i = 0; i < _numOfColumns; ++i)
               // _matCoeff[i, 0] = arrCoefficients[i];
        }
        // add an L2-Regularization
        public void AddL2Regularization(double penaltyParameter)
        {
            _L2PenaltyParameter = penaltyParameter;
        }
        public void RunRegression(double[,] X, double[] Y)
        {
            // coeff = (XT.X)-1.XT.Y
            int numOfObs = Y.Length;
            _numOfColumns = X.GetLength(1);

            double[,] Y2 = new double[numOfObs,1];

            // create a 2 dimential array from Y
            for (int i = 0; i < numOfObs; i++)
                Y2[i, 0] = Y[i];
            
            Matrix _matX = new Matrix(X);
            Matrix _matY = new Matrix(Y2);
            _matXTX = Matrix.Transpose(_matX) * _matX;

            // add L2 regularization penalty
            if (_L2PenaltyParameter > 0)
            {
                Matrix I = new Matrix(Matrix.Identity(_numOfColumns));
                I[0, 0] = 0; // no penalty for the intercept
                _matXTX = _matXTX + Matrix.ScalarMultiply(_L2PenaltyParameter, I);
            }

            // calculate coefficients                
            Matrix B = Matrix.Inverse(_matXTX);
            _matCoeff = B * Matrix.Transpose(_matX) * _matY;

            // update the coefficient array
            _arrCoefficients = new double[_numOfColumns];
            for (int i = 0; i < _numOfColumns; ++i)
                _arrCoefficients[i] = _matCoeff[i, 0];
        }

        public double yValue(double[] x)
        {
            double sum = 0;
            for (int i = 0; i < x.Length; ++i)
                sum += _arrCoefficients[i] * x[i];
            return sum;
        }

        public void SetupTraining(int numOfColumns, int multiplyNumOfColumnsByThisFactorToBeginTraining = 1)//, double perturbationRatio)
        {
            //_RNG = new Random(2);
            _numOfColumns = numOfColumns;
            _minObsToInitializeTraining = Math.Max(_numOfColumns, _numOfColumns * multiplyNumOfColumnsByThisFactorToBeginTraining);
            _itrNumber = 0;

            _matX = new Matrix(_minObsToInitializeTraining, numOfColumns);
            _matXTX = new Matrix(_numOfColumns, _numOfColumns);
            _matY = new Matrix(_minObsToInitializeTraining, 1);
            _matH = new Matrix(numOfColumns, numOfColumns);
            _matB = new Matrix(numOfColumns, numOfColumns);
            _matCoeff = new Matrix(numOfColumns, 1);

            _arrCoefficients = new double[numOfColumns];
        }

        //Page 249 from Powel (2007) - Page 350 from Powel (2011)
        public void Update(double[] x, double y, double observationDiscountRate)
        {
            Matrix tempMatrix, tempMatrix2;
            double gamma;

            // increatement the iteration number
            _itrNumber += 1;

            if (_itrNumber < _minObsToInitializeTraining)
            {
                // add x to the last line of _arrX0
                for (int i = 0; i < _numOfColumns; i++)
                    _matX[_itrNumber - 1, i] = x[i];
                _matY[(int)_itrNumber - 1, 0] = y;

            }
            else if (_itrNumber == _minObsToInitializeTraining)
            {
                // add x to the last line of _arrX0
                for (int i = 0; i < _numOfColumns; i++)
                    _matX[_itrNumber - 1, i] = x[i];
                _matY[(int)_itrNumber - 1, 0] = y;

                // XTX
                _matXTX = Matrix.Transpose(_matX) * _matX;

                // add L2 regularization penalty
                if (_L2PenaltyParameter > 0)
                {
                    Matrix I = new Matrix(Matrix.Identity(_numOfColumns));
                    I[0, 0] = 0; // no penalty for the intercept
                    _matXTX = _matXTX + Matrix.ScalarMultiply(_L2PenaltyParameter, I);
                }

                // check if the rank of _X is zero       
                double absDet = Math.Abs(Matrix.Det(_matXTX));

               //// double epsilon;                       
                while (absDet < double.Epsilon)
                {
                    // make matrix XTX identity matrix
                    for (int i = 0; i < _numOfColumns; ++i)
                        for (int j = 0; j < _numOfColumns; ++j)
                            if (i == j)
                                _matXTX[i, j] = 1;
                            else
                                _matXTX[i, j] = 0;

                    absDet = Math.Abs(Matrix.Det(_matXTX));
                }
                
                // calculate B                
                _matB = Matrix.Inverse(_matXTX);                
                // calculate coefficient
                _matCoeff = _matB * Matrix.Transpose(_matX) * _matY;
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
                    MessageBox.Show("Error in OLS updating: gamma is either not a number or infinity.");
                    // make the matrix X identity matrix
                    for (int i = 0; i < _numOfColumns; ++i)
                        for (int j = 0; j < _numOfColumns; ++j)
                            if (i == j)
                                _matB[i, j] = 1;
                            else
                                _matB[i, j] = 0;
                    tempMatrix = Matrix.Transpose(mat_x) * _matB * mat_x;
                }
                gamma = observationDiscountRate + tempMatrix[0, 0];

                // adjust gamma for stability if needed
                if (Math.Abs(tempMatrix[0, 0]) <= 0.00001)
                    gamma += 0.00001;

                //// check if gamma is close to zero
                //if (Math.Abs(gamma) < double.Epsilon ) //_minDeterminant
                //{
                //    //MessageBox.Show("Gamma close to 0 occurs!");
                //    gamma += double.Epsilon;//_perturbationRatio
                //}
                //if (double.IsInfinity(gamma) || double.IsNaN(gamma))
                //    gamma = double.MaxValue;

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

                _matB = Matrix.ScalarDivide(observationDiscountRate, tempMatrix2);
            }

            // update the coefficient array
            for (int i = 0; i < _numOfColumns; ++i)
                _arrCoefficients[i] = _matCoeff[i, 0];
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
