using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Windows.Forms;

namespace ComputationLib
{
    // One Dimention Linear Splines class
    public class LinearSplines1Dimension
    {
        private string _name;
        private enumBracketType _bracketType;
        private double _minX, _maxX;
        private double[] _breakPoints;
        private double[] _fBreakPoints;
        private int[] _numOfObservables;
        private double _bracketSize;
        private int _numOfBrackets;
        private LeastSquares _OLS;
        private double _observationDiscountRate;
        private double _perturbationRatio;

        private enum enumBracketType : int
        {
            FixedBracketSize = 0,
            VariableBracketSize = 1,
        }

        public LinearSplines1Dimension(string name)
        {
            _name = name;
        }

        public int NumOfBreakPoints
        {
            get
            {
                return _breakPoints.Length;
            }
        }
        public double[] fBreakPoints
        {
            get
            {
                return _fBreakPoints;
            }
        }
        public double MaxBreakPoint
        {
            get
            {
                return _maxX;
            }
        }

        public void Reset()
        {
            // reset OLS  
            SetupTraining(_observationDiscountRate, _perturbationRatio);

            for (int i = 0; i < _numOfBrackets + 1; ++i)
            {
                _fBreakPoints[i] = 0;
            }
            for (int i = 0; i < _numOfBrackets; ++i)
            {
                _numOfObservables[i] = 0;
            }
        }
        public void AssignBreakPoints(double minX, double maxX, double bracketSize)
        {
            if (Math.Truncate((maxX - minX) / bracketSize) != (maxX - minX) / bracketSize)
            {
                MessageBox.Show("Error in the choice of bracket size!");
            }            

            _bracketType = enumBracketType.FixedBracketSize;
            _minX = minX;
            _maxX = maxX;

            _bracketSize = bracketSize;
            // find the number of brackets
            _numOfBrackets = (int)((_maxX - _minX) / _bracketSize);
            // redim arrays
            _breakPoints = new double[_numOfBrackets + 1];
            _fBreakPoints = new double[_numOfBrackets + 1];
            _numOfObservables = new int[_numOfBrackets];

            // find the break points
            double x = minX;
            int i = 0;
            while (x <= maxX)
            {
                _breakPoints[i] = x;
                _fBreakPoints[i] = 0;
                i += 1;
                x += _bracketSize;
            }
        }
        public void AssignBreakPoints(double[] breakPoints)
        {
            _bracketType = enumBracketType.VariableBracketSize;
            _breakPoints = breakPoints;            
            _numOfBrackets = breakPoints.Length - 1;
            _numOfObservables = new int[_numOfBrackets];
            _fBreakPoints = new double[_numOfBrackets + 1];
            _minX = breakPoints[0];
            _maxX = breakPoints[_numOfBrackets];
        }        

        public void SetupTraining(double observationDiscountRate, double perturbationRatio)
        {
            _observationDiscountRate = observationDiscountRate;
            _perturbationRatio = perturbationRatio;

            _OLS = new LeastSquares(_name);

            _OLS.SetupTraining(_numOfBrackets + 1);
        }

        public void Update(double x, double fValue)
        {
            double[] OLS_x = new double[_numOfBrackets+1];

            // find the index of the lower break point
            int bracketIndex = BracketIndex(x);
            
            // increment the number of observation in this bracket
            _numOfObservables[bracketIndex] += 1;

            // find lambda
            double lambda = Lambda(x, bracketIndex);
            
            OLS_x[bracketIndex]= 1 - lambda;
            OLS_x[bracketIndex + 1] = lambda;
            // update OLS
            _OLS.Update(OLS_x, fValue, _observationDiscountRate);
            // update value of break points
            _fBreakPoints = _OLS.Coefficients;
        }
        public double f(double x)
        {
            // find the index of the lower break point
            int bracketIndex = BracketIndex(x);

            // find lambda
            double lambda = Lambda(x, bracketIndex);

            // find the value of the function
            return (1 - lambda) * _fBreakPoints[bracketIndex] + lambda * _fBreakPoints[bracketIndex + 1];
        }

        public int BracketIndex(double x)
        {
            int racketIndex = 0; ;

            if (x >= _maxX)
            {
               racketIndex = _numOfBrackets - 1;
            }
            else if (x<= _minX)
            {
                racketIndex = 0;
            }
            else
            {
                switch (_bracketType)
                {
                    case enumBracketType.FixedBracketSize:
                        {
                            racketIndex = (int)Math.Truncate((x - _minX) / _bracketSize);                            
                        }
                        break;
                        
                    case enumBracketType.VariableBracketSize:
                        {
                            for (int bracketIndex = 0; bracketIndex < _numOfBrackets; ++bracketIndex)
                                {
                                    if (_breakPoints[bracketIndex+1] > x)
                                    {
                                        racketIndex = bracketIndex;
                                        break;
                                    }
                                }
                        }
                        break;
                }
            }
            return racketIndex;
        }
        public double Lambda(double x, int bracketIndex)
        {
            if (x >= _maxX)
            {
                return 1;
            }
            else if (x <= _minX)
            {
                return 0;
            }
            else
            {
                return (x - _breakPoints[bracketIndex]) / (_breakPoints[bracketIndex + 1] - _breakPoints[bracketIndex]);
            }
        }
        public double Lambda(double x)
        {
            int bracketIndex = BracketIndex(x);

            return Lambda(x, bracketIndex);
        }

    }// end of LinearSplines class

    // Two Dimention Linear Splines class
    public class LinearSplines2Dimensions
    {
        public enum enumDimention : int
        {
            X1 = 0,
            X2 = 1,            
        }

        double _observationDiscountRate;
        string _name;
        double _minX1, _maxX1;
        double _minX2, _maxX2;
        double[] _breakPointsX1;
        double[] _breakPointsX2;
        double[,] _fBreakPoints;
        double _bracketSizeX1, _bracketSizeX2;
        int _numOfX1Brackets, _numOfX2Brackets;
        LeastSquares _OLS;

        public LinearSplines2Dimensions(string name)
        {
            _name = name;
        }

        public double[,] fBreakPoints
        {
            get
            {
                return _fBreakPoints;
            }
        }

        public int NumOfBreakPoints(enumDimention dim)
        {
            int numOfBreakPoints = 0;
            switch (dim)
            {
                case enumDimention.X1:
                    numOfBreakPoints = _breakPointsX1.Length;
                    break;
                case enumDimention.X2:
                    numOfBreakPoints = _breakPointsX2.Length;
                    break;
            }
            return numOfBreakPoints;
        }

        public void AssignBreakPoints(double minX1, double maxX1, double bracketSizeX1,
                                      double minX2, double maxX2, double bracketSizeX2)
        {
            if ((Math.Truncate((maxX1 - minX1) / bracketSizeX1) != (maxX1 - minX1) / bracketSizeX1)
                || (Math.Truncate((maxX2 - minX2) / bracketSizeX2) != (maxX2 - minX2) / bracketSizeX2))
            {
                MessageBox.Show("Error in the choice of bracket size!");
            }
            _minX1 = minX1;
            _maxX1 = maxX1;
            _minX2 = minX2;
            _maxX2 = maxX2;

            _bracketSizeX1 = bracketSizeX1;
            _bracketSizeX2 = bracketSizeX2;

            // find the number of brackets
            _numOfX1Brackets = (int)((_maxX1 - _minX1) / _bracketSizeX1);
            _numOfX2Brackets = (int)((_maxX2 - _minX2) / _bracketSizeX2);

            // redim arrays
            _breakPointsX1 = new double[_numOfX1Brackets + 1];
            _breakPointsX2 = new double[_numOfX2Brackets + 1];
            _fBreakPoints = new double[_numOfX1Brackets + 1, _numOfX2Brackets + 1];

            // find the break points
            double x1 = minX1;
            double x2 = minX2;
            int i1 = 0, i2 = 0;
            while (x1 <= maxX1)
            {
                _breakPointsX1[i1] = x1;
                i1 += 1;
                x1 += _bracketSizeX1;
            }
            while (x2 <= maxX2)
            {
                _breakPointsX2[i2] = x2;
                i2 += 1;
                x2 += _bracketSizeX2;
            }            
            
        }
        public void AssignBreakPoints(double[] breakPointsX1, double[] breakPointsX2)
        {
            _breakPointsX1 = breakPointsX1;
            _breakPointsX2 = breakPointsX2;
            _numOfX1Brackets = _breakPointsX1.Length - 1;
            _numOfX2Brackets = _breakPointsX2.Length - 1;
            _minX1 = _breakPointsX1[0];
            _minX2 = _breakPointsX2[0];
            _maxX1 = breakPointsX1[breakPointsX1.Length - 1];
            _maxX2 = breakPointsX2[breakPointsX2.Length - 1];
            _fBreakPoints = new double[_numOfX1Brackets + 1, _numOfX2Brackets + 1];
        }

        public void Reset()
        {

        }

        public void SetupTraining(double observationDiscountRate, double perturbationRatio)
        {
            _observationDiscountRate = observationDiscountRate;
            _OLS = new LeastSquares(_name);
            int numOfCols = (_numOfX1Brackets + 1) * (_numOfX2Brackets + 1) - 1;            
            // set up OLS
            _OLS.SetupTraining(numOfCols);
        }

        public double f(double x1, double x2)
        {
            double lambda1, lambda2;
            int bracketIndexX1, bracketIndexX2;
            // find lambdas
            double[] lambda = ConvertXtoLambda(x1, x2);
            // find bracket indeces
            int[] bracketIndeces = BracketIndices(x1, x2);            

            lambda1 = lambda[0];
            lambda2 = lambda[1];
            bracketIndexX1 = bracketIndeces[0];
            bracketIndexX2 = bracketIndeces[1];

            // calcualte function value
            double fValue = (1 - lambda1 - lambda2)* _fBreakPoints[bracketIndexX1, bracketIndexX2]
                            + lambda1 * _fBreakPoints[bracketIndexX1 + 1, bracketIndexX2]
                            + lambda2 * _fBreakPoints[bracketIndexX1, bracketIndexX2 + 1];
            return fValue;
        }        

        public void Update(double x1, double x2, double fValue)
        {
            // create a design row
            double[] OLS_x = new double[(_numOfX1Brackets + 1) * (_numOfX2Brackets + 1)-1];            
            // find bracket indeces
            int[] bracketIndeces = BracketIndices(x1, x2);
            // find lambdas
            double[] lambda = ConvertXtoLambda(x1, x2);

            int bracketIndexX1, bracketIndexX2;
            double lambda1, lambda2;

            bracketIndexX1 = bracketIndeces[0];
            bracketIndexX2 = bracketIndeces[1];
            lambda1 = lambda[0];
            lambda2 = lambda[1];            
            
            // populate design row      
            OLS_x[OLSCoeffIndex(bracketIndexX1, bracketIndexX2)] = 1 - lambda1 - lambda2;
            OLS_x[OLSCoeffIndex(bracketIndexX1 + 1, bracketIndexX2)] = lambda1 ;
            OLS_x[OLSCoeffIndex(bracketIndexX1, bracketIndexX2 + 1)] = lambda2;
            
            // update OLS
            _OLS.Update(OLS_x, fValue, _observationDiscountRate);

            // update fbreakpoints
            _fBreakPoints = ConvertOLSCoeffTofBreakPoints(_OLS.Coefficients);

        }

        private int OLSCoeffIndex(int bracketIndexX1, int bracketIndexX2)
        {
            return bracketIndexX1 * (_numOfX2Brackets+1) + bracketIndexX2;
        }
        private int[] BracketIndices(double x1, double x2)
        {
            int[] bracketIdices = new int[2];

            int bracketIndexX1, bracketIndexX2;
            if (x1 == _maxX1)
            {
                bracketIndexX1 = (int)Math.Truncate((x1 - _minX1) / _bracketSizeX1) - 1;
            }
            else
            {
                bracketIndexX1 = (int)Math.Truncate((x1 - _minX1) / _bracketSizeX1);
            }
            if (x2 == _maxX2)
            {
                bracketIndexX2 = (int)Math.Truncate((x2 - _minX2) / _bracketSizeX2) - 1;
            }
            else
            {
                bracketIndexX2 = (int)Math.Truncate((x2 - _minX2) / _bracketSizeX2);
            }

            bracketIdices[0] = bracketIndexX1;
            bracketIdices[1] = bracketIndexX2;

            return bracketIdices;
        }
        private double[,] ConvertOLSCoeffTofBreakPoints(double[] OLSCoeff)
        {
            double[,] fBreakPoints = new double[_numOfX1Brackets + 1, _numOfX2Brackets + 1];
            int thisOLSCoeffIndex;

            for (int x1Index = 0; x1Index <= _numOfX1Brackets; ++x1Index)
            {
                for (int x2Index = 0; x2Index <= _numOfX2Brackets; ++x2Index)
                {
                    thisOLSCoeffIndex=OLSCoeffIndex(x1Index,x2Index);
                    if (thisOLSCoeffIndex < OLSCoeff.Length)
                    {
                        fBreakPoints[x1Index, x2Index] = OLSCoeff[thisOLSCoeffIndex];
                    }
                }
            }
            return fBreakPoints;
        }
        private double[,] ConvertfBreakPointsToOLSCoeff(double[,] fBreakPoints)
        {
            double[] OLSCoeff = new double[(_numOfX1Brackets + 1)*(_numOfX2Brackets + 1)-1];
            int thisOLSCoeffIndex;

            for (int x1Index = 0; x1Index <= _numOfX1Brackets; ++x1Index)
            {
                for (int x2Index = 0; x2Index <= _numOfX2Brackets; ++x2Index)
                {
                    thisOLSCoeffIndex = OLSCoeffIndex(x1Index, x2Index);
                    if (thisOLSCoeffIndex < OLSCoeff.Length)
                    {
                        OLSCoeff[thisOLSCoeffIndex] = fBreakPoints[x1Index, x2Index];
                    }
                }
            }
            return fBreakPoints;
        }
        private double[] ConvertXtoLambda(double x1, double x2)
        {
            int bracketIndexX1, bracketIndexX2;
            double[] lambda = new double[2];

            // find the index of the lower break point
            
            int[] bracketIndeces = BracketIndices(x1, x2);
            bracketIndexX1 = bracketIndeces[0];
            bracketIndexX2 = bracketIndeces[1];

            // find lambda
            double lambda1 = 0, lambda2 = 0;
            lambda1 = (x1 - _breakPointsX1[bracketIndexX1]) / (_breakPointsX1[bracketIndexX1 + 1] - _breakPointsX1[bracketIndexX1]);
            lambda2 = (x2 - _breakPointsX2[bracketIndexX2]) / (_breakPointsX2[bracketIndexX2 + 1] - _breakPointsX2[bracketIndexX2]);

            lambda[0] = lambda1;
            lambda[1] = lambda2;

            return lambda;
        }
        


    }// end of LinearSplines class

}
