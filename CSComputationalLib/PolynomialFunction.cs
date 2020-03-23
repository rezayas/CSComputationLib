using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ComputationLib
{
    public class PolynomialFunction
    {
        string _name;
        int _degree;        
        LeastSquares _leastSquares;
        double[] _coefficients;

        public PolynomialFunction(string name, int degree)
        {
            _name = name;
            _degree = degree;
        }

        public int Degree
        {
            get{ return _degree; }
        }
        public double[] Coefficient
        {
            get{ return _coefficients; }
        }

        public void AssignCoefficient(double[] coefficients)
        {
            _coefficients = coefficients;
        }
        public void SetupTraining(int multiplyNumOfColumnsByThisFactorToBeginTraining = 1)
        {
            _leastSquares = new LeastSquares();
            _leastSquares.SetupTraining(_degree+1);
        }

        // update the parameters using least squares
        public void Update(double x, double y)
        {
            double[] xs = new double[_degree +1];
            // build a row
            for (int i =0; i <= _degree;++i)
                xs[i] = Math.Pow(x,i);
            // update
            _leastSquares.Update(xs, y, 1);
            // get the coefficient
            _coefficients = _leastSquares.Coeff.ToArray();
        }

        // find a the minimum
        public double FindaMinimizer(double initialx, double minX, double maxX, double precisionError)
        {
            bool optimalFound = false;
            double x=0, nextX;
            double fxPrime, fxDoublePrime;
            double step;
            int counter = 0;
            int boundryCounters = 0;
            Random rnd = new Random(1);

            nextX = initialx;
            while (!optimalFound)
            {
                x = nextX;
                fxPrime = fPrimeValue(x);
                fxDoublePrime = fDoublePrimeValue(x);
                step = -fxPrime / fxDoublePrime;

                // consider the boundries in finding the next x
                if (step >= 0)
                {
                    //nextX = StatisticalFunctions.Min(maxX, x + step);
                }
                else
                {
                    //nextX = StatisticalFunctions.Max(minX, x + step);
                }
                counter += 1;

                // check if this is a minimum
                if (counter > 100) 
                {
                    
                    nextX = minX + rnd.NextDouble() * (maxX - minX);
                    counter = 0;
                }
                else if (nextX == minX && fxPrime >= 0)
                {
                    if (boundryCounters < 5)
                    {                        
                        nextX = minX + rnd.NextDouble() * (maxX - minX);
                        boundryCounters += 1;
                    }
                    else
                    {
                        optimalFound = true;
                    }
                }
                else if (nextX == maxX && fxPrime <= 0)
                {
                    if (boundryCounters < 5)
                    {
                        nextX = minX + rnd.NextDouble() * (maxX - minX);
                        boundryCounters += 1;
                    }
                    else
                    {
                        optimalFound = true;
                    }
                }
                else if (Math.Abs(nextX - x) < precisionError)
                {
                    if (fxDoublePrime >= 0)
                    {
                        optimalFound = true;
                    }
                    else
                    {
                        nextX = minX + rnd.NextDouble() * (maxX - minX);                                               
                    }
                }

                 

            }

            return nextX;
            

        }

        // find zero of the function using Newton method
        public double FindAZero(double minRange, double maxRange, double maxError)
        {
            double x_n, x_nPlus1 = 0;
            double f_xn, fPrime_xn;
            double error = double.MaxValue;
            
            // initialize x_n
            x_n = (minRange + maxError)/2;

            // do while error is not less thatn the maxError
            while (error >= maxError)
            {
                f_xn = fValue(x_n);
                fPrime_xn = fPrimeValue(x_n);

                x_nPlus1 = x_n - f_xn / fPrime_xn;

                error = Math.Abs(x_nPlus1 - x_n);
            }
            return x_nPlus1;
        }

        public double fValue(double x)
        {
            double sum = 0;
            for (int i = 0; i <= _degree; ++i)
            {
                sum += _coefficients[i] * Math.Pow(x, i);
            }
            return sum;
        }
        public double fPrimeValue(double x)
        {
            double sum = 0;
            for (int i = 1; i <= _degree; ++i)
            {
                sum += i * _coefficients[i] * Math.Pow(x, i-1);
            }
            return sum;
        }
        public double fDoublePrimeValue(double x)
        {
            double sum = 0;
            for (int i = 2; i <= _degree; ++i)
            {
                sum += i * (i-1) * _coefficients[i] * Math.Pow(x, i - 2);
            }
            return sum;
        }

    }
}
