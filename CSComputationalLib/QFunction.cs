using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ComputationLib
{
    public abstract class QFunction
    {
       
        // Fields
        public string Name { get; }

        // Instantiation
        public QFunction(string name)
        {
            Name = name;
        }

        // update
        public virtual void Update(double[] continuousVar, double fValue, double discountRate=1) { }
        public virtual void Update(int[] indicatorVar, double[] continuousVar, double fValue, double discountRate = 1) { }
        // function value
        public virtual double fValue(double[] continuousVar) { return 0; }
        public virtual double fValue(int[] indicatorVar, double[] continuousVar) { return 0; }
        // reset Q-function
        public virtual void Reset() {  }
    }

    public class PolynomialQFunction : QFunction
    {
        // Fields
        LeastSquares _leastSquares;
        public int NumOfContinuousVars { get; }
        public int NumOfIndicatorVariables { get; } // f(i1, i2, x1, x2) = g0(x1, x2) + i1*g1(x1, x2) + i2*g2(x1, x2)
        public int NumOfRegressionColumns { get; }
        public int PolynomialDegree { get; }
        // each row represents a term and columns represent power of each variable for each term
        public int[][] DegreesOfContinuousVarsInPolynomialTerms { get; } 
                
        /// <summary>
        /// Instantiation of polynomial function of the form:
        /// f(i1, i2, x1, x2) = g0(x1, x2) + i1*g1(x1, x2) + i2*g2(x1, x2)
        /// </summary>
        /// <param name="name"></param>
        /// <param name="numOfIndicatorVariables"> set 0 for a regular regression</param>
        /// <param name="numOfContinuousVariables"></param>
        /// <param name="polynomialDegree"></param>
        public PolynomialQFunction(string name, int numOfIndicatorVariables, 
            int numOfContinuousVariables, int polynomialDegree, double l2Penalty=0)
            : base(name)
        {
            NumOfContinuousVars = numOfContinuousVariables;
            NumOfIndicatorVariables = numOfIndicatorVariables;
            PolynomialDegree = polynomialDegree;

            // find the number of regression columns
            DegreesOfContinuousVarsInPolynomialTerms = DegreesOfVariablesInPolynomialTerms(numOfContinuousVariables, polynomialDegree);
            NumOfRegressionColumns = DegreesOfContinuousVarsInPolynomialTerms.GetLength(0) * (numOfIndicatorVariables + 1);

            // setup the least squares
            _leastSquares = new LeastSquares(l2Penalty);
            _leastSquares.SetupTraining(NumOfRegressionColumns);
        }        

        public int[,] RegressionTermDegrees
        {
            get { return SupportFunctions.ConvertJaggedArrayToRegularArray(DegreesOfContinuousVarsInPolynomialTerms, NumOfContinuousVars); }
        }
        public double[] Coefficients
        {
            get { return _leastSquares.Coeff.ToArray(); }
        }

        // update
        public override void Update(double[] continuousVar, double fValue, double discountRate=1)
        {
            //_leastSquares.Update(ConvertToRowDesign(var), fValue, _stepSizeRule.ObservationDiscountRate(itr));
            _leastSquares.Update(ConvertToRowDesign(continuousVar), fValue, discountRate);
        }        
        public override void Update(int[] indicatorVar, double[] continuousVar, double fValue, double discountRate = 1)
        {
            //_leastSquares.Update(ConvertToRowDesign(var), fValue, _stepSizeRule.ObservationDiscountRate(itr));
            _leastSquares.Update(ConvertToRowDesign(indicatorVar, continuousVar), fValue, discountRate);
        }
                        
        /// <summary>
        /// Function value
        /// </summary>
        /// <param name="continuousVar"> </param>
        /// <returns></returns>
        public override double fValue(double[] continuousVar)
        {
            return _leastSquares.yValue(ConvertToRowDesign(continuousVar));
        }
        /// <summary>
        /// Function value
        /// </summary>
        /// <param name="indicatorVar"></param>
        /// <param name="continuousVar"> </param>
        /// <returns></returns>
        public override double fValue(int[] indicatorVar, double[] continuousVar)
        {
            return _leastSquares.yValue(ConvertToRowDesign(indicatorVar, continuousVar));
        }
                
        // update estimates
        public void UpdateCoefficients(double[] coefficients)
        {
            // update the estimates
            _leastSquares.UpdateCoefficients(coefficients);            
        }       
        
        // reset
        public override void Reset()
        {
            _leastSquares.Reset();
        }

        // Private Functions
        // convert the variables into a row design
        private double[] ConvertToRowDesign(double[] continuousVariableValues)
        {
            double[] rowDesign = new double[NumOfRegressionColumns];
            SupportFunctions.MakeArrayEqualTo(ref rowDesign, 1);

            for (int polynomialTermIndex = 0; polynomialTermIndex < NumOfRegressionColumns; polynomialTermIndex++)
                for (int varIndex = 0; varIndex < NumOfContinuousVars; varIndex++)
                    rowDesign[polynomialTermIndex] *= Math.Pow(continuousVariableValues[varIndex], DegreesOfContinuousVarsInPolynomialTerms[polynomialTermIndex][varIndex]);

            return rowDesign;
        }
        private double[] ConvertToRowDesign(int[] indicatorVarValues, double[] continuousVarValues)
        {
            int regColIndex = 0;
            double[] rowDesign = new double[NumOfRegressionColumns];
            SupportFunctions.MakeArrayEqualTo(ref rowDesign, 1);

            // for the first regression model 
            for (int polynomialTermIndex = 0; polynomialTermIndex < DegreesOfContinuousVarsInPolynomialTerms.GetLength(0); polynomialTermIndex++)
                for (int continuousVarIndex = 0; continuousVarIndex < NumOfContinuousVars; continuousVarIndex++)
                    rowDesign[polynomialTermIndex] *= Math.Pow(continuousVarValues[continuousVarIndex], DegreesOfContinuousVarsInPolynomialTerms[polynomialTermIndex][continuousVarIndex]);

            // for the remaining regression models
            regColIndex += DegreesOfContinuousVarsInPolynomialTerms.GetLength(0);
            for (int indicatorVarIndex = 0; indicatorVarIndex < NumOfIndicatorVariables; ++indicatorVarIndex)
            {
                for (int polynomialTermIndex = 0; polynomialTermIndex < DegreesOfContinuousVarsInPolynomialTerms.GetLength(0); polynomialTermIndex++)
                {
                    rowDesign[regColIndex] *= indicatorVarValues[indicatorVarIndex];
                    if (rowDesign[regColIndex] > 0)
                        for (int continuousVarIndex = 0; continuousVarIndex < NumOfContinuousVars; continuousVarIndex++)
                            rowDesign[regColIndex] *= Math.Pow(continuousVarValues[continuousVarIndex], DegreesOfContinuousVarsInPolynomialTerms[polynomialTermIndex][continuousVarIndex]);
                    ++regColIndex;
                }                
            }
            return rowDesign;
        }
                        
        /// <summary>
        /// Find the degree (exponent) of each variable in terms of a polynomial function of a certain degree
        /// Return the matrix where rows represent the terms of a polynomial functions and columns represent the variables. Each entity is the degree of a variable in a polynomial function term.
        /// </summary>
        /// <param name="numOfVariables"></param>
        /// <param name="degree"></param>
        /// <returns></returns>
        private int[][] DegreesOfVariablesInPolynomialTerms(int numOfVariables, int degree)
        {
            int[][] result = new int[0][];
            // if number of variables is 1 the result is a column vector
            if (numOfVariables == 1)
            {
                result = new int[degree + 1][];
                for (int deg = 0; deg <= degree; deg++)
                {
                    result[deg] = new int[1];
                    result[deg][0] = deg;
                }
            }
            else // if number of variables is greater than 1
            {
                // iterate over first variable from degree 0 to max degree 
                for (int firstVarDeg = 0; firstVarDeg <= degree; firstVarDeg++)
                {
                    // find the possible permutations of remaining variables if the degree of the first variable is set
                    int[][] nextPermutation = DegreesOfVariablesInPolynomialTerms(numOfVariables - 1, degree - firstVarDeg);                    
                    // go over the permutations of the remaining variables                    
                    for (int nextVarPermultationRowIndex = 0; nextVarPermultationRowIndex < nextPermutation.GetLength(0); nextVarPermultationRowIndex++)
                    {
                        int[][] thisRow = new int[1][];
                        thisRow[0] = new int[numOfVariables];
                        thisRow[0][0] = firstVarDeg;
                        for (int varIndex = 0; varIndex < numOfVariables - 1; varIndex++)
                            thisRow[0][varIndex + 1] = nextPermutation[nextVarPermultationRowIndex][varIndex];
                        result = SupportFunctions.ConcatJaggedArray(result, thisRow);
                    }
                }
            }
            return result;
        }    
    }

}
