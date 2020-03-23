using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ComputationLib
{
    public abstract class OldQFunction
    {

        // Fields
        public string Name { get; }

        public double[] MinimumX { get; protected set; }
        public double FuncValueAtMinimumX { get; protected set; }
        public int NumOfOptimizationItrs { get; protected set; }

        // Instantiation
        public OldQFunction(string name)
        {
            Name = name;
        }

        // update
        public virtual void Update(double[] continuousVar, double fValue, int itr) { }
        public virtual void Update(double[] continuousVar, double fValue) { }
        public virtual void Update(int[] indicatorVar, double[] continuousVar, double fValue, int itr) { }
        public virtual void Update(int[] indicatorVar, double[] continuousVar, double fValue) { }
        // function value
        public virtual double fValue(double[] continuousVar) { return 0; }
        public virtual double fValue(int[] indicatorVar, double[] continuousVar) { return 0; }
        // gradient 
        public virtual double[] fGradientValue(double[] var) { return null; }
        // jacobian
        public virtual double[,] fJacobianValue(double[] var) { return null; }
        // reset Q-function
        public virtual void Reset() { }
    }

    public class OldPolynomialQFunction : OldQFunction
    {
        // Fields
        LeastSquares _leastSquares;
        int _numOfContinuousVariables;
        int _numOfIndicatorVariables; // f(i1, i2, x1, x2) = g0(x1, x2) + i1*g1(x1, x2) + i2*g2(x1, x2)
        int _polynomialDegree;
        int _numOfRegColumns;
        int[][] _degreesOfContinuousVariablesInPolynomialTerms; // each row represents a term and columns represent power of each variable for each term

        /// <summary>
        /// Instantiation of polynomial function of the form:
        /// f(i1, i2, x1, x2) = g0(x1, x2) + i1*g1(x1, x2) + i2*g2(x1, x2)
        /// </summary>
        /// <param name="name"></param>
        /// <param name="numOfIndicatorVariables"> set 0 for a regular regression</param>
        /// <param name="numOfContinuousVariables"></param>
        /// <param name="polynomialDegree"></param>
        public OldPolynomialQFunction(string name, int numOfIndicatorVariables, int numOfContinuousVariables, int polynomialDegree, int multiplyNumOfColumnsByThisFactorToBeginTraining = 1)
            : base(name)
        {
            _numOfContinuousVariables = numOfContinuousVariables;
            _numOfIndicatorVariables = numOfIndicatorVariables;
            _polynomialDegree = polynomialDegree;

            // find the number of regression columns
            _degreesOfContinuousVariablesInPolynomialTerms = DegreesOfVariablesInPolynomialTerms(numOfContinuousVariables, polynomialDegree);
            _numOfRegColumns = _degreesOfContinuousVariablesInPolynomialTerms.GetLength(0) * (numOfIndicatorVariables + 1);

            // setup the least squares
            _leastSquares = new LeastSquares();
            _leastSquares.SetupTraining(_numOfRegColumns);
        }

        // Properties
        public int NumberOfRegressionColumns
        {
            get { return _numOfRegColumns; }
        }
        public int[,] RegressionTermDegrees
        {
            get { return SupportFunctions.ConvertJaggedArrayToRegularArray(_degreesOfContinuousVariablesInPolynomialTerms, _numOfContinuousVariables); }
        }
        public double[] CoeffientEstimates
        {
            get { return _leastSquares.Coeff.ToArray(); }
        }

        // update
        public override void Update(double[] continuousVar, double fValue, int itr)
        {

            //_leastSquares.Update(ConvertToRowDesign(var), fValue, _stepSizeRule.ObservationDiscountRate(itr));
            _leastSquares.Update(ConvertToRowDesign(continuousVar), fValue, 1);
        }
        public override void Update(double[] continuousVar, double fValue)
        {
            _leastSquares.Update(ConvertToRowDesign(continuousVar), fValue, 1);
        }
        public override void Update(int[] indicatorVar, double[] continuousVar, double fValue, int itr)
        {
            //_leastSquares.Update(ConvertToRowDesign(var), fValue, _stepSizeRule.ObservationDiscountRate(itr));
            _leastSquares.Update(ConvertToRowDesign(indicatorVar, continuousVar), fValue, 1);
        }
        public override void Update(int[] indicatorVar, double[] continuousVar, double fValue)
        {
            _leastSquares.Update(ConvertToRowDesign(indicatorVar, continuousVar), fValue, 1);
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
        // calculate the gradient
        public override double[] fGradientValue(double[] var)
        {
            double[] gradient = new double[_numOfContinuousVariables];
            // find the gradient for each variable            
            for (int varIndex = 0; varIndex < _numOfContinuousVariables; ++varIndex)
                gradient[varIndex] = _leastSquares.yValue(ConvertToGradientRowDesignWithRespectToAVariable(var, varIndex));

            return gradient;
        }
        // calculate the jacobian
        public override double[,] fJacobianValue(double[] var)
        {
            double[,] jacobian = new double[_numOfContinuousVariables, _numOfContinuousVariables];

            // find the jacobian
            for (int var1Index = 0; var1Index < _numOfContinuousVariables; ++var1Index)
                for (int var2Index = var1Index; var2Index < _numOfContinuousVariables; ++var2Index)
                {
                    jacobian[var1Index, var2Index] = _leastSquares.yValue(ConvertToJacobianRowDesignWithRespectToAVariable(var, var1Index, var2Index));
                    if (var2Index > var1Index)
                        jacobian[var2Index, var1Index] = jacobian[var1Index, var2Index];
                }

            return jacobian;
        }

        // update estimates
        public void UpdateCoefficients(double[] estimates)
        {
            // update the estimates
            _leastSquares.UpdateCoefficients(estimates);
        }

        // minimize using gradient descent method
        //public void MinimizeUsingGradientDescent(EnumLineSearchMethod lineSearchMethod, double[] initialVariableValues, double stepSize, double normOfGradientToStop)
        //{
        //    int numOfIterations = 0;
        //    double[] currentVar = initialVariableValues;
        //    double[] currentGradient; // = new double[_numOfContinuousVariables];
        //    // = new double[_numOfContinuousVariables, _numOfContinuousVariables];
        //    double[] direction = new double[_numOfContinuousVariables];
        //    double[] newVar = new double[_numOfContinuousVariables];
        //    double norm = double.MaxValue;

        //    // calculate the gradient
        //    currentGradient = fGradientValue(initialVariableValues);
        //    // calculate the norm of gradient
        //    norm = LinearAlgebraFunctions.Norm(currentGradient, LinearAlgebraFunctions.enumVectorNorm.L_inf);

        //    // do while error is still too big
        //    while (norm >= normOfGradientToStop)
        //    {
        //        // calculate the new direction
        //        switch (lineSearchMethod)
        //        {
        //            case EnumLineSearchMethod.SteepestDescent:
        //                {
        //                    for (int varIndex = 0; varIndex < _numOfContinuousVariables; varIndex++)
        //                        direction[varIndex] = -stepSize * currentGradient[varIndex];
        //                }
        //                break;
        //            case EnumLineSearchMethod.Newton:
        //                {
        //                    double[,] currentJacobian = fJacobianValue(currentVar);
        //                    while (!LinearAlgebraFunctions.Matrix_IfInvertable(currentJacobian))
        //                    {
        //                        currentVar[0] += 0.001;
        //                        currentJacobian = fJacobianValue(currentVar);
        //                    }
        //                    double[,] invOfJacubian = LinearAlgebraFunctions.Matrix_Inverse(currentJacobian);
        //                    double[] invJTimesG = LinearAlgebraFunctions.Matrix_Multiply(invOfJacubian, currentGradient);

        //                    for (int varIndex = 0; varIndex < _numOfContinuousVariables; varIndex++)
        //                        direction[varIndex] = -stepSize * invJTimesG[varIndex];
        //                }
        //                break;
        //        }
        //        // find a new variable
        //        for (int varIndex = 0; varIndex < _numOfContinuousVariables; varIndex++)
        //            newVar[varIndex] = currentVar[varIndex] + direction[varIndex];

        //        // update current variable
        //        currentVar = (double[])newVar.Clone();

        //        // calculate the gradient
        //        currentGradient = fGradientValue(currentVar);

        //        // calculate the norm of the gradient
        //        norm = LinearAlgebraFunctions.Norm(currentGradient, LinearAlgebraFunctions.enumVectorNorm.L_inf);

        //        // increment the number of iterations
        //        ++numOfIterations;
        //    }

        //    MinimumX = (double[])newVar.Clone();
        //    FuncValueAtMinimumX = fValue(MinimumX);
        //    NumOfOptimizationItrs = numOfIterations;
        //}

        //// minimize using Nelder-Mead solver
        //public void MinimizeUsingNelderMeadSolver(double[] initialVariableValues)
        //{
        //    var solution = NelderMeadSolver.Solve(
        //       x => fValue(x), initialVariableValues);

        //    _minimumX = new double[initialVariableValues.Length];
        //    for (int varIndex = 0; varIndex < initialVariableValues.Length; ++ varIndex)
        //        _minimumX[varIndex] = solution.GetValue(varIndex + 1);

        //    _fAtMinimumX = fValue(_minimumX);            
        //}
        //public void MinimizeUsingNelderMeadSolver(double[] initialVariableValues, double[] xLBounds, double[] xUBounds)
        //{
        //    var solution = NelderMeadSolver.Solve(
        //       x => fValue(x), initialVariableValues, xLBounds, xUBounds);

        //    _minimumX = new double[initialVariableValues.Length];
        //    for (int varIndex = 0; varIndex < initialVariableValues.Length; ++varIndex)
        //        _minimumX[varIndex] = solution.GetValue(varIndex + 1);

        //    _fAtMinimumX = fValue(_minimumX);
        //}

        // reset
        public override void Reset()
        {
            _leastSquares.Reset();
        }

        // Private Functions
        // convert the variables into a row design
        private double[] ConvertToRowDesign(double[] continuousVariableValues)
        {
            double[] rowDesign = new double[_numOfRegColumns];
            SupportFunctions.MakeArrayEqualTo(ref rowDesign, 1);

            for (int polynomialTermIndex = 0; polynomialTermIndex < _numOfRegColumns; polynomialTermIndex++)
                for (int varIndex = 0; varIndex < _numOfContinuousVariables; varIndex++)
                    rowDesign[polynomialTermIndex] *= Math.Pow(continuousVariableValues[varIndex], _degreesOfContinuousVariablesInPolynomialTerms[polynomialTermIndex][varIndex]);

            return rowDesign;
        }
        private double[] ConvertToRowDesign(int[] indicatorVariableValues, double[] continuousVariableValues)
        {
            int regColIndex = 0;
            double[] rowDesign = new double[_numOfRegColumns];
            SupportFunctions.MakeArrayEqualTo(ref rowDesign, 1);

            // for the first regression model 
            for (int polynomialTermIndex = 0; polynomialTermIndex < _degreesOfContinuousVariablesInPolynomialTerms.GetLength(0); polynomialTermIndex++)
                for (int continuousVarIndex = 0; continuousVarIndex < _numOfContinuousVariables; continuousVarIndex++)
                    rowDesign[polynomialTermIndex] *= Math.Pow(continuousVariableValues[continuousVarIndex], _degreesOfContinuousVariablesInPolynomialTerms[polynomialTermIndex][continuousVarIndex]);

            // for the remaining regression models
            regColIndex += _degreesOfContinuousVariablesInPolynomialTerms.GetLength(0);
            for (int indicatorVarIndex = 0; indicatorVarIndex < _numOfIndicatorVariables; ++indicatorVarIndex)
            {
                for (int polynomialTermIndex = 0; polynomialTermIndex < _degreesOfContinuousVariablesInPolynomialTerms.GetLength(0); polynomialTermIndex++)
                {
                    rowDesign[regColIndex] *= indicatorVariableValues[indicatorVarIndex];
                    if (rowDesign[regColIndex] > 0)
                        for (int continuousVarIndex = 0; continuousVarIndex < _numOfContinuousVariables; continuousVarIndex++)
                            rowDesign[regColIndex] *= Math.Pow(continuousVariableValues[continuousVarIndex], _degreesOfContinuousVariablesInPolynomialTerms[polynomialTermIndex][continuousVarIndex]);
                    ++regColIndex;
                }
            }
            return rowDesign;
        }

        // convert the variables into a row design for gradient with respect to a variable
        private double[] ConvertToGradientRowDesignWithRespectToAVariable(double[] continuousVariableValues, int gradientVarIndex)
        {
            double[] rowDesign = new double[_numOfRegColumns];
            SupportFunctions.MakeArrayEqualTo(ref rowDesign, 1);

            for (int polynomialTermIndex = 0; polynomialTermIndex < _numOfRegColumns; polynomialTermIndex++)
                for (int varIndex = 0; varIndex < _numOfContinuousVariables; varIndex++)
                {
                    int degreeOfThisVariable = _degreesOfContinuousVariablesInPolynomialTerms[polynomialTermIndex][varIndex];
                    if (varIndex == gradientVarIndex)
                    {
                        // calculate the derivative if gradient is calculated with respect to this variable
                        if (degreeOfThisVariable > 0)
                            rowDesign[polynomialTermIndex] *= degreeOfThisVariable * Math.Pow(continuousVariableValues[gradientVarIndex], degreeOfThisVariable - 1);
                        else // if degree of this variable in this term is less than 0, then the derivative is 0
                            rowDesign[polynomialTermIndex] = 0;
                    }
                    else
                        rowDesign[polynomialTermIndex] *= Math.Pow(continuousVariableValues[varIndex], degreeOfThisVariable);
                }
            return rowDesign;
        }
        // convert the variables into a row design for calculating jacobian with respect to 2 variables
        private double[] ConvertToJacobianRowDesignWithRespectToAVariable(double[] continuousVariableValues, int firstVarIndex, int secondVarIndex)
        {
            double[] rowDesign = new double[_numOfRegColumns];
            SupportFunctions.MakeArrayEqualTo(ref rowDesign, 1);

            for (int polynomialTermIndex = 0; polynomialTermIndex < _numOfRegColumns; polynomialTermIndex++)
                for (int varIndex = 0; varIndex < _numOfContinuousVariables; varIndex++)
                {
                    int degreeOfThisVariable = _degreesOfContinuousVariablesInPolynomialTerms[polynomialTermIndex][varIndex];
                    if (varIndex == firstVarIndex || varIndex == secondVarIndex)
                    {
                        // if both variables which jabobian is calculated with respect to are the same (diagonal elements of jacobian)
                        if (firstVarIndex == secondVarIndex)
                        {
                            if (degreeOfThisVariable > 1)
                                rowDesign[polynomialTermIndex] *= degreeOfThisVariable * (degreeOfThisVariable - 1) * Math.Pow(continuousVariableValues[firstVarIndex], degreeOfThisVariable - 2);
                            else
                                rowDesign[polynomialTermIndex] = 0;
                        }
                        else // off diagonal elements of jacobian
                        {
                            if (degreeOfThisVariable > 0)
                                rowDesign[polynomialTermIndex] *= degreeOfThisVariable * Math.Pow(continuousVariableValues[firstVarIndex], degreeOfThisVariable - 1);
                            else
                                rowDesign[polynomialTermIndex] = 0;
                        }
                    }
                    else
                        rowDesign[polynomialTermIndex] *= Math.Pow(continuousVariableValues[varIndex], degreeOfThisVariable);
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
            // if number of variables is 1 the result is just a vector
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
