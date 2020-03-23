using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearRegression;
using MathNet.Numerics.Distributions;
using ComputationLib;

namespace TestQFunctions
{
    class Program
    {
        static void Main(string[] args)
        {
            PolynomialQFunction Q = new PolynomialQFunction(
                name: "Testing Q function",
                numOfIndicatorVariables: 0,
                numOfContinuousVariables: 2,
                polynomialDegree: 2,
                l2Penalty: 0);

            // z = 1 + 2x + 3y + 4xy + 5x^2 + 6y^2
            double[] beta = new double[6] { 1, 2, 3, 4, 5, 6 };
            double y;

            for (int i = 0; i < 1000; i++)
            {
                double[] var = new double[2];
                var[0] = Normal.Sample(mean: 0, stddev: 1);
                var[1] = Normal.Sample(mean: 1, stddev: 1);

                y = beta[0]
                    + beta[1] * var[0]
                    + beta[2] * var[1]
                    + beta[3] * var[0] * var[1]
                    + beta[4] * Math.Pow(var[0], 2)
                    + beta[5] * Math.Pow(var[1], 2);

                y += Normal.Sample(mean: 0, stddev: 1);

                Q.Update(continuousVar: var, fValue: y, discountRate: 1);
            }

            for (int i = 0; i < Q.Coefficients.Length; i ++)
            {
                Array.ForEach(Q.DegreesOfContinuousVarsInPolynomialTerms[i], Console.WriteLine);
                Console.WriteLine(Q.Coefficients[i]);
            }

            Console.ReadKey();

        }
    }
}
