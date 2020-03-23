using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearRegression;
using MathNet.Numerics.Distributions;
using ComputationLib;

namespace TestLeastSquares
{
    class Program
    {
        static void Main(string[] args)
        {
                        
            // y = 1 + 2x1 + 3x2
            Vector<double> beta = Vector<double>.Build.Dense(new double[3] { 1, 2, 3 });
            Matrix<double> X= Matrix<double>.Build.Dense(1000, 3);
            Vector<double> y = Vector<double>.Build.Dense(1000);
            
            // pupulate non-collinear
            PopulateNonCollinear(beta, X, y);

            // run tests
            Console.WriteLine("---- Testing on non-collinear data:");
            RunTests(X, y);

            // populate collinear
            Console.WriteLine("---- Testing on collinear data:");
            PopulateCollinear(beta, X, y);

            // run tests
            RunTests(X, y);

            Console.ReadKey();
        }

        static void PopulateNonCollinear(Vector<double> beta, Matrix<double> X, Vector<double> y)
        {
            for (int i = 0; i < 1000; i++)
            {
                X[i, 0] = 1;
                X[i, 1] = Normal.Sample(mean: 0, stddev: 1);
                X[i, 2] = Normal.Sample(mean: 1, stddev: 1);

                for (int j = 0; j < 3; j++)
                    y[i] += beta[j] * X[i, j];
                y[i] += Normal.Sample(mean: 0, stddev: 1);
            }
        }

        static void PopulateCollinear(Vector<double> beta, Matrix<double> X, Vector<double> y)
        {
            for (int i = 0; i < 1000; i++)
            {
                X[i, 0] = 1;
                X[i, 1] = Normal.Sample(mean: 0, stddev: 1);
                X[i, 2] = X[i, 1];

                for (int j = 0; j < 3; j++)
                    y[i] += beta[j] * X[i, j];
                y[i] += Normal.Sample(mean: 0, stddev: 1);
            }
        }

        static void RunTests(Matrix<double> X, Vector<double> y)
        {
            // Math.Net algorithm 
            try
            {
                Vector<double> p = MultipleRegression.NormalEquations(X, y);
                Console.WriteLine("Math.Net algorithm:");
                Console.WriteLine(p);
            }
            catch { }

            // my algorithm
            try
            {
                LeastSquares LS = new LeastSquares();
                LS.RunRegression(X: X.ToArray(), y: y.ToArray());
                Console.WriteLine("My algorithm:");
                Console.WriteLine(LS.Coeff);
            }
            catch { }

            // my algorithm with L2 regularization
            LeastSquares LS2 = new LeastSquares(l2Penalty: 0.01);
            LS2.RunRegression(X: X.ToArray(), y: y.ToArray());
            Console.WriteLine("My algorithm with L2 regularilization:");
            Console.WriteLine(LS2.Coeff);

            // recursive version with l2 regularization 
            LeastSquares LS_Recursive = new LeastSquares(l2Penalty: 0.01);
            LS_Recursive.SetupTraining(numOfColumns: 3);

            for (int i = 0; i < 1000; i++)
            {
                double[] x = new double[3];
                for (int j = 0; j < 3; j++)
                    x[j] = X[i, j];

                LS_Recursive.Update(x: x, y: y[i], discountRate: 1);
            }
            Console.WriteLine("My recursive algorithm with L2 regularilization:");
            Console.WriteLine(LS_Recursive.Coeff);
        }
    }
}
