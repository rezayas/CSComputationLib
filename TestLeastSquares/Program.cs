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

            for (int i = 0; i < 1000; i++)
            {
                X[i, 0] = 1;
                X[i, 1] = Normal.Sample(mean: 0, stddev: 1);
                X[i, 2] = Normal.Sample(mean: 1, stddev: 1);

                for (int j = 0; j < 3; j++)
                    y[i] += beta[j] * X[i, j];
                y[i] += Normal.Sample(mean: 0, stddev: 1);
            }

            Vector<double> p = MultipleRegression.NormalEquations(X, y);

            Console.WriteLine("Math.Net algorithm:");
            Console.WriteLine(p);


            // my version 
            LeastSquares LS = new LeastSquares();
            LS.RunRegression(X: X.ToArray(), y: y.ToArray());
            Console.WriteLine("My algorithm:");
            Console.WriteLine(LS.Coeff);

            // my version with L2 regularization
            LeastSquares LS2 = new LeastSquares();
            LS2.AddL2Regularization(0.01);
            LS2.RunRegression(X: X.ToArray(), y: y.ToArray());
            Console.WriteLine("My algorithm with L2 regularilization:");
            Console.WriteLine(LS2.Coeff);


            Console.ReadKey();

        }
    }
}
