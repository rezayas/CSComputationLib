using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ComputationLib;
using MathNet.Numerics.LinearAlgebra;

namespace TestParallelStochasticApprox
{
    class Program
    {
        static void Main(string[] args)
        {
            double[] stepSizeGH_a0s = new double[2] { 5, 10};
            double[] stepSizeGH_bs = new double[2] { 50, 100};
            double[] stepSizeDf_cs = new double[2] { 5, 10};

            // build models
            List<SimModel> models = new List<SimModel>();
            foreach (double a in stepSizeGH_a0s)
                foreach (double b in stepSizeGH_bs)
                    foreach (double c in stepSizeDf_cs)
                        models.Add(
                            new TestBedX2Y2XY(errorVar: 100)
                            );

            // build a parallel optimizer
            MultipleStochasticApproximation optimization = new MultipleStochasticApproximation(
                simModels: models,
                stepSizeGH_a0s: stepSizeGH_a0s,
                stepSizeGH_bs: stepSizeGH_bs,
                stepSizeDf_cs: stepSizeDf_cs
                ); 

            // initial value of x
            double[] x0 = new double[2] { -100, 200 };
            double[] xScale = new double[2] { 100, 0.1 };

            // minimize
            optimization.Minimize(
                nItrs: 5000,
                nLastItrsToAve: 500,
                x0: Vector<double>.Build.DenseOfArray(x0),
                xScale: Vector<double>.Build.DenseOfArray(xScale),
                modelProvidesDerivatives: false,
                ifTwoSidedDerivative: true,
                ifParallel: true               
                );

            Console.WriteLine("Optimal x: " + optimization.xStar);
            Console.WriteLine("Optimal dx/x: " + optimization.dxOverXAveStar);
            Console.WriteLine("Optimal f: " + optimization.fStar);
            Console.WriteLine("Optimal a0: " + optimization.a0Star);
            Console.WriteLine("Optimal b: " + optimization.bStar);
            Console.WriteLine("Optimal c0: " + optimization.c0Star);
            optimization.ExportResultsToCSV("TestX2Y2");

            // print summary
            Console.WriteLine("");
            Print2DArray(optimization.GetSummary(f_digits:1, x_digits:3));

            Console.ReadKey();
        }

        public static void Print2DArray<T>(T[,] matrix)
        {
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                for (int j = 0; j < matrix.GetLength(1); j++)
                {
                    Console.Write(matrix[i, j] + "\t");
                }
                Console.WriteLine();
            }
        }
    }
}
