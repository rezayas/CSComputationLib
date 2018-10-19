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
            double[] stepSizeGH_a0s = new double[2] { 10, 2};
            double[] stepSizeGH_bs = new double[2] { 1, 10};
            double[] stepSizeDf_cs = new double[2] { 1, 5};

            // build models
            List<SimModel> models = new List<SimModel>();
            foreach (double a in stepSizeGH_a0s)
                foreach (double b in stepSizeGH_bs)
                    foreach (double c in stepSizeDf_cs)
                        models.Add(
                            new TestBedX2Y2XY(errorVar: 10)
                            );

            // build a parallel optimizer
            ParallelStochasticApproximation optimization = new ParallelStochasticApproximation(
                simModels: models,
                stepSizeGH_a0s: stepSizeGH_a0s,
                stepSizeGH_bs: stepSizeGH_bs,
                stepSizeDf_cs: stepSizeDf_cs
                ); 

            // initial value of x
            double[] x0 = new double[2] { -100, 200 };

            // minimize
            optimization.Minimize(
                maxItrs: 5000,
                nLastItrsToAve: 500,
                x0: Vector<double>.Build.DenseOfArray(x0),
                ifTwoSidedDerivative: true,
                ifParallel: true, 
                modelProvidesDerivatives: false
                );

            Console.WriteLine("Optimal x: " + optimization.xStar);
            Console.WriteLine("Optimal f: " + optimization.fStar);
            Console.WriteLine("Optimal a0: " + optimization.a0Star);
            Console.WriteLine("Optimal b: " + optimization.bStar);
            Console.WriteLine("Optimal c0: " + optimization.c0Star);
            optimization.ExportResultsToCSV("TestX2Y2");

            Console.ReadKey();
        }
    }
}
