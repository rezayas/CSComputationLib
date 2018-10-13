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
            double[] stepSize_as = new double[5] { 10, 25, 50, 75, 100 };
            double[] stepSize_cs = new double[4] { 1, 5, 10, 25};

            // build models
            List<SimModel> models = new List<SimModel>();
            foreach (double a in stepSize_as)
                foreach (double c in stepSize_cs)
                    models.Add(
                        new TestBedX2Y2XY(errorVar: 100)
                        );

            // build a parallel optimizer
            ParallelStochasticApproximation optimization = new ParallelStochasticApproximation(
                simModels: models,
                stepSize_as: stepSize_as,
                stepSize_cs: stepSize_cs
                ); 

            // initial value of x
            double[] x0 = new double[2] { -10, 20 };

            // minimize
            optimization.Minimize(
                maxItrs: 500,
                nLastItrsToAve: 500,
                x0: Vector<double>.Build.DenseOfArray(x0),
                ifTwoSidedDerivative: true,
                ifParallel: true, 
                modelProvidesDerivatives: false
                );

            Console.WriteLine("Optimal x: " + optimization.xStar);
            Console.WriteLine("Optimal f: " + optimization.fStar);
            Console.WriteLine("Optimal a0: " + optimization.aStar);
            Console.WriteLine("Optimal c0: " + optimization.cStar);
            optimization.ExportResultsToCSV("TestX2Y2");

            Console.ReadKey();
        }
    }
}
