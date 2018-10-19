using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ComputationLib;
using MathNet.Numerics.LinearAlgebra;


namespace TestStochApproximation
{
    class Program
    {
        static void Main(string[] args)
        {

            StochasticApproximation optimization = new StochasticApproximation(
                simModel: new TestBedX2Y2XY(errorVar: 100),
                stepSize_a: new StepSize_GH(a0: 20, b: 10),
                stepSize_Df: new StepSize_Df(c0: 5)
                );

            // initial value of x
            double[] x0 = new double[2] { -100, 200 }; 
            
            // minimize
            optimization.Minimize(
                maxItrs: 5000,
                nLastItrsToAve: 500,
                x0: Vector<double>.Build.DenseOfArray(x0), 
                ifTwoSidedDerivative: true);

            Console.WriteLine("Optimal x: " + optimization.xStar);
            Console.WriteLine("Optimal f: " + optimization.fStar);
            optimization.ExportResultsToCSV("TestX2Y2.csv");

            Console.ReadKey();
        }
    }
}
