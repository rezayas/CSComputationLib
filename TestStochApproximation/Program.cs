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

            StochasticApproximation optProb = new StochasticApproximation(
                simModel: new TestBedX2Y2XY(errorVar: 10),
                stepSize_a: new StepSize_a(a0: 20),
                stepSize_Df: new StepSize_Df(c0: 1)
                );

            // initial value of x
            double[] x0 = new double[2] { -10, 20 }; 
            
            // minimize
            optProb.Minimize(
                maxItrs: 5000,
                nLastItrsToAve: 500,
                x0: Vector<double>.Build.DenseOfArray(x0), 
                ifTwoSidedDerivative: true);

            Console.WriteLine("Optimal x:" + optProb.xStar);
            Console.WriteLine("Optimal f:" + optProb.fStar);
            optProb.ExportResultsToCSV("TestX2Y2.csv");
            
        }
    }
}
