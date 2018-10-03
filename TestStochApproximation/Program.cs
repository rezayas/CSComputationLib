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
                derivativeStep: 1,
                stepSize: new StepSize(a: 1)
                );

            // minimize
            double[] x0 = new double[2] { -10, 20 };
            optProb.Minimize(
                maxItrs: 5000, 
                x0: Vector<double>.Build.DenseOfArray(x0), 
                ifTwoSidedDerivative: true);

            Console.WriteLine(optProb.XStar);
            optProb.ExportResultsToCSV("TestX2Y2.csv");
            
        }
    }
}
