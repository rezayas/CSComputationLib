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
                stepSize_GH: new StepSize_GH(a0: 5, b: 100),
                stepSize_Df: new StepSize_Df(c0: 5)
                );

            // initial value of x
            double[] x0 = new double[2] { -50, 50 };
            double[] xScale = new double[2] { 1000, 0.1 };
            
            // minimize
            optimization.Minimize(
                nItrs: 5000,
                nLastItrsToAve: 500,
                x0: Vector<double>.Build.DenseOfArray(x0),
                xScale: Vector<double>.Build.DenseOfArray(xScale),
                modelProvidesDerivatives: false,
                ifTwoSidedDerivative: true);

            Console.WriteLine("Optimal x: " + optimization.xStar);
            Console.WriteLine("Average dx over x: " + optimization.dx_over_x_ave);
            Console.WriteLine("Optimal f: " + optimization.fStar);
            optimization.ExportResultsToCSV("TestX2Y2.csv");

            Console.ReadKey();
        }
    }
}
