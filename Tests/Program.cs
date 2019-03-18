using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ComputationLib;
using RandomVariateLib;
using MathNet.Numerics.LinearAlgebra;

namespace TestComputationalLib
{
    class Program
    {
        static void Main(string[] args)
        {
            RNG myRND = new RNG(1);

            var array = new double[] { 1.0, 2.0, 3.0, 4.0, 5.0 };

            SupportFunctions.AddToEndOfArray(ref array, 10.0);

            SupportFunctions.Shuffle(myRND, array);

            Console.Write(string.Join(" ", array));

            double[] array1 = new double[] { 0.01231423, 1.13412341 };
            Vector<double> vec = Vector<double>.Build.DenseOfArray(array1);

            Console.WriteLine("Testing vector to string: ");
            Console.WriteLine(vec.ToString());
            Console.WriteLine(vec.ToVectorString("G3"));



            Console.WriteLine(string.Join(",", array1));
            string formated = string.Join(",", array1.Select(x => Math.Round(x,3)).ToArray());
            Console.WriteLine(formated);

            Console.ReadKey();

        }
    }
}
