using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ComputationLib;
using RandomNumberGeneratorLib;

namespace TestComputationalLib
{
    class Program
    {
        static void Main(string[] args)
        {
            ThreadSpecificRNG myRND = new ThreadSpecificRNG();
            myRND.Reset(1);

            var array = new double[] { 1.0, 2.0, 3.0, 4.0, 5.0 };

            SupportFunctions.AddToEndOfArray(ref array, 10.0);

            SupportFunctions.Shuffle(myRND, array);

            Console.Write(string.Join(" ", array));
        }
    }
}
