using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MathNet.Numerics.Statistics;

namespace ComputationLib
{
    public static class StatisticalFunctions
    {        
        public static double Variance(double[] data)
        {
            return Statistics.Variance(data);
        }
        public static double StDev(double[] data)
        {
            return Statistics.StandardDeviation(data);
        }
        public static double SumSq(double[] data)
        {
            double sum = 0;
            for (int i = 0; i < data.Length; i++)
                sum += Math.Pow(data[i], 2);    
            return sum;
        }               
    }
}
