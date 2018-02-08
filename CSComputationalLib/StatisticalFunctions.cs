using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Windows.Forms;
using System.Web.UI.DataVisualization.Charting;

namespace ComputationLib
{
    public static class StatisticalFunctions
    {
        
        public static double Variance(double[] data)
        {
            int count = data.Length;
            //double sum2 = 0;            
            //for (int i = 0; i < count; ++i)
            //{
            //    sum2 += data[i]*data[i];
            //}
            //double ave2 = Math.Pow(Average(data), 2);
            //double var = count * (sum2 / count - ave2)/ (count - 1);

            //return Math.Sqrt(var);

            double aveSofar = 0, varSofar = 0;
            for (int i = 0; i < count; ++i)
            {
                if (i == 0)
                {
                    varSofar = 0;
                    aveSofar = data[0];
                }
                else
                {
                    varSofar = varSofar + i * Math.Pow(data[i] - aveSofar, 2) / (i + 1);
                    aveSofar = i * aveSofar / (i + 1) + data[i] / (i + 1);
                }
            }
            return varSofar / (count - 1);
        }
        public static double StDev(double[] data)
        {
            return Math.Sqrt(Variance(data));
        }
        public static double SumSq(double[] data)
        {
            double sum = 0;
            for (int i = 0; i < data.Length; i++)
                sum += Math.Pow(data[i], 2);

            return sum;
        }
        
        /// <summary>
        ///  returns two-tail inverse of t ( return t such that P(X less than -t or X greater than t). = probability )
        /// </summary>
        /// <param name="probability"></param>
        /// <param name="degreeOfFreedom"></param>
        /// <returns></returns>
        public static double TInv_TwoTails(double probability, int degreeOfFreedom)
        {
            Chart thisChart = new Chart();
            return thisChart.DataManipulator.Statistics.InverseTDistribution(probability, degreeOfFreedom);
        }
        
    }
}
