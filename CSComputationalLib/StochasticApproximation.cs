using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.Data.Text;

namespace ComputationLib
{

    public abstract class SimModel
    {
        public abstract double GetAReplication(Vector<double> x, bool ifResampleSeeds = true);
    }

    public class TestBedX2Y2XY : SimModel
    {
        RandomVariateLib.Normal _err; // a normally distributed error term
        RandomVariateLib.RNG _rnd;
        int _currentSeed = 0;

        public TestBedX2Y2XY(double errorVar)
        {
            _err = new RandomVariateLib.Normal("Error term", 0, errorVar);
            _rnd = new RandomVariateLib.RNG(_currentSeed);
        }

        public override double GetAReplication(Vector<double> x, bool ifResampleSeeds)
        {
            if (ifResampleSeeds)
            {
                ++_currentSeed;
                _rnd = new RandomVariateLib.RNG(_currentSeed);
            }

            return Math.Pow(x[0], 2) + Math.Pow(x[1], 2) + x[0]*x[1] + _err.SampleContinuous(_rnd);
        }
    }

    public class StepSize
    {
        private double _a;

        public StepSize(double a)
        {
            _a = a;
        }
        public double GetValue(int itr)
        {
            return _a / (itr + 1);
        }
    }

    public class StochasticApproximation
    {
        private double _derivativeStep = 0;
        private StepSize _stepSize = null;
        private SimModel _simModel = null;

        public List<int> Itr_i { get; private set; } = new List<int>();
        public List<Vector<double>> Itr_x { get; private set; } = new List<Vector<double>>();
        public List<double> Itr_f { get; private set; } = new List<double>();
        public List<Vector<double>> Itr_Df { get; private set; } = new List<Vector<double>>();

        public Vector<double> xStar { get; private set; }
        public double fStar { get; private set; }

        public StochasticApproximation(SimModel simModel, double derivativeStep, StepSize stepSize)
        {
            _simModel = simModel;
            _derivativeStep = derivativeStep;
            _stepSize = stepSize;
        }

        public void Minimize(int maxItrs, int nLastItrsToAve, Vector<double> x0, bool ifTwoSidedDerivative = true)
        {
            // iteration 0
            Vector<double> x = x0;
            double f = _simModel.GetAReplication(x, ifResampleSeeds: true);

            // store information at iteration 0
            Itr_i.Add(0);
            Itr_x.Add(x);
            Itr_f.Add(f);            

            // build epsilon matrix
            Matrix<double> epsilonMatrix = Matrix<double>.Build.DenseDiagonal(x0.Count(), _derivativeStep);

            // iterations of the algorithm
            for (int itr = 1; itr <= maxItrs; itr++)
            {
                // estimate the derivative of f at x
                Vector<double> Df = Vector<double>.Build.Dense(x0.Count());
                for (int i = 0; i < x0.Count(); i++){
                    if (ifTwoSidedDerivative)
                    {
                        Df[i] = 
                            (
                                _simModel.GetAReplication(x + epsilonMatrix.Row(i), ifResampleSeeds: false) - 
                                _simModel.GetAReplication(x - epsilonMatrix.Row(i), ifResampleSeeds: false)
                            ) 
                            / (2*_derivativeStep);
                    }
                    else
                    {
                        Df[i] = 
                            (_simModel.GetAReplication(x + epsilonMatrix.Row(i), ifResampleSeeds: false) - f) 
                            / _derivativeStep;
                    }
                }

                // normalize derivative
                Vector<double> nDf = Df.Normalize(p: 2);

                // find a new x: x_new = x - stepSize*f'(x)
                x = x - _stepSize.GetValue(itr) * nDf;

                // get f(x)
                f = _simModel.GetAReplication(x, ifResampleSeeds: true);
                //Console.WriteLine(f+x.ToString()+derivative.ToString());

                // store information at iteration 0
                Itr_i.Add(itr);
                Itr_x.Add(x);
                Itr_f.Add(f);
                Itr_Df.Add(nDf);
            }

            // store the optimal x and optimal objective value 
            double fSum = 0;
            Vector<double> xSum = Vector<double>.Build.Dense(x0.Count);   
            for (int itr = maxItrs; itr > maxItrs - nLastItrsToAve; itr--)
            {
                fSum += Itr_f[itr - 1];
                xSum += Itr_x[itr - 1];
            }
            xStar = xSum / nLastItrsToAve;
            fStar = fSum / nLastItrsToAve;

            // assumed 0 for the derivative at xStar
            Itr_Df.Add(Vector<double>.Build.DenseOfArray(new double[x.Count]));
        }

        public double[,] GetResultsInAMatrix()
        {
            double[,] result = new double[Itr_i.Count, 1 + 1 + 2*xStar.Count];

            for (int itr = 0; itr<Itr_i.Count; itr++)
            {
                int j = 0;
                result[itr, j++] = Itr_i[itr];
                result[itr, j++] = Itr_f[itr];

                for (int i = 0; i < xStar.Count; i++)
                    result[itr, j++] = Itr_x[itr][i];
                for (int i = 0; i < xStar.Count; i++)
                    result[itr, j++] = Itr_Df[itr][i];
            }
            return result;
        }

        public void ExportResultsToCSV(string filename)
        {
            Matrix<double> matrix = Matrix<double>.Build.DenseOfArray(GetResultsInAMatrix());

            List<string> colHeader = new List<string>();
            colHeader.Add("Iteration");
            colHeader.Add("f");
            for (int i = 0; i < xStar.Count; i++)
                colHeader.Add("x"+i);
            for (int i = 0; i < xStar.Count; i++)
                colHeader.Add("Df" + i);

            DelimitedWriter.Write(filename, matrix, ",", columnHeaders: colHeader);
        }
    }
}
