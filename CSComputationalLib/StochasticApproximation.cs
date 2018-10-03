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
        public abstract double GetAReplication(Vector<double> x);
    }

    public class TestBedX2Y2XY : SimModel
    {
        RandomVariateLib.Normal _err; // a normally distributed error term
        RandomVariateLib.RNG _rnd = new RandomVariateLib.RNG(0);

        public TestBedX2Y2XY(double errorVar)
        {
            _err = new RandomVariateLib.Normal("Error term", 0, errorVar);
        }

        public override double GetAReplication(Vector<double> x)
        {
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
        public Vector<double> XStar { get; private set; }

        public StochasticApproximation(SimModel simModel, double derivativeStep, StepSize stepSize)
        {
            _simModel = simModel;
            _derivativeStep = derivativeStep;
            _stepSize = stepSize;
        }

        public void Minimize(int maxItrs, Vector<double> x0, bool ifTwoSidedDerivative = true)
        {
            // iteration 0
            Vector<double> x = x0;
            double f = _simModel.GetAReplication(x);

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
                Vector<double> derivative = Vector<double>.Build.Dense(x0.Count());
                for (int i = 0; i < x0.Count(); i++){
                    if (ifTwoSidedDerivative)
                    {
                        derivative[i] = 
                            (
                                _simModel.GetAReplication(x + epsilonMatrix.Row(i)) - 
                                _simModel.GetAReplication(x - epsilonMatrix.Row(i))
                            ) 
                            / (2*_derivativeStep);
                    }
                    else
                    {
                        derivative[i] = 
                            (_simModel.GetAReplication(x + epsilonMatrix.Row(i)) - f) 
                            / _derivativeStep;
                    }
                }

                // find a new x: x_new = x - stepSize*f'(x)
                x = x - _stepSize.GetValue(itr) * derivative;

                // get f(x)
                f = _simModel.GetAReplication(x);
                //Console.WriteLine(f+x.ToString()+derivative.ToString());

                // store information at iteration 0
                Itr_i.Add(itr);
                Itr_x.Add(x);
                Itr_f.Add(f);
            }

            // store the optimal x and optimal objective value
            XStar = x;
        }

        public double[,] GetResultsInAMatrix()
        {
            double[,] result = new double[Itr_i.Count, 1+1+XStar.Count];

            for (int i = 0; i<Itr_i.Count; i++)
            {
                result[i, 0] = Itr_i[i];
                result[i, 1] = Itr_f[i];
                for (int j = 0; j < XStar.Count; j++)
                    result[i, 2 + j] = Itr_x[i][j];
            }
            return result;
        }

        public void ExportResultsToCSV(string filename)
        {
            Matrix<double> matrix = Matrix<double>.Build.DenseOfArray(GetResultsInAMatrix());

            DelimitedWriter.Write(filename, matrix, ",");
        }
    }
}
