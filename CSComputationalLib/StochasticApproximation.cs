using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.LinearAlgebra;

namespace CSComputationalLib
{

    abstract class SimModel
    {
        public abstract double GetAReplication(Vector<double> x);
    }

    class StepSize
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

    class StochasticApproximation
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

        public void Minimize(int maxItrs, Vector<double> x0)
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

            // irerations of the algorithm
            for (int itr = 0; itr < maxItrs; itr++)
            {
                // estimate the derivative of f at x
                Vector<double> derivative = Vector<double>.Build.Dense(x0.Count());
                for (int i = 0; i < x0.Count(); i++){
                    derivative[i] = (_simModel.GetAReplication(x + epsilonMatrix.Row(i)) - f) / _derivativeStep;
                }

                // find a new x: x_new = x - stepSize*f'(x)
                x = x - _stepSize.GetValue(itr) * derivative;

                // get f(x)
                f = _simModel.GetAReplication(x);

                // store information at iteration 0
                Itr_i.Add(itr);
                Itr_x.Add(x);
                Itr_f.Add(f);
            }

            // store the optimal x and optimal objective value
            XStar = x;
        }
    }
}
