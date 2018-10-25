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
        public virtual double GetAReplication(Vector<double> x, bool ifResampleSeeds = true) { return 0; }
        public virtual Vector<double> GetDerivativeEstimate(Vector<double> x, double derivative_step) { return null; }
        public virtual void ResetSeedAtItr0() { }

        public virtual void Sample_f_and_Df(
            Vector<double> x, 
            double derivative_step, 
            Vector<double> xScale = null,
            bool ifResampleSeeds = true) { }

        public virtual double Get_f() { return 0; }
        public virtual Vector<double> Get_Df() { return null; }
        
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
        public override void ResetSeedAtItr0()
        {
            _currentSeed = 0;
        }
    }

    public class StepSize_GH
    {
        // generalized harmonic (GH) stepsize
        // step_n = a0 * b / (b + n) for n >= 0, a0 > 0, and b >= 1

        public double a0 { get; }
        public double b { get; }

        public StepSize_GH(double a0, double b=1)
        {
            this.a0 = a0;
            this.b = b;
        }
        public double GetValue(int itr)
        {
            return a0 * b/ (itr + b);
        }
    }

    public class StepSize_Df
    {
        public double c0 { get; }

        public StepSize_Df(double c0)
        {
            this.c0 = c0;
        }
        public double GetValue(int itr)
        {
            return c0 * Math.Pow(itr+1, -0.25);
        }
    }

    public class StochasticApproximation
    {
        private StepSize_GH _stepSize_GH = null;
        private StepSize_Df _stepSize_Df = null;
        private SimModel _simModel = null;

        public List<int> Itr_i { get; private set; } = new List<int>();                
        public List<Vector<double>> Itr_x { get; private set; } = new List<Vector<double>>();
        public List<double> Itr_f { get; private set; } = new List<double>();
        public List<Vector<double>> Itr_Df { get; private set; } = new List<Vector<double>>();
        public List<double> Itr_step_Df { get; private set; } = new List<double>();
        public List<double> Itr_step_GH { get; private set; } = new List<double>();

        public Vector<double> xStar { get; private set; }
        public double fStar { get; private set; }

        public StochasticApproximation(SimModel simModel, StepSize_GH stepSize_a, StepSize_Df stepSize_Df)
        {
            _simModel = simModel;
            _stepSize_GH = stepSize_a;
            _stepSize_Df = stepSize_Df;
        }

        public void Minimize(int maxItrs, int nLastItrsToAve, Vector<double> x0, Vector<double> xScale = null,
            bool ifTwoSidedDerivative = true, bool modelProvidesDerivatives = false)
        {

            if (xScale is null)
                xScale = Vector<double>.Build.Dense(length: x0.Count, value: 1);

            // reset seed of the simulation model at iteration 0
            // note that this method could be empty if there is no need to reset the seed 
            _simModel.ResetSeedAtItr0();

            // iteration 0
            Vector<double> x = x0;
            double f;         

            // iterations of the algorithm
            for (int itr = 0; itr < maxItrs; itr++)
            {
                // current derivative step size
                double step_Df = _stepSize_Df.GetValue(itr);                

                // estimate the derivative of f at x
                Vector<double> Df = Vector<double>.Build.Dense(x0.Count());

                // calcualte derivative 
                if (modelProvidesDerivatives)
                {
                    // calcualte f and Df 
                    _simModel.Sample_f_and_Df(x, step_Df, xScale: xScale, ifResampleSeeds: true);

                    // get f(x)
                    f = _simModel.Get_f();

                    // get the derivative from the model
                    Df = _simModel.Get_Df();
                    
                }
                else
                {
                    // get f(x)
                    f = _simModel.GetAReplication(x, ifResampleSeeds: true);

                    // build epsilon matrix
                    Matrix<double> epsilonMatrix = Matrix<double>.Build.DenseDiagonal(x0.Count(), step_Df);

                    for (int i = 0; i < x0.Count(); i++)
                    {
                        if (ifTwoSidedDerivative)
                        {                         
                            // estimate the derivative here
                            Df[i] =
                                (
                                _simModel.GetAReplication(x + epsilonMatrix.Row(i) * xScale[i], ifResampleSeeds: false) -
                                _simModel.GetAReplication(x - epsilonMatrix.Row(i) * xScale[i], ifResampleSeeds: false)
                                ) / (2 * step_Df * xScale[i]);
                        }
                        else
                        {
                            Df[i] =
                                (_simModel.GetAReplication(x + epsilonMatrix.Row(i) * xScale[i], ifResampleSeeds: false) - f)
                                / (step_Df * xScale[i]);
                        }
                    }
                }

                // normalize derivative
                Vector<double> nDf = Df.Normalize(p: 2);
                
                // find step size
                double step_GH = _stepSize_GH.GetValue(itr);

                // store information of this iteration 
                Itr_i.Add(itr);
                Itr_x.Add(x);
                Itr_f.Add(f);
                Itr_Df.Add(nDf);
                Itr_step_Df.Add(step_Df);
                Itr_step_GH.Add(step_GH);

                // find a new x: x_new = x - stepSize*f'(x)
                for (int i = 0; i < x.Count; i++)
                    x[i] = x[i] - step_GH * nDf[i] * xScale[i];              
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
        }

        public double[,] GetResultsInAMatrix()
        {
            double[,] result = new double[Itr_i.Count, 2 + 2*xStar.Count + 2];

            for (int itr = 0; itr<Itr_i.Count; itr++)
            {
                int j = 0;

                // iteration
                result[itr, j++] = Itr_i[itr];
                // f
                result[itr, j++] = Itr_f[itr];

                // x
                for (int i = 0; i < xStar.Count; i++)
                    result[itr, j++] = Itr_x[itr][i];
                // Df(x)
                for (int i = 0; i < xStar.Count; i++)
                    result[itr, j++] = Itr_Df[itr][i];
                
                // steps
                result[itr, j++] = Itr_step_Df[itr];
                result[itr, j++] = Itr_step_GH[itr];
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
            colHeader.Add("Step_Df");
            colHeader.Add("Step_GH");

            DelimitedWriter.Write(filename, matrix, ",", columnHeaders: colHeader);
        }

        public string Get_a0_b_c0()
        {
            return "a0-" + _stepSize_GH.a0 + "- b-" + _stepSize_GH.b + "- c0-" + _stepSize_Df.c0; //.ToString("F2")
        }
        public double Get_a0()
        {
            return _stepSize_GH.a0;
        }
        public double Get_b()
        {
            return _stepSize_GH.b;
        }
        public double Get_c0()
        {
            return _stepSize_Df.c0;
        }
    }

    public class ParallelStochasticApproximation
    {
        List<StochasticApproximation> stochasticApproximations = new List<StochasticApproximation>();
        public double fStar { get; private set; } = double.MaxValue;
        public Vector<double> xStar { get; private set; }
        public double a0Star { get; private set; } = double.NaN;
        public double bStar { get; private set; } = double.NaN;
        public double c0Star { get; private set; } = double.NaN;

        public ParallelStochasticApproximation(List<SimModel> simModels, double[] stepSizeGH_a0s, double[] stepSizeGH_bs, double[] stepSizeDf_cs)
        {

            // build the stochastic approximations
            int i = 0;
            foreach (double a0 in stepSizeGH_a0s)
            {
                foreach (double b in stepSizeGH_bs)
                {
                    foreach (double c in stepSizeDf_cs)
                    {
                        stochasticApproximations.Add(
                            new StochasticApproximation(
                                simModel: simModels[i++],
                                stepSize_a: new StepSize_GH(a0, b),
                                stepSize_Df: new StepSize_Df(c))
                                );
                    }
                }
            }
        }

        public void Minimize(int maxItrs, int nLastItrsToAve, Vector<double> x0, Vector<double> xScale = null,
            bool ifTwoSidedDerivative = true, bool ifParallel = true, bool modelProvidesDerivatives = false)
        {
            if (ifParallel && stochasticApproximations.Count > 1)
            {
                Parallel.ForEach(stochasticApproximations, stocApprx =>
                {
                    stocApprx.Minimize(maxItrs, nLastItrsToAve, x0, xScale, ifTwoSidedDerivative, modelProvidesDerivatives);
                });
            }
            else
            {
                foreach (StochasticApproximation stocApprx in stochasticApproximations)
                {
                    stocApprx.Minimize(maxItrs, nLastItrsToAve, x0, xScale, ifTwoSidedDerivative, modelProvidesDerivatives);
                }
            }

            // find the optimizer
            // find the a value that minimizes f
            xStar = Vector<double>.Build.Dense(x0.Count);
            fStar = double.MaxValue;
            foreach (StochasticApproximation stocApprx in stochasticApproximations)
            {
                // if this a led to the minimum f
                if (stocApprx.fStar < fStar)
                {
                    fStar = stocApprx.fStar;
                    xStar = stocApprx.xStar;
                    a0Star = stocApprx.Get_a0();
                    bStar = stocApprx.Get_b();
                    c0Star = stocApprx.Get_c0();
                }
            }
        }

        public void ExportResultsToCSV(string filename)
        {
            foreach (StochasticApproximation stocApprx in stochasticApproximations)
            {
                stocApprx.ExportResultsToCSV(filename + stocApprx.Get_a0_b_c0() + ".csv");
            }
        }
    }
}
