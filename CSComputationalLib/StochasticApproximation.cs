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

        /// <param name="x"></param>
        /// <param name="ifResampleSeeds"> set this to true if a new seed should be chosen and used for this replication </param>
        /// <returns> a sample for f at this x </returns>
        public virtual double Sample_f(Vector<double> x, bool ifResampleSeeds = true) { return 0; }

        /// <param name="x"></param>
        /// <param name="derivative_step"></param>
        /// <returns> a sample for derivative at this x </returns>
        public virtual Vector<double> Sample_Df(Vector<double> x, double derivative_step) { return null; }

        /// <summary>
        /// reset the seed of the simulation at the iteration 0 of the stochastic approximation algorithm 
        /// </summary>
        public virtual void ResetSeedAtItr0() { }

        /// <summary>
        /// Samples f and Df at x. 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="derivative_step"></param>
        /// <param name="xScale"> detremines the relative movement toward each direction </param>
        /// <param name="ifResampleSeeds">set this to true if a new seed should be chosen and used for this replication</param>
        public virtual void Sample_f_and_Df(
            Vector<double> x, 
            double derivative_step, 
            Vector<double> xScale = null,
            bool ifResampleSeeds = true) { }

        /// <returns> the current sampled value of f </returns>
        public virtual double Get_f() { return 0; }
        /// <returns> the current sampled value of Df </returns>
        public virtual Vector<double> Get_Df() { return null; }
        
    }

    public class TestBedX2Y2XY : SimModel
    {
        // x^2 + (1000y)^2 + x*(1000*y)

        RandomVariateLib.Normal _err; // a normally distributed error term
        RandomVariateLib.RNG _rnd;
        int _currentSeed = 0;

        public TestBedX2Y2XY(double errorVar)
        {
            _err = new RandomVariateLib.Normal("Error term", 0, errorVar);
            _rnd = new RandomVariateLib.RNG(_currentSeed);
        }

        public override double Sample_f(Vector<double> x, bool ifResampleSeeds)
        {
            if (ifResampleSeeds)
            {
                ++_currentSeed;
                _rnd = new RandomVariateLib.RNG(_currentSeed);
            }

            return Math.Pow(x[0], 2) + Math.Pow(x[1]*1000, 2) + x[0]*x[1]*1000 + _err.SampleContinuous(_rnd);
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
        // stepsize for derivative step
        // step_n = c0 * (n+1)^(-1/4) for n >= 0, a0 > 0, and b >= 1

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
        public List<Vector<double>> Itr_dx_over_x { get; private set; } = new List<Vector<double>>();
        public List<double> Itr_f { get; private set; } = new List<double>();
        public List<Vector<double>> Itr_Df { get; private set; } = new List<Vector<double>>();
        public List<double> Itr_step_Df { get; private set; } = new List<double>();
        public List<double> Itr_step_GH { get; private set; } = new List<double>();

        public Vector<double> xStar { get; private set; }
        public Vector<double> dx_over_x_ave { get; private set; }
        public double fStar { get; private set; }

        /// <param name="simModel"> the simulation model </param>
        /// <param name="stepSize_GH"> generalized harmonic step size </param>
        /// <param name="stepSize_Df"> step size for derivative steps </param>
        public StochasticApproximation(SimModel simModel, StepSize_GH stepSize_GH, StepSize_Df stepSize_Df)
        {
            _simModel = simModel;
            _stepSize_GH = stepSize_GH;
            _stepSize_Df = stepSize_Df;
        }

        /// <summary>
        /// finds a minimum of a stochastic function 
        /// </summary>
        /// <param name="nItrs"> number of iterations </param>
        /// <param name="nLastItrsToAve"> number of last iterations to calcualte the average of f </param>
        /// <param name="x0"> initial value of x </param>
        /// <param name="xScale"> relative scale of x variables (the defauls is [1, 1, 1 ...]) </param>
        /// <param name="modelProvidesDerivatives"> set to true if the derivative at each point is provided by the simulation model
        /// instead of being calculated by the algorithm </param>
        /// <param name="ifTwoSidedDerivative"> set to true if two-sided derivative (instead of one-sided) should be used </param>
        public void Minimize(int nItrs, int nLastItrsToAve, Vector<double> x0, Vector<double> xScale = null,
            bool modelProvidesDerivatives = false, bool ifTwoSidedDerivative = true)
        {
            double f;
            Vector<double> dx_over_x = Vector<double>.Build.Dense(x0.Count());

            // if x-scale is not provided, use [1, 1, 1, ...] as scale
            if (xScale is null)
                xScale = Vector<double>.Build.Dense(length: x0.Count, value: 1);

            // reset seed of the simulation model at iteration 0
            // note that this method could be empty if there is no need to reset the seed 
            _simModel.ResetSeedAtItr0();

            // set current x to x0
            Vector<double> x = x0;                     

            // iterations of the algorithm
            for (int itr = 0; itr < nItrs; itr++)
            {
                // get the current size of derivative step 
                double step_Df = _stepSize_Df.GetValue(itr);
                
                // initialize Df to a 0 vector
                Vector<double> Df = Vector<double>.Build.Dense(x0.Count());

                // calcualte derivative of f at x
                // if the model provides the derivative 
                if (modelProvidesDerivatives)
                {
                    // sample f and Df 
                    _simModel.Sample_f_and_Df(
                        x: x, 
                        derivative_step: step_Df, 
                        xScale: xScale, 
                        ifResampleSeeds: true);

                    // get f(x)
                    f = _simModel.Get_f();

                    // get the derivative from the model
                    Df = _simModel.Get_Df();                    
                }
                else 
                {
                    // if derivatives are not provided by the model we calculate it here
                    // for one sided Df(x) = (f(x+ε) - f(x))/ε
                    // for two sided Df(x) = (f(x+ε) - f(x-ε))/2ε

                    // build epsilon matrix
                    Matrix<double> epsilonMatrix = Matrix<double>.Build.DenseDiagonal(x0.Count(), step_Df);

                    // sample f(x)
                    f = _simModel.Sample_f(x, ifResampleSeeds: true);

                    // calcualte Df
                    for (int i = 0; i < x0.Count(); i++)
                    {
                        if (ifTwoSidedDerivative)
                        {                         
                            Df[i] =(
                                _simModel.Sample_f(x + epsilonMatrix.Row(i) * xScale[i], ifResampleSeeds: false) -
                                _simModel.Sample_f(x - epsilonMatrix.Row(i) * xScale[i], ifResampleSeeds: false)
                                ) / (2 * step_Df * xScale[i]);
                        }
                        else
                        {
                            Df[i] =
                                (_simModel.Sample_f(x + epsilonMatrix.Row(i) * xScale[i], ifResampleSeeds: false) - f)
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
                // dx over x                
                if (itr > 0)
                    dx_over_x = (x - Itr_x.Last()) / x;
                Itr_dx_over_x.Add(dx_over_x.PointwiseAbs());
                Itr_x.Add(x);
                Itr_f.Add(f);
                Itr_Df.Add(nDf);
                Itr_step_Df.Add(step_Df);
                Itr_step_GH.Add(step_GH);                

                // find a new x: x_new = x - stepSize*f'(x)*scale
                Vector<double> tempX = Vector<double>.Build.Dense(x.Count);
                for (int i = 0; i < x.Count; i++)
                    tempX[i] = x[i] - step_GH * nDf[i] * xScale[i];

                x = tempX;
            }

            // optimal x and optimal value of f is the average of last nLastItrsToAve iterations  
            double fSum = 0;
            Vector<double> xSum = Vector<double>.Build.Dense(x0.Count);
            Vector<double> dx_over_x_Sum = Vector<double>.Build.Dense(x0.Count);
            for (int itr = nItrs; itr > nItrs - nLastItrsToAve; itr--)
            {
                fSum += Itr_f[itr - 1];
                xSum += Itr_x[itr - 1];
                dx_over_x_Sum += Itr_dx_over_x[itr - 1];
            }
            xStar = xSum / nLastItrsToAve;
            dx_over_x_ave = dx_over_x_Sum / nLastItrsToAve;
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

    public class MultipleStochasticApproximation
    {
        List<StochasticApproximation> listStochApprox = new List<StochasticApproximation>();
        public double fStar { get; private set; } = double.MaxValue;
        public Vector<double> xStar { get; private set; }
        public Vector<double> dxOverXAveStar { get; private set; }
        public double a0Star { get; private set; } = double.NaN;
        public double bStar { get; private set; } = double.NaN;
        public double c0Star { get; private set; } = double.NaN;

        /// <summary>
        /// creates multiple stochastic approximation algorithsm
        /// </summary>
        /// <param name="simModels"> list of simulation models, one for each combination of optimization parameters (a0, b, c) </param>
        /// <param name="stepSizeGH_a0s">list of a0's for the generalized harmonic (GH) stepsize </param>
        /// <param name="stepSizeGH_bs">list of bs for generalized harmonic (GH) stepsize </param>
        /// <param name="stepSizeDf_cs">list of cs for the derivative step size </param>
        public MultipleStochasticApproximation(List<SimModel> simModels, double[] stepSizeGH_a0s, double[] stepSizeGH_bs, double[] stepSizeDf_cs)
        {
            // build the algorithms of stochastic approximation
            int i = 0;
            foreach (double a0 in stepSizeGH_a0s)
            {
                foreach (double b in stepSizeGH_bs)
                {
                    foreach (double c in stepSizeDf_cs)
                    {
                        listStochApprox.Add(
                            new StochasticApproximation(
                                simModel: simModels[i++],
                                stepSize_GH: new StepSize_GH(a0, b),
                                stepSize_Df: new StepSize_Df(c))
                                );
                    }
                }
            }
        }

        /// <summary>
        /// minimize the collection of f functions
        /// </summary>
        /// <param name="nItrs"></param>
        /// <param name="nLastItrsToAve"></param>
        /// <param name="x0"></param>
        /// <param name="xScale"></param>
        /// <param name="modelProvidesDerivatives"></param>
        /// <param name="ifTwoSidedDerivative"></param>
        /// <param name="ifParallel"> if the collection of stochastic approximation algorithms should be run in parallel </param>
        public void Minimize(int nItrs, int nLastItrsToAve, Vector<double> x0, Vector<double> xScale = null,
            bool modelProvidesDerivatives = false, bool ifTwoSidedDerivative = true, bool ifParallel = true)
        {
            // run all stochastic approximation algorithsm 
            if (ifParallel && listStochApprox.Count > 1)
            {
                Parallel.ForEach(listStochApprox, stocApprx =>
                {
                    stocApprx.Minimize(
                        nItrs: nItrs, 
                        nLastItrsToAve: nLastItrsToAve, 
                        x0: x0, 
                        xScale: xScale, 
                        modelProvidesDerivatives: modelProvidesDerivatives, 
                        ifTwoSidedDerivative: ifTwoSidedDerivative);
                });
            }
            else
            {
                foreach (StochasticApproximation stocApprx in listStochApprox)
                {
                    stocApprx.Minimize(
                        nItrs: nItrs,
                        nLastItrsToAve: nLastItrsToAve,
                        x0: x0,
                        xScale: xScale,
                        modelProvidesDerivatives: modelProvidesDerivatives,
                        ifTwoSidedDerivative: ifTwoSidedDerivative);
                }
            }

            // find the optimizer
            // find the a value that minimizes f
            xStar = Vector<double>.Build.Dense(x0.Count);
            fStar = double.MaxValue;
            foreach (StochasticApproximation stocApprx in listStochApprox)
            {
                // if this led to the minimum f
                if (stocApprx.fStar < fStar)
                {
                    fStar = stocApprx.fStar;
                    xStar = stocApprx.xStar;
                    dxOverXAveStar = stocApprx.dx_over_x_ave;
                    a0Star = stocApprx.Get_a0();
                    bStar = stocApprx.Get_b();
                    c0Star = stocApprx.Get_c0();
                }
            }
        }

        public void ExportResultsToCSV(string filename)
        {
            foreach (StochasticApproximation stocApprx in listStochApprox)
            {
                stocApprx.ExportResultsToCSV(filename + stocApprx.Get_a0_b_c0() + ".csv");
            }
        }
    }
}
