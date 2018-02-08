using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ComputationLib
{
    public abstract class cStepSizeRule
    {
        protected string _name;
        // Instantiation
        public cStepSizeRule(string name)
        {
            _name = name;
        }
        // step size
        public virtual double StepSize(long iteration)
        {
            return 0;
        }
        // discount rate corresponding to this step size
        public double ObservationDiscountRate(long iteration)
        {
            return 1 / (1 + StepSize(iteration));

            //double alpha_nMinus1, alpha_n;
            //alpha_nMinus1 = StepSize(iteration);
            //alpha_n = StepSize(iteration + 1);

            //return alpha_nMinus1 * (1-alpha_n)/alpha_n;            
        }
    }

    public class cConstantStepSize : cStepSizeRule
    {
        double _stepSize;
        // Instantiation
        public cConstantStepSize(string name, double stepSize)
            : base(name)
        {
            _stepSize = stepSize;
        }
        // step size
        public override double StepSize(long iteration)
        {            
            return _stepSize;
        }        
    }

    public class cHarmonicStepSize : cStepSizeRule
    {
        double _a;
        // Instantiation
        public cHarmonicStepSize(string name, double a)
            : base(name)
        {
            _a = a;
        }
        // step size
        public override double StepSize(long iteration)
        {
            if (iteration <= 0)
                return 1;

            if (_a <= 0) return 0;

            return _a/(_a + iteration - 1);
        }
    }

    public class cPolynomialStepSize : cStepSizeRule
    {
        private double _beta;
        // Instantiation
        public cPolynomialStepSize(string name, double beta)
            : base(name)
        {
            _beta = beta;
        }
        // step size
        public override double StepSize(long iteration)
        {
            if (iteration <= 0)
                return 1;

            return 1/Math.Pow(iteration, _beta);
        }
    }
}
