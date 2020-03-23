using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ComputationLib
{
    public abstract class StepSizeRule
    {
        public string Name { get; }

        // Instantiation
        public StepSizeRule(string name)
        {
            Name = name;
        }
        // step size
        public virtual double GetStepSize(long itr)
        {
            return 0;
        }
        // discount rate corresponding to this step size
        public double GetDiscountRate(long itr)
        {
            return 1 / (1 + GetStepSize(itr));         
        }
    }

    public class ConstantStepSize : StepSizeRule
    {
        private readonly double _stepSize;
        // Instantiation
        public ConstantStepSize(string name, double stepSize)
            : base(name)
        {
            _stepSize = stepSize;
        }
        // step size
        public override double GetStepSize(long itr)
        {            
            return _stepSize;
        }        
    }

    public class HarmonicStepSize : StepSizeRule
    {
        public double a { get; }

        // Instantiation
        public HarmonicStepSize(string name, double a)
            : base(name)
        {
            this.a = a;
        }
        // step size
        public override double GetStepSize(long itr)
        {
            if (itr <= 0)
                return 1;

            if (a <= 0) return 0;

            return a/(a + itr - 1);
        }
    }

    public class PolynomialStepSize : StepSizeRule
    {
        public double beta { get; }
        // Instantiation
        public PolynomialStepSize(string name, double beta)
            : base(name)
        {
            this.beta = beta;
        }
        // step size
        public override double GetStepSize(long itr)
        {
            if (itr <= 0)
                return 1;

            return 1/Math.Pow(itr, beta);
        }
    }
}
