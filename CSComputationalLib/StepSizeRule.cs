using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ComputationLib
{

    public abstract class ExplorationRule
    {
        // Instantiation
        public ExplorationRule()
        {
        }
        public abstract double GetEpsilon(long itr);
    }

    public class EpsilonGreedy : ExplorationRule
    {
        // select the greedy option with probability 1-epsilon, where 
        // epsilon = 1/n^beta, beta \in (0.5, 1]

        public double Beta { get; }
        public EpsilonGreedy(double beta): base()
        {
            Beta = beta;
        }
        // update epsilon greedy
        public override double GetEpsilon(long itr)
        {
            return Math.Pow(itr, -Beta);
        }
    }


    public abstract class LearningRule
    {
        // Instantiation
        public LearningRule()
        {
        }
        // step size
        public abstract double GetStepSize(long itr);

        // discount rate corresponding to this step size
        public double GetDiscountRate(long itr)
        {
            return 1 / (1 + GetStepSize(itr));         
        }
    }

    public class ConstantStepSize : LearningRule
    {
        private readonly double _stepSize;
        // Instantiation
        public ConstantStepSize(double stepSize)
            : base()
        {
            _stepSize = stepSize;
        }
        // step size
        public override double GetStepSize(long itr)
        {            
            return _stepSize;
        }        
    }

    public class HarmonicStepSize : LearningRule
    {
        public double a { get; }

        // Instantiation
        public HarmonicStepSize(double a)
            : base()
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

    public class PolynomialStepSize : LearningRule
    {
        public double beta { get; }
        // Instantiation
        public PolynomialStepSize(double beta)
            : base()
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
