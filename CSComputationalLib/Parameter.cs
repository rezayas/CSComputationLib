using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using RandomVariateLib;

namespace ComputationLib
{
    public abstract class Parameter
    {
        public enum EnumType
        {
            Independet = 1,     // random variate generator
            Correlated = 2,
            LinearCombination = 3,
            Product = 4,
            Multiplicative = 5,
            TimeDependentLinear = 6,
            TimeDependentOscillating = 7,
            TimeDependentExponential = 8,
            TimeDependentSigmoid = 9,
            ComorbidityDisutility = 10,
        }

        public int ID { get; }
        public string Name { get; }
        public double Value { get; protected set; } = 0;
        public bool IncludedInCalibration { get ; set ; }
        public bool ShouldBeUpdatedByTime { get ; set ; }
        public EnumType Type { get; protected set; }
        public abstract double Sample(double time, RNG rng);

        public Parameter(int id, string name)
        {
            ID = id;
            Name = name;
        }

        public static EnumType FindParameterType(string type)
        {
            EnumType thisEnum;
            switch (type.ToLower())
            {
                case "correlated":
                    thisEnum = EnumType.Correlated;
                    break;
                case "linear combination":
                    thisEnum = EnumType.LinearCombination;
                    break;
                case "product":
                    thisEnum = EnumType.Product;
                    break;
                case "multiplicative":
                    thisEnum = EnumType.Multiplicative;
                    break;
                case "time-dependent linear":
                    thisEnum = EnumType.TimeDependentLinear;
                    break;
                case "time-dependent oscillating":
                    thisEnum = EnumType.TimeDependentOscillating;
                    break;
                case "time-dependent exponential":
                    thisEnum = EnumType.TimeDependentExponential;
                    break;
                case "time-dependent sigmoid":
                    thisEnum = EnumType.TimeDependentSigmoid;
                    break;
                case "comorbidity disutility":
                    thisEnum = EnumType.ComorbidityDisutility;
                    break;
                default:
                    thisEnum = EnumType.Independet;
                    break;
            }
            return thisEnum;
        }
    } 

    public class IndependetParameter : Parameter
    {
        // Properties
        public EnumRandomVariates EnumRandomVariate { get; }

        private double _par1, _par2, _par3, _par4;
        private RVG _RVG = null;
               
        public IndependetParameter(int ID, string name, EnumRandomVariates enumRandomVariateGenerator, double par1, double par2, double par3, double par4)
            : base(ID, name)
        {
            EnumRandomVariate = enumRandomVariateGenerator;
            _par1 = par1;
            _par2 = par2;
            _par3 = par3;
            _par4 = par4;
            Type = EnumType.Independet;

            _RVG = SupportProcedures.ReturnARandomVariateGenerator(enumRandomVariateGenerator, name, par1, par2, par3, par4);
        }

        // sample this parameter
        public override double Sample(double time, RNG rng)
        {
            Value = _RVG.SampleContinuous(rng);
            return Value;
        }

        public double[] GetEquallyDistancedPoints(int nOfPoints)
        {
            return _RVG.GetEquallyDistributedPoints(nOfPoints);
        }

    } 

    public class CorrelatedParameter : Parameter
    {
        private Parameter _depedentPar;
        private double _slope, _intercept;

        // Instantiation
        public CorrelatedParameter(int ID, string name, Parameter depedentPar, double slope, double intercept) 
            : base(ID, name)
        {
            Type = EnumType.Correlated;
            _depedentPar = depedentPar;
            _slope = slope;
            _intercept = intercept;
        }

        // sample this parameter
        public override double Sample(double time, RNG rng)
        {
            Value = _slope * _depedentPar.Value + _intercept;
            return Value;
        }
    } 

    public class LinearCombination : Parameter
    {
        private List<Parameter> _parameters;
        private double[] _arrCoefficients;

        // Instantiation
        public LinearCombination(int ID, string name, List<Parameter> parameters, double[] coefficients) 
            : base(ID,name)
        {
            Type = EnumType.LinearCombination;
            _parameters = parameters;
            _arrCoefficients = (double[])coefficients.Clone();
        }

        // sample this parameters
        public override double Sample(double time, RNG rng)
        {
            Value =0;
            for (int i = 0; i < _arrCoefficients.Length; i++)
                Value += _arrCoefficients[i] * _parameters[i].Value;
            return Value;
        }

    }

    public class ProductParameter : Parameter
    {
        private List<Parameter> _parameters;

        // Instantiation
        public ProductParameter(int ID, string name, List<Parameter> parameters)
            : base(ID, name)
        {
            Type = EnumType.Product;
            _parameters = parameters;
        }

        // sample this parameters
        public override double Sample(double time, RNG rng)
        {
            Value = 1;
            foreach (Parameter par in _parameters)
                Value *= par.Value;
            return Value;
        }
    }

    public class MultiplicativeParameter : Parameter
    {
        bool _inverseFirstParameter;
        private Parameter _par1;
        private Parameter _par2;

        // Instantiation
        public MultiplicativeParameter(int ID, string name, Parameter par1, Parameter par2, bool ifInversePar1)
            : base(ID, name)
        {
            Type = EnumType.Multiplicative;
            _par1 = par1;
            _par2 = par2;
            _inverseFirstParameter = ifInversePar1;
        }

        // sample this parameter
        public override double Sample(double time, RNG rng)
        {
            if (_inverseFirstParameter)
                Value = _par2.Value / _par1.Value;                
            else
                Value = _par1.Value * _par2.Value;                
            return Value; ;
        }
    }

    public class ComorbidityDisutility : Parameter
    {
        private Parameter _par1;
        private Parameter _par2;

        public ComorbidityDisutility(int id, string name, Parameter par1, Parameter par2)
            : base(id, name)
        {
            Type = EnumType.ComorbidityDisutility;
            _par1 = par1;
            _par2 = par2;
        }

        public override double Sample(double time, RNG rng)
        {
            Value = 1 - (1 - _par1.Value) * (1 - _par2.Value);
            return Value;
        }
    }

    public class TimeDependetLinear : Parameter
    {
        private Parameter _interceptPar { get; }
        private Parameter _slopePar { get; }
        private double _timeOn { get; }
        private double _timeOff { get; }

        // Instantiation 
        public TimeDependetLinear(int ID, string name, Parameter interceptPar, Parameter slopePar, double timeOn, double timeOff)
            : base(ID, name)
        {
            Type = EnumType.TimeDependentLinear;
            ShouldBeUpdatedByTime = true;

            _interceptPar = interceptPar;
            _slopePar = slopePar;
            _timeOn = timeOn;
            _timeOff = timeOff;
        }

        // sample this parameter
        public override double Sample(double time, RNG rng)
        {
            if (time < _timeOn || time >= _timeOff)
                Value = 0;
            else
                Value = _interceptPar.Value + _slopePar.Value * time;

            return Value;
        }
    }

    public class TimeDependetOscillating : Parameter
    {
        private Parameter _a0Par;
        private Parameter _a1Par;
        private Parameter _a2Par;
        private Parameter _a3Par;

        // Instantiation 
        public TimeDependetOscillating(int ID, string name, Parameter a0Par, Parameter a1Par, Parameter a2Par, Parameter a3Par)
            : base(ID, name)
        {
            Type = EnumType.TimeDependentOscillating;
            ShouldBeUpdatedByTime = true;

            _a0Par = a0Par;
            _a1Par = a1Par;
            _a2Par = a2Par;
            _a3Par = a3Par;
        }

        // sample this parameter
        public override double Sample(double time, RNG rng)
        {
            Value = _a0Par.Value + _a1Par.Value * Math.Cos((time+_a2Par.Value)*2*Math.PI/_a3Par.Value);
            return Value; 
        }
    }

    public class TimeDependentExponential : Parameter
    {
        private Parameter _minPar;
        private Parameter _maxPar;
        private Parameter _bPar;
        private Parameter _tStartPar;

        // Instantiation 
        public TimeDependentExponential(int ID, string name, Parameter bPar, Parameter minPar, Parameter maxPar, Parameter tStartPar)
            : base(ID, name)
        {
            Type = EnumType.TimeDependentExponential;
            ShouldBeUpdatedByTime = true;

            _minPar = minPar;
            _maxPar = maxPar;
            _bPar = bPar;
            _tStartPar = tStartPar;
        }

        // sample this parameter
        public override double Sample(double time, RNG rng)
        {
            if (time >= _tStartPar.Value)
                Value = _maxPar.Value - (_maxPar.Value- _minPar.Value) * Math.Exp(-_bPar.Value * (time - _tStartPar.Value));
            else
                Value = _minPar.Value;
            return Value;
        }
    }

    public class TimeDependentSigmoid : Parameter
    {
        private Parameter _minPar;
        private Parameter _bPar;
        private Parameter _tStartPar;

        // Instantiation 
        public TimeDependentSigmoid(int ID, string name, Parameter bPar, Parameter minPar, Parameter tStartPar)
            : base(ID, name)
        {
            Type = EnumType.TimeDependentSigmoid;
            ShouldBeUpdatedByTime = true;

            _minPar = minPar;
            _bPar = bPar;
            _tStartPar = tStartPar;
        }

        // sample this parameter
        public override double Sample(double time, RNG rng)
        {
            Value = _minPar.Value + (1-_minPar.Value)/(1+Math.Exp(-_bPar.Value*(time-_tStartPar.Value)));
            return Value;
        }
    }
}
