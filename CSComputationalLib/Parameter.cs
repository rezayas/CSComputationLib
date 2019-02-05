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
            ComorbidityDisutility = 9,
        }

        public int ID { get; }
        public string Name { get; }
        public double Value { get; protected set; } = 0;
        public bool IncludedInCalibration { get ; set ; }
        public bool ShouldBeUpdatedByTime { get ; set ; }
        public EnumType Type { get; protected set; }

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
        public double Sample(RNG rng)
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
        public double Sample()
        {
            Value = _slope * _depedentPar.Value + _intercept;
            return Value;
        }
    } 

    public class LinearCombination : Parameter
    {
        private double[] _arrCoefficients;

        // Instantiation
        public LinearCombination(int ID, string name, int[] parIDs, double[] coefficients) 
            : base(ID,name)
        {
            Type = EnumType.LinearCombination;
            this.ParIDs = (int[])parIDs.Clone();
            _arrCoefficients = (double[])coefficients.Clone();
        }

        // Properties 
        public int[] ParIDs { get; }

        // sample this parameters
        public double Sample(double[] valuesOfParams)
        {
            Value =0;
            for (int i = 0; i < valuesOfParams.Length; i++)
                Value += _arrCoefficients[i] * valuesOfParams[i];
            return Value;
        }

    }

    public class ProductParameter : Parameter
    {
        // Instantiation
        public ProductParameter(int ID, string name, int[] parIDs)
            : base(ID, name)
        {
            Type = EnumType.Product;
            this.ParIDs = (int[])parIDs.Clone();
        }

        // Properties 
        public int[] ParIDs { get; }

        // sample this parameters
        public double Sample(double[] valueOfParams)
        {
            Value = 1;
            for (int i = 0; i < valueOfParams.Length; i++)
                Value *= valueOfParams[i];
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
        public double Sample()
        {
            if (_inverseFirstParameter)
                Value = _par2.Value / _par1.Value;                
            else
                Value = _par1.Value * _par2.Value;                
            return Value; ;
        }
    } 

    public class TimeDependetLinear : Parameter
    {
        public int InterceptParID { get; }
        public int SlopeParID { get; }
        public double TimeOn { get; }
        public double TimeOff { get; }

        // Instantiation 
        public TimeDependetLinear(int ID, string name, int interceptParID, int slopeParID, double timeOn, double timeOff)
            : base(ID, name)
        {
            Type = EnumType.TimeDependentLinear;
            ShouldBeUpdatedByTime = true;

            InterceptParID = interceptParID;
            SlopeParID = slopeParID;
            TimeOn = timeOn;
            TimeOff = timeOff;
        }

        // sample this parameter
        public double Sample(double time, double intercept, double slope, double timeOn, double timeOff)
        {
            if (time < timeOn || time >= timeOff)
                Value = 0;
            else
                Value = intercept + slope * time;

            return Value;
        }
    }

    public class TimeDependetOscillating : Parameter
    {
        public int a0ParID { get; }
        public int a1ParID { get; }
        public int a2ParID { get; }
        public int a3ParID { get; }

        // Instantiation 
        public TimeDependetOscillating(int ID, string name, int a0ParID, int a1ParID, int a2ParID, int a3ParID)
            : base(ID, name)
        {
            Type = EnumType.TimeDependentOscillating;
            ShouldBeUpdatedByTime = true;

            this.a0ParID = a0ParID;
            this.a1ParID = a1ParID;
            this.a2ParID = a2ParID;
            this.a3ParID = a3ParID;
        }

        // sample this parameter
        public double Sample(double time, double a0, double a1, double a2, double a3)
        {
            Value = a0 + a1 * Math.Cos((time+a2)*2*Math.PI/a3);
            return Value; 
        }
    }

    public class TimeDependentExponential : Parameter
    {
        public int minParID { get; }
        public int maxParID { get; }
        public int bParID { get; }
        public int tStartParID { get; }

        // Instantiation 
        public TimeDependentExponential(int ID, string name, int bParID, int minParID, int maxParID, int tStartParID)
            : base(ID, name)
        {
            Type = EnumType.TimeDependentExponential;
            ShouldBeUpdatedByTime = true;

            this.minParID = minParID;
            this.maxParID = maxParID;
            this.bParID = bParID;
            this.tStartParID = tStartParID;
        }

        // sample this parameter
        public double Sample(double time, double b, double min, double max, double tStart)
        {
            if (time >= tStart)
                Value = max - (max - min) * Math.Exp(-b * (time - tStart));
            else
                Value = min;
            return Value;
        }
    }

    public class ComorbidityDisutility : Parameter
    {
        public int Par1ID { get; }
        public int Par2ID { get; }

        public ComorbidityDisutility(int id, string name, int par1ID, int par2ID): base(id, name)
        {
            Type = EnumType.ComorbidityDisutility;
            Par1ID = par1ID;
            Par2ID = par2ID;
        }

        public double Sample(double valDistulity1, double valDistulity2)
        {
            Value = 1 - (1 - valDistulity1) * (1 - valDistulity2);
            return Value;
        }
    }

}
