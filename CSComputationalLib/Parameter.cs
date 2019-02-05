using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using RandomVariateLib;

namespace ComputationLib
{
    public abstract class Parameter
    {
        protected double _value = 0;
        protected EnumType _type;
        protected bool _includedInCalibration;
        protected bool _shouldBeUpdatedByTime;

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
        public double Value { get => _value; }
        public bool IncludedInCalibration { get => _includedInCalibration; set => _includedInCalibration = value; }
        public bool ShouldBeUpdatedByTime { get => _shouldBeUpdatedByTime; set => _shouldBeUpdatedByTime = value; }
        public EnumType Type { get => _type; }

        public Parameter(int id, string name)
        {
            ID = id;
            Name = name;
            _value = 0;
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
            _type = EnumType.Independet;

            _RVG = SupportProcedures.ReturnARandomVariateGenerator(enumRandomVariateGenerator, name, par1, par2, par3, par4);
        }

        // sample this parameter
        public double Sample(RNG rng)
        {
            _value = _RVG.SampleContinuous(rng);
            return _value;
        }

        public double[] GetEquallyDistancedPoints(int nOfPoints)
        {
            return _RVG.GetEquallyDistributedPoints(nOfPoints);
        }

    } 

    public class CorrelatedParameter : Parameter
    {
        public int IDOfDepedentPar { get; }
        private double _slope, _intercept;

        // Instantiation
        public CorrelatedParameter(int ID, string name, int idOfDepedentPar, double slope, double intercept) 
            : base(ID, name)
        {
            _type = EnumType.Correlated;
            IDOfDepedentPar = idOfDepedentPar;
            _slope = slope;
            _intercept = intercept;
        }

        // sample this parameter
        public double Sample(double valueOfDependentPar)
        {
            _value = _slope * valueOfDependentPar + _intercept;
            return _value;
        }

    } 

    public class LinearCombination : Parameter
    {
        private double[] _arrCoefficients;

        // Instantiation
        public LinearCombination(int ID, string name, int[] parIDs, double[] coefficients) 
            : base(ID,name)
        {
            _type = EnumType.LinearCombination;
            this.ParIDs = (int[])parIDs.Clone();
            _arrCoefficients = (double[])coefficients.Clone();
        }

        // Properties 
        public int[] ParIDs { get; }

        // sample this parameters
        public double Sample(double[] valuesOfParams)
        {
            _value =0;
            for (int i = 0; i < valuesOfParams.Length; i++)
                _value += _arrCoefficients[i] * valuesOfParams[i];
            return _value;
        }

    }

    public class ProductParameter : Parameter
    {
        // Instantiation
        public ProductParameter(int ID, string name, int[] parIDs)
            : base(ID, name)
        {
            _type = EnumType.Product;
            this.ParIDs = (int[])parIDs.Clone();
        }

        // Properties 
        public int[] ParIDs { get; }

        // sample this parameters
        public double Sample(double[] valueOfParams)
        {
            _value = 1;
            for (int i = 0; i < valueOfParams.Length; i++)
                _value *= valueOfParams[i];
            return _value;
        }
    }

    public class MultiplicativeParameter : Parameter
    {
        bool _inverseFirstParameter;

        // Instantiation
        public MultiplicativeParameter(int ID, string name, int par1ID, int par2ID, bool ifInversePar1)
            : base(ID, name)
        {
            _type = EnumType.Multiplicative;
            FirstParameterID = par1ID;
            SecondParameterID = par2ID;
            _inverseFirstParameter = ifInversePar1;
        }
        // Properties 
        public int FirstParameterID { get; }
        public int SecondParameterID { get; }

        // sample this parameter
        public double Sample(double valueOfPar1, double valueOfPar2)
        {
            if (_inverseFirstParameter)
                _value = valueOfPar2 / valueOfPar1;                
            else
                _value = valueOfPar1 * valueOfPar2;                
            return _value; ;
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
            _type = EnumType.TimeDependentLinear;
            _shouldBeUpdatedByTime = true;

            InterceptParID = interceptParID;
            SlopeParID = slopeParID;
            TimeOn = timeOn;
            TimeOff = timeOff;
        }

        // sample this parameter
        public double Sample(double time, double intercept, double slope, double timeOn, double timeOff)
        {
            if (time < timeOn || time >= timeOff)
                _value = 0;
            else
                _value = intercept + slope * time;

            return _value;
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
            _type = EnumType.TimeDependentOscillating;
            _shouldBeUpdatedByTime = true;

            this.a0ParID = a0ParID;
            this.a1ParID = a1ParID;
            this.a2ParID = a2ParID;
            this.a3ParID = a3ParID;
        }

        // sample this parameter
        public double Sample(double time, double a0, double a1, double a2, double a3)
        {
            _value = a0 + a1 * Math.Cos((time+a2)*2*Math.PI/a3);
            return _value; 
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
            _type = EnumType.TimeDependentExponential;
            _shouldBeUpdatedByTime = true;

            this.minParID = minParID;
            this.maxParID = maxParID;
            this.bParID = bParID;
            this.tStartParID = tStartParID;
        }

        // sample this parameter
        public double Sample(double time, double b, double min, double max, double tStart)
        {
            if (time >= tStart)
                _value = max - (max - min) * Math.Exp(-b * (time - tStart));
            else
                _value = min;
            return _value;
        }
    }

    public class ComorbidityDisutility : Parameter
    {
        public int Par1ID { get; }
        public int Par2ID { get; }

        public ComorbidityDisutility(int id, string name, int par1ID, int par2ID): base(id, name)
        {
            _type = EnumType.ComorbidityDisutility;
            Par1ID = par1ID;
            Par2ID = par2ID;
        }

        public double Sample(double valDistulity1, double valDistulity2)
        {
            _value = 1 - (1 - valDistulity1) * (1 - valDistulity2);
            return _value;
        }
    }

}
