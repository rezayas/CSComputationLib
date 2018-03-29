using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using RandomVariateLib;

namespace SimulationLib
{
    public abstract class Parameter
    {
        protected int _ID;
        protected string _name;
        protected double _defaultValue;
        protected double _value;
        protected bool _includedInCalibration;
        protected bool _shouldBeUpdatedByTime = false;
        protected EnumType _type; 
        
        public enum EnumType
        {
            Independet = 1,
            Correlated = 2,
            LinearCombination = 3,
            MultipleCombination = 4,
            Multiplicative = 5,
            TimeDependetLinear = 6,
            TimeDependetOscillating = 7,
        }
 
        public Parameter(int ID, string name, double defaultValue)
        {
            _ID = ID;
            _name = name;
            _defaultValue = defaultValue;
            _value = _defaultValue;
        }

        public int ID
        { get { return _ID; } }

        public string Name
        { get { return _name; } }

        public bool IncludedInCalibration
        { 
            get { return _includedInCalibration; }
            set { _includedInCalibration = value; }        
        }
        public double Value
        {get {return _value; } } 
        public EnumType Type
        { get { return _type; } }
        public bool ShouldBeUpdatedByTime
        { 
            get { return _shouldBeUpdatedByTime; }
            set { _shouldBeUpdatedByTime = value; }
        }
    } // end of Parameter class

    public class IndependetParameter : Parameter
    {
        private double _par1, _par2, _par3, _par4;        
        private EnumRandomVariates _enumRandomVariateGenerator;
        private RVG _RVG = null;
               
        /// <summary>
        /// Return a random variate generators
        /// </summary>
        /// <param name="par1">Constant -> constant value; Uniform -> Min; Normal -> Mean.</param>
        /// <param name="par2">Uniform -> Max ; Normal -> StDev.</param>
        /// <param name="par3">.</param>
        /// <param name="par4">.</param>
        /// <returns> A token class representing the unit of work.</returns>
        public IndependetParameter(int ID, string name, double defaultValue, EnumRandomVariates enumRandomVariateGenerator, double par1, double par2, double par3, double par4)
            : base(ID, name, defaultValue)
        {
            _enumRandomVariateGenerator = enumRandomVariateGenerator;
            _par1 = par1;
            _par2 = par2;
            _par3 = par3;
            _par4 = par4;
            _type = EnumType.Independet;

            _RVG = SupportProcedures.ReturnARandomVariateGenerator(enumRandomVariateGenerator, name, par1, par2, par3, par4);
        }        
        
        // Properties
        public EnumRandomVariates Distrbution
        { get { return _enumRandomVariateGenerator; } }

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

    } // end of IndependetParameter

    public class CorrelatedParameter : Parameter
    {
        int _idOfDepedentPar;
        double _slope, _intercept;

        // Instantiation
        public CorrelatedParameter(int ID, string name, double defaultValue, int idOfDepedentPar, double slope, double intercept) 
            : base(ID, name,defaultValue)
        {
            _type = EnumType.Correlated;
            _idOfDepedentPar = idOfDepedentPar;
            _slope = slope;
            _intercept = intercept;
        }
        // Properties
        public int IDOfDepedentPar
        { get { return _idOfDepedentPar; } }

        // sample this parameter
        public double Sample(double valueOfDependentPar)
        {
            _value = _slope * valueOfDependentPar + _intercept;
            return _value;
        }

    } // end of IndependetParameter

    public class LinearCombination : Parameter
    {
        int[] _arrParIDs;
        double[] _arrCoefficients;

        // Instantiation
        public LinearCombination(int ID, string name, double defaultValue, int[] arrParIDs, double[] arrCoefficients) 
            : base(ID,name,defaultValue)
        {
            _type = EnumType.LinearCombination;
            _arrParIDs = (int[])arrParIDs.Clone();
            _arrCoefficients = (double[])arrCoefficients.Clone();
        }

        // Properties 
        public int[] arrParIDs
        { get { return _arrParIDs; } }

        // sample this parameters
        public double Sample(double[] arrValuesOfParameters)
        {
            _value =0;
            for (int i = 0; i < arrValuesOfParameters.Length; i++)
                _value += _arrCoefficients[i] * arrValuesOfParameters[i];
            return _value;
        }

    }

    public class MultipleCombination : Parameter
    {
        int[] _arrParIDs;

        // Instantiation
        public MultipleCombination(int ID, string name, double defaultValue, int[] arrParIDs)
            : base(ID, name, defaultValue)
        {
            _type = EnumType.MultipleCombination;
            _arrParIDs = (int[])arrParIDs.Clone();
        }

        // Properties 
        public int[] arrParIDs
        { get { return _arrParIDs; } }

        // sample this parameters
        public double Sample(double[] arrValueOfParameters)
        {
            _value = 1;
            for (int i = 0; i < arrValueOfParameters.Length; i++)
                _value *= arrValueOfParameters[i];
            return _value;
        }
    }

    public class MultiplicativeParameter : Parameter
    {
        int _firstParameterID, _secondParameterID;
        bool _inverseFirstParameter;

        // Instantiation
        public MultiplicativeParameter(int ID, string name, double defaultValue, int firstParameterID, int secondParameterID, bool ifInverseFirstParameter)
            : base(ID, name, defaultValue)
        {
            _type = EnumType.Multiplicative;
            _firstParameterID = firstParameterID;
            _secondParameterID = secondParameterID;
            _inverseFirstParameter = ifInverseFirstParameter;
        }
        // Properties 
        public int FirstParameterID
        { get { return _firstParameterID; } }
        public int SecondParameterID
        { get { return _secondParameterID; } }

        // sample this parameter
        public double Sample(double valueOfFirstParameter, double valueOfSecondParameter)
        {
            if (_inverseFirstParameter)
                _value = valueOfSecondParameter / valueOfFirstParameter;                
            else
                _value = valueOfFirstParameter * valueOfSecondParameter;                
            return _value; ;
        }

    } // end of MultiplicativeParameter

    public class TimeDependetLinear : Parameter
    {
        int _interceptParID, _slopeParID;
        double _timeOn, _timeOff;

        public int InterceptParID
        { get { return _interceptParID; } }
        public int SlopeParID
        { get { return _slopeParID; } }
        public double TimeOn
        { get { return _timeOn; } }
        public double TimeOff
        { get { return _timeOff; } }

        // Instantiation 
        public TimeDependetLinear(int ID, string name, double defaultValue, int interceptParID, int slopeParID, double timeOn, double timeOff)
            : base(ID, name, defaultValue)
        {
            _type = EnumType.TimeDependetLinear;
            _shouldBeUpdatedByTime = true;

            _interceptParID = interceptParID;
            _slopeParID = slopeParID;
            _timeOn = timeOn;
            _timeOff = timeOff;
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

    }// end of TimeDependetLinear

    public class TimeDependetOscillating : Parameter
    {
        int _a0ParID, _a1ParID, _a2ParID, _a3ParID;

        public int a0ParID
        { get { return _a0ParID; } }
        public int a1ParID
        { get { return _a1ParID; } }
        public int a2ParID
        { get { return _a2ParID; } }
        public int a3ParID
        { get { return _a3ParID; } }

        // Instantiation 
        public TimeDependetOscillating(int ID, string name, double defaultValue, int a0ParID, int a1ParID, int a2ParID, int a3ParID)
            : base(ID, name, defaultValue)
        {
            _type = EnumType.TimeDependetOscillating;
            _shouldBeUpdatedByTime = true;

            _a0ParID = a0ParID;
            _a1ParID = a1ParID;
            _a2ParID = a2ParID;
            _a3ParID = a3ParID;
        }

        // sample this parameter
        public double Sample(double time, double a0, double a1, double a2, double a3)
        {
            _value = a0 + a1 * Math.Cos((time+a2)*2*Math.PI/a3);
            return _value; 
        }

    }// end of TimeDependetOscillating
}
