using System;
using System.Collections.Generic;
using System.Text;
using MathNet.Numerics.Distributions;

namespace ComputationLib
{
    public class ObsBasedStat
    {
        double _sumOfVarNuminator;
        int _numOfObservationsToStore;
        bool _ifStoreObservations = false;
        
        public ObsBasedStat(string name, int numOfObservationsToStore = 0)
        {
            Name = name;
            NumOfRecordings = 0;
            Mean = 0;
            Variance = 0;
            Min = double.MaxValue;
            Max = double.MinValue;
            Total = 0;
            _sumOfVarNuminator = 0;

            if (numOfObservationsToStore>0)
                _ifStoreObservations = true;
            _numOfObservationsToStore = numOfObservationsToStore;
            Observations = new double[numOfObservationsToStore];
        }

        public string Name { get; }
        public double[] Observations { get; private set; }
        public long NumOfRecordings { get; private set; }
        public double Mean { get; private set; }
        public double Max { get; private set; }
        public double Min { get; private set; }
        public double Variance { get; private set; }
        public double StDev
        {
            get { return Math.Sqrt(Variance); }
        }
        public double Total { get; private set; }
        public double StErr
        {
            get
            {
                if (NumOfRecordings > 1)
                    return Math.Sqrt(Variance / NumOfRecordings);
                else
                    return 0;
            }
        }
        public double[] GetMeanStDevStErr()
        {
            double[] results = new double[3];
            results[0] = Mean;
            results[1] = StDev;
            results[2] = StErr;

            return results;
        }
        public double HalfWidth(double significanceLevel)
        {
            // 95% confidence interval reflects a significance level of 0.05

            if (NumOfRecordings <= 1)
                return 0;

            double coeff = StudentT.InvCDF(0, 1, NumOfRecordings - 1, significanceLevel);
            return coeff * StErr;
        }
        public double UBoundConfidenceInterval(double significanceLevel)
        {
            if (NumOfRecordings <= 1)
                return 0;
            return Mean + HalfWidth(significanceLevel);
        }
        public double LBoundConfidenceInterval(double significanceLevel)
        {
            if (NumOfRecordings <= 1)
                return 0;

            return Mean - HalfWidth(significanceLevel);
        }
        public bool IfStatisticallyDifferent(double nullHypotheis, double significanceLevel)
        {
            bool result = false;

            if (nullHypotheis < LBoundConfidenceInterval(significanceLevel) ||
                nullHypotheis > UBoundConfidenceInterval(significanceLevel))
                result = true;

            return result;
        }
        public void Record(double obs)
        {
            Record(obs, NumOfRecordings - 1);
        }
        public void Record(double obs, long obsLocationIndex)
        {
            if (double.IsNaN(obs))
                return;

            double inc = obs - Mean;
            ++NumOfRecordings;
            Mean += inc / NumOfRecordings; // incremental change in mean
            _sumOfVarNuminator += (NumOfRecordings - 1) * inc * (inc / NumOfRecordings); // running variance numerator
            Total += obs;

            if (NumOfRecordings > 1) Variance = _sumOfVarNuminator / (NumOfRecordings - 1);
            if (obs < Min) Min = obs;
            if (obs > Max) Max = obs;

            if (_ifStoreObservations)
                Observations[obsLocationIndex] = obs;
        }

        public void Reset()
        {
            NumOfRecordings = 0;
            Mean = 0;
            Variance = 0;
            Min = double.MaxValue;
            Max = double.MinValue;
            Total = 0;
            _sumOfVarNuminator = 0;

            if (_ifStoreObservations)
                Observations = new double[_numOfObservationsToStore];
        }
    }

    public class ContinuousTimeStat
    {
        double _baseTime;
        double _obsPreviousValue, _obsPreviousTime, _obsValue;
        double tot;

        // Instantiation
        public ContinuousTimeStat()
        {
            NumOfObservations = 0;
            _baseTime = 0;
            Mean = 0.0;
            Variance = 0.0;
            Min = Double.MaxValue;
            Max = Double.MinValue;
            _obsPreviousValue = 0;
            _obsPreviousTime = 0;
            _obsValue = 0;
        }

        // Properties
        public long NumOfObservations { get; private set; }
        public double Mean { get; private set; }
        public double Max { get; private set; }
        public double Min { get; private set; }
        public double Variance { get; private set; }
        public double StDev
        {
            get { return Math.Sqrt(Variance); }
        }

        // Methods
        public void Reset(double time)
        {
            NumOfObservations = 0;
            _baseTime = time;
            Mean = 0.0;
            Variance = 0.0;
            Min = Double.MaxValue;
            Max = Double.MinValue;
            _obsPreviousValue = 0;
            _obsPreviousTime = _baseTime;
            _obsValue = 0;
        }

        public void Record(double obsTime, double obsValue)//, double interestRate, int numOfDiscountingPeriods)
        {
            if (double.IsNaN(obsValue))
                return;

            ++NumOfObservations;
            if (obsTime == _baseTime)
            {
                _obsPreviousTime = obsTime;
                _obsPreviousValue = obsValue;
            }
            else if (obsTime > _baseTime)
            {
                // **** check the calculation of variance *****
                tot += (_obsPreviousTime - _baseTime) * (obsTime - _obsPreviousTime) * Math.Pow(_obsValue - Mean, 2) / (obsTime - _baseTime);
                Mean = (Mean * (_obsPreviousTime - _baseTime) + (obsTime - _obsPreviousTime) * _obsPreviousValue) / (obsTime - _baseTime);
                _obsPreviousTime = obsTime;
                _obsPreviousValue = obsValue;
                if (_obsValue < Min) Min = _obsValue;
                if (_obsValue > Max) Max = _obsValue;
                if (NumOfObservations > 1) Variance = tot / (_obsPreviousTime - _baseTime);
            }
        }
    }
}
