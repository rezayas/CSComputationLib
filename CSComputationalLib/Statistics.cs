using System;
using System.Collections.Generic;
using System.Text;
using MathNet.Numerics.Distributions;

namespace ComputationLib
{
    public class ObservationBasedStatistics
    {
        string _name;
        double _mean, _variance, _min, _max, _total;
        long _count;
        double _sumOfVarNominator;

        double[] _observations;
        int _numOfObservationsToStore;
        bool _ifStoreObservations = false;

        public ObservationBasedStatistics(string name)
        {
            _name = name;
            _count = 0;
            _mean = 0;
            _variance = 0;
            _min = double.MaxValue;
            _max = double.MinValue;
            _total = 0;
            _sumOfVarNominator = 0;
        }
        public ObservationBasedStatistics(string name, int numOfObservationsToStore)
        {
            _name = name;
            _count = 0;
            _mean = 0;
            _variance = 0;
            _min = double.MaxValue;
            _max = double.MinValue;
            _total = 0;
            _sumOfVarNominator = 0;

            _ifStoreObservations = true;
            _numOfObservationsToStore = numOfObservationsToStore;
            _observations = new double[numOfObservationsToStore];
        }

        public string Name
        {
            get { return _name; }
        }
        public double[] Observations
        {
            get { return _observations; }
        }
        public long NumOfObservations
        {
            get { return _count; }
        }
        public double Mean
        {
            get { return _mean; }
        }
        public double Max
        {
            get { return _max; }
        }
        public double Min
        {
            get { return _min; }
        }
        public double Variance
        {
            get { return _variance; }
        }
        public double StDev
        {
            get { return Math.Sqrt(_variance); }
        }
        public double Total
        {
            get { return _total; }
        }
        public double StErr
        {
            get
            {
                if (_count > 1)
                    return Math.Sqrt(_variance / _count);
                else
                    return 0;
            }
        }
        /// <summary>
        /// 95% confidence interval reflects a significance level of 0.05
        /// </summary>
        /// <param name="significanceLevel"></param>
        /// <returns></returns>
        public double HalfWidth(double significanceLevel)
        {
            if (_count <= 1)
                return 0;

            double coeff = StudentT.InvCDF(0, 1, _count - 1, significanceLevel);
            //double coefficient = StatisticalFunctions.TInv_TwoTails(significanceLevel, (int)_count - 1);
            return coeff * this.StErr;
        }
        /// <summary>
        /// 95% confidence interval reflects a significance level of 0.05
        /// </summary>
        /// <param name="significanceLevel"></param>
        /// <returns></returns>
        public double UBoundConfidenceInterval(double significanceLevel)
        {
            if (_count <= 1)
                return 0;
            return _mean + this.HalfWidth(significanceLevel);
        }
        /// <summary>
        /// 95% confidence interval reflects a significance level of 0.05
        /// </summary>
        /// <param name="significanceLevel"></param>
        /// <returns></returns>
        public double LBoundConfidenceInterval(double significanceLevel)
        {
            if (_count <= 1)
                return 0;

            return _mean - this.HalfWidth(significanceLevel);
        }
        /// <summary>
        /// 95% confidence interval reflects a significance level of 0.05
        /// </summary>
        /// <param name="value"></param>
        /// <param name="significanceLevel"></param>
        /// <returns></returns>
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
            this.Record(obs, _count - 1);
        }
        public void Record(double obs, long locationIndex)
        {
            if (double.IsNaN(obs))
                return;

            double inc = obs - _mean;
            ++_count;
            _mean += inc / _count; // incremental change in mean
            _sumOfVarNominator += (_count - 1) * inc * (inc / _count); // running variance numerator
            _total += obs;

            if (_count > 1) _variance = _sumOfVarNominator / (_count - 1);
            if (obs < _min) _min = obs;
            if (obs > _max) _max = obs;

            if (_ifStoreObservations)
                _observations[locationIndex] = obs;
        }

        public void Reset()
        {
            _count = 0;
            _mean = 0;
            _variance = 0;
            _min = double.MaxValue;
            _max = double.MinValue;
            _total = 0;
            _sumOfVarNominator = 0;

            if (_ifStoreObservations)
                _observations = new double[_numOfObservationsToStore];
        }
    }

    public class TimePersistentStatistics
    {
        double _mean, _variance, _minimum, _maximum;
        long _count;
        double _baseTime;
        double _obsPreviousValue, _obsPreviousTime, _obsValue;
        double tot;

        // Instantiation
        public TimePersistentStatistics()
        {
            _count = 0;
            _baseTime = 0;
            _mean = 0.0;
            _variance = 0.0;
            _minimum = Double.MaxValue;
            _maximum = Double.MinValue;
            _obsPreviousValue = 0;
            _obsPreviousTime = 0;
            _obsValue = 0;
        }

        // Properties
        public long NumOfObservations
        {
            get { return _count; }
        }
        public double Mean
        {
            get { return _mean; }
        }
        public double Max
        {
            get { return _maximum; }
        }
        public double Min
        {
            get { return _minimum; }
        }
        public double Variance
        {
            get { return _variance; }
        }
        public double StDev
        {
            get { return Math.Sqrt(_variance); }
        }

        // Methods
        public void Reset(double time)
        {
            _count = 0;
            _baseTime = time;
            _mean = 0.0;
            _variance = 0.0;
            _minimum = Double.MaxValue;
            _maximum = Double.MinValue;
            _obsPreviousValue = 0;
            _obsPreviousTime = _baseTime;
            _obsValue = 0;
        }

        public void Record(double obsTime, double obsValue)//, double interestRate, int numOfDiscountingPeriods)
        {
            if (double.IsNaN(obsValue))
                return;

            ++_count;
            if (obsTime == _baseTime)
            {
                _obsPreviousTime = obsTime;
                _obsPreviousValue = obsValue;
            }
            else if (obsTime > _baseTime)
            {
                // **** check the calculation of variance *****
                tot += (_obsPreviousTime - _baseTime) * (obsTime - _obsPreviousTime) * Math.Pow(_obsValue - _mean, 2) / (obsTime - _baseTime);
                _mean = (_mean * (_obsPreviousTime - _baseTime) + (obsTime - _obsPreviousTime) * _obsPreviousValue) / (obsTime - _baseTime);
                _obsPreviousTime = obsTime;
                _obsPreviousValue = obsValue;
                if (_obsValue < _minimum) _minimum = _obsValue;
                if (_obsValue > _maximum) _maximum = _obsValue;
                if (_count > 1) _variance = tot / (_obsPreviousTime - _baseTime);
            }
        }
    }
}
