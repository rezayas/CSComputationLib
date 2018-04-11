using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ComputationLib
{
    
    public class OldTimeSeries
    {
        // Variables
        #region Variables

        public enum enumPredictionModel
        {
            Nothing = 0,
            Linear = 1,
            Quadratic = 2,
        }

        // Fields
        string _name;
        double _currentAggregatedObsInThisObsPeriod;
        int _currentNumOfRecordingsInThisObsPeriod = 0;
        double[] _arrAggregatedObs;                
        // prediction
        enumPredictionModel _predictionModel;
        int _numOfObsPeriods;
        int _numOfRecodingsInEachObsPeriod;
        double[,] _matX;
        double[] _arrRow;
        #endregion

        // Instantiation       
       
        /// <summary>
        /// Creates time series (observations over observation periods will NOT be aggregated)
        /// </summary>
        /// <param name="numOfPastObsPeriodsToStore"> Specify the number of past observations to store </param>        
        public OldTimeSeries(string name, int numOfObsPeriods, enumPredictionModel predictionModel)
        {
            _name = name;
            _predictionModel = predictionModel;
            _numOfObsPeriods = numOfObsPeriods;
            _numOfRecodingsInEachObsPeriod = 1;

            // array to store observations
            _arrAggregatedObs = new double[numOfObsPeriods];
            // setup prediction
            SetupPrediction();
        }        
        /// <summary>
        /// Creates time series (observations over observation periods will be aggregated)
        /// </summary>
        /// <param name="name"></param>
        /// <param name="numOfPastObsPeriodsToStore">Specify the number of past observation periods to store</param>
        /// <param name="numOfObsInEachObsPeriod">Specify the number of data to be aggregated in each observation period</param>
        /// <param name="predictionModel"></param>
        public OldTimeSeries(string name, int numOfObsPeriods, int numOfRecodingsInEachObsPeriod, enumPredictionModel predictionModel)
        {
            _name = name;
            _predictionModel = predictionModel;
            _numOfObsPeriods = numOfObsPeriods;
            _numOfRecodingsInEachObsPeriod = numOfRecodingsInEachObsPeriod;

            // array to store aggregated observations
            _arrAggregatedObs = new double[numOfObsPeriods];
            // setup prediction
            SetupPrediction();
        }

        // Properties
        public string Name
        { 
            get{return _name;} 
        }
        public enumPredictionModel PredictionModel
        { 
            get{return _predictionModel;}
        }
        public int NumberOfObsPeriods
        { 
            get{return _numOfObsPeriods;}
        }
        public int NumOfRecodingsInEachObsPeriod
        {
            get{return _numOfRecodingsInEachObsPeriod;}
        }
        public double CurrentAggregatedObsInLastObsPeriod
        {
            get{return _currentAggregatedObsInThisObsPeriod;}
        }
        public double SumOfObservations
        {
            get { return _arrAggregatedObs.Sum(); }
        }

        // Methods
        // add an Obs
        public void AddAnObs(double obs)
        {
            // if there will be a single observation in each observation period
            if (_numOfRecodingsInEachObsPeriod == 1)
            {
                // store the current obs
                _currentAggregatedObsInThisObsPeriod = obs;
                if (_arrAggregatedObs != null)
                    SupportFunctions.AddToEndOfArrayFixedSize(ref _arrAggregatedObs, obs);
            }
            else // if (_numOfObsInEachObsPeriod > 1)
            {                 
                if (_currentNumOfRecordingsInThisObsPeriod == 0)
                    _currentAggregatedObsInThisObsPeriod = 0;

                // increment the counter for the number of observations in this observation period
                ++_currentNumOfRecordingsInThisObsPeriod;
                // increment the aggregate observations
                _currentAggregatedObsInThisObsPeriod += obs;

                // check aggregated observations should be stored
                if (_currentNumOfRecordingsInThisObsPeriod == _numOfRecodingsInEachObsPeriod)
                {                    
                    if (_arrAggregatedObs != null)
                        SupportFunctions.AddToEndOfArrayFixedSize(ref _arrAggregatedObs, _currentAggregatedObsInThisObsPeriod);
                    _currentNumOfRecordingsInThisObsPeriod = 0;
                }
            }
        }
        // prediction
        public double Prediction(int numOfObsPeriodsInFuture, bool integrateOverFutureObsPeriods)
        {
            double prediction = 0;

            // Select the prediction model
            switch (_predictionModel)
            {
                case enumPredictionModel.Linear:
                    #region enumPredictionModel.Linear
                    {
                        double beta0, beta1;
                        double xBar = _numOfObsPeriods / 2;
                        double yBar;
                        double nomin, denom;

                        yBar = _arrAggregatedObs.Average();
                        nomin = 0;
                        denom = 0;

                        for (int t = 0; t < _numOfObsPeriods; ++t)
                        {
                            nomin += (t - xBar) * (_arrAggregatedObs[t] - yBar);
                            denom += Math.Pow(t - xBar, 2);
                        }

                        // trend
                        beta1 = nomin / denom;
                        // intercept
                        beta0 = yBar - beta1 * xBar;

                        prediction = beta0 + beta1 * (_numOfObsPeriods + numOfObsPeriodsInFuture);
                    }
                    break;
                    #endregion
                case enumPredictionModel.Quadratic:
                    #region enumPredictionModel.Quadratic
                    {
                        // define a regression model
                        LeastSquares thisLS = new LeastSquares("Prediction");
                        thisLS.RunRegression(_matX, _arrAggregatedObs);

                        if (integrateOverFutureObsPeriods == false) // prediction should not be integrated over the prediction period
                        {
                            // create a new design row                        
                            _arrRow[0] = 1;
                            _arrRow[1] = numOfObsPeriodsInFuture + _numOfObsPeriods - 1;
                            _arrRow[2] = Math.Pow(_arrRow[1], 2);
                            // get the prediction
                            prediction = thisLS.yValue(_arrRow);
                        }
                        else // prediction integrated over the prediction period
                        {
                            for (int pointInFuture = 1; pointInFuture <= numOfObsPeriodsInFuture; ++pointInFuture)
                            {
                                // create a new design row                        
                                _arrRow[0] = 1;
                                _arrRow[1] = pointInFuture + _numOfObsPeriods - 1;
                                _arrRow[2] = Math.Pow(_arrRow[1], 2);
                                // get the prediction
                                prediction += thisLS.yValue(_arrRow);
                            }
                        }
                    }
                    break;
                    #endregion
            }
            return prediction;
        }
        // return trend
        public double Trend()
        {
            double beta0, beta1 = 0;
            double xBar = _numOfObsPeriods / 2;
            double yBar;
            double nomin, denom;

            yBar = _arrAggregatedObs.Average();
            nomin = 0;
            denom = 0;

            for (int t = 0; t < _numOfObsPeriods; ++t)
            {
                nomin += (t - xBar) * (_arrAggregatedObs[t] - yBar);
                denom += Math.Pow(t - xBar, 2);
            }

            // trend
            beta1 = nomin / denom;
            // intercept
            beta0 = yBar - beta1 * xBar;
                    
            return beta1;
        }
        // return sum
        public double Sum(int firstObsPeriod, int lastObsPeriod)
        {
            double sum = 0;
            for (int i = firstObsPeriod; i <= lastObsPeriod; i++)
            {
                sum += _arrAggregatedObs[i];
            }
            return sum;
        }
        // return average
        public double Mean()
        {
            return _arrAggregatedObs.Average();
        }
        // return average
        public double Mean(int firstObsPeriod, int lastObsPeriod)
        {
            double sum = 0;
            for (int i = firstObsPeriod; i <= lastObsPeriod; i++)
            {
                sum += _arrAggregatedObs[i];
            }
            return sum/(lastObsPeriod - firstObsPeriod + 1);
        }

        // reset current aggregated observation
        public void ResetCurrentAggregatedObsInThisObsPeriodervation()
        {
            _currentAggregatedObsInThisObsPeriod = 0;
        }
        // reset
        public void Reset()
        {
            if (_arrAggregatedObs != null)
                _arrAggregatedObs = new double[_numOfObsPeriods];
            //_currentObs = 0;
            _currentAggregatedObsInThisObsPeriod = 0;
            _currentNumOfRecordingsInThisObsPeriod = 0;
        }

        // setup prediction
        private void SetupPrediction()
        {   
            // setup a regression for if prediction uses quadratic function
            if (_predictionModel == enumPredictionModel.Quadratic)
            {
                // create X
                _matX = new double[_numOfObsPeriods, 3];
                _arrRow = new double[3];
                for (int t = 0; t < _numOfObsPeriods; ++t)
                {
                    _matX[t, 0] = 1; // intercept
                    _matX[t, 1] = t; // main effect
                    _matX[t, 2] = t * t; // quadratic effect
                }
            }
        }
    }
}
