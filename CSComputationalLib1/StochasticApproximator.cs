using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ADPLibrary
{
    public class StochasticApproximator
    {
        private enumMethod _method;
        private int _iteration =0;
        private int _numOfParametorsToEstimate;
        private double _constantStepSize;
        private double _estimate;
        private double[] _estimates;

        public enum enumMethod : int
        {
            ConstantStepsize =0,
            OneOverN=1,

        }        
        // Properties
        public double Estimate
        {
            get
            {
                return _estimate;
            }
        }
        public double[] Estimates
        {
            get
            {
                return _estimates;
            }
        }

        // Methods
        public StochasticApproximator(enumMethod method, int numOfParametorsToEstimate)
        {
            _method = method;
            _numOfParametorsToEstimate = numOfParametorsToEstimate;
            // check if a vector of parameters should be estimated
            if (_numOfParametorsToEstimate > 1)
                _estimates = new double[_numOfParametorsToEstimate];
            else
                _estimate = 0;

        }
        public void SetupConstantStepSize(double stepSize)
        {
            _constantStepSize = stepSize;
        }

        public void Update(double data)
        {
            _iteration += 1; 

            // which stochastic approximation method to use
            switch (_method)
            {
                case enumMethod.ConstantStepsize:
                    _estimate = (1 - _constantStepSize) * _estimate + _constantStepSize * data;
                    break;
                case enumMethod.OneOverN:
                    {
                        if (_iteration == 1)
                        {
                            _estimate = data;
                        }
                        else
                        {
                            double stepSize = (double)1 / _iteration;
                            _estimate = (1 - stepSize) * _estimate + stepSize * data;
                        }
                    }
                    break;
            }
        }
        public void Update(double[] data)
        {
            _iteration += 1; 
            // which stochastic approximation method to use
            switch (_method)
            {
                case enumMethod.ConstantStepsize:
                    {
                        for (int i = 0; i < _numOfParametorsToEstimate; ++i)
                        {
                            _estimates[i] = (1 - _constantStepSize) * _estimates[i] + _constantStepSize * data[i];
                        }
                    }
                    break;
                case enumMethod.OneOverN:
                    {
                        if (_iteration == 1)
                        {
                            for (int i = 0; i < _numOfParametorsToEstimate; ++i)
                            {
                                _estimates[i] = data[i];
                            }
                        }
                        else
                        {
                            double stepSize = (double)1/_iteration;
                            for (int i = 0; i < _numOfParametorsToEstimate; ++i)
                            {
                                _estimates[i] = (1 - stepSize) * _estimates[i] + stepSize * data[i];
                            }
                        }
                    }
                    break;

            }

        }

    }
}
