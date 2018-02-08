using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Collections;

namespace ComputationLib
{
    public class cSingleLayerNeuralNet
    {
        int _numOfInputVariables;
        int _numOfHiddenNodes;
        double _regulationPenalty;
        Random thisRandom = new Random();

        ArrayList _hiddenNodes = new ArrayList();
        double[] _arrOutputWeights;

        double _currentOutput;

        // instantiation
        public cSingleLayerNeuralNet(int numOfInputVariables, int numOfHiddenNodes)
        {
            _numOfInputVariables = numOfInputVariables;
            _numOfHiddenNodes = numOfHiddenNodes;

            // build hidden nodes
            for (int i = 0; i < numOfHiddenNodes; i++)
            {
                // create a hidden node
                HiddenNode thisHiddenNode = new HiddenNode(numOfInputVariables, 1);
                // initialize the weights
                double[] thisWeights = new double[numOfInputVariables + 1];
                for (int varIndex = 0; varIndex <= numOfInputVariables; ++varIndex)
                    thisWeights[varIndex] = thisRandom.NextDouble() - 0.5;  // (double)(varIndex/ numOfHiddenNodes) - 0.5;
                thisHiddenNode.InitializeWeights(thisWeights);
                // add the hidden node
                _hiddenNodes.Add(thisHiddenNode);
            }

            _arrOutputWeights = new double[numOfHiddenNodes + 1]; // the first position is for the intercept            
            for (int i = 0; i < _arrOutputWeights.Length; i++)
                _arrOutputWeights[i] = thisRandom.NextDouble() - 0.5;
        }

        // properties
        public double CurrentOutput
        { get { return _currentOutput; } }

        // customize the training settings
        public void CustomizeTrainingSetting(double activationFunctionScaleParameter, double regulationPenalty)
        {
            foreach (HiddenNode thisHiddenNode in _hiddenNodes)
            {
                thisHiddenNode.ActivationFunctionScaleParameter = activationFunctionScaleParameter;
                thisHiddenNode.RegulationPenalty = regulationPenalty;
            }
            _regulationPenalty = regulationPenalty;
        }
        
        // return output
        public double CalculateOutput(double[] inputVariables)
        {
            RunFeedForward(inputVariables);
            return _currentOutput;
        }

        // update the neural network parameters
        public void Update(double[] x, double y, double learningRate)
        {
            double[,] xMatrix = new double[1, _numOfInputVariables];
            double[] yVector = new double[1];

            // find this observation
            for (int i = 0; i < _numOfInputVariables; ++i)
                xMatrix[0, i] = x[i];
            yVector[0] = y;

            // train
            Train(xMatrix, yVector, learningRate);
        }

        // train with a fixed data set
        public void Train(double[,] xMatrix, double[] y, double learningRate)
        {
            int numOfObservations = y.Length;
            double[] inputVariables = new double[_numOfInputVariables];
            double currentError; 
            //
            double[] arrDeltaOutputWeights = new double[_numOfHiddenNodes + 1]; // first spot for the intercept

            for (int obsIndex = 0; obsIndex < numOfObservations; obsIndex++)
            {
                // find this input variable
                for (int j = 0; j < _numOfInputVariables; j++)
                    inputVariables[j] = xMatrix[obsIndex, j];
                
                // run feedforward
                RunFeedForward(inputVariables);

                // calculate the current error
                currentError = _currentOutput - y[obsIndex];

                // update the intercept            
                arrDeltaOutputWeights[0] += -learningRate * (currentError + _regulationPenalty * _arrOutputWeights[0]);

                // calculate the main output delta weights
                foreach (HiddenNode thisHiddenNode in _hiddenNodes)
                    arrDeltaOutputWeights[thisHiddenNode.ID + 1] += -learningRate *
                        (thisHiddenNode.CurrentOutputValue * currentError + _regulationPenalty * _arrOutputWeights[thisHiddenNode.ID + 1]);

                // calculate delta weights for each hidden nodes
                foreach (HiddenNode thisHiddenNode in _hiddenNodes)
                    thisHiddenNode.UpdateDeltaWeightsUsingThisObservation(inputVariables, currentError, _arrOutputWeights[thisHiddenNode.ID + 1], learningRate);
            }

            // update output weights - first index is for the interface
            for (int hiddenNodeIndex = 0; hiddenNodeIndex <= _numOfHiddenNodes; ++hiddenNodeIndex)
                _arrOutputWeights[hiddenNodeIndex] += arrDeltaOutputWeights[hiddenNodeIndex];
            // update output weights in hidden node
            foreach (HiddenNode thisHiddenNode in _hiddenNodes)
                thisHiddenNode.UpdateWeigthsUsingDeltaWeights();
        }

        // calculate the sum of squared errors
        public double CalculateSumOfSequaredErrors(double[,] xMatrix, double[] y)
        {
            int numOfObservations = y.Length;
            double[] inputVariables = new double[_numOfInputVariables];
            double sumOfSqrErrs = 0;

            for (int obsIndex = 0; obsIndex < numOfObservations; obsIndex++)
            {
                // find this input variable
                for (int j = 0; j < _numOfInputVariables; j++)
                    inputVariables[j] = xMatrix[obsIndex, j];

                // run feedforward
                RunFeedForward(inputVariables);

                // update the sum of squared errors
                sumOfSqrErrs += Math.Pow(_currentOutput - y[obsIndex], 2);
            }

            return sumOfSqrErrs;
        }

        // feedforward
        private void RunFeedForward(double[] inputVariables)
        {
            // update input and output of each node and the outpu of the NN
            double value = _arrOutputWeights[0]; // intercept

            foreach (HiddenNode thisHiddenNode in _hiddenNodes)
            {
                // update the input-output of each hidden node
                thisHiddenNode.UpdateInputOutput(inputVariables);
                // update the NN output
                value += thisHiddenNode.CurrentOutputValue * _arrOutputWeights[thisHiddenNode.ID + 1];
            }
            // update current output
            _currentOutput = value;
        }

        class HiddenNode
        {
            int _ID;
            static int _nextID = 0;
            int _numOfInputVariables;
            double _currentInputValue;
            double _currentOutputValue;
            double[] _arrWeights;
            double[] _arrDeltaWeights;
            double _regulationPenalty;
            double _activationFunctionScaleParameter;

            public double ActivationFunctionScaleParameter
            {
                get { return _activationFunctionScaleParameter; }
                set { _activationFunctionScaleParameter = value; }
            }

            // instantiate
            public HiddenNode(int numOfInputVariables, double activationFunctionScaleParameter)
            {
                _ID = _nextID ++;                
                _numOfInputVariables = numOfInputVariables;
                _arrWeights = new double[numOfInputVariables + 1];
                _arrDeltaWeights = new double[numOfInputVariables + 1];
                _activationFunctionScaleParameter = activationFunctionScaleParameter;
                _regulationPenalty = 0;
            }

            // Properties
            public int ID
            { get { return _ID; } }
            public double RegulationPenalty
            {
                get { return _regulationPenalty; }
                set { _regulationPenalty = value; }
            }
            public double CurrentInputValue
            { get { return _currentInputValue; } }
            public double CurrentOutputValue
            { get { return _currentOutputValue; } }

            // initialize weigths
            public void InitializeWeights(double[] arrWeights)
            {
                _arrWeights = arrWeights;
            }

            // update the input and output of this hidden node
            public double UpdateInputOutput(double[] inputVariables)
            {
                double value = _arrWeights[0]; // intercept
                for (int varIndex = 0; varIndex < _numOfInputVariables; varIndex++)
                    value += inputVariables[varIndex] * _arrWeights[varIndex + 1];

                _currentInputValue = value;
                _currentOutputValue = ActivationFunctionValue(value);
                return _currentOutputValue;
            }

            // update weights
            public void UpdateWeights(double[] x, double outputWeightOfThisNode, double outputError, double learningRate)
            {
                double hiddenNodeError = 0;

                // find the hidden node error
                double sum = _arrWeights[0]; // intercept
                for (int varIndex = 0; varIndex < _numOfInputVariables; varIndex++)
                    sum += _arrWeights[varIndex + 1] * x[varIndex];
                hiddenNodeError = ActivationFunctionDerivative(sum) * outputWeightOfThisNode * outputError;

                // update the weights
                // the intercept
                _arrWeights[0] += -learningRate * (hiddenNodeError + _regulationPenalty * _arrWeights[0]);
                // the main effects
                for (int varIndex = 0; varIndex < _numOfInputVariables; varIndex++)
                    _arrWeights[varIndex + 1] += -learningRate * 
                        (x[varIndex] * hiddenNodeError + _regulationPenalty * _arrWeights[varIndex + 1]);
            }
            
            // update delta weights using this observation
            public void UpdateDeltaWeightsUsingThisObservation(double[] inputVariables, double currentNNOutputError, double outputWeightOfThisNode, double learningRate)
            {
                // the intercept
                _arrDeltaWeights[0] += -learningRate * 
                    (currentNNOutputError * outputWeightOfThisNode * ActivationFunctionDerivative(_currentInputValue) + _regulationPenalty * _arrWeights[0]);
                // the main effects
                for (int varIndex = 0; varIndex < _numOfInputVariables; varIndex++)
                    _arrDeltaWeights[varIndex + 1] += -learningRate *
                     (currentNNOutputError * outputWeightOfThisNode * ActivationFunctionDerivative(_currentInputValue) * inputVariables[varIndex] + _regulationPenalty * _arrWeights[varIndex + 1]);
            }

            // update weights using delta weights
            public void UpdateWeigthsUsingDeltaWeights()
            {
                // the main effects
                for (int varIndex = 0; varIndex <= _numOfInputVariables; varIndex++)
                    _arrWeights[varIndex] += _arrDeltaWeights[varIndex];

                _arrDeltaWeights = new double[_numOfInputVariables + 1];
            }

            #region Private Functions
            // activation function value
            private double ActivationFunctionValue(double v)
            {
                return 1 / (1 + Math.Exp(-_activationFunctionScaleParameter * v));
            }

            // activation function derivative
            private double ActivationFunctionDerivative(double v)
            {
                double exp = Math.Exp(-_activationFunctionScaleParameter * v);
                double denominator = Math.Pow(1 + exp, 2);
                return _activationFunctionScaleParameter * exp / denominator;
            }
            #endregion            
        } // end of class HiddenNode
    }
}
