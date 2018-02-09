using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ComputationLib;

namespace ComputationLib
{
    public class ADP_State
    {
        // Fields
        private double[] _observationFeatureValues;
        private int[] _selectedNextPeriodActionCombination;
        private int[] _previousPeriodActionCombination;        
        
        private bool _validStateToUpdateQFunctions = true;
        private double _rewardToGo;
        private double _decisionIntervalReward;

        // Instantiation
        public ADP_State(double[] observationFeatureValues, int[] selectedNextPeriodActionCombination)
        {
            _observationFeatureValues = (double[])observationFeatureValues.Clone();
            _selectedNextPeriodActionCombination = (int[])selectedNextPeriodActionCombination.Clone();
        }
        public ADP_State(double[] observationFeatureValues, int[] selectedNextPeriodActionCombination, int[] previousPeriodActionCombination)
        {
            _observationFeatureValues = (double[])observationFeatureValues.Clone();
            _selectedNextPeriodActionCombination = (int[])selectedNextPeriodActionCombination.Clone();
            _previousPeriodActionCombination = previousPeriodActionCombination;                        
        }
        
        // Properties
        public double[] ObservationFeatureValues
        {
            get { return _observationFeatureValues; }
        }
        public int[] SelectedNextPeriodActionCombination
        {
            get { return _selectedNextPeriodActionCombination; }
        }
        public int[] PreviousPeriodActionCombination
        {
            get { return _previousPeriodActionCombination; }
        }
        public bool ValidStateToUpdateQFunction
        {
            get { return _validStateToUpdateQFunctions; }
            set { _validStateToUpdateQFunctions = value; }
        }
        public double DecisoinIntervalReward
        {
            get{ return _decisionIntervalReward; }
        }
        public double RewardToGo
        {
            get { return _rewardToGo; }
            set { _rewardToGo = value; }
        }
        
        // add reward
        public void AddToDecisionIntervalReward(double reward)
        {
            _decisionIntervalReward += reward;
        }
    }
}
