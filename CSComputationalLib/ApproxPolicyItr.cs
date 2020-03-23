﻿using System;
using System.Collections.Generic;
using System.Collections;
using System.Linq;
using System.Text;
using RandomVariateLib;
using ComputationLib;

namespace ComputationLib
{

    public enum EnumTransformMethod : int
    {
        None = 0,
        NaturalLog_PositiveArgument = 1,
        NaturalLog_NegativeArgument = 2,
        SquaredRoot_PositiveArgument = 3,
        SquaredRoot_NegativeArgument = 4,
    }
    public enum EnumQFunctionApproximationMethod
    {
        Q_Approximation = 0,
        A_Approximation = 1,
        H_Approximation = 2,
    }

    /// <summary>
    /// Dynamic programing state 
    /// </summary>
    public class DPState
    {
        public double[] FeatureValues { get; }
        public int[] NextPeriodActionCombination { get; }
        public bool ValidStateToUpdateQFunctions { get; set; } = true;
        public double CostToGo { get; set; }
        public double DecisionIntervalCost { get; set; }

        // Instantiation
        public DPState(double[] featureValues, int[] nextPeriodActionCombination)
        {
            FeatureValues = (double[])featureValues.Clone();
            NextPeriodActionCombination = (int[])nextPeriodActionCombination.Clone();
        }
    }

    public class ApproximatePolicyIteration
    {
        public List<DPState> DPStates { get; private set; } = new List<DPState>();   
        public int[][] ActionCombinations { get; private set; } // list of action combinations 
        private readonly int _nOfActions; // n of action combinations = (n of actions)^2        

        // approximation models  
        public EnumQFunctionApproximationMethod QFunctionApproxMethod { get; private set; }
        public List<QFunction> QFunctions { get; private set; }
        public List<QFunction> HFunctions_On { get; private set; }
        public List<QFunction> HFunctions_Off { get; private set; }
        public PolynomialQFunction AFunction { get; private set; } // additive approximation of Q-functions
        private EnumTransformMethod transformMethod;

        // algorithm
        private int _itr;
        private double _discountFactor;
        private ExplorationRule _explorationRule;
        private LearningRule _learningRule;
        public double[] Errors { get; private set; }

        // debugging
        private double[][] _matOf_IfEliggible_Actions_FeatureValues_Responses;

        public ApproximatePolicyIteration(LearningRule learningRule, ExplorationRule explorationRule, double discountFactor)
        {
            DPStates = new List<DPState>();
            _discountFactor = discountFactor;
            _learningRule = learningRule;
            _explorationRule = explorationRule;
        }

        // set up approximation model
        public void SetUpQFunctionApproximationModel(
            EnumQFunctionApproximationMethod qFunctionApproximationMethod, EnumTransformMethod transformMethod,
            int nOfFeatures, int polynomialDegree, double l2Penalty = 0)
        {
            QFunctionApproxMethod = qFunctionApproximationMethod;
            this.transformMethod = transformMethod;

            switch (QFunctionApproxMethod)
            {
                case EnumQFunctionApproximationMethod.Q_Approximation:
                    {
                        QFunctions = new List<QFunction>();
                    }
                    break;
                case EnumQFunctionApproximationMethod.A_Approximation:
                    {
                    }
                    break;
                case EnumQFunctionApproximationMethod.H_Approximation:
                    {
                        HFunctions_Off = new List<QFunction>();
                        HFunctions_On = new List<QFunction>();
                    }
                    break;
            }

            switch (QFunctionApproxMethod)
            {
                case EnumQFunctionApproximationMethod.Q_Approximation:
                    {
                        for (int i = 0; i < Math.Pow(2, _nOfActions); i++)
                        {
                            QFunctions.Add(new PolynomialQFunction(
                                name: "Action combination index for dynamically controlled actions: " + i,
                                numOfIndicatorVariables: 0,
                                numOfContinuousVariables: nOfFeatures,
                                polynomialDegree: polynomialDegree,
                                l2Penalty: l2Penalty));
                        }
                    }
                    break;
                case EnumQFunctionApproximationMethod.A_Approximation:
                    {
                        AFunction = new PolynomialQFunction(
                            name: "Additive approximation function",
                            numOfIndicatorVariables: _nOfActions,
                            numOfContinuousVariables: nOfFeatures,
                            polynomialDegree: polynomialDegree,
                            l2Penalty: l2Penalty);
                    }
                    break;
                case EnumQFunctionApproximationMethod.H_Approximation:
                    {
                        // set up H functions for actions that are guided by the policy
                        for (int a = 0; a < _nOfActions; a++)
                        {
                            HFunctions_Off.Add(new PolynomialQFunction(
                                name: a + "- Off",
                                numOfIndicatorVariables: 0,
                                numOfContinuousVariables: nOfFeatures,
                                polynomialDegree: polynomialDegree,
                                l2Penalty: l2Penalty));
                        }
                    }
                    break;
            }
        }

        // add an ADP state
        public void AddAnADPState(DPState DPState)
        {
            DPStates.Add(DPState);
        }
        // add to the cost of current decision interval
        public void AddToCostOfCurrentDecisionInterval(double cost)
        {
            DPStates.Last().DecisionIntervalCost += cost;
        }
        // add to a decision interval cost
        public void AddToDecisionIntervald(int DPStateIndexInCollection, double cost)
        {
            DPStates[DPStateIndexInCollection].DecisionIntervalCost += cost;
        }

        // find the optimal action combination
        public int[] FindOptimalActionCombination(double[] featureValues)
        {
            int[] result = null;
            int optActionCombIndex = 0;

            switch (QFunctionApproxMethod)
            {
                case EnumQFunctionApproximationMethod.Q_Approximation:
                    #region Q-Approximation
                    {
                        int actionCombIndex = 0;
                        double min = double.MaxValue, qValue = 0;
                        // for each available action
                        foreach (QFunction thisQFunction in QFunctions)
                        {
                            // find the value of this action
                            qValue = TransformBack(thisQFunction.fValue(featureValues), transformMethod);
                            if (qValue < min)
                            {
                                min = qValue;
                                optActionCombIndex = actionCombIndex;
                            }
                            ++actionCombIndex;
                        }
                        result = (int[])ActionCombinations[optActionCombIndex].Clone();
                    }
                    break;
                #endregion
                case EnumQFunctionApproximationMethod.A_Approximation:
                    #region Additive-Approximation
                    {
                        double min = double.MaxValue, qValue = 0;
                        // for each possible combination                        
                        for (int i = 0; i < Math.Pow(2, _nOfActions); ++i)
                        {
                            // find the value of this action
                            qValue = TransformBack(
                                AFunction.fValue(ActionCombinations[i], featureValues), transformMethod);
                            if (qValue < min)
                            {
                                min = qValue;
                                optActionCombIndex = i;
                            }
                        }
                        result = (int[])ActionCombinations[optActionCombIndex].Clone();
                    }
                    break;
                #endregion
                case EnumQFunctionApproximationMethod.H_Approximation:
                    #region H-Approximation
                    {
                        // start with assuming that all vertices are in the optimal set (1: if yes, 0: if no)
                        int[] statusOfActionCombinationVertices = new int[(int)Math.Pow(2, _nOfActions)];
                        SupportFunctions.MakeArrayEqualTo(ref statusOfActionCombinationVertices, 1);


                        // eliminate vertices that are not in the optimal set 
                        for (int a = 0; a < _nOfActions; a++)
                        {
                            // if turning off this action lead to better outcome
                            if (((QFunction)HFunctions_Off[a]).fValue(featureValues) >=
                                ((QFunction)HFunctions_On[a]).fValue(featureValues))
                            {
                                // find the vertices that have this action on
                                for (int actionCombIndex = 0; actionCombIndex < Math.Pow(2, _nOfActions); ++actionCombIndex)
                                {
                                    if (ActionCombinations[actionCombIndex][a] == 1)
                                        // then remove this vertix
                                        statusOfActionCombinationVertices[actionCombIndex] = 0;
                                }
                            }
                            else
                            {
                                // find the vertices that have this action off
                                for (int actionCombIndex = 0; actionCombIndex < Math.Pow(2, _nOfActions); ++actionCombIndex)
                                {
                                    if (ActionCombinations[actionCombIndex][a] == 0)
                                        // then remove this vertix
                                        statusOfActionCombinationVertices[actionCombIndex] = 0;
                                }
                            }
                        }
                        // return the optimal vertix
                        for (int actionCombIndex = 0; actionCombIndex < Math.Pow(2, _nOfActions); ++actionCombIndex)
                        {
                            if (statusOfActionCombinationVertices[actionCombIndex] == 1)
                                result = (int[])ActionCombinations[actionCombIndex].Clone();
                        }
                    }
                    break;
                    #endregion
            }
            return result;
        }

        // find an epsilon greedy decision        
        public int[] FindEpsilonGreedyActionCombination (RNG rng, double[] featureValues)
        {
            int[] anActionCombination;
            // with probability of epsilon, make a random decision
            if (rng.NextDouble() <= _explorationRule.GetEpsilon(_itr))
                anActionCombination = GetARandomActionCombination(rng);//, timeIndex, arrAvailableResources);
            else // make a greedy decision
                anActionCombination = FindOptimalActionCombination(featureValues);//, timeIndex, arrAvailableResources);

            return anActionCombination;
        }
        
        // backpropagation 
        public bool DoBackpropagation(int itr, bool stoppedDueToEradication, bool useDecisionsAsFeature)
        {

            // is there any DP state to process?
            if (DPStates.Count == 0)
                return false;

            _backPropogationResults[dim] = true;
            // array of errors
            //_ADPPredictionErrors[dim] = new double[numOfADPStates];
            double[] thisPredictionErrorsForEligibleStates = new double[0];

            // get the last ADP state-decision
            #region last ADP state
            DPState lastADPStateDecision = (DPState)thisCollectionOfADPStates[numOfADPStates - 1];
            // update the reward to go of the last ADP state-decision
            if (stoppedDueToEradication == true)
                lastADPStateDecision.CostToGo = lastADPStateDecision.DecisoinIntervalReward;
            else // not eradicated
            {
                lastADPStateDecision.CostToGo =
                    EstimatedTransformedRewardToGo(lastADPStateDecision.SelectedNextPeriodActionCombination, lastADPStateDecision.ObservationFeatureValues);
            }
            #endregion

            // do back propagation for the rest of the state-decisions
            #region calculate other state's reward to go
            for (int i = thisCollectionOfADPStates.Count - 1; i >= 1; --i)
            {
                // get the last ADP state-decision
                DPState thisADPState = (DPState)thisCollectionOfADPStates[i - 1];
                DPState nextADPState = (DPState)thisCollectionOfADPStates[i];

                // update the reward to go of this ADP state-decision
                thisADPState.CostToGo = thisADPState.DecisoinIntervalReward
                    + discountFactor * nextADPState.CostToGo;
            }
            #endregion

            // debugging information
            #region debugging
            _matOf_IfEliggible_Actions_FeatureValues_Responses = new double[0][];
            foreach (DPState thisADPState in thisCollectionOfADPStates)
            {
                if (thisADPState.ValidStateToUpdateQFunction)
                {
                    double[] thisRowOf_Actions_FeatureValues_Responses = new double[0];
                    // action code
                    SupportFunctions.AddToEndOfArray(ref thisRowOf_Actions_FeatureValues_Responses, SupportFunctions.ConvertToBase10FromBase2(thisADPState.SelectedNextPeriodActionCombination));
                    // feature values
                    double[] featureValues = thisADPState.ObservationFeatureValues;
                    for (int i = 0; i < featureValues.Length; i++)
                        SupportFunctions.AddToEndOfArray(ref thisRowOf_Actions_FeatureValues_Responses, featureValues[i]);
                    // response
                    SupportFunctions.AddToEndOfArray(ref thisRowOf_Actions_FeatureValues_Responses, thisADPState.CostToGo);
                    // concatinate
                    _matOf_IfEliggible_Actions_FeatureValues_Responses =
                        SupportFunctions.ConcatJaggedArray(_matOf_IfEliggible_Actions_FeatureValues_Responses, thisRowOf_Actions_FeatureValues_Responses);
                }
            }
            #endregion

            // find errors
            #region errors
            int ADPStateIndex = 0;
            foreach (DPState thisADPState in thisCollectionOfADPStates)
            {
                if (thisADPState.ValidStateToUpdateQFunction)
                {
                    double transformedObservedRewardToGo = Transform(thisADPState.CostToGo, transformMethod);

                    // find the estimate reward to go
                    double transformedEstimatedRewardToGo =
                        EstimatedTransformedRewardToGo(thisADPState.SelectedNextPeriodActionCombination, thisADPState.ObservationFeatureValues);

                    SupportFunctions.AddToEndOfArray(ref thisPredictionErrorsForEligibleStates, transformedObservedRewardToGo - transformedEstimatedRewardToGo);
                    ++ADPStateIndex;
                }
            }
            Errors[dim] = thisPredictionErrorsForEligibleStates;
            #endregion

            // update Q-functions
            #region update q-functions
            for (int i = thisCollectionOfADPStates.Count - 1; i >= 0; i--)
            {
                // get the ADP State
                DPState thisADPStateDecision = (DPState)thisCollectionOfADPStates[i];

                // first check if this ADP state can be used to update Q-functions
                if (thisADPStateDecision.ValidStateToUpdateQFunction)
                {
                    int[] switchStatusOfActionsDynamicallyControlled = GetSwitchStatusOfActionsControlledDynamically(thisADPStateDecision.SelectedNextPeriodActionCombination);

                    switch (_qFunctionApproximationMethod)
                    {
                        case enumQFunctionApproximationMethod.Q_Approximation:
                            {
                                ((QFunction)_colOfQFunctions[SupportFunctions.ConvertToBase10FromBase2(switchStatusOfActionsDynamicallyControlled)])
                                    .Update(thisADPStateDecision.ObservationFeatureValues, Transform(thisADPStateDecision.CostToGo, transformMethod), itr);
                            }
                            break;
                        case enumQFunctionApproximationMethod.A_Approximation:
                            _qFunctionApproximationModel_Additive.Update(
                                switchStatusOfActionsDynamicallyControlled,
                                thisADPStateDecision.ObservationFeatureValues,
                                Transform(thisADPStateDecision.CostToGo, transformMethod), itr);
                            break;
                        case enumQFunctionApproximationMethod.H_Approximation:
                            {
                                for (int a = 0; a < _numOfActionsControlledDynamically; a++)
                                {
                                    // if this vertex is on
                                    if (switchStatusOfActionsDynamicallyControlled[a] == 1)
                                        ((PolynomialQFunction)_colOfHFunctions_On[a]).Update(
                                            thisADPStateDecision.ObservationFeatureValues,
                                            Transform(thisADPStateDecision.CostToGo, transformMethod), itr);
                                    else // if this vertex is off
                                        ((PolynomialQFunction)_colOfHFunctions_Off[a]).Update(
                                            thisADPStateDecision.ObservationFeatureValues,
                                            Transform(thisADPStateDecision.CostToGo, transformMethod), itr);
                                }
                            }
                            break;
                    }
                }
            }
            #endregion

            return true;
        }
        
        // get selected next period action combination of a dp state
        public int[] GetSelectedNextPeriodActionCombinationOfAnADPState(int ADPStateIndexInCollection)
        {
            return (int[])((DPState)DPStates[ADPStateIndexInCollection]).NextPeriodActionCombination.Clone();
        }

        // get cost-to-go
        public double GetCostToGo(int DPStateIndex)
        {
            return DPStates[DPStateIndex].CostToGo;
        }
        // prediction errors
        public double ADPPredictionErrors(int DPStateIndex)
        {
            if (Errors.Length == 0)
                return 0;
            else
                return Errors[DPStateIndex];
        }

        // reset for another simulation run
        public void ResetForAnotherSimulationRun()//(ref double initialCost)
        {
            if (DPStates != null)
                DPStates.Clear();
        }

        // get approximating q-function polynomial terms
        public int[,] GetQFunctionPolynomialTerms()
        {
            int[,] result = new int[0, 0];
            switch (QFunctionApproxMethod)
            {
                case EnumQFunctionApproximationMethod.Q_Approximation:
                    //TODO: update GetQFunctionPolynomialTerms function for Q-Approximation
                    result = new int[1, 1];// _qFunctionApproximationModel_Additive.RegressionTermDegrees;
                    break;
                case EnumQFunctionApproximationMethod.A_Approximation:
                    result = (int[,])AFunction.RegressionTermDegrees.Clone();
                    break;
                case EnumQFunctionApproximationMethod.H_Approximation:
                    //TODO: update GetQFunctionPolynomialTerms function for H-Approximation
                    result = new int[1, 1];// _qFunctionApproximationModel_Additive.RegressionTermDegrees;
                    break;
            }
            return result;
        }
        // get the Q-function estimates
        public double[] GetQFunctionCoefficientEstimates()
        {
            double[] result = new double[0];
            switch (QFunctionApproxMethod)
            {
                case EnumQFunctionApproximationMethod.Q_Approximation:
                    //TODO: update GetQFunctionCoefficientEstimates function for Q-Approximation
                    result = new double[1];// _qFunctionApproximationModel_Additive.CoeffientEstimates;
                    break;
                case EnumQFunctionApproximationMethod.A_Approximation:
                    result = (double[])AFunction.Coefficients.Clone();
                    break;
                case EnumQFunctionApproximationMethod.H_Approximation:
                    //TODO: update GetQFunctionCoefficientEstimates function for H-Approximation
                    result = new double[1];// _qFunctionApproximationModel_Additive.CoeffientEstimates;
                    break;
            }
            return result;
        }
        // update the estimate of the Q-function
        public void UpdateQFunctionCoefficients(double[] coefficients)
        {
            switch (QFunctionApproxMethod)
            {
                case EnumQFunctionApproximationMethod.Q_Approximation:

                    //TODO: update UpdateQFunctionCoefficients function for Q-Approximation
                    break;
                case EnumQFunctionApproximationMethod.A_Approximation:
                    AFunction.UpdateCoefficients(coefficients);
                    break;
                case EnumQFunctionApproximationMethod.H_Approximation:
                    //TODO: update UpdateQFunctionCoefficients function for H-Approximation
                    //_qFunctionApproximationModel_Additive.UpdateCoefficients(coefficients);
                    break;
            }
        }

        // PRIVATE SUBS
        #region Private Subs
        // get a random action combination 
        private int[] GetARandomActionCombination(RNG rng)
        {
            int index = rng.Next(_nOfActions);
            return ActionCombinations[index];
        }

        // transformation
        private double Transform(double value, EnumTransformMethod transformMethod)
        {
            double result = 0;
            switch (transformMethod)
            {
                case EnumTransformMethod.None:
                    result = value;
                    break;
                case EnumTransformMethod.NaturalLog_PositiveArgument:
                    result = Math.Log(Math.Max(0, value) + 1);
                    break;
                case EnumTransformMethod.NaturalLog_NegativeArgument:
                    result = Math.Log(Math.Max(0, -value) + 1);
                    break;
                case EnumTransformMethod.SquaredRoot_PositiveArgument:
                    result = Math.Sqrt(Math.Max(0, value));
                    break;
                case EnumTransformMethod.SquaredRoot_NegativeArgument:
                    result = Math.Sqrt(Math.Max(0, -value));
                    break;
            }

            return result;
        }
        private double TransformBack(double transformedValue, EnumTransformMethod transformMethod)
        {
            double result = 0;
            switch (transformMethod)
            {
                case EnumTransformMethod.None:
                    result = transformedValue;
                    break;
                case EnumTransformMethod.NaturalLog_PositiveArgument:
                    result = Math.Max(0, Math.Exp(transformedValue) - 1);
                    break;
                case EnumTransformMethod.NaturalLog_NegativeArgument:
                    result = Math.Min(0, 1 - Math.Exp(transformedValue));
                    break;
                case EnumTransformMethod.SquaredRoot_PositiveArgument:
                    result = Math.Pow(transformedValue, 2);
                    break;
                case EnumTransformMethod.SquaredRoot_NegativeArgument:
                    result = -Math.Pow(transformedValue, 2);
                    break;
            }

            return result;
        }

        // return the estimated reward to go
        private double EstimatedTransformedRewardToGo(int[] nextPeriodActionCombination, double[] featureValues)
        {
            double estimatedTransformedRewardToGo = 0;

            // find the estimate reward to go
            switch (QFunctionApproxMethod)
            {
                case EnumQFunctionApproximationMethod.Q_Approximation:
                    {
                        int actionCombIndex = SupportFunctions.ConvertToBase10FromBase2(nextPeriodActionCombination);
                        estimatedTransformedRewardToGo = QFunctions[actionCombIndex].fValue(featureValues);
                    }
                    break;
                case EnumQFunctionApproximationMethod.A_Approximation:
                    {
                        estimatedTransformedRewardToGo = AFunction.fValue(nextPeriodActionCombination, featureValues);
                    }
                    break;
                case EnumQFunctionApproximationMethod.H_Approximation:
                    {
                        for (int i = 0; i < _nOfActions; i++)
                        {
                            // if this vertex is on
                            if (nextPeriodActionCombination[i] == 1)
                                estimatedTransformedRewardToGo += HFunctions_On[i].fValue(featureValues);
                            else // if this vertex is off
                                estimatedTransformedRewardToGo += HFunctions_Off[i].fValue(featureValues);
                        }
                        estimatedTransformedRewardToGo = estimatedTransformedRewardToGo / _nOfActions;
                    }
                    break;
            }

            return estimatedTransformedRewardToGo;

        }

        #endregion
    }

    
}