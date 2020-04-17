using System;
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
        Ln_PosArgument = 1,
        Ln_NegArgument = 2,
        Sqrt_PosArgument = 3,
        Sqrt_NegArgument = 4,
    }
    public enum EnumQFuncApproximationMethod
    {
        Q_Approx = 0,
        A_Approx = 1,
        H_Approx = 2,
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
        public EnumQFuncApproximationMethod QFunctionApproxMethod { get; private set; }
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
        public double[] Errors { get => _errors; private set => _errors = value; }

        // debugging
        private double[][] _matOf_IfEliggible_Actions_FeatureValues_Responses;
        private double[] _errors;

        public ApproximatePolicyIteration(LearningRule learningRule, ExplorationRule explorationRule, double discountFactor)
        {
            DPStates = new List<DPState>();
            _discountFactor = discountFactor;
            _learningRule = learningRule;
            _explorationRule = explorationRule;
        }

        // set up approximation model
        public void SetUpQFunctionApproximationModel(
            EnumQFuncApproximationMethod qFunctionApproximationMethod, EnumTransformMethod transformMethod,
            int nOfFeatures, int polynomialDegree, double l2Penalty = 0)
        {
            QFunctionApproxMethod = qFunctionApproximationMethod;
            this.transformMethod = transformMethod;

            switch (QFunctionApproxMethod)
            {
                case EnumQFuncApproximationMethod.Q_Approx:
                    {
                        QFunctions = new List<QFunction>();
                    }
                    break;
                case EnumQFuncApproximationMethod.A_Approx:
                    {
                    }
                    break;
                case EnumQFuncApproximationMethod.H_Approx:
                    {
                        HFunctions_Off = new List<QFunction>();
                        HFunctions_On = new List<QFunction>();
                    }
                    break;
            }

            switch (QFunctionApproxMethod)
            {
                case EnumQFuncApproximationMethod.Q_Approx:
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
                case EnumQFuncApproximationMethod.A_Approx:
                    {
                        AFunction = new PolynomialQFunction(
                            name: "Additive approximation function",
                            numOfIndicatorVariables: _nOfActions,
                            numOfContinuousVariables: nOfFeatures,
                            polynomialDegree: polynomialDegree,
                            l2Penalty: l2Penalty);
                    }
                    break;
                case EnumQFuncApproximationMethod.H_Approx:
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
                case EnumQFuncApproximationMethod.Q_Approx:
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
                case EnumQFuncApproximationMethod.A_Approx:
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
                case EnumQFuncApproximationMethod.H_Approx:
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
        public int[] FindEpsilonGreedyActionCombination(RNG rng, double[] featureValues)
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

            // array of errors
            Errors = new double[DPStates.Count];

            // get the last ADP state-decision
            #region last ADP state
            DPState lastADPStateDecision = DPStates.Last();
            // update the reward to go of the last ADP state-decision
            if (stoppedDueToEradication == true)
                lastADPStateDecision.CostToGo = lastADPStateDecision.DecisionIntervalCost;
            else // not eradicated
            {
                lastADPStateDecision.CostToGo =
                    EstimatedTransformedRewardToGo(
                        lastADPStateDecision.NextPeriodActionCombination,
                        lastADPStateDecision.FeatureValues);
            }
            #endregion

            // do back propagation for the rest of the state-decisions
            #region calculate other state's reward to go
            for (int i = DPStates.Count - 1; i >= 1; --i)
            {
                // update the reward to go of this ADP state-decision
                DPStates[i - 1].CostToGo = DPStates[i - 1].DecisionIntervalCost + _discountFactor * DPStates[i].CostToGo;
            }
            #endregion

            // debugging information
            #region debugging
            _matOf_IfEliggible_Actions_FeatureValues_Responses = new double[0][];
            foreach (DPState thisDPState in DPStates)
            {
                double[] thisRowOf_Actions_FeatureValues_Responses = new double[0];
                // action code
                SupportFunctions.AddToEndOfArray(
                    ref thisRowOf_Actions_FeatureValues_Responses, 
                    SupportFunctions.ConvertToBase10FromBase2(thisDPState.NextPeriodActionCombination));
                // feature values
                double[] featureValues = thisDPState.FeatureValues;
                for (int i = 0; i < featureValues.Length; i++)
                    SupportFunctions.AddToEndOfArray(ref thisRowOf_Actions_FeatureValues_Responses, featureValues[i]);
                // response
                SupportFunctions.AddToEndOfArray(ref thisRowOf_Actions_FeatureValues_Responses, thisDPState.CostToGo);
                // concatinate
                _matOf_IfEliggible_Actions_FeatureValues_Responses =
                    SupportFunctions.ConcatJaggedArray(_matOf_IfEliggible_Actions_FeatureValues_Responses, thisRowOf_Actions_FeatureValues_Responses);
            }
            #endregion

            // find errors
            #region errors
            int ADPStateIndex = 0;
            foreach (DPState thisDPState in DPStates)
            {
                double transformedObservedRewardToGo = Transform(thisDPState.CostToGo, transformMethod);

                // find the estimate reward to go
                double transformedEstimatedRewardToGo =
                    EstimatedTransformedRewardToGo(
                        thisDPState.NextPeriodActionCombination,
                        thisDPState.FeatureValues);

                SupportFunctions.AddToEndOfArray(ref _errors, transformedObservedRewardToGo - transformedEstimatedRewardToGo);
                ++ADPStateIndex;
            }
            #endregion

            // update Q-functions
            #region update q-functions
            for (int i = DPStates.Count - 1; i >= 0; i--)
            {
                // get the ADP State
                DPState thisADPState = DPStates[i];

                // first check if this ADP state can be used to update Q-functions
                switch (QFunctionApproxMethod)
                {
                    case EnumQFuncApproximationMethod.Q_Approx:
                        {
                            QFunctions[SupportFunctions.ConvertToBase10FromBase2(thisADPState.NextPeriodActionCombination)]
                                .Update(thisADPState.FeatureValues, Transform(thisADPState.CostToGo, transformMethod), itr);
                        }
                        break;
                    case EnumQFuncApproximationMethod.A_Approx:
                        AFunction.Update(
                            thisADPState.NextPeriodActionCombination,
                            thisADPState.FeatureValues,
                            Transform(thisADPState.CostToGo, transformMethod), itr);
                        break;
                    case EnumQFuncApproximationMethod.H_Approx:
                        {
                            for (int a = 0; a < _nOfActions; a++)
                            {
                                // if this vertex is on
                                if (thisADPState.NextPeriodActionCombination[a] == 1)
                                    ((PolynomialQFunction)HFunctions_On[a]).Update(
                                        thisADPState.FeatureValues,
                                        Transform(thisADPState.CostToGo, transformMethod), itr);
                                else // if this vertex is off
                                    ((PolynomialQFunction)HFunctions_Off[a]).Update(
                                        thisADPState.FeatureValues,
                                        Transform(thisADPState.CostToGo, transformMethod), itr);
                            }
                        }
                        break;
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
                case EnumQFuncApproximationMethod.Q_Approx:
                    //TODO: update GetQFunctionPolynomialTerms function for Q-Approximation
                    result = new int[1, 1];// _qFunctionApproximationModel_Additive.RegressionTermDegrees;
                    break;
                case EnumQFuncApproximationMethod.A_Approx:
                    result = (int[,])AFunction.RegressionTermDegrees.Clone();
                    break;
                case EnumQFuncApproximationMethod.H_Approx:
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
                case EnumQFuncApproximationMethod.Q_Approx:
                    //TODO: update GetQFunctionCoefficientEstimates function for Q-Approximation
                    result = new double[1];// _qFunctionApproximationModel_Additive.CoeffientEstimates;
                    break;
                case EnumQFuncApproximationMethod.A_Approx:
                    result = (double[])AFunction.Coefficients.Clone();
                    break;
                case EnumQFuncApproximationMethod.H_Approx:
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
                case EnumQFuncApproximationMethod.Q_Approx:

                    //TODO: update UpdateQFunctionCoefficients function for Q-Approximation
                    break;
                case EnumQFuncApproximationMethod.A_Approx:
                    AFunction.UpdateCoefficients(coefficients);
                    break;
                case EnumQFuncApproximationMethod.H_Approx:
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
                case EnumTransformMethod.Ln_PosArgument:
                    result = Math.Log(Math.Max(0, value) + 1);
                    break;
                case EnumTransformMethod.Ln_NegArgument:
                    result = Math.Log(Math.Max(0, -value) + 1);
                    break;
                case EnumTransformMethod.Sqrt_PosArgument:
                    result = Math.Sqrt(Math.Max(0, value));
                    break;
                case EnumTransformMethod.Sqrt_NegArgument:
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
                case EnumTransformMethod.Ln_PosArgument:
                    result = Math.Max(0, Math.Exp(transformedValue) - 1);
                    break;
                case EnumTransformMethod.Ln_NegArgument:
                    result = Math.Min(0, 1 - Math.Exp(transformedValue));
                    break;
                case EnumTransformMethod.Sqrt_PosArgument:
                    result = Math.Pow(transformedValue, 2);
                    break;
                case EnumTransformMethod.Sqrt_NegArgument:
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
                case EnumQFuncApproximationMethod.Q_Approx:
                    {
                        int actionCombIndex = SupportFunctions.ConvertToBase10FromBase2(nextPeriodActionCombination);
                        estimatedTransformedRewardToGo = QFunctions[actionCombIndex].fValue(featureValues);
                    }
                    break;
                case EnumQFuncApproximationMethod.A_Approx:
                    {
                        estimatedTransformedRewardToGo = AFunction.fValue(nextPeriodActionCombination, featureValues);
                    }
                    break;
                case EnumQFuncApproximationMethod.H_Approx:
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
