using System;
using System.Collections.Generic;
using System.Collections;
using System.Linq;
using System.Text;
using CSRandomVariateLib;
using ComputationLib;
using SimulationLib;

namespace OptimizationLib
{
    public class POMDP_ADP
    {
        public enum enumResponseTransformation : int
        {
            None = 0,
            NaturalLog_PositiveArgument = 1,
            NaturalLog_NegativeArgument = 2,
            SquaredRoot_PositiveArgument = 3,
            SquaredRoot_NegativeArgument = 4,
        }
        public enum enumQFunctionApproximationMethod
        {
            Q_Approximation = 0,
            A_Approximation = 1,
            H_Approximation = 2,
        }

        // Fields
        private ArrayList _actions = new ArrayList();

        private int _numOfActions = 0;
        private int[] _currentActionCombination;
        private int[] _defaultActionCombination;

        private int _numOfActionsControlledDynamically;
        private int[] _indicesOfActionsControlledDynamically;
        private int[][] _dynamicallyControlledActionCombinations;
        private int[] _availabilityOfDynamicallyControlledActionCombinations;
        private int[] _currentDynamicallyControlledActionCombinations;

        private double _currentCostPerUnitOfTime;
        private double _currentFixedAnPenaltyCost;

        // adp algorithm
        private int _numOfADPIterations;
        private int _numOfSimRunsToBackPropogate = 1;
        private ArrayList _ADPStates;
        private ArrayList _ADPStateCollections; // for when multiple simulation runs are obtained before back-propogation is done
        private double[][] _ADPPredictionErrors;
        private bool[] _backPropogationResults;
        // approximation models        
        private enumQFunctionApproximationMethod _qFunctionApproximationMethod;
        private ArrayList _colOfQFunctions;
        private ArrayList _colOfHFunctions_On;
        private ArrayList _colOfHFunctions_Off;
        private PolynomialQFunction _qFunctionApproximationModel_Additive;
        private enumResponseTransformation _responseTransformation;

        //private int[] _indecesOfInitiallyAvailableActionCombinations;
        //private int[] _inverseOf_indecesOfInitiallyAvailableActionCombinations;

        // exploration and exploitation                
        private double _epsilonGreedy, _epsilonGreedy_beta, _epsilonGreedy_delta;
        private double _BoltzmannTemperature, _Boltzmann_initialTemperature, _Boltzmann_beta; // T*exp(-beta*i)

        // debugging
        private double[][] _matOf_IfEliggible_Actions_FeatureValues_Responses;

        // Instantiation
        public POMDP_ADP()
        {
        }

        // Properties
        #region properties
        public ArrayList Actions
        {
            get { return _actions; }
        }
        public int NumberOfActions
        {
            get { return _numOfActions; }
        }
        public int NumOfActionsControlledDynamically
        {
            get { return _numOfActionsControlledDynamically; }
        }
        public int[] IndicesOfActionsControlledDynamically
        {
            get { return _indicesOfActionsControlledDynamically; }
        }
        public int[][] DynamicallyControlledActionCombinations
        {
            get { return _dynamicallyControlledActionCombinations; }
        }
        public int[] DefaultActionCombination
        {
            get { return _defaultActionCombination; }
        }
        public int[] CurrentActionCombination
        {
            get { return _currentActionCombination; }
        }
        public enumQFunctionApproximationMethod QFunctionApproximationMethod
        {
            get { return _qFunctionApproximationMethod; }
        }
        public ArrayList ADPStates
        {
            get { return _ADPStates; }
        }
        public int NumberOfADPIterations
        {
            get { return _numOfADPIterations; }
            set { _numOfADPIterations = value; }
        }
        public int NumOfSimRunsToBackPropogate
        {
            get { return _numOfSimRunsToBackPropogate; }
        }
        public double CurrentCostPerUnitOfTime
        {
            get { return _currentCostPerUnitOfTime; }
        }
        public double CurrentFixedAnPenaltyCost
        {
            get { return _currentFixedAnPenaltyCost; }
        }
        #endregion

        // Methods
        // setup
        public void SetupTraining()
        {
            _ADPStates = new ArrayList();
            _ADPPredictionErrors = new double[1][];
            _backPropogationResults = new bool[1];
        }
        public void SetupTraining(int numOfSimRunsToBackPropogate)
        {
            _numOfSimRunsToBackPropogate = numOfSimRunsToBackPropogate;

            // set up adp state collection
            _ADPStateCollections = new ArrayList();
            for (int dim = 0; dim < _numOfSimRunsToBackPropogate; ++dim)
                _ADPStateCollections.Add(new ArrayList());

            _ADPPredictionErrors = new double[_numOfSimRunsToBackPropogate][];
            _backPropogationResults = new bool[_numOfSimRunsToBackPropogate];
        }
        // add a decision
        public void AddAnAction(SimulationAction action)
        {
            // add the actions
            _actions.Add(action);
            ++_numOfActions;
        }
        // add an ADP state
        public void AddAnADPState(ADP_State ADPState)
        {
            _ADPStates.Add(ADPState);
        }
        public void AddAnADPState(int dimension, ADP_State ADPState)
        {
            ((ArrayList)_ADPStateCollections[dimension]).Add(ADPState);
        }
        // gather information after all actions are entered
        public void GatherInformationAfterAllActionsAreEntered()
        {
            _currentActionCombination = new int[_numOfActions];
            _defaultActionCombination = new int[_numOfActions];
            foreach (SimulationAction thisAction in _actions)
            {
                // check if this action can be off
                if (thisAction.Type == SimulationAction.enumActionType.Default)
                    _defaultActionCombination[thisAction.Index] = 1;
                else
                    _defaultActionCombination[thisAction.Index] = 0;
                // check if this action is controlled dynamically
                if (thisAction.OnOffSwitchSetting == SimulationAction.enumOnOffSwitchSetting.Dynamic)
                    ++_numOfActionsControlledDynamically;
            }

            // find the indeces of interventions that are controlled dynamically
            _indicesOfActionsControlledDynamically = new int[_numOfActionsControlledDynamically];
            int a = 0;
            foreach (SimulationAction thisAction in _actions)
            {
                // check if this action is controlled dynamically
                if (thisAction.OnOffSwitchSetting == SimulationAction.enumOnOffSwitchSetting.Dynamic)
                    _indicesOfActionsControlledDynamically[a++] = thisAction.Index;
            }

            if (_numOfActionsControlledDynamically > 0)
            {
                // all action combinations controlled dynamically are assumed to be available
                _availabilityOfDynamicallyControlledActionCombinations = new int[(int)Math.Pow(2, _numOfActionsControlledDynamically)];
                SupportFunctions.MakeArrayEqualTo(ref _availabilityOfDynamicallyControlledActionCombinations, 1);

                // define the action combinations controlled dynamically
                _dynamicallyControlledActionCombinations = new int[0][];
                for (int i = 0; i < Math.Pow(2, _numOfActions); _numOfActionsControlledDynamically++)
                {
                    int[] thisActionCombination = ComputationLib.SupportFunctions.ConvertToBase2FromBase10(i, _numOfActionsControlledDynamically);
                    _dynamicallyControlledActionCombinations = ComputationLib.SupportFunctions.ConcatJaggedArray(_dynamicallyControlledActionCombinations, thisActionCombination);
                }
            }

        }
        // specify the available action combinations controlled dynamically
        public void SpecifyAvailabilityOfDynamicallyControlledActionCombinations(int[] availabilityOfDynamicallyControlledActionCombinations)
        {
            _availabilityOfDynamicallyControlledActionCombinations = (int[])availabilityOfDynamicallyControlledActionCombinations.Clone();
        }
        //// reset the set of available action combinations to the initial setting
        //public void ResetAvailableActionCombintions()
        //{
        //    _ifThisActionCombinationAvailable = (int[])_ifThisActionCombinationAvailable_default.Clone();
        //}
        public void MakeAllDynamicallyControlledActionsAvailable()
        {
            SupportFunctions.MakeArrayEqualTo(ref _availabilityOfDynamicallyControlledActionCombinations, 1);
        }

        // get the current action combination 
        public int[] GetCurrentActionCombination()
        {
            return _currentActionCombination;
        }
        // get the code of the current dynamically controled action combination 
        public int GetCodeOfCurrentDynamicallyControlledActionCombinations()
        {
            return ComputationLib.SupportFunctions.ConvertToBase10FromBase2(_currentDynamicallyControlledActionCombinations);
        }
        // update the currect action combination       
        public void ChangeCurrentActionCombination(int[] newActionCombination, ref double resultingCost)
        {
            int thisActionIndex = 0;
            _currentCostPerUnitOfTime = 0;
            foreach (SimulationAction thisAction in _actions)
            {
                thisActionIndex = thisAction.Index;
                // calculate the number of switches        
                if (_currentActionCombination[thisActionIndex] != newActionCombination[thisActionIndex])
                    thisAction.NumOfSwitchesOccured += 1;
                // calculate the fixed cost
                if (_currentActionCombination[thisActionIndex] == 0 && newActionCombination[thisActionIndex] == 1)
                {
                    resultingCost += thisAction.FixedCost;
                    thisAction.NumOfSwitchesOccured += 1;
                }
                // calculate the penalty cost for switching from on to off
                if (thisAction.PenaltyForSwitchingFromOnToOff > 0)
                {
                    if (_currentActionCombination[thisActionIndex] == 1 && newActionCombination[thisActionIndex] == 0)
                        resultingCost += thisAction.PenaltyForSwitchingFromOnToOff;
                }
                // update the new action
                _currentActionCombination[thisActionIndex] = newActionCombination[thisActionIndex];
                // update the cost per unit of time for this action combination
                if (_currentActionCombination[thisActionIndex] == 1)
                {
                    _currentCostPerUnitOfTime += thisAction.CostPerUnitOfTime;
                    ++thisAction.NumOfDecisionPeriodsOverWhichThisInterventionWasUsed;
                }
            }
        }
        public void ChangeCurrentActionCombination(int[] newActionCombination)
        {
            int thisActionIndex = 0;
            foreach (SimulationAction thisAction in _actions)
            {
                thisActionIndex = thisAction.Index;
                // calculate the number of switches        
                if (_currentActionCombination[thisActionIndex] != newActionCombination[thisActionIndex])
                    thisAction.NumOfSwitchesOccured += 1;
                // update the new action
                _currentActionCombination[thisActionIndex] = newActionCombination[thisActionIndex];
                // update the number of decision periods over which this intervention was used
                if (_currentActionCombination[thisActionIndex] == 1)
                    ++thisAction.NumOfDecisionPeriodsOverWhichThisInterventionWasUsed;
            }
        }
        // update the fixed and penalty costs
        public void UpdateCurrentFixedAndPenaltyCosts(int[] newActionCombination)
        {
            int thisActionIndex = 0;
            foreach (SimulationAction thisAction in _actions)
            {
                thisActionIndex = thisAction.Index;
                // calculate the fixed cost
                if (_currentActionCombination[thisActionIndex] == 0 && newActionCombination[thisActionIndex] == 1)
                {
                    _currentFixedAnPenaltyCost += thisAction.FixedCost;
                }
                // calculate the penalty cost for switching from on to off
                if (thisAction.PenaltyForSwitchingFromOnToOff > 0)
                {
                    if (_currentActionCombination[thisActionIndex] == 1 && newActionCombination[thisActionIndex] == 0)
                        _currentFixedAnPenaltyCost += thisAction.PenaltyForSwitchingFromOnToOff;
                }
            }
        }
        // update the variable cost
        public void UpdateCurrentCostPerUnitOfTime(int[] newActionCombination)
        {
            _currentCostPerUnitOfTime = 0;
            int thisActionIndex = 0;
            foreach (SimulationAction thisAction in _actions)
            {
                thisActionIndex = thisAction.Index;
                // update the cost per unit of time for this action combination
                if (newActionCombination[thisActionIndex] == 1)
                    _currentCostPerUnitOfTime += thisAction.CostPerUnitOfTime;
            }
        }

        // get approximating q-function polynomial terms
        public int[,] GetQFunctionPolynomialTerms()
        {
            int[,] result = new int[0, 0];
            switch (_qFunctionApproximationMethod)
            {
                case enumQFunctionApproximationMethod.Q_Approximation:
                    //TODO: update GetQFunctionPolynomialTerms function for Q-Approximation
                    result = new int[1, 1];// _qFunctionApproximationModel_Additive.RegressionTermDegrees;
                    break;
                case enumQFunctionApproximationMethod.A_Approximation:
                    result = (int[,])_qFunctionApproximationModel_Additive.RegressionTermDegrees.Clone();
                    break;
                case enumQFunctionApproximationMethod.H_Approximation:
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
            switch (_qFunctionApproximationMethod)
            {
                case enumQFunctionApproximationMethod.Q_Approximation:
                    //TODO: update GetQFunctionCoefficientEstimates function for Q-Approximation
                    result = new double[1];// _qFunctionApproximationModel_Additive.CoeffientEstimates;
                    break;
                case enumQFunctionApproximationMethod.A_Approximation:
                    result = (double[])_qFunctionApproximationModel_Additive.CoeffientEstimates.Clone();
                    break;
                case enumQFunctionApproximationMethod.H_Approximation:
                    //TODO: update GetQFunctionCoefficientEstimates function for H-Approximation
                    result = new double[1];// _qFunctionApproximationModel_Additive.CoeffientEstimates;
                    break;
            }
            return result;
        }
        // update the estimate of the Q-function
        public void UpdateQFunctionCoefficients(double[] coefficients)
        {
            switch (_qFunctionApproximationMethod)
            {
                case enumQFunctionApproximationMethod.Q_Approximation:
                    //TODO: update UpdateQFunctionCoefficients function for Q-Approximation
                    break;
                case enumQFunctionApproximationMethod.A_Approximation:
                    _qFunctionApproximationModel_Additive.UpdateCoefficients(coefficients);
                    break;
                case enumQFunctionApproximationMethod.H_Approximation:
                    //TODO: update UpdateQFunctionCoefficients function for H-Approximation
                    //_qFunctionApproximationModel_Additive.UpdateCoefficients(coefficients);
                    break;
            }
        }

        // find the greedy action combinations dynamically controlled
        public int[] FindTheGreedyActionCombinationsDynamicallyControlled(double[] arrObservationFeatureValues)//, ref double resultingCost) //, long timeIndex, long[] arrAvailableResources
        {
            return FindTheOptimalDynamicallyControlledActionCombination(arrObservationFeatureValues);//, timeIndex, arrAvailableResources);
            //// announce the decision
            //ChangeCurrentActionCombination(newActionCombination);//, ref resultingCost);
        }
        // find an epsilon greedy action combinations dynamically controlled
        public int[] FindAnEpsilongGreedyActionCombinationsDynamicallyControlled(RNG rng, double[] arrObservationFeatureValues)//, ref double resultingCost) //long timeIndex, long[] arrAvailableResources,
        {
            return EpsilonGreedyActionCombination(rng, arrObservationFeatureValues);//, timeIndex, arrAvailableResources);
            //// announce the decision
            //ChangeCurrentActionCombination(newActionCombination);//, ref resultingCost);
        }

        //-------------------------------        
        // find the optimal action combination
        public int[] FindTheOptimalDynamicallyControlledActionCombination(double[] arrObservationFeatureValues)//, long timeIndex, long[] arrAvailableResources = null)
        {
            int[] optimalDynamicallyControlledActionCombination = null;
            int optActionCombIndex = 0;

            switch (_qFunctionApproximationMethod)
            {
                case enumQFunctionApproximationMethod.Q_Approximation:
                    #region Q-Approximation
                    {
                        int actionCombIndex = 0;
                        double max = double.MinValue, qValue = 0;
                        // for each available action
                        foreach (QFunction thisQFunction in _colOfQFunctions)
                        {
                            // check if this combination is feasible
                            if (_availabilityOfDynamicallyControlledActionCombinations[actionCombIndex] == 1)
                            {
                                // find the value of this action
                                qValue = TransformBack(thisQFunction.fValue(arrObservationFeatureValues), _responseTransformation);
                                if (qValue > max)
                                {
                                    max = qValue;
                                    optActionCombIndex = actionCombIndex;
                                }
                            }
                            ++actionCombIndex;
                        }
                        optimalDynamicallyControlledActionCombination = (int[])_dynamicallyControlledActionCombinations[optActionCombIndex].Clone();
                    }
                    break;
                #endregion
                case enumQFunctionApproximationMethod.A_Approximation:
                    #region Additive-Approximation
                    {
                        double max = double.MinValue, qValue = 0;
                        // for each possible combination                        
                        for (int i = 0; i < Math.Pow(2, _numOfActionsControlledDynamically); ++i)
                        {
                            // build an action combiniation
                            int[] thisActionCombination = _dynamicallyControlledActionCombinations[i];
                            // check if this combination is feasible
                            if (_availabilityOfDynamicallyControlledActionCombinations[i] == 1)
                            {
                                // find the value of this action
                                qValue = TransformBack(
                                    _qFunctionApproximationModel_Additive.fValue(thisActionCombination, arrObservationFeatureValues), _responseTransformation);
                                if (qValue > max)
                                {
                                    max = qValue;
                                    optActionCombIndex = i;
                                }
                            }
                        }
                        optimalDynamicallyControlledActionCombination = (int[])_dynamicallyControlledActionCombinations[optActionCombIndex].Clone();
                    }
                    break;
                #endregion
                case enumQFunctionApproximationMethod.H_Approximation:
                    #region H-Approximation
                    {
                        // start with assuming that all vertices are in the optimal set (1: if yes, 0: if no)
                        int[] statusOfActionCombinationVertices = new int[(int)Math.Pow(2, _numOfActionsControlledDynamically)];
                        SupportFunctions.MakeArrayEqualTo(ref statusOfActionCombinationVertices, 1);

                        // remove the vertices that are not feasible
                        for (int i = 0; i < Math.Pow(2, _numOfActionsControlledDynamically); ++i)
                        {
                            // check if this combination is feasible
                            if (_availabilityOfDynamicallyControlledActionCombinations[i] == 0)
                                // eliminate this action combination
                                statusOfActionCombinationVertices[i] = 0;
                        }

                        // eliminate vertices that are not in the optimal set 
                        if (_availabilityOfDynamicallyControlledActionCombinations.Sum() > 0)
                        {
                            for (int a = 0; a < _numOfActionsControlledDynamically; a++)
                            {
                                // if turning off this action lead to better outcome
                                if (((QFunction)_colOfHFunctions_Off[a]).fValue(arrObservationFeatureValues) >=
                                    ((QFunction)_colOfHFunctions_On[a]).fValue(arrObservationFeatureValues))
                                {
                                    // find the vertices that have this action on
                                    for (int actionCombIndex = 0; actionCombIndex < Math.Pow(2, _numOfActionsControlledDynamically); ++actionCombIndex)
                                    {
                                        if (_dynamicallyControlledActionCombinations[actionCombIndex][a] == 1)
                                            // then remove this vertix
                                            statusOfActionCombinationVertices[actionCombIndex] = 0;
                                    }
                                }
                                else
                                {
                                    // find the vertices that have this action off
                                    for (int actionCombIndex = 0; actionCombIndex < Math.Pow(2, _numOfActionsControlledDynamically); ++actionCombIndex)
                                    {
                                        if (_dynamicallyControlledActionCombinations[actionCombIndex][a] == 0)
                                            // then remove this vertix
                                            statusOfActionCombinationVertices[actionCombIndex] = 0;
                                    }
                                }
                            }
                            // return the optimal vertix
                            if (statusOfActionCombinationVertices.Sum() > 0)
                            {
                                for (int actionCombIndex = 0; actionCombIndex < Math.Pow(2, _numOfActionsControlledDynamically); ++actionCombIndex)
                                {
                                    if (statusOfActionCombinationVertices[actionCombIndex] == 1)
                                        optimalDynamicallyControlledActionCombination = (int[])_dynamicallyControlledActionCombinations[actionCombIndex].Clone();
                                }
                            }
                            else
                            {
                                optimalDynamicallyControlledActionCombination = (int[])_defaultActionCombination.Clone();
                                // Error
                            }
                        }
                        else // if (_availabilityOfDynamicallyControlledActionCombinations.Sum() == 0)
                            optimalDynamicallyControlledActionCombination = (int[])_dynamicallyControlledActionCombinations[0].Clone();

                    }
                    break;
                    #endregion
            }

            return optimalDynamicallyControlledActionCombination;
        }

        // find the a feasible action combination using this action combination
        //public int[] FindAFeasibleActionCombinationBasedOnThisActionCombination(int[] thisActionCombination, long timeIndex, long[] arrAvailableResources = null)
        //{
        //    //todo: remove this function 
        //    // check if the initial action is feasible
        //    if (arrAvailableResources == null)
        //    {
        //        if (IfThisActionCombinationIsFeasible(thisActionCombination, timeIndex) == true)
        //            return thisActionCombination;
        //    }
        //    else
        //    {
        //        if (IfThisActionCombinationIsFeasible(thisActionCombination, timeIndex, arrAvailableResources) == true)
        //            return thisActionCombination;
        //    }

        //    // if the initial action combination is not feasible then make it feasible
        //    int[] tempActionCombination = (int[])thisActionCombination.Clone();

        //    foreach (SimulationAction thisAction in _actions)
        //    {
        //        if (tempActionCombination[thisAction.ID] == 0)
        //        {
        //            if (thisAction.Type == SimulationAction.enumActionType.Default)
        //                tempActionCombination[thisAction.ID] = 1;
        //            if (thisAction.RemainsOnOnceSwitchedOn == true && thisAction.NumOfDecisionPeriodsOverWhichThisInterventionWasUsed > 0)
        //                tempActionCombination[thisAction.ID] = 1;
        //        }
        //        else if (tempActionCombination[thisAction.ID] == 1)
        //        {
        //            if (thisAction.Type == SimulationAction.enumActionType.Default)
        //            {
        //                tempActionCombination[thisAction.ID] = 1;
        //            }
        //            else if (thisAction.OnOffSwitchSetting == SimulationAction.enumOnOffSwitchSetting.Dynamic ||
        //                thisAction.OnOffSwitchSetting == SimulationAction.enumOnOffSwitchSetting.Predetermined)
        //            {
        //                if (timeIndex < thisAction.TimeIndexBecomeAvailable || timeIndex >= thisAction.TimeIndexBecomeUnavailable)
        //                    tempActionCombination[thisAction.ID] = 0;

        //                //if (arrAvailableResources != null)
        //                //    if (thisAction.HasResourceCondition && arrAvailableResources[thisAction.IDOfTheResourceRequiredToBeAvailable] <= 0)
        //                //        tempActionCombination[thisAction.ID] = 0;
        //            }
        //        }
        //    }
        //    return (int[])tempActionCombination.Clone();
        //}
        // find the optimal action combination
        //public int[] FindTheOptimalActionCombination(double[] arrObservationFeatureValues)//, long timeIndex, long[] arrAvailableResources = null)
        //{
        //    int[] optimalActionCombination = new int[_numOfActions];
        //    _epsilonGreedyDecisionSelectedAmongThisManyAlternatives = 0;

        //    switch (_qFunctionApproximationMethod)
        //    {
        //        case enumQFunctionApproximationMethod.Q_Approximation:
        //            #region Q-Approximation
        //            {
        //                double max = double.MinValue, qValue = 0;
        //                // for each available action
        //                int i = 0;
        //                foreach (QFunction thisQFunction in _colOfQFunctions)
        //                {
        //                    // build an action combiniation
        //                    int[] thisActionCombination = _actionCombinations[_actionCombinationIndecesOfAvailableActionCombinations[i]];
        //                    // check if this combination is feasible
        //                    if (IfThisActionCombinationIsFeasible(thisActionCombination, timeIndex, arrAvailableResources))
        //                    {
        //                        // increment the number of feasible alternatives
        //                        ++_epsilonGreedyDecisionSelectedAmongThisManyAlternatives;
        //                        // find the value of this action
        //                        qValue = TransformBack(thisQFunction.fValue(arrObservationFeatureValues), _responseTransformation);
        //                        if (qValue > max)
        //                        {
        //                            max = qValue;
        //                            optimalActionCombination = (int[])thisActionCombination.Clone();
        //                        }
        //                    }
        //                    ++i;
        //                }                        
        //            }
        //            break;
        //            #endregion
        //        case enumQFunctionApproximationMethod.Additive_Approximation:
        //            #region Additive-Approximation
        //            {
        //                double max = double.MinValue, qValue = 0;
        //                // for each possibel combination                        
        //                for (int actionCombIndex = 0; actionCombIndex < Math.Pow(2, _numOfActions); ++actionCombIndex)
        //                {
        //                    // build an action combiniation
        //                    int[] thisActionCombination = _actionCombinations[actionCombIndex];
        //                    // check if this combination is feasible
        //                    if (IfThisActionCombinationIsFeasible(thisActionCombination, timeIndex, arrAvailableResources))
        //                    {
        //                        // increment the number of feasible alternatives
        //                        ++_epsilonGreedyDecisionSelectedAmongThisManyAlternatives;
        //                        // find the value of this action
        //                        qValue = TransformBack(_qFunctionApproximationModel_Additive.fValue(
        //                            GetSwitchStatusOfActionsControlledDynamically(thisActionCombination), arrObservationFeatureValues), _responseTransformation);                                
        //                        if (qValue > max)
        //                        {
        //                            max = qValue;
        //                            optimalActionCombination = (int[])thisActionCombination.Clone();
        //                        }
        //                    }
        //                }
        //            }                    
        //            break;
        //            #endregion
        //        case enumQFunctionApproximationMethod.H_Approximation:
        //            #region H-Approximation
        //            {
        //                // start with assuming that all vertices are in the optimal set (1: if yes, 0: if no)
        //                int[] statusOfActionCombinationVertices= new int[(int)Math.Pow(2, _numOfActions)];
        //                SupportFunctions.MakeArrayEqualTo(ref statusOfActionCombinationVertices, 1);

        //                // remove the vertices that are not feasible
        //                for (int actionCombIndex = 0; actionCombIndex < Math.Pow(2, _numOfActions); ++actionCombIndex)
        //                {
        //                    // build an action combiniation
        //                    int[] thisActionCombination = _actionCombinations[actionCombIndex];
        //                    // check if this combination is feasible
        //                    if (IfThisActionCombinationIsFeasible(thisActionCombination, timeIndex, arrAvailableResources))
        //                        // increment the number of feasible alternatives
        //                        ++_epsilonGreedyDecisionSelectedAmongThisManyAlternatives;
        //                    else
        //                        // eliminate this action combination
        //                        statusOfActionCombinationVertices[actionCombIndex] = 0;
        //                }

        //                // eliminate vertices that are not in the optimal set 
        //                if (_epsilonGreedyDecisionSelectedAmongThisManyAlternatives > 1)
        //                {                            
        //                    for (int a = 0; a < _numOfActionsControlledDynamically; a++)
        //                    {
        //                        // if turning off this action lead to better outcome
        //                        if (((QFunction)_colOfHFunctions_Off[a]).fValue(arrObservationFeatureValues) >=
        //                            ((QFunction)_colOfHFunctions_On[a]).fValue(arrObservationFeatureValues))
        //                        {
        //                            // find the vertices that have this action on
        //                            for (int actionCombIndex = 0; actionCombIndex < Math.Pow(2, _numOfActions); ++actionCombIndex)
        //                            {
        //                                if (_actionCombinations[actionCombIndex][_indicesOfActionsControlledDynamically[a]] == 1)
        //                                    // then remove this vertix
        //                                    statusOfActionCombinationVertices[actionCombIndex] = 0;
        //                            }
        //                        }
        //                        else
        //                        {
        //                            // find the vertices that have this action off
        //                            for (int actionCombIndex = 0; actionCombIndex < Math.Pow(2, _numOfActions); ++actionCombIndex)
        //                            {
        //                                if (_actionCombinations[actionCombIndex][_indicesOfActionsControlledDynamically[a]] == 0)
        //                                    // then remove this vertix
        //                                    statusOfActionCombinationVertices[actionCombIndex] = 0;
        //                            }
        //                        }
        //                    }
        //                }                    

        //                // return the optimal vertix
        //                if (statusOfActionCombinationVertices.Sum() > 0)
        //                {
        //                    for (int actionCombIndex = 0; actionCombIndex < Math.Pow(2, _numOfActions); ++actionCombIndex)
        //                    {
        //                        if (statusOfActionCombinationVertices[actionCombIndex] == 1)
        //                            optimalActionCombination = (int[])_actionCombinations[actionCombIndex].Clone();
        //                    }
        //                }
        //                else
        //                {
        //                    optimalActionCombination = FindAFeasibleActionCombinationBasedOnThisActionCombination
        //                        (_defaultActionCombination, timeIndex);
        //                    // Error
        //                }
        //            }
        //            break;
        //            #endregion
        //    }

        //    // if no feasible intervention combination found
        //    if (_epsilonGreedyDecisionSelectedAmongThisManyAlternatives == 0)
        //        optimalActionCombination = FindAFeasibleActionCombinationBasedOnThisActionCombination(_defaultActionCombination, timeIndex, arrAvailableResources);

        //    return optimalActionCombination;
        //}                
        // find an epsilon greedy decision        
        public int[] EpsilonGreedyActionCombination
            (RNG rng, double[] arrObservationFeatureValues)//, long timeIndex, long[] arrAvailableResources = null)
        {
            int[] anActionCombination;
            // with probability of epsilon, make a random decision
            if (rng.NextDouble() <= _epsilonGreedy)
                anActionCombination = GetARandomActionCombinationAmongAvailableDynamicallyControlledActionCombinations(rng);//, timeIndex, arrAvailableResources);
            else // make a greedy decision
                anActionCombination = FindTheOptimalDynamicallyControlledActionCombination(arrObservationFeatureValues);//, timeIndex, arrAvailableResources);

            return anActionCombination;
        }
        // find a Boltzmann decision        
        public int[] BoltzmannActionCombination
            (RNG rng, double[] arrObservationFeatureValues, long timeIndex, long[] arrAvailableResources = null)
        {
            int[] anActionCombination;
            double[] expQOverT = new double[(int)Math.Pow(2, _numOfActionsControlledDynamically)];

            // with probability of epsilon, make a random decision
            if (rng.NextDouble() <= _epsilonGreedy)
                anActionCombination = GetARandomActionCombinationAmongAvailableDynamicallyControlledActionCombinations(rng);//, timeIndex, arrAvailableResources);
            else // make a greedy decision
                anActionCombination = FindTheOptimalDynamicallyControlledActionCombination(arrObservationFeatureValues);//, timeIndex, arrAvailableResources);

            return anActionCombination;
        }

        //--------------------------
        // set up approximation model
        public void SetUpQFunctionApproximationModel(
            enumQFunctionApproximationMethod qFunctionApproximationMethod, enumResponseTransformation responseTransformation,
            int numOfObservationFeatures, int polynomialDegree, int multiplyNumOfColumnsByThisFactorToBeginTraining)
        {
            _qFunctionApproximationMethod = qFunctionApproximationMethod;
            _responseTransformation = responseTransformation;

            switch (_qFunctionApproximationMethod)
            {
                case enumQFunctionApproximationMethod.Q_Approximation:
                    {
                        _colOfQFunctions = new ArrayList();
                    }
                    break;
                case enumQFunctionApproximationMethod.A_Approximation:
                    {
                    }
                    break;
                case enumQFunctionApproximationMethod.H_Approximation:
                    {
                        _colOfHFunctions_Off = new ArrayList();
                        _colOfHFunctions_On = new ArrayList();
                    }
                    break;
            }

            switch (_qFunctionApproximationMethod)
            {
                case enumQFunctionApproximationMethod.Q_Approximation:
                    {
                        //_indecesOfInitiallyAvailableActionCombinations = new int[0];
                        //_inverseOf_indecesOfInitiallyAvailableActionCombinations = new int[(int)Math.Pow(2, _numOfActions)];
                        for (int i = 0; i < Math.Pow(2, _numOfActionsControlledDynamically); i++)
                        {
                            ////if (IfThisActionCombinationIsFeasible(_actionCombinations[actionCombinationIndex]))
                            //if (_ifThisActionCombinationAvailable[actionCombinationIndex] == 1)
                            //{
                            //    SupportFunctions.AddToEndOfArray(ref _indecesOfInitiallyAvailableActionCombinations, actionCombinationIndex);
                            //    _inverseOf_indecesOfInitiallyAvailableActionCombinations[actionCombinationIndex] = i++;
                            // define the Q-function
                            _colOfQFunctions.Add(new PolynomialQFunction("Action combination index for dynamically controlled actions: " + i, 0, numOfObservationFeatures, polynomialDegree, multiplyNumOfColumnsByThisFactorToBeginTraining));
                            //}
                        }
                    }
                    break;
                case enumQFunctionApproximationMethod.A_Approximation:
                    {
                        _qFunctionApproximationModel_Additive = new PolynomialQFunction
                            ("Additive approximation function", _numOfActionsControlledDynamically, numOfObservationFeatures, polynomialDegree, multiplyNumOfColumnsByThisFactorToBeginTraining);
                    }
                    break;
                case enumQFunctionApproximationMethod.H_Approximation:
                    {
                        // set up H functions for actions that are guided by the policy
                        foreach (SimulationAction thisAction in _actions)
                        {
                            if (thisAction.OnOffSwitchSetting == SimulationAction.enumOnOffSwitchSetting.Dynamic)
                            {
                                _colOfHFunctions_Off.Add(new PolynomialQFunction(thisAction.Name + " - Off", 0, numOfObservationFeatures, polynomialDegree, multiplyNumOfColumnsByThisFactorToBeginTraining));
                                _colOfHFunctions_On.Add(new PolynomialQFunction(thisAction.Name + " - On", 0, numOfObservationFeatures, polynomialDegree, multiplyNumOfColumnsByThisFactorToBeginTraining));
                            }
                        }
                    }
                    break;
            }
        }
        // add L2 regularization
        public void AddL2Regularization(double penaltyParameter)
        {
            switch (_qFunctionApproximationMethod)
            {
                case enumQFunctionApproximationMethod.Q_Approximation:
                    {
                        foreach (QFunction thisQFunction in _colOfQFunctions)
                            thisQFunction.AddL2Regularization(penaltyParameter);
                    }
                    break;
                case enumQFunctionApproximationMethod.A_Approximation:
                    _qFunctionApproximationModel_Additive.AddL2Regularization(penaltyParameter);
                    break;
                case enumQFunctionApproximationMethod.H_Approximation:
                    {
                        foreach (QFunction thisQFunction in _colOfHFunctions_On)
                            thisQFunction.AddL2Regularization(penaltyParameter);
                        foreach (QFunction thisQFunction in _colOfHFunctions_Off)
                            thisQFunction.AddL2Regularization(penaltyParameter);
                    }
                    break;
            }
        }
        // setup constant learning rule for all decisions
        public void SetupAConstantStepSizeRule(double stepSize)
        {
            switch (_qFunctionApproximationMethod)
            {
                case enumQFunctionApproximationMethod.Q_Approximation:
                    {
                        foreach (QFunction thisQFunction in _colOfQFunctions)
                            thisQFunction.SetupAConstantStepSizeRule(stepSize);
                    }
                    break;
                case enumQFunctionApproximationMethod.A_Approximation:
                    _qFunctionApproximationModel_Additive.SetupAConstantStepSizeRule(stepSize);
                    break;
                case enumQFunctionApproximationMethod.H_Approximation:
                    {
                        foreach (QFunction thisQFunction in _colOfHFunctions_On)
                            thisQFunction.SetupAConstantStepSizeRule(stepSize);
                        foreach (QFunction thisQFunction in _colOfHFunctions_Off)
                            thisQFunction.SetupAConstantStepSizeRule(stepSize);
                    }
                    break;
            }
        }
        // setup harmonic learning rule for all decisions
        public void SetupAHarmonicStepSizeRule(double harmonic_a)
        {
            switch (_qFunctionApproximationMethod)
            {
                case enumQFunctionApproximationMethod.Q_Approximation:
                    {
                        foreach (QFunction thisQFunction in _colOfQFunctions)
                            thisQFunction.SetupAHarmonicStepSizeRule(harmonic_a);
                    }
                    break;
                case enumQFunctionApproximationMethod.A_Approximation:
                    _qFunctionApproximationModel_Additive.SetupAHarmonicStepSizeRule(harmonic_a);
                    break;
                case enumQFunctionApproximationMethod.H_Approximation:
                    {
                        foreach (QFunction thisQFunction in _colOfHFunctions_On)
                            thisQFunction.SetupAHarmonicStepSizeRule(harmonic_a);
                        foreach (QFunction thisQFunction in _colOfHFunctions_Off)
                            thisQFunction.SetupAHarmonicStepSizeRule(harmonic_a);
                    }
                    break;
            }
        }
        // setup polynomial learning rule for all decisions
        public void SetupAPolynomialStepSizeRule(double polynomial_beta)
        {
            switch (_qFunctionApproximationMethod)
            {
                case enumQFunctionApproximationMethod.Q_Approximation:
                    {
                        foreach (QFunction thisQFunction in _colOfQFunctions)
                            thisQFunction.SetupAPolynomialStepSizeRule(polynomial_beta);
                    }
                    break;
                case enumQFunctionApproximationMethod.A_Approximation:
                    _qFunctionApproximationModel_Additive.SetupAPolynomialStepSizeRule(polynomial_beta);
                    break;
                case enumQFunctionApproximationMethod.H_Approximation:
                    {
                        foreach (QFunction thisQFunction in _colOfHFunctions_On)
                            thisQFunction.SetupAPolynomialStepSizeRule(polynomial_beta);
                        foreach (QFunction thisQFunction in _colOfHFunctions_Off)
                            thisQFunction.SetupAPolynomialStepSizeRule(polynomial_beta);
                    }
                    break;
            }
        }
        // setup epsilon greedy exploration rule
        public void SetupAnEpsilonGreedyExplorationRule(double epsilonGreedy_beta, double epsilonGreedy_delta)
        {
            _epsilonGreedy_beta = epsilonGreedy_beta;
            _epsilonGreedy_delta = epsilonGreedy_delta;
        }
        // update epsilon greedy
        public void UpdateEpsilonGreedy(long itr)
        {
            _epsilonGreedy = (1 - _epsilonGreedy_delta) * Math.Pow(itr, -_epsilonGreedy_beta) + _epsilonGreedy_delta;
        }
        // set up a Boltzmann exploration rule
        public void SetupABoltzmannExplorationRule(double boltzmann_initialTemperature, double boltzmann_beta)
        {
            _Boltzmann_initialTemperature = boltzmann_initialTemperature;
            _Boltzmann_beta = boltzmann_beta;
        }
        // update Boltzmann exploration temperature
        public void UpdateBoltzmannTemperature(long itr)
        {
            _BoltzmannTemperature = _Boltzmann_initialTemperature * Math.Exp(-_Boltzmann_beta * itr);
        }

        // backpropagation 
        public void DoBackpropagation(int itr, double discountFactor, bool stoppedDueToEradication,
            bool useDecisionsAsFeature)
        {
            // do back propagation for each simulation iterations
            for (int dim = 0; dim < _numOfSimRunsToBackPropogate; ++dim)
            {
                ArrayList thisCollectionOfADPStates;
                if (_ADPStates != null)
                    thisCollectionOfADPStates = _ADPStates;
                else
                    thisCollectionOfADPStates = (ArrayList)_ADPStateCollections[dim];

                int numOfADPStates = thisCollectionOfADPStates.Count;
                // is there any ADP state to process?
                if (numOfADPStates == 0)
                    _backPropogationResults[dim] = false;
                else
                {
                    _backPropogationResults[dim] = true;
                    // array of errors
                    //_ADPPredictionErrors[dim] = new double[numOfADPStates];
                    double[] thisPredictionErrorsForEligibleStates = new double[0];

                    // get the last ADP state-decision
                    #region last ADP state
                    ADP_State lastADPStateDecision = (ADP_State)thisCollectionOfADPStates[numOfADPStates - 1];
                    // update the reward to go of the last ADP state-decision
                    if (stoppedDueToEradication == true)
                        lastADPStateDecision.RewardToGo = lastADPStateDecision.DecisoinIntervalReward;
                    else // not eradicated
                    {
                        lastADPStateDecision.RewardToGo =
                            EstimatedTransformedRewardToGo(lastADPStateDecision.SelectedNextPeriodActionCombination, lastADPStateDecision.ObservationFeatureValues);
                    }
                    #endregion

                    // do back propagation for the rest of the state-decisions
                    #region calculate other state's reward to go
                    for (int i = thisCollectionOfADPStates.Count - 1; i >= 1; --i)
                    {
                        // get the last ADP state-decision
                        ADP_State thisADPState = (ADP_State)thisCollectionOfADPStates[i - 1];
                        ADP_State nextADPState = (ADP_State)thisCollectionOfADPStates[i];

                        // update the reward to go of this ADP state-decision
                        thisADPState.RewardToGo = thisADPState.DecisoinIntervalReward
                            + discountFactor * nextADPState.RewardToGo;
                    }
                    #endregion

                    // debugging information
                    #region debugging
                    _matOf_IfEliggible_Actions_FeatureValues_Responses = new double[0][];
                    foreach (ADP_State thisADPState in thisCollectionOfADPStates)
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
                            SupportFunctions.AddToEndOfArray(ref thisRowOf_Actions_FeatureValues_Responses, thisADPState.RewardToGo);
                            // concatinate
                            _matOf_IfEliggible_Actions_FeatureValues_Responses =
                                SupportFunctions.ConcatJaggedArray(_matOf_IfEliggible_Actions_FeatureValues_Responses, thisRowOf_Actions_FeatureValues_Responses);
                        }
                    }
                    #endregion

                    // find errors
                    #region errors
                    int ADPStateIndex = 0;
                    foreach (ADP_State thisADPState in thisCollectionOfADPStates)
                    {
                        if (thisADPState.ValidStateToUpdateQFunction)
                        {
                            double transformedObservedRewardToGo = Transform(thisADPState.RewardToGo, _responseTransformation);

                            // find the estimate reward to go
                            double transformedEstimatedRewardToGo =
                                EstimatedTransformedRewardToGo(thisADPState.SelectedNextPeriodActionCombination, thisADPState.ObservationFeatureValues);

                            SupportFunctions.AddToEndOfArray(ref thisPredictionErrorsForEligibleStates, transformedObservedRewardToGo - transformedEstimatedRewardToGo);
                            ++ADPStateIndex;
                        }
                    }
                    _ADPPredictionErrors[dim] = thisPredictionErrorsForEligibleStates;
                    #endregion

                    // update Q-functions
                    #region update q-functions
                    for (int i = thisCollectionOfADPStates.Count - 1; i >= 0; i--)
                    {
                        // get the ADP State
                        ADP_State thisADPStateDecision = (ADP_State)thisCollectionOfADPStates[i];

                        // first check if this ADP state can be used to update Q-functions
                        if (thisADPStateDecision.ValidStateToUpdateQFunction)
                        {
                            int[] switchStatusOfActionsDynamicallyControlled = GetSwitchStatusOfActionsControlledDynamically(thisADPStateDecision.SelectedNextPeriodActionCombination);

                            switch (_qFunctionApproximationMethod)
                            {
                                case enumQFunctionApproximationMethod.Q_Approximation:
                                    {
                                        ((QFunction)_colOfQFunctions[SupportFunctions.ConvertToBase10FromBase2(switchStatusOfActionsDynamicallyControlled)])
                                            .Update(thisADPStateDecision.ObservationFeatureValues, Transform(thisADPStateDecision.RewardToGo, _responseTransformation), itr);
                                    }
                                    break;
                                case enumQFunctionApproximationMethod.A_Approximation:
                                    _qFunctionApproximationModel_Additive.Update(
                                        switchStatusOfActionsDynamicallyControlled,
                                        thisADPStateDecision.ObservationFeatureValues,
                                        Transform(thisADPStateDecision.RewardToGo, _responseTransformation), itr);
                                    break;
                                case enumQFunctionApproximationMethod.H_Approximation:
                                    {
                                        for (int a = 0; a < _numOfActionsControlledDynamically; a++)
                                        {
                                            // if this vertex is on
                                            if (switchStatusOfActionsDynamicallyControlled[a] == 1)
                                                ((PolynomialQFunction)_colOfHFunctions_On[a]).Update(
                                                    thisADPStateDecision.ObservationFeatureValues,
                                                    Transform(thisADPStateDecision.RewardToGo, _responseTransformation), itr);
                                            else // if this vertex is off
                                                ((PolynomialQFunction)_colOfHFunctions_Off[a]).Update(
                                                    thisADPStateDecision.ObservationFeatureValues,
                                                    Transform(thisADPStateDecision.RewardToGo, _responseTransformation), itr);
                                        }
                                    }
                                    break;
                            }
                        }
                    }
                    #endregion

                } // end if (numOfADPStates == 0)                
            } // end for (int dim = 0; dim < _numOfSimRunsToBackPropogate; ++dim)
        }
        // return back-propagation result
        public bool BackPropagationResult()
        {
            return _backPropogationResults[0];
        }
        public bool BackPropagationResult(int dimension)
        {
            return _backPropogationResults[dimension];
        }

        // num of ADP states
        public int NumberOfADPStates()
        {
            return _ADPStates.Count;
        }
        public int NumberOfADPStates(int dim)
        {
            return ((ArrayList)_ADPStateCollections[dim]).Count;
        }
        // num of valid ADP states
        public int NumberOfValidADPStates()
        {
            return _ADPPredictionErrors[0].Length;
        }
        public int NumberOfValidADPStates(int dim)
        {
            return _ADPPredictionErrors[dim].Length;
        }
        // get selected next period action combination of an adp state
        public int[] GetSelectedNextPeriodActionCombinationOfAnADPState(int ADPStateIndexInCollection)
        {
            return (int[])((ADP_State)_ADPStates[ADPStateIndexInCollection]).SelectedNextPeriodActionCombination.Clone();
        }
        public int[] GetSelectedNextPeriodActionCombinationOfAnADPState(int dimension, int ADPStateIndexInCollection)
        {
            return (int[])((ADP_State)((ArrayList)_ADPStateCollections[dimension])[ADPStateIndexInCollection]).SelectedNextPeriodActionCombination.Clone();
        }

        // add to a decision interval reward
        public void AddToDecisionIntervalReward(int ADPStateIndexInCollection, double reward)
        {
            ((ADP_State)_ADPStates[ADPStateIndexInCollection]).AddToDecisionIntervalReward(reward);
        }
        public void AddToDecisionIntervalReward(int dimension, int ADPStateIndexInCollection, double reward)
        {
            ((ADP_State)((ArrayList)_ADPStateCollections[dimension])[ADPStateIndexInCollection]).AddToDecisionIntervalReward(reward);
        }
        // get reward-to-go
        public double GetRewardToGo(int ADPStateIndexInCollection)
        {
            return ((ADP_State)_ADPStates[ADPStateIndexInCollection]).RewardToGo;
        }
        public double GetRewardToGo(int dimension, int ADPStateIndexInCollection)
        {
            return ((ADP_State)((ArrayList)_ADPStateCollections[dimension])[ADPStateIndexInCollection]).RewardToGo;
        }
        // prediction errors
        public double ADPPredictionErrors(int ADPStateIndex)
        {
            if (_ADPPredictionErrors[0].Length == 0)
                return 0;
            else
                return _ADPPredictionErrors[0][ADPStateIndex];
        }
        public double ADPPredictionErrors(int dimension, int ADPStateIndex)
        {
            if (_ADPPredictionErrors[dimension].Length == 0)
                return 0;
            else
                return _ADPPredictionErrors[dimension][ADPStateIndex];
        }
        public double PredictionErrorForTheFirstEligibleADPState(int dimension)
        {
            double result = 0;

            ArrayList thisCollectionOfADPStates = (ArrayList)_ADPStateCollections[dimension];
            for (int ADPStateIndex = 0; ADPStateIndex < thisCollectionOfADPStates.Count; ++ADPStateIndex)
            {
                if (((ADP_State)thisCollectionOfADPStates[ADPStateIndex]).ValidStateToUpdateQFunction)
                {
                    result = _ADPPredictionErrors[dimension][ADPStateIndex];
                    return result;
                }
            }
            return result;
        }
        public double PredictionErrorForTheLastEligibleADPState(int dimension)
        {
            double result = 0;

            ArrayList thisCollectionOfADPStates = (ArrayList)_ADPStateCollections[dimension];
            for (int ADPStateIndex = thisCollectionOfADPStates.Count - 1; ADPStateIndex >= 0; --ADPStateIndex)
            {
                if (((ADP_State)thisCollectionOfADPStates[ADPStateIndex]).ValidStateToUpdateQFunction)
                {
                    result = _ADPPredictionErrors[dimension][ADPStateIndex];
                    return result;
                }
            }
            return result;
        }

        // reset for another simulation run
        public void ResetForAnotherSimulationRun()//(ref double initialCost)
        {
            if (_ADPStates != null)
                _ADPStates.Clear();
            if (_ADPStateCollections != null)
                foreach (ArrayList thisArrayList in _ADPStateCollections)
                    thisArrayList.Clear();

            // reset the number of times each intervention is used
            foreach (SimulationAction thisAction in _actions)
            {
                thisAction.NumOfSwitchesOccured = 0;
                thisAction.NumOfDecisionPeriodsOverWhichThisInterventionWasUsed = 0;
            }

            // find a feasible action
            _currentActionCombination = (int[])_defaultActionCombination.Clone();
            //// update costs
            //_currentCostPerUnitOfTime = 0;
            //foreach (SimulationAction thisAction in _actions)
            //{
            //    // find costs
            //    if (_currentActionCombination[thisAction.Index] == 1)
            //    {
            //        // fixed cost
            //        initialCost += thisAction.FixedCost;
            //        // variable cost
            //        _currentCostPerUnitOfTime += thisAction.CostPerUnitOfTime;
            //    }
            //}
        }
        // clear all
        public void Clear()
        {
            _actions.Clear();
            if (_ADPStates != null)
                _ADPStates.Clear();
            if (_ADPStateCollections != null)
                foreach (ArrayList thisArrayList in _ADPStateCollections)
                    thisArrayList.Clear();
            _numOfActions = 0;
            _numOfActionsControlledDynamically = 0;
        }

        // PRIVATE SUBS
        #region Private Subs
        // get a random action combination among available dynamically controlled action combinations
        private int[] GetARandomActionCombinationAmongAvailableDynamicallyControlledActionCombinations(RNG rng)
        {
            int index = rng.Next(_availabilityOfDynamicallyControlledActionCombinations.Sum());
            //int sum = 0;
            //int i = 0;
            //while (sum <= rnd)
            //{
            //    if (_ifThisActionCombinationAvailable[i] == 1)
            //        ++sum;
            //    ++i;
            //}
            return _dynamicallyControlledActionCombinations[index];
        }

        //// get if this is the default action combination 
        //private bool IfThisIsTheDefaultActionCombination(int[] actionCombination)
        //{
        //    bool isFeasible = true;
        //    int actionIndex = 0;
        //    foreach (SimulationAction thisAction in _actions)
        //    {
        //        if (thisAction.Type == SimulationAction.enumActionType.Default && actionCombination[actionIndex] == 0)
        //            return false;
        //        ++actionIndex;
        //    }
        //    return isFeasible;
        //}
        //private bool IfThisActionCombinationIsFeasible(int[] actionCombination, long timeIndex)
        //{
        //    bool isFeasible = true;

        //    // check if this combination is feasible regardless of resource availability
        //    isFeasible = IfThisActionCombinationIsFeasible(actionCombination);
        //    if (isFeasible == false) return false;

        //    int actionIndex = 0;
        //    foreach (SimulationAction thisAction in _actions)
        //    {
        //        switch (thisAction.OnOffSwitchSetting)
        //        {
        //            case SimulationAction.enumOnOffSwitchSetting.Predetermined:
        //                {
        //                    if (actionCombination[actionIndex] == 1)
        //                    {
        //                        if (timeIndex < thisAction.TimeIndexBecomeAvailable || timeIndex >= thisAction.TimeIndexBecomeUnavailable)
        //                            return false;
        //                    }
        //                }
        //                break;
        //            case SimulationAction.enumOnOffSwitchSetting.Dynamic:
        //                {
        //                    if (actionCombination[actionIndex] == 1)
        //                    {
        //                        if (timeIndex < thisAction.TimeIndexBecomeAvailable || timeIndex >= thisAction.TimeIndexBecomeUnavailable)
        //                            return false;
        //                    }                            
        //                }
        //                break;
        //        }

        //        ++actionIndex;
        //    }
        //    return isFeasible;
        //}
        //private bool IfThisActionCombinationIsFeasible(int[] actionCombination, long timeIndex, long[] arrAvailableResources)
        //{
        //    bool isFeasible = true;

        //    // check if this combination is feasible regardless of resource availability
        //    isFeasible = IfThisActionCombinationIsFeasible(actionCombination, timeIndex);
        //    if (isFeasible == false) return false;

        //    // if no resource information is provided
        //    if (arrAvailableResources == null || arrAvailableResources.Length == 0)
        //        return isFeasible;
        //    else 
        //    {
        //        // for now, let's assume that the availability of resources does not affect the feasibility of actions
        //        //foreach (SimulationAction thisAction in _actions)
        //        //{
        //        //    if (thisAction.HasResourceCondition && actionCombination[thisAction.ID] == 1 &&
        //        //        arrAvailableResources[thisAction.IDOfTheResourceRequiredToBeAvailable] <= 0)
        //        //        return false;                   
        //        //}
        //    }            
        //    return isFeasible;
        //}

        // get the switch status of actions that are controlled dynamically
        private int[] GetSwitchStatusOfActionsControlledDynamically(int[] actionCombination)
        {
            int[] result = new int[_numOfActionsControlledDynamically];

            for (int i = 0; i < _numOfActionsControlledDynamically; i++)
                result[i] = actionCombination[_indicesOfActionsControlledDynamically[i]];

            //int i = 0;
            //foreach (SimulationAction thisAction in _actions)
            //{
            //    if (thisAction.OnOffSwitchSetting == SimulationAction.enumOnOffSwitchSetting.Dynamic)
            //        result[i++] = actionCombination[thisAction.Index];
            //}

            return (int[])result.Clone();
        }

        // transformation
        private double Transform(double value, enumResponseTransformation responseTransformation)
        {
            double result = 0;
            switch (responseTransformation)
            {
                case enumResponseTransformation.None:
                    result = value;
                    break;
                case enumResponseTransformation.NaturalLog_PositiveArgument:
                    result = Math.Log(Math.Max(0, value) + 1);
                    break;
                case enumResponseTransformation.NaturalLog_NegativeArgument:
                    result = Math.Log(Math.Max(0, -value) + 1);
                    break;
                case enumResponseTransformation.SquaredRoot_PositiveArgument:
                    result = Math.Sqrt(Math.Max(0, value));
                    break;
                case enumResponseTransformation.SquaredRoot_NegativeArgument:
                    result = Math.Sqrt(Math.Max(0, -value));
                    break;
            }

            return result;
        }
        private double TransformBack(double transformedValue, enumResponseTransformation responseTransformation)
        {
            double result = 0;
            switch (responseTransformation)
            {
                case enumResponseTransformation.None:
                    result = transformedValue;
                    break;
                case enumResponseTransformation.NaturalLog_PositiveArgument:
                    result = Math.Max(0, Math.Exp(transformedValue) - 1);
                    break;
                case enumResponseTransformation.NaturalLog_NegativeArgument:
                    result = Math.Min(0, 1 - Math.Exp(transformedValue));
                    break;
                case enumResponseTransformation.SquaredRoot_PositiveArgument:
                    result = Math.Pow(transformedValue, 2);
                    break;
                case enumResponseTransformation.SquaredRoot_NegativeArgument:
                    result = -Math.Pow(transformedValue, 2);
                    break;
            }

            return result;
        }

        // return the estimated reward to go
        private double EstimatedTransformedRewardToGo(int[] selectedNextPeriodActionCombination, double[] observationFeatureValues)
        {
            double estimatedTransformedRewardToGo = 0;
            int[] switchStatusOfActionsControlledDynamically = GetSwitchStatusOfActionsControlledDynamically(selectedNextPeriodActionCombination);

            // find the estimate reward to go
            switch (_qFunctionApproximationMethod)
            {
                case enumQFunctionApproximationMethod.Q_Approximation:
                    {
                        int actionCombIndex = SupportFunctions.ConvertToBase10FromBase2(selectedNextPeriodActionCombination);
                        estimatedTransformedRewardToGo =
                            ((QFunction)_colOfQFunctions[SupportFunctions.ConvertToBase10FromBase2(switchStatusOfActionsControlledDynamically)])
                            .fValue(observationFeatureValues);
                    }
                    break;
                case enumQFunctionApproximationMethod.A_Approximation:
                    estimatedTransformedRewardToGo = _qFunctionApproximationModel_Additive.fValue(switchStatusOfActionsControlledDynamically, observationFeatureValues);
                    break;
                case enumQFunctionApproximationMethod.H_Approximation:
                    {
                        for (int i = 0; i < _numOfActionsControlledDynamically; i++)
                        {
                            // if this vertex is on
                            if (switchStatusOfActionsControlledDynamically[i] == 1)
                                estimatedTransformedRewardToGo += ((PolynomialQFunction)_colOfHFunctions_On[i]).fValue(observationFeatureValues);
                            else // if this vertex is off
                                estimatedTransformedRewardToGo += ((PolynomialQFunction)_colOfHFunctions_Off[i]).fValue(observationFeatureValues);
                        }
                        estimatedTransformedRewardToGo = estimatedTransformedRewardToGo / _numOfActionsControlledDynamically;
                    }
                    break;
            }

            return estimatedTransformedRewardToGo;

        }
        #endregion
    }
}
