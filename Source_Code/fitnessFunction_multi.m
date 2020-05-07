function fitReport = fitnessFunction_multi (x,strExperimentalData, sequenceParameters, plottingCondition, folderName, conditionNoSkip, optimizationMethod, Model_type, Nrep, N_sims, modelReduction)
global NoEvaluations EvaluationCondition1 EvaluationCondition2 bestFitValue
NoEvaluations=NoEvaluations+1;

%% This function performs calculates the -Loglikelihood functions used for the parameter optimization.

if optimizationMethod ==3 % this section makes the EvaluationConditions variables equal to 100. To
    EvaluationCondition1 =0; % condition to start evaluating HT experiments after 100 evaluations.
    EvaluationCondition2 =0; % condition to start evaluating NaAs experiments after 100 evaluations.
end
percetageTolerance = 0.1; % percentage to accept as valid parameter and perform additonal simulations (these are HT and NaAs)
param.reportPromoter =0;
param.ki_CAP = x(1);
param.ki_IRES = x(2);
param.ke_avg = x(3); % average elongation rate to be optimized
param.k12 = x(4);
param.k21 = x(5);
param.k13 = x(6);
param.k31 = x(7);
param.k34 = x(8);  % k_34 = k_12;
param.k43 = x(9);  % k_43 = k_21;
param.k24 = x(10); % k_24 = k_13;
param.k42 = x(11);  % k_42 = k_31;
param.crossover = x(12);

% This section associates parameters according to each model.

if Model_type ==5 || Model_type ==6
    param.k34 = param.k12;  % k_34 = k_12;
    param.k43 = param.k21;  % k_43 = k_21;
    param.k24 = param.k13; % k_24 = k_13;
    param.k42 = param.k31;  % k_42 = k_31;
end

if Model_type == 7 ||  Model_type == 8
    param.k34 = param.k12;  % k_34 = k_12;
    param.k43 = param.k21;  % k_43 = k_21;
end

if  Model_type == 9 || Model_type == 10
    param.k24 = param.k13; % k_24 = k_13;
    param.k42 = param.k31;  % k_42 = k_31;
end

if Model_type == 11  || Model_type == 12  % k31  ~= k42;
    param.k34 = param.k12;  % k_34 = k_12;
    param.k43 = param.k21;  % k_43 = k_21;
    param.k24 = param.k13; % k_24 = k_13;
end

if Model_type == 13 || Model_type == 14 % k13 ~= k24;
    param.k34 = param.k12;  % k_34 = k_12;
    param.k43 = param.k21;  % k_43 = k_21;
    param.k42 = param.k31;  % k_42 = k_31;
end

%% Simulation to generate intiensities
delta_t = 10; t_max = 50*60; % maximum simulation time. The burnin time is automatically added in the SSA.
[CAP_int,IRES_int] = fun_initialSimulations(delta_t,t_max,N_sims,param,sequenceParameters,modelReduction);

%% Calculating percentage of spots per frame
[LL_Percent,~,~,~,~] = percentagePerFrame(N_sims, folderName, CAP_int,IRES_int,plottingCondition,strExperimentalData);
%% Calculating Intensity Distributions
[fit_LL_Distribution,~, ~] = comparing_IntensityDistributions(CAP_int,IRES_int,strExperimentalData,folderName,plottingCondition);
%% Calculate ratio
if NoEvaluations>=1 && fit_LL_Distribution< 1e4
    if plottingCondition==1 || conditionNoSkip ==1 || (LL_Percent < bestFitValue(1)*(1+percetageTolerance)) && (fit_LL_Distribution < bestFitValue(2)*(1+percetageTolerance)) % NoEvaluations > EvaluationCondition1
        %% Harringtonine assays.
        [LL_HT, ~, ~]  = fun_Harringtonine(Nrep,N_sims,param,sequenceParameters,folderName, strExperimentalData,plottingCondition,modelReduction);
    else
        LL_HT_1G=1e4; LL_HT_2G=1e4; LL_HT=1e4;
    end
end
LL_NA = 0; % Notice that the NaAs is not part of the optimization

%% Deffining the Objective Function

fitValue(1)=LL_Percent;
fitValue(2)=fit_LL_Distribution;
fitValue(3)=LL_HT;
fitValue(4)=LL_NA;

if optimizationMethod ==1
    fitReport = fitValue;
else
    fitReport =  sum(fitValue);
end

% Saves the result as new bestFitValue if the evaluated parameters are
% better than the previous evaluation
if all(fitValue<bestFitValue)
    bestFitValue  = fitValue;
elseif all(fitValue<bestFitValue) && optimizationMethod==1
    fitReport = fitValue.*Inf;
elseif all(fitValue<bestFitValue) && optimizationMethod~=1
    fitReport = inf;
end

end
