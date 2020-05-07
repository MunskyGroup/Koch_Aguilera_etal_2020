function [fit_Percent,fit_Distribution,fit_HT] = fun_parameterUncertainty(x,strExperimentalData, sequenceParameters,Nrep,N_sims,evaluationType,modelReduction)
plottingCondition =0; folderName =[];
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
param.reportPromoter =0;
if evaluationType ==1
    %% Simulation to generate intiensities
    delta_t = 10; t_max = 50*60; % maximum simulation time. The burnin time is automatically added in the SSA.
    [CAP_int,IRES_int] = fun_initialSimulations(delta_t,t_max,N_sims,param,sequenceParameters,modelReduction);
    %% Calculating percentage of spots per frame
    [fit_Percent,~,~,~,~] = percentagePerFrame(N_sims, folderName, CAP_int,IRES_int,plottingCondition,strExperimentalData);
    %% Calculating Intensity Distributions
    [fit_Distribution,~, ~] = comparing_IntensityDistributions(CAP_int,IRES_int,strExperimentalData,folderName,plottingCondition);
    fit_HT=0;
elseif evaluationType ==2
    %% Harringtonine assays.
    [fit_HT, ~, ~]  = fun_Harringtonine(Nrep,N_sims,param,sequenceParameters,folderName, strExperimentalData,plottingCondition,modelReduction);
    fit_Percent=0;
    fit_Distribution=0;
end
