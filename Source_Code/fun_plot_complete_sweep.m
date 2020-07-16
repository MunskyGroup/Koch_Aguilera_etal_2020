function fitReport = fun_plot_complete_sweep (x,strExperimentalData, sequenceParameters,plottingCondition,optimizationMethod,Nrep,N_sims,modelReduction)
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
param.crossover = x(12); % cross-over
param.reportPromoter =0;
param.k_inhi1 = 0;   % inhibitor type I
param.k_inhi2 = 0;   % inhibitor type II

% nameFile = 'z_parameters.mat';
% save (nameFile, 'param')
% movefile (nameFile, folderName)
folderName =[];

%% Plotting a Time course
%delta_t =1; t_max= 6e4;
%timeCourse(delta_t,t_max,param,sequenceParameters,plottingCondition,folderName,'');
%% Simulation to generate intiensities
delta_t = 10; t_max = 1e5; % maximum simulation time. The burnin time is automatically added in the SSA.
[CAP_int,IRES_int] = fun_initialSimulations(delta_t,t_max,N_sims,param,sequenceParameters,modelReduction);
%% Calculating percentage of spots per frame
[LL_Percent,~,~,~,~] = percentagePerFrame(N_sims, folderName, CAP_int,IRES_int,plottingCondition,strExperimentalData);
%% Calculating Intensity Distributions
[fit_LL_Distribution,~, ~] = comparing_IntensityDistributions(CAP_int,IRES_int,strExperimentalData,folderName,plottingCondition);
%% Harringtonine assays.
[LL_HT, ~, ~]  = fun_Harringtonine(Nrep,N_sims,param,sequenceParameters,folderName, strExperimentalData,plottingCondition,modelReduction);
%% NaAs stress
% NaAs stress type I
N_parameters_to_Test = 3;
value_sweep = [0:1/N_parameters_to_Test:1];
Nsweep = length(value_sweep);

for k=1:Nsweep
    param.k_inhi1 = value_sweep(k);   % inhibitor type I
    param.k_inhi2 = value_sweep(k);   % inhibitor type II
    param.InhibitorType =1;
    [LL_NA_1(1,k), ~, ~, ~]  = fun_NaAS(Nrep,N_sims,param,sequenceParameters,folderName, strExperimentalData,plottingCondition,modelReduction);
    % NaAs stress type II
    param.InhibitorType =2;
    [LL_NA_2(1,k), ~, ~, ~]  = fun_NaAS(Nrep,N_sims,param,sequenceParameters,folderName, strExperimentalData,plottingCondition,modelReduction);
    % NaAs stress type III
    param.InhibitorType =3;
    [LL_NA_3(1,k), ~, ~, ~]  = fun_NaAS(Nrep,N_sims,param,sequenceParameters,folderName, strExperimentalData,plottingCondition,modelReduction);
    
    %% DTT stress
    % DTT stress type I
    param.InhibitorType =1;
    [LL_DTT_1(1,k), ~, ~, ~]  = fun_DTT(Nrep,N_sims,param,sequenceParameters,folderName, strExperimentalData,plottingCondition,modelReduction);
    % DTT stress type II
    param.InhibitorType =2;
    [LL_DTT_2(1,k), ~, ~, ~]  = fun_DTT(Nrep,N_sims,param,sequenceParameters,folderName, strExperimentalData,plottingCondition,modelReduction);
    % DTT stress type III
    param.InhibitorType =3;
    [LL_DTT_3(1,k), ~, ~, ~]  = fun_DTT(Nrep,N_sims,param,sequenceParameters,folderName, strExperimentalData,plottingCondition,modelReduction);
end

%% Deffining the Objective Function
fitValue(1)=LL_Percent;
fitValue(2)=fit_LL_Distribution;
fitValue(3)=LL_HT;

fitValue=[fitValue,LL_NA_1]; % 4 - 7  % type I. k_on_c = [0,33, 66, 100]
fitValue=[fitValue,LL_NA_2]; % 8 - 11 % type II. k_init_c = [0,33, 66, 100]
fitValue=[fitValue,LL_NA_3]; % 12 - 15 % type III. k_init_c, k_init_i = [0,33, 66, 100]

fitValue=[fitValue,LL_DTT_1]; % 16 - 19 % type I. k_on_c = [0,33, 66, 100]
fitValue=[fitValue,LL_DTT_2]; % 20 - 23 % type II. k_init_c = [0,33, 66, 100]
fitValue=[fitValue,LL_DTT_3]; % 24 - 27 % type III. k_init_c, k_init_i = [0,33, 66, 100]

if optimizationMethod ==1
    fitReport = fitValue;
else
    fitReport =  sum(fitValue);
end

end


