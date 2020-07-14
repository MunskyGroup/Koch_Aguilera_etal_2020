clear all; close all; clc;
Selected_Model = [13];  % 4S_Im2 %   Model 13. 4 States - IRES activation rates dependent (KI+ ~= KCI+) and inactivation rates independent (KI- = KCI- ).% k13 ~= k24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This function is intended to perform stochastic simulations for the selected model.
% INPUT:
%- Model_type. An integer between 1 and 4. Where the number indicates one of the following models:
% 3S      %  Model 1. 3 States.
% 3S_C    %  Model 2. 3 States + Cross-over.
% 4S_DD   %  Model 3. 4 States.
% 4S_DDC  %  Model 4. 4 States + Cross-over.
% 4S_II   %  Model 5. 4 States - Independent. k12=k34; k21=k43; k13=k24;  k31=% k42;
% 4S_IIC  %  Model 6. 4 States - Independent. k12=k34; k21=k43; k13=k24; k31= k42 +CO
% 4S_DI   %  Model 7. 4 States - CAP is independent of IRES, but IRES depends on CAP. k12=k34; k21=k43;
% 4S_DIC  %  Model 8. 4 States - CAP is independent of IRES, but IRES depends on CAP.  k12=k34; k21=k43; + CO
% 4S_ID   %  Model 9. 4 States - IRES is independent of CAP, but CAP depends on IRES.   k13=k24;  k31= k42;
% 4S_IDC  %  Model 10. 4 States - IRES is independent of CAP, but CAP depends on IRES. k13=k24;  k31= k42; + CO
% 4S_Im1  %  Model 11. 4 States - IRES activation rates dependent (KI+ = KCI+) and inactivation rates dependent (KI- ~= KCI- ).% k31  ~= k42;
% 4S_Im1C %  Model 12. 4 States - IRES activation rates dependent (KI+ = KCI+) and inactivation rates dependent (KI- ~= KCI- ) + CO.% k31  ~= k42;
% 4S_Im2  %  Model 13. 4 States - IRES activation rates dependent (KI+ ~= KCI+) and inactivation rates independent (KI- = KCI- ).% k13 ~= k24
% 4S_Im2C %  Model 14. 4 States - IRES activation rates dependent (KI+ ~= KCI+) and inactivation rates independent (KI- = KCI-) + CO % k13 ~= k24;
% OUTPUT:
% The code generates the folder Selected_Model_Selected_Model. Where _Selected_Model is a placeholder that represents the number of the selected model. In the folder, the simulation results are saved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd Source_Code
modelReduction =1; % 1 = Uses the fast approximaiton. 0 = Uses a complete model (the complete model is was not used in this publication).
for kk=1:length(Selected_Model)
    %% Selecting the model type
    Model_type = Selected_Model(kk);
    if Model_type ==4 % Model 4 is the most complex model.
        folderName = horzcat('Complex_Model_',num2str(Model_type)); mkdir(folderName);
    else
        folderName = horzcat('Selected_Model_',num2str(Model_type)); mkdir(folderName);
    end
    Nrep =10; % 5
    N_sims = 4000;    % number of repetitions. 4000
    plottingCondition =1;
    strExperimentalData = extractingExperimentalData; % loading experimental data
        
    sequenceParameters = fun_sequenceParameters; % loading sequence parameters
    [initPop] = initialPopulation (Model_type);
    param.ki_CAP = initPop(1);
    param.ki_IRES = initPop(2);
    param.ke_avg = initPop(3); % average elongation rate to be optimized
    param.k12 = initPop(4);  % k_on_c
    param.k21 = initPop(5);  % k_off_c
    param.k13 = initPop(6);  % k_on_i
    param.k31 = initPop(7);  % k_off_i
    param.k34 = initPop(8);  % k'_on_c
    param.k43 = initPop(9);  % k'_off_c
    param.k24 = initPop(10); % k'_on_i
    param.k42 = initPop(11); % k'_off_i
    param.crossover = initPop(12);
    param.reportPromoter=1;
    
    % Display parameters in command window
    k_INIT_C= param.ki_CAP
    k_INIT_I= param.ki_IRES
    k_e = param.ke_avg
    k_ON_C = param.k12
    k_OFF_C = param.k21
    k_ON_I = param.k13
    k_OFF_I = param.k31
    k_p_ON_I = param.k24
    
    %% Calculating model biophysical parameters
    mean_elongation = param.ke_avg
    ki_cap = 1/param.ki_CAP
    ki_ires = 1/param.ki_IRES
    cap_duration = (1/param.k21)/60
    cap_refractory = (1/param.k12)/60
    burst_size_CAP = (param.ki_CAP/param.k21)
    ires_duration =  (1/param.k31)/60
    ires_refractory = (1/param.k13)/60
    IRES_on_frac = (param.k13 /(param.k13+param.k31))
    burst_size_IRES = (param.ki_IRES/param.k31)
    IRES_ref_time_dec = (1/ param.k24)/60
    CAP_IRES_on_frac = param.k24/(param.k24+param.k42)
    Cap_enhancement = CAP_IRES_on_frac/IRES_on_frac
        
    %% Simulation to generate intiensities
    delta_t = 10; t_max = 50*60; % maximum simulation time. The burnin time is automatically added in the SSA.
    [CAP_int,IRES_int] = fun_initialSimulations(delta_t,t_max,N_sims,param,sequenceParameters,modelReduction);
    %% Plotting a Time course
    %     %% Calculating percentage of spots per frame
    [LL_Percent,conditionPercentage,fit_Percent_1G,fit_Percent_2G,fit_Percent_BG] = percentagePerFrame(N_sims, folderName, CAP_int,IRES_int,plottingCondition,strExperimentalData);
    %     %% Calculating Intensity Distributions
    [fit_LL_Distribution,fit_LL_Distribution_CAP,fit_LL_Distribution_IRES] = comparing_IntensityDistributions(CAP_int,IRES_int,strExperimentalData,folderName,plottingCondition);
        %%     %% Harringtonine assays.
    [LL_HT, LL_HT_1G, LL_HT_2G]  = fun_Harringtonine(Nrep,N_sims,param,sequenceParameters,folderName, strExperimentalData,plottingCondition,modelReduction);
    
    %% NaAs stress
    % NaAs stress type I - Blocking k_on_c
    param.InhibitorType =1;
    fun_NaAS_sweep(N_sims,param,sequenceParameters,folderName, strExperimentalData,plottingCondition,modelReduction);
    % NaAs stress type II - Blocking k_init_c
    param.InhibitorType =2;
    fun_NaAS_sweep(N_sims,param,sequenceParameters,folderName, strExperimentalData,plottingCondition,modelReduction);
    % NaAs stress type III - Blocking k_init_c and - Blocking k_init_i
    %param.InhibitorType =3;
    %fun_NaAS_sweep(N_sims,param,sequenceParameters,folderName, strExperimentalData,plottingCondition,modelReduction);
    %% DTT stress
    % DTT stress type I - Blocking k_on_c
    param.InhibitorType =1;
    fun_DTT_sweep(N_sims,param,sequenceParameters,folderName, strExperimentalData,plottingCondition,modelReduction);
    % DTT stress type II - Blocking k_init_c
    param.InhibitorType =2;
    fun_DTT_sweep(N_sims,param,sequenceParameters,folderName, strExperimentalData,plottingCondition,modelReduction);
    % DTT stress type III - Blocking k_init_c and - Blocking k_init_i
    %param.InhibitorType =3;
    %fun_DTT_sweep(N_sims,param,sequenceParameters,folderName, strExperimentalData,plottingCondition,modelReduction);
end
cd ..