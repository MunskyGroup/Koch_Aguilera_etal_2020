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
modelReduction =1; % 1 = Uses the fast approximaiton. 0 = Uses the complete model.
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
    % mean ke
    param.ke_avg
    % initiation rates were found to be K_i = 1/30 sec^{-1}
    ki_cap = 1/param.ki_CAP
    ki_cap = 1/param.ki_IRES
    % Cap translation would have a duration of 1/K_{C-}
    cap_duration = (1/param.k21)/60
    % a refractory time of 1/K_{C-} = XXX min?
    cap_refractory = (1/param.k12)/60
    % CAP state ON fraction of f_C = K_{C+}/(K_{C+} + K_{C-}) = XXX
    %  CAP_on_frac = param.k12 /(param.k12+param.k21)
    %burst size CAP
    burst_size_CAP = (param.ki_CAP/param.k21)
    % Because single proteins takes XXX minutes to translate
    % mean_elongation_Time = mean ([sequenceParameters.codon_time_CAP, sequenceParameters.codon_time_IRES] );
    % mean_elongation_Time = mean_elongation_Time/60
    % IRES has a much smaller duration of 1/K_{I-} = 22 min
    ires_duration =  (1/param.k31)/60
    % IRES refractory time of 1/K_{I-}
    ires_refractory = (1/param.k13)/60
    % IRES steady state ON fraction of f_I = K_{I+}/(K_{I+} + K_{I-}) = XXX
    IRES_on_frac = (param.k13 /(param.k13+param.k31))
    % burst size IRES
    burst_size_IRES = (param.ki_IRES/param.k31)
    % In the presence of Cap the burst duration increases to 1/K_{IC+} = XXX min
    %IRES_burst_Increase_with_CAP = (1/ param.k42)/60
    % refractory time decreases to 1/K_{CI-} =XXX min
    IRES_ref_time_dec = (1/ param.k24)/60
    % ON fraction increases to \tilde{f}_I = XXX,
    CAP_IRES_on_frac = param.k24/(param.k24+param.k42)
    % Cap-mediated enhancement of IRES activation
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
    % NaAs stress type I
    param.InhibitorType =1;
    fun_NaAS_sweep(N_sims,param,sequenceParameters,folderName, strExperimentalData,plottingCondition,modelReduction);
    % NaAs stress type II
    param.InhibitorType =2;
    fun_NaAS_sweep(N_sims,param,sequenceParameters,folderName, strExperimentalData,plottingCondition,modelReduction);
    % NaAs stress type III
    %param.InhibitorType =3;
    %fun_NaAS_sweep(N_sims,param,sequenceParameters,folderName, strExperimentalData,plottingCondition,modelReduction);
    %% DTT stress
    % DTT stress type I
    param.InhibitorType =1;
    fun_DTT_sweep(N_sims,param,sequenceParameters,folderName, strExperimentalData,plottingCondition,modelReduction);
    % DTT stress type II
    param.InhibitorType =2;
    fun_DTT_sweep(N_sims,param,sequenceParameters,folderName, strExperimentalData,plottingCondition,modelReduction);
    % DTT stress type III
    %param.InhibitorType =3;
    %fun_DTT_sweep(N_sims,param,sequenceParameters,folderName, strExperimentalData,plottingCondition,modelReduction);
end
cd ..