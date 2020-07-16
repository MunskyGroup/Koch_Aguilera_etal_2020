function [LL_HT, LL_HT_1G, LL_HT_2G]  = fun_Harringtonine(Nrep,N_sims,param,sequenceParameters,folderName, strExperimentalData,plottingCondition,modelReduction)
delta_t = 60; % time space
n_timePoints =35;
t_max = n_timePoints*delta_t; % maximum simulation time. The burnin time is automatically added in the SSA.
param.using_harringtonine = 1; param.using_stress = 0;
param.tim_inhibitor = 5*delta_t; % burnin is added in the SSA function
pointOfNormalization =5;
% Pre-alocating memory
sim_CAP_ONLY_HT = zeros(Nrep, n_timePoints);
sim_IRES_ONLY_HT = zeros(Nrep, n_timePoints);
error_CAP_ONLY_HT = zeros(Nrep, n_timePoints);
error_IRES_ONLY_HT = zeros(Nrep, n_timePoints);

for k=1:Nrep
    if modelReduction==1
        parfor i =1:N_sims
            [~,~,CAP_int_HT(i,:),IRES_int_HT(i,:)] =compilingTimes_Intensities (param,sequenceParameters,t_max,delta_t);
        end
    else
        parfor i =1:N_sims
            [RibsomePositions,~] = ssa_completeModel(param,sequenceParameters,t_max,delta_t);
            [CAP_int_HT(i,:),IRES_int_HT(i,:)] = SSAtoIntensity(RibsomePositions,sequenceParameters);
        end
    end
    [sim_CAP_ONLY_HT(k,:), sim_IRES_ONLY_HT(k,:), error_CAP_ONLY_HT(k,:), error_IRES_ONLY_HT(k,:)] = intesities_selection_HT(CAP_int_HT, IRES_int_HT,pointOfNormalization);
end
[LL_HT, LL_HT_1G, LL_HT_2G] = compareRunOffs(folderName, strExperimentalData, sim_CAP_ONLY_HT, sim_IRES_ONLY_HT,error_CAP_ONLY_HT, error_IRES_ONLY_HT, plottingCondition);
end