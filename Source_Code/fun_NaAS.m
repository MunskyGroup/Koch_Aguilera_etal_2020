function [LL_NA, LL_NA_1G, LL_NA_2G, LL_NA_CI]  = fun_NaAS(Nrep,N_sims,param,sequenceParameters,folderName, strExperimentalData,plottingCondition,modelReduction)
delta_t = 180; % time space
N_TimePoints =25;
t_max = N_TimePoints*delta_t; % maximum simulation time. The burnin time is automatically added in the SSA.
param.using_harringtonine = 0; param.using_stress = 1;
param.tim_inhibitor = 900; % burnin is added in the SSA function
pointOfNormalization =6; % number of points used to normalize.
% Pre-alocating memory
sim_CAP_ONLY = zeros(Nrep, N_TimePoints);
sim_IRES_ONLY = zeros(Nrep, N_TimePoints);
sim_CI = zeros(Nrep, 25);
error_CAP_ONLY = zeros(Nrep, N_TimePoints);
error_IRES_ONLY = zeros(Nrep, N_TimePoints);
error_CI = zeros(Nrep, N_TimePoints);
for k=1:Nrep
    if modelReduction ==1
        parfor i =1:N_sims
            [~,~,CAP_int_NA(i,:),IRES_int_NA(i,:)] =compilingTimes_Intensities (param,sequenceParameters,t_max,delta_t);
        end
    else
        parfor i =1:N_sims
            [RibsomePositions,~] = ssa_completeModel(param,sequenceParameters,t_max,delta_t);
            [CAP_int_NA(i,:),IRES_int_NA(i,:)] = SSAtoIntensity(RibsomePositions,sequenceParameters);
        end
    end
    [sim_CAP_ONLY(k,:), sim_IRES_ONLY(k,:),sim_CI(k,:), error_CAP_ONLY(k,:), error_IRES_ONLY(k,:),error_CI(k,:)] = intesities_selection_3int(CAP_int_NA, IRES_int_NA,pointOfNormalization);
end
[LL_NA, LL_NA_1G, LL_NA_2G, LL_NA_CI ] = compareNaAs_3int(param, folderName, strExperimentalData, sim_CAP_ONLY, sim_IRES_ONLY, sim_CI,error_CAP_ONLY, error_IRES_ONLY,error_CI, plottingCondition);
end