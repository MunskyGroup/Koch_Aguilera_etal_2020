function [CAP_int,IRES_int] = fun_initialSimulations(delta_t,t_max,N_sims,param,sequenceParameters,modelReduction)
param.using_harringtonine = 0; param.using_stress = 0; param.tim_inhibitor = 0;
%% Running the SSA
t_array = 1:delta_t:t_max; N_steps = length(t_array);
CAP_int = zeros (N_sims,N_steps);IRES_int = zeros (N_sims,N_steps);
if modelReduction ==1
    parfor i =1:N_sims
        [~,~,CAP_int(i,:),IRES_int(i,:)] =compilingTimes_Intensities (param,sequenceParameters,t_max,delta_t);
    end
else
    parfor i =1:N_sims
        [RibsomePositions,~] = ssa_completeModel(param,sequenceParameters,t_max,delta_t);
        [CAP_int(i,:),IRES_int(i,:)] = SSAtoIntensity(RibsomePositions,sequenceParameters);
    end
end
end