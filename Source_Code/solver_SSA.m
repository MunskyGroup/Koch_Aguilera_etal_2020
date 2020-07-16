function[t_array,S_arr,CAP_int,IRES_int] =solver_SSA (param,sequenceParameters,t_max,delta_t,N_sims)

%% Running the SSA
t_array = 1:delta_t:t_max; N_steps = length(t_array);
CAP_int = zeros (N_sims,N_steps);IRES_int = zeros (N_sims,N_steps);
S_arr = zeros (N_sims,N_steps);
parfor i =1:N_sims
    [RibsomePositions,S_arr(i,:)] = ssa_completeModel(param,sequenceParameters,t_max,delta_t);
    [CAP_int(i,:),IRES_int(i,:)] = SSAtoIntensity(RibsomePositions,sequenceParameters);
end

end