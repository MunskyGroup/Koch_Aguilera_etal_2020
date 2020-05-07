function [t_array,S_arr,CAP_int,IRES_int] =compilingTimes_Intensities (param,sequenceParameters,t_max,delta_t)
%% Running the SSA
[S_arr,cap_times_IN,ires_times_IN] = ssa_fast_implementation(param,sequenceParameters,t_max);
%% Converting initiation times to Intensities.
[CAP_int,IRES_int,t_array ] = initiationTimesToIntensity (param,sequenceParameters,t_max,delta_t,cap_times_IN,ires_times_IN);
S_arr = S_arr(1:delta_t:end-1);
end
