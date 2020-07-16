function [CAP_int,IRES_int,t_array ] = initiationTimesToIntensity (param,sequenceParameters,t_max,delta_t,cap_times_IN,ires_times_IN)

Elongation_time_CAP = sequenceParameters.codon_time_CAP/param.ke_avg;
Elongation_time_IRES = sequenceParameters.codon_time_IRES/param.ke_avg;
% time of probe addition after starting translation
t_vs_position_cap = cumsum(1./sequenceParameters.ratio_ke_cap*param.ke_avg);
t_vs_intensity_cap = t_vs_position_cap(sequenceParameters.probePosition_CAP);
t_vs_position_ires = cumsum(1./sequenceParameters.ratio_ke_ires*param.ke_avg);
t_vs_intensity_ires = t_vs_position_ires(sequenceParameters.probePosition_IRES);
% deffine the time to terminate translation.
times_OUT_cap = cap_times_IN + Elongation_time_CAP;
times_OUT_ires = ires_times_IN + Elongation_time_IRES;
% deffining minimum and maximum times
tmin =1;
%% Generating output
% generating time array
t_array = tmin:delta_t:t_max;
Nt = length(t_array);
% Prealocating memory for CAP and IRES intensity
CAP_int = zeros(Nt,1);
IRES_int = zeros(Nt,1);
% Selecting the initiation events that are completed within the simulation time
G_cap = times_OUT_cap;
J_cap = (G_cap>=tmin);
cap_times_IN = cap_times_IN(J_cap);
G_ires = times_OUT_ires;
J_ires = (G_ires>=tmin);
ires_times_IN = ires_times_IN(J_ires);
% Time needed to translate the genes
dtmax_cap = Elongation_time_CAP;% sum(1./Ke_cap);
dtmax_ires = Elongation_time_IRES; % sum(1./Ke_ires);
% Convert event times to molecule number at all times.
for it = 1:Nt    %put this into coder of its own function
    t = t_array(it);
    dt_cap = t-cap_times_IN;  % time since initiation.
    CAP_int(it) = get_intens(dt_cap,dtmax_cap,t_vs_intensity_cap);
    dt_ires = t-ires_times_IN;  % time since initiation.
    IRES_int(it) = get_intens(dt_ires,dtmax_ires,t_vs_intensity_ires);
end
t_array = (t_array - tmin)';
% converting to UMP
CAP_int = CAP_int/length(sequenceParameters.probePosition_CAP);
IRES_int = IRES_int/length(sequenceParameters.probePosition_IRES);
end
