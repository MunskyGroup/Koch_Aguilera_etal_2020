function [S_arr,cap_times_IN,ires_times_IN] = ssa_fast_implementation(param,sequenceParameters,t_max)
if param.using_harringtonine==1  &&  param.using_stress==1
    error('using harringtonine and stress at the same time')
end
% Approximate time to reach steady state.
burnin = (max(sum(1./sequenceParameters.ke_cap*param.ke_avg),sum(1./sequenceParameters.ke_ires*param.ke_avg))*10);
% increase time to have multiple stant transitions.
if param.k24 ~= 0 % detecting when the system has only 3 states
    burnin = round(burnin+10*max(1./[param.k12,param.k21,param.k13,param.k31,param.k34,param.k43,param.k24,param.k42 ]));
else
    burnin = round(burnin+10*max(1./[param.k12,param.k21,param.k13,param.k31]));
end
burnin(burnin>1e4)=1e4;
% parameters affected by inhibitor
val_k12 =  param.k12 ;
val_k34 =   param.k34;
param.Elongation_time_CAP = sequenceParameters.codon_time_CAP/param.ke_avg;
param.Elongation_time_IRES = sequenceParameters.codon_time_IRES/param.ke_avg;
param.crossover = param.crossover/(param.crossover+mean(sequenceParameters.ke_cap*param.ke_avg)); % adjusting the cross-over rate to consider competition with termination reaction.

if param.using_harringtonine==1  ||  param.using_stress==1
    param.tim_inhibitor = normrnd(param.tim_inhibitor,10); % This line adds =/-10 seconds of varibility for the inhibitor application
    param.tim_inhibitor = burnin+ param.tim_inhibitor;
end
t_max = t_max+burnin ;

% if the model has 3 states just make k5 to k8 equal to zero.
t=0;
x=1;
S = [1 -1 2 -2 1 -1 2 -2];
t_arr = [1:1:t_max]';
S_arr = 0*t_arr;
j=1;
tstop = max(t_arr);
Nj = length(t_arr);
cap_times_IN = zeros(1,1e4);
ires_times_IN = zeros(1,1e4);
cap_times_OUT = zeros(1,1e4);
counter_kcap = 1;
counter_kires = 1;

W = zeros(1,11);

while t<tstop
    % states in the system
    % S1 = non translating
    % S2 = CAP translation
    % S3 = IRES only translation.
    % S4 = CAP and IRES translation
    % State transition
    W(1) = val_k12*(x==1); % K12
    W(2) = param.k21*(x==2); % K21
    W(3) = param.k13*(x==1); % K13
    W(4) = param.k31*(x==3); % K31
    W(5) = val_k34*(x==3); % K34
    W(6) = param.k43*(x==4); % K43
    W(7) = param.k24*(x==2); % K24
    W(8) = param.k42*(x==4); % K42
    % inititation and cross-over events
    % harringtonine
    if (param.using_harringtonine==1 && t>= param.tim_inhibitor)
        W(9) = 0; % cap init
        W(10) = 0; % ires init
        W(11) = param.crossover*(any(cap_times_OUT>t-1 & cap_times_OUT<t+1)); % Cross-over is allowd if a ribosome is finishing +/- a second of the termination time. (any(cap_times_OUT>=t-1 & cap_times_OUT<=t+1))
    elseif (param.using_stress==1 && t>= param.tim_inhibitor)
        if  param.InhibitorType ==1
            W(9) = param.ki_CAP*(x==2 || x==4); % cap init
            W(10) = param.ki_IRES*(x==3 || x==4); % ires init
            W(11) = param.crossover*(any(cap_times_OUT>t-1 & cap_times_OUT<t+1)); %(any(cap_times_OUT>t-1 & cap_times_OUT<t+1))   Cross-over is allowd if a ribosome is finishing +/- a second of the termination time. (any(cap_times_OUT>=t-1 & cap_times_OUT<=t+1))
            val_k12 =  param.k_inhi1* param.k12 ;
            val_k34 =  param.k_inhi1* param.k34;
        elseif  param.InhibitorType ==2
            W(9) = param.k_inhi1 * param.ki_CAP*(x==2 || x==4); % cap init
            W(10) = param.ki_IRES*(x==3 || x==4); % ires init
            W(11) = param.crossover*(any(cap_times_OUT>t-1 & cap_times_OUT<t+1)); % Cross-over is allowd if a ribosome is finishing +/- a second of the termination time. (any(cap_times_OUT>=t-1 & cap_times_OUT<=t+1))
        elseif  param.InhibitorType ==3
            W(9) = param.k_inhi1 * param.ki_CAP*(x==2 || x==4); % cap init
            W(10) = param.k_inhi2 * param.ki_IRES*(x==3 || x==4); % ires init
            W(11) = param.crossover*(any(cap_times_OUT>t-1 & cap_times_OUT<t+1)); % Cross-over is allowd if a ribosome is finishing +/- a second of the termination time. (any(cap_times_OUT>=t-1 & cap_times_OUT<=t+1))
        end
        
    elseif (param.using_harringtonine==0  &&  param.using_stress==0) || (param.using_harringtonine==1 && t<= param.tim_inhibitor) || (param.using_stress==1 && t<= param.tim_inhibitor)
        W(9) = param.ki_CAP*(x==2 || x==4); % cap init
        W(10) = param.ki_IRES*(x==3 || x==4); % ires init
        W(11) = param.crossover*(any(cap_times_OUT>t-1 & cap_times_OUT<t+1)); %(any(cap_times_OUT>t-1 & cap_times_OUT<t+1))   Cross-over is allowd if a ribosome is finishing +/- a second of the termination time. (any(cap_times_OUT>=t-1 & cap_times_OUT<=t+1))
        
    end
    % SSA
    a0 = sum(W);
    t = t+1/a0*log(1/rand);
    if param.reportPromoter ==1
        while j<=Nj&&t>t_arr(j)
            S_arr(j) = x;
            j=j+1;
        end
    end
    if t<=tstop
        r2=rand;
        i=1;
        w = W(1)/a0;
        while r2>w
            i=i+1;
            w=w+W(i)/a0;
        end
        if i<=8
            x = x+S(i);
        elseif i==9 % cap init
            cap_times_IN(counter_kcap) = t;
            cap_times_OUT(counter_kcap) = t+ param.Elongation_time_CAP;
            counter_kcap = counter_kcap+1;
        elseif i==10 % ires init
            ires_times_IN(counter_kires) = t;
            %ires_times_OUT(counter_kcap) = t + param.Elongation_IRES;
            counter_kires =counter_kires+1;
        elseif i==11 % cross-over
            ires_times_IN(counter_kires) = t;
            counter_kires =counter_kires+1;
        end
    end
    
end

% Removing zeros
cap_times_IN=cap_times_IN(cap_times_IN~=0);
ires_times_IN=ires_times_IN(ires_times_IN~=0);
% Removing initiation events before the burning time
cap_times_IN=cap_times_IN(cap_times_IN>burnin-param.Elongation_time_CAP );
ires_times_IN=ires_times_IN(ires_times_IN>burnin-param.Elongation_time_IRES);
% removing the burnin time
cap_times_IN = cap_times_IN - burnin;
ires_times_IN = ires_times_IN - burnin;
% Removing burning elements
S_arr = S_arr (burnin:end);
end
