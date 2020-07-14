function [RibsomePositions,S_arr] = ssa_completeModel(param,sequenceParameters,t_max,delta_t)

param.crossover = param.crossover/(param.crossover+mean(sequenceParameters.ratio_ke_cap*param.ke_avg)); % adjusting the cross-over rate to consider competition with termination reaction. This is just to make simple and complete models equally deffined.
if param.using_harringtonine==1  &&  param.using_stress==1
    error('using harringtonine and stress at the same time')
end
% Approximate time to reach steady state
burnin =1e4;
if param.using_harringtonine==1  ||  param.using_stress==1
    param.tim_inhibitor = normrnd(param.tim_inhibitor,10); % This line adds =/-10 seconds of varibility for the inhibitor application
    param.tim_inhibitor = burnin + param.tim_inhibitor;
end
t_final = t_max + burnin ;
t_array = [1:delta_t:t_final];

%% Gene states.
s1 =1; % non translating mRNA
s2 =0; % cap translating mRNA
s3 =0; % ires translating mRNA
s4 =0; % cap and ires translating mRNA

%% Deffining parameters
exclusion = 9;
maxNoRibosomes =1e3;
maximum_Number_Ribosomes = maxNoRibosomes; % maximum number of ribosomes that can be loaded on the mRNA
ki_CAP = param.ki_CAP;  % Initiation rate in CAP
ki_IRES = param.ki_IRES;  % Initiation rate in IRES ************
k_cap_ires =param.crossover; % cross-over parameter
k_cap_ires_o =param.crossover; % cross-over parameter
ki_CAP_o = param.ki_CAP;  % Initiation rate in CAP
ki_IRES_o = param.ki_IRES;  % Initiation rate in IRES ************
ki_CO_NaAs_o = param.ki_CO_NaAs;
k_term = param.ke_avg;
k_elongation_1G = sequenceParameters.ratio_ke_cap.*param.ke_avg; % elongation constants for the 0 frame
k_elongation_2G = sequenceParameters.ratio_ke_ires.*param.ke_avg;
geneLength_1G = length(k_elongation_1G);  % number of codons first gene.
geneLength_2G = length(k_elongation_2G);  % number of codons second gene.
gene_total= geneLength_1G+geneLength_2G;
k_elongation_total = [k_elongation_1G, k_elongation_2G];
X_State = zeros(maximum_Number_Ribosomes,1);                       % Initial state (length of the longest protein = 1)
t = t_array(1); % time array
%% Prealocating memory
number_TimePoints = round(t_max/delta_t);
X_output = zeros(maximum_Number_Ribosomes,number_TimePoints);
st = zeros(number_TimePoints,4);

%% Run SSA
iteration = 1;
while t < t_final

    if  param.using_stress==1 && t>= param.tim_inhibitor
        ki_CAP = ki_CAP_o * exp(-param.gamma_1*(max(0,t-param.tim_inhibitor)));          % Block of CAP initiation
        ki_IRES = ki_IRES_o + param.ki_NaAs*(1-exp(-param.gamma_2*(max(0,t-param.tim_inhibitor)))); % Increase of IRES initiation
        k_cap_ires = k_cap_ires_o + ki_CO_NaAs_o*(1-exp(-param.gamma_2*(max(0,t-param.tim_inhibitor))));
    end

    %% Compute propensity functions
    % updating state
    X_State = sort(X_State,'descend');
    Nribosomes = sum(X_State > 0); % counts the current number of active ribosomes in the mRNA
    % Deffining the stoichiometry matrix for the elongation reactions
    temp_Stoichiometry = eye(Nribosomes);
    % Deffining the propensity function
    temp_Propensities = zeros(Nribosomes + 12,1);% Nribosomes = elongations/termination. 4 additional rows for: 2 initiation, 1 circularization, 1 re-use

    if Nribosomes>=1 && exclusion>1
        selEl=X_State(X_State > 0);
        temp_Propensities(X_State > 0) = k_elongation_total(selEl); % assigning values to the elongations
    end
    % initiation in CAP
    if s2 >0 || s4 ~=0
        if Nribosomes == 0 || X_State(Nribosomes) > exclusion   %
            temp_Propensities(Nribosomes+1) = ki_CAP ;
        else
            temp_Propensities(Nribosomes+1) = 0 ;
        end
    else
        temp_Propensities(Nribosomes+1) =0;
    end

    % initiation in IRES
    if s3 >0 || s4 ~=0
        if Nribosomes == 0 || nnz(X_State(1:Nribosomes) >= geneLength_1G+1 & X_State(1:Nribosomes) <= geneLength_1G + exclusion)==0
            temp_Propensities(Nribosomes+2) = ki_IRES;
        else
            temp_Propensities(Nribosomes+2) = 0;
        end
    else
        temp_Propensities(Nribosomes+2) = 0 ;
    end

    % intiation under Inhibitory conditions. Harringtonine
    if t >= param.tim_inhibitor  && param.using_harringtonine==1
        temp_Propensities(Nribosomes+1) =0;
        temp_Propensities(Nribosomes+2) =0;
    end

    % elongation
    if Nribosomes > 1
        elongation_Condition = ones (Nribosomes,1);
        if (X_State(1)> 1 && X_State(1)< gene_total)
            elongation_Condition(1) = 1;
        else
            elongation_Condition(1) =0;
        end
        elongation_Condition(2:Nribosomes) = ((X_State(2:Nribosomes) + exclusion) < X_State(1:Nribosomes-1));
        temp_Propensities(1:Nribosomes) = temp_Propensities(1:Nribosomes) .* elongation_Condition;
    elseif  Nribosomes == 1 && (X_State(1) <=10 && X_State(1)<geneLength_1G)
        temp_Propensities(1) = k_elongation_1G(10);
    elseif Nribosomes == 1 && (X_State(1) >1 && X_State(1) <gene_total)
        temp_Propensities(1) = k_elongation_1G(10);
    end

    % Propensity to go from CAP to IRES
    c_i_condition = (X_State(1:Nribosomes)) ==  geneLength_1G;  % X(i)=1 if that ribosome is at the end of gene 1.  Otherwise zero.
    temp_Propensities(1:Nribosomes) = temp_Propensities(1:Nribosomes).* ~c_i_condition + (c_i_condition .* k_cap_ires);


    % Propensity to terminate CAP
    if sum(X_State==geneLength_1G) >= 1
        Index_CAP_termination = find(X_State==geneLength_1G,1);
        temp_Propensities(Nribosomes+3) = k_term;
    else
        Index_CAP_termination = maximum_Number_Ribosomes-10;
    end

    % Propensity to terminate IRES
    if sum(X_State==gene_total) >= 1
        Index_IRES_termination = find(X_State==gene_total,1);
        temp_Propensities(Nribosomes+4) = k_term;
    else
        Index_IRES_termination = maximum_Number_Ribosomes-10;
    end

    %% Propensities to change states.
    temp_Propensities(Nribosomes+5) = param.k12 * s1; % ... propensity to make s1 = 0 and s2=1;
    temp_Propensities(Nribosomes+6) = param.k21 * s2; % ... propensity to make s1 = 1 and s2 = 0;

    temp_Propensities(Nribosomes+7) = param.k13 * s1; % ... propensity to make s1 = 0 and s3 =1;
    temp_Propensities(Nribosomes+8) = param.k31 * s3; % ... propensity to make s3 = 0 and s1=1;

    temp_Propensities(Nribosomes+9) =  param.k34 * s3; % ... propensity to make s4 = 1 and s3=0;
    temp_Propensities(Nribosomes+10) = param.k43 * s4; % ... propensity to make s3 = 1 and s4 =0;

    temp_Propensities(Nribosomes+11) = param.k24 * s2; % ... propensity to make s2 = 0 and s4=1;
    temp_Propensities(Nribosomes+12) = param.k42 * s4; % ... propensity to make s4 = 0 and s2 =1;

    %% Updating sum of propensities
    sum_Propensities = sum(temp_Propensities);
    %% Update time
    t = t - log(rand) / sum_Propensities;
    %% Generate output
    while (t >= burnin) && (iteration <= number_TimePoints) && t > t_array(iteration)+burnin
        X_output(1:size(X_State,1),iteration) = X_State;
        st (iteration,:) = [s1,s2,s3,s4];
        iteration = iteration + 1;
    end

    %% Update state
    if t < t_final
        %% Select reaction
        selectedReaction = sum_Propensities * rand;
        i = 1;
        tmp = temp_Propensities(i);
        while tmp < selectedReaction  % iterating all propensities and update the system.
            i = i + 1;
            tmp = tmp + temp_Propensities(i);
        end
        %% update states
        if i <= length(temp_Propensities) - 12              % UPDATING STATE BY ELONGATION REACTIONS
            X_State(1:Nribosomes,1) = X_State(1:Nribosomes,1) + temp_Stoichiometry(:,i);
        elseif i == length(temp_Propensities) - 11  && (s2 ~=0 || s4 ~=0)         % CAP INITIATION
            X_State (Nribosomes+1,1) = 1;
        elseif i == length(temp_Propensities) - 10  && (s3 ~=0 || s4 ~=0)       % IRES INITIATION
            X_State (Nribosomes+1,1) = geneLength_1G+1;
        elseif i == length(temp_Propensities)  -9            % CAP TERMINATION
            X_State (Index_CAP_termination) = 0;
        elseif i == length(temp_Propensities) -8             % IRES TERMINATION
            X_State (Index_IRES_termination) = 0;
        elseif i == length(temp_Propensities)-7              % propensity to make s1 = 0 and s2=1;
            s1 = s1-1;
            s2 = s2+1;
        elseif i == length(temp_Propensities)-6              % propensity to make s1 = 1 and s2 = 0;
            s1 = s1+1;
            s2 = s2-1;
        elseif i == length(temp_Propensities)-5              % propensity to make s1 = 0 and s3 =1;
            s1 = s1-1;
            s3 = s3+1;
        elseif i == length(temp_Propensities)-4              % propensity to make s3 = 0 and s1=1;
            s1 = s1+1;
            s3 = s3-1;
        elseif i == length(temp_Propensities)-3              % propensity to make s4 = 1 and s3=0;
            s4 = s4+1;
            s3 = s3-1;
        elseif i == length(temp_Propensities)-2              % propensity to make s3 = 1 and s4 =0;
            s4 = s4-1;
            s3 = s3+1;
        elseif i == length(temp_Propensities)-1              % propensity to make s2 = 0 and s4=1;
            s4 = s4+1;
            s2 = s2-1;
        elseif i == length(temp_Propensities)                % propensity to make s4 = 0 and s2 =1;
            s4 = s4-1;
            s2 = s2+1;
        end

    end
    % Returning the system to Starting position in case of RNA
    % circularization.
    X_State(X_State>gene_total) = 0;
end
RibsomePositions = X_output(:,1:end);
S_arr = st(1:end,:);

end
