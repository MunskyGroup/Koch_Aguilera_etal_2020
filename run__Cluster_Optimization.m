function run__Cluster_Optimization(Model_type,randInit, optimizationMethod, PopulationSize, generations,maxSimulationTime)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This function is intended to perform a parameter optimization routine.
% INPUT:
%- Model_type. An integer between 1 and 4. Where the number indicates one of the following models:
    % 3S      %  Model 1. 3 States .
    % 3S_C    %  Model 2. 3 States + Cross-over.
    % 4S_DD   %  Model 3. 4 States 
    % 4S_DDC  %  Model 4. 4 States + Cross-over
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
%  -  randInit. An integer with values 0 or 1, where 1 uses a random initial condition. and 0 uses as an inital condition the best fit stored values as initial conditon. 
%  -  optimizationMethod. an Integer with values 1,2, or 3.  1 = Multi-onjective Optimization. 2 = GA. 3 = Pattern Search.
%  -  PopulationSize. If this integer is larger than 1. Indicates the number of individuals used for the genetic algorithm.
%  -  generations. If this integer is larger than 1. Indicates the number of generations used for the genetic algorithm.
%  -  maxSimulationTime. Time in hours indicating the maximum simulation time. It can be used only for method 3 (Pattern search).
% OUTPUT:
% The optimization generates the folder Fit_OP_MT_usedModel. Where _usedModel is a placeholder that represents the model number used. In the folder, the simulation results are saved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd Source_Code
%% This section creates the folders to save the fits.
if isnumeric(Model_type) ==0
    Model_type = str2double(Model_type);
    PopulationSize = str2double(PopulationSize);
    generations = str2double(generations);
    randInit = str2double (randInit);
    optimizationMethod = str2double(optimizationMethod); %  1 = Multi-onjective Optimization. 2 = GA. 3 = Pattern Search.
    maxSimulationTime =   str2double(maxSimulationTime);
end
modelReduction =1; % 1 = Uses the fast approximaiton. 0 = Uses the complete model.
Nrep =4;
N_sims = 4000;    % number of repetitions.
global NoEvaluations EvaluationCondition1 EvaluationCondition2 bestFitValue
NoEvaluations = 0;
conditionNoSkip =0;
EvaluationCondition1 = 0;
EvaluationCondition2 = 0;
bestFitValue = [Inf,Inf,Inf,Inf];

mainFolder = 'Best_Fits';
mainName = 'Fit_OP_MT_';
folderName =  [mainFolder,'/' ,mainName,num2str(Model_type)];
if ispc
    folderName =  [mainFolder,'\' ,mainName,num2str(Model_type)];
end
mkdir(folderName);
if optimizationMethod==3
    [initPop] = initialPopulation (Model_type);
    nParameters = size(initPop,2);
    if randInit ==1
        randSign =  2*randi([0 1],[1,nParameters])-1;
        variability =0.05;% percentage of variability in intial popultion
        initPop = (randSign*rand(nParameters).* initPop*variability) +  initPop;
    end
else
    [initPop] = initialPopulation_Matrix (Model_type);
end
initialMatrix =initPop;
initialMatrix(:,13:16) = 0;
strExperimentalData = extractingExperimentalData; % loading experimental data
sequenceParameters = fun_sequenceParameters; % loading sequence parameters

%% Parameter ranges
ki_cap_min = 1/70; ki_cap_max = 1/20;
ki_ires_min = 1/70; ki_ires_max = 1/20;
ke_min= 1.4; ke_max = 1.8;
k_12_min =1e-5; k_12_max =0.01;
k_21_min =1e-5; k_21_max =0.01;
k_13_min =1e-5; k_13_max =0.01;
k_31_min =1e-5; k_31_max =0.01;

if Model_type ==3 || Model_type ==4
    k_34_min =1e-5; k_34_max =0.01;
    k_43_min =1e-5; k_43_max =0.01;
    k_24_min =1e-5; k_24_max =0.01;
    k_42_min =1e-5; k_42_max =0.01;
    
elseif Model_type ==1 || Model_type ==2 || Model_type ==5 || Model_type ==6
    k_34_min =0; k_34_max =0;
    k_43_min =0; k_43_max =0;
    k_24_min =0; k_24_max =0;
    k_42_min =0; k_42_max =0;
    initialMatrix(:,8:11) = 0;
elseif Model_type == 7 || Model_type == 8
    k_24_min =1e-5; k_24_max =0.01;
    k_42_min =1e-5; k_42_max =0.01;
    k_34_min =0; k_34_max =0;
    k_43_min =0; k_43_max =0;
    
elseif Model_type == 9 || Model_type == 10
    k_34_min =1e-5; k_34_max =0.01;
    k_43_min =1e-5; k_43_max =0.01;
    k_24_min =0; k_24_max =0;
    k_42_min =0; k_42_max =0;
    
elseif Model_type == 11  || Model_type == 12 % k31  ~= k42;
    k_34_min =0; k_34_max =0;
    k_43_min =0; k_43_max =0;
    k_24_min =0; k_24_max =0;
    k_42_min =1e-5; k_42_max =0.01;
    
elseif Model_type == 13 || Model_type == 14  % k13  ~=  k24;
    k_34_min =0; k_34_max =0;
    k_43_min =0; k_43_max =0;
    k_24_min =1e-5; k_24_max =0.01;
    k_42_min =0; k_42_max =0;
end

if Model_type ==2 || Model_type ==4 || Model_type ==6 || Model_type ==8 || Model_type ==10 || Model_type ==12 || Model_type ==14
    k_crossover_min =1e-3; k_crossover_max=1;
    ki_CO_NaAs_min =0; ki_CO_NaAs_max=0;
else
    k_crossover_min =0; k_crossover_max=0;
    ki_CO_NaAs_min =0; ki_CO_NaAs_max =0;
    initPop(12)=0;
    initPop(16)=0;
end

k_gamma_1_min =0; k_gamma_1_max=0;
k_gamma_2_min =0; k_gamma_2_max=0;
ki_NaAs_min =0; ki_NaAs_max=0;

%% Parameter matrix
range =[ ki_cap_min ,ki_cap_max; ki_ires_min , ki_ires_max; ke_min , ke_max; k_12_min , k_12_max; k_21_min , k_21_max;k_13_min , k_13_max;k_31_min , k_31_max;k_34_min , k_34_max;k_43_min , k_43_max;k_24_min , k_24_max;k_42_min , k_42_max;k_crossover_min , k_crossover_max;k_gamma_1_min , k_gamma_1_max;k_gamma_2_min , k_gamma_2_max;ki_NaAs_min , ki_NaAs_max; ki_CO_NaAs_min , ki_CO_NaAs_max];
NoOptimizedParameters=size(range,1);

%% deffining parameter ranges
lb =[ki_cap_min ,ki_ires_min, ke_min, k_12_min, k_21_min, k_13_min, k_31_min, k_34_min, k_43_min, k_24_min, k_42_min, k_crossover_min, k_gamma_1_min, k_gamma_2_min, ki_NaAs_min, ki_CO_NaAs_min];
ub = [ki_cap_max ,ki_ires_max, ke_max, k_12_max, k_21_max, k_13_max, k_31_max, k_34_max, k_43_max, k_24_max, k_42_max, k_crossover_max, k_gamma_1_max, k_gamma_2_max, ki_NaAs_max, ki_CO_NaAs_max];
nvars= length(lb);

%% Running the optimization
fitFnc = @(x)fitnessFunction_multi(x,strExperimentalData, sequenceParameters,0,folderName,conditionNoSkip,optimizationMethod,Model_type,Nrep,N_sims,modelReduction);
if optimizationMethod ==1
    options = optimoptions('gamultiobj','MaxGenerations',generations,'PopulationSize',PopulationSize,'Display', 'iter','InitialPopulationMatrix',initialMatrix);
    [x,fval,~,~] = gamultiobj(fitFnc,nvars,[],[],[],[],lb,ub,options);
elseif optimizationMethod ==2
    optionsGA = optimoptions('ga','MaxGenerations',generations,'PopulationSize',PopulationSize,'ConstraintTolerance',1e-6,'InitialPopulationMatrix',initialMatrix,'Display', 'iter');
    [x,fval] = ga(fitFnc,nvars,[],[],[],[],lb,ub,[],[],optionsGA);
elseif optimizationMethod ==3
    options_PS = optimoptions('patternsearch','Display','iter','MaxIterations',1e6,'MeshTolerance',1e-8,'MaxTime',3600*maxSimulationTime);
    [x,fval] = patternsearch(fitFnc,initialMatrix,[],[],[],[],lb,ub,options_PS);
end

%% Plotting Results
if optimizationMethod ==1
    mean_fit = nanmean(fval);
    norm_fval = fval./mean_fit;
    sum_fit = sum(norm_fval,2);
    [~,Ind]=min(sum_fit);
else
    Ind=1;
end
x = x(Ind,:);
% re-defining parameters for models 5 and 6
if Model_type ==5 || Model_type ==6
    x(:,8) = x(:,4);  % k_34 = k_12;
    x(:,9) = x(:,5);  % k_43 = k_21;
    x(:,10) = x(:,6); % k_24 = k_13;
    x(:,11) = x(:,7);  % k_42 = k_31;
end

if  Model_type == 7 || Model_type == 8
    x(:,8) = x(:,4);  % k_34 = k_12;
    x(:,9) = x(:,5);  % k_43 = k_21;
end

if  Model_type == 9 || Model_type == 10
    x(:,10) = x(:,6); % k_24 = k_13;
    x(:,11) = x(:,7);  % k_42 = k_31;
end

if Model_type == 11  || Model_type == 12  % k31  ~= k42;
    x(:,8) = x(:,4);  % k_34 = k_12;
    x(:,9) = x(:,5);  % k_43 = k_21;
    x(:,10) = x(:,6); % k_24 = k_13;
end

if Model_type == 13 || Model_type == 14  % k13  ~=  k24;
    x(:,8) = x(:,4);  % k_34 = k_12;
    x(:,9) = x(:,5);  % k_43 = k_21;
    x(:,11) = x(:,7);  % k_42 = k_31;
end

% Plotting results for the best fit.
fitValue = fitnessFunction_multi (x(Ind,:),strExperimentalData, sequenceParameters,1,folderName,conditionNoSkip,optimizationMethod,Model_type,Nrep,N_sims,modelReduction);
%% Saving and moving the parameters
dlmwrite('zplotted_Parameters.txt',x(1,:),'delimiter','\t','precision',3)
dlmwrite('zplotted_Fit.txt',fval(1,:),'delimiter','\t','precision',3)
movefile ('zplotted_Parameters.txt', folderName)
movefile ('zplotted_Fit.txt', folderName)
save op_val.mat x fval
movefile ('op_val.mat', folderName)

cd ..
end
