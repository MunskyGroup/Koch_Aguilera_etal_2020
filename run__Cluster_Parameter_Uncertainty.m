%% This function is intended to perform a parameter uncertainty after an optimization search.
function run__Cluster_Parameter_Uncertainty(Model_type,numberTotalEvaluations,Nrep,N_sims,justMakePlot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This function is intended to perform a parameter uncertainty routine for the selected model.
% INPUT:
%- Model_type. An integer between 1 and 4. Where the number indicates one of the following models:
% 3S      %  Model 1. 3 States .
% 3S_C    %  Model 2. 3 States + Cross-over.
% 4S_DD   %  Model 3. 4 States .
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
%  - Nrep. Integer larger than 1. Indicates the number of independent repetitions for Harringtonine.
%  - N_sims. If Integer is larger than 1. Indicates the number of simulated mRNA spots.
%  - justMakePlot. Integer with values 0 or 1, where 1 indicates that the code only generates plots with the save data. and 0 performs simulations and make plots.
% OUTPUT:
% The optimization generates the folder ParameterUncertainty_usedModel. Where _usedModel is a placeholder that represents the model number used. In the folder, the simulation results are saved.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd Source_Code
if isnumeric(Model_type) ==0
    Model_type = str2double(Model_type);
    numberTotalEvaluations = str2double(numberTotalEvaluations);
    Nrep = str2double(Nrep);
    N_sims = str2double(N_sims);
    justMakePlot= str2double(justMakePlot);
end
modelReduction =1; % 1 = Uses the fast approximaiton. 0 = Uses the complete model.

folderName = horzcat('ParameterUncertainty_',num2str(Model_type)); if exist (folderName, 'dir') ~= 7; mkdir(folderName); end

%% Provided model parameters and experimental data.
paramNamesComplete = {'ki_{cap}','ki_{ires}','k_{e}','k_{12}','k_{21}','k_{13}','k_{31}','k_{34}','k_{43}','k_{24}','k_{42}','k_{CO}'};

% This section selects the specific parameters that are used for each model.
if Model_type ==1
    selPar=[1:7];
    paramNames =paramNamesComplete(selPar);
    [initPop] = initialPopulation (Model_type);
elseif Model_type ==2
    selPar=[1:7,12];
    paramNames =paramNamesComplete(selPar);
    [initPop] = initialPopulation (Model_type);
elseif Model_type ==3
    selPar=[1:11];
    paramNames =paramNamesComplete(selPar);
    [initPop] = initialPopulation (Model_type);
elseif Model_type ==4
    selPar=[1:12];
    paramNames =paramNamesComplete(selPar);
    [initPop] = initialPopulation (Model_type);
elseif Model_type ==5
    selPar=[1:7];
    paramNames =paramNamesComplete(selPar);
    [initPop] = initialPopulation (Model_type);
elseif Model_type ==6
    selPar=[1:7,12];
    paramNames =paramNamesComplete(selPar);
    [initPop] = initialPopulation (Model_type);
elseif Model_type ==7
    selPar=[1:7,10:11];
    paramNames =paramNamesComplete(selPar);
    [initPop] = initialPopulation (Model_type);
elseif Model_type ==8
    selPar=[1:7,10:12];
    paramNames =paramNamesComplete(selPar);
    [initPop] = initialPopulation (Model_type);
elseif Model_type ==9
    selPar=[1:9];
    paramNames =paramNamesComplete(selPar);
    [initPop] = initialPopulation (Model_type);
elseif Model_type ==10
    selPar=[1:9,12];
    paramNames =paramNamesComplete(selPar);
    [initPop] = initialPopulation (Model_type);
elseif Model_type ==11
    selPar=[1:7,11];
    paramNames =paramNamesComplete(selPar);
    [initPop] = initialPopulation (Model_type);
elseif Model_type ==12
    selPar=[1:7,11,12];
    paramNames =paramNamesComplete(selPar);
    [initPop] = initialPopulation (Model_type);
elseif Model_type ==13
    selPar = [1:11]; % displaying all parameter values
    paramNames =paramNamesComplete(selPar);
    [initPop] = initialPopulation (Model_type);
elseif Model_type ==14
    selPar=[1:7,10,12];
    paramNames =paramNamesComplete(selPar);
    [initPop] = initialPopulation (Model_type);
end
strExperimentalData = extractingExperimentalData; % loading experimental data
sequenceParameters = fun_sequenceParameters; % loading sequence parameters

% Provide parameters for the Parameter Scan.
numberOfObjectiveFunctions = 3;                       %  Percentage, Intensities, HT.
ftvals_best = ones(1,numberOfObjectiveFunctions).*1e4; % Prealocating matrix with the values for the best objective function
Ftvals = ones(1,numberOfObjectiveFunctions);          % Prealocating matrix with the values for the objective function.
pars = initPop;
pars= pars(1:12);
numberOfFreeParameters = length(pars);
pars_real= pars;
pars_sv = [];
inorout = [];
burninParameters = 50; %number of initial parameters to not consider.
fig_format = 'tif';
namePlot = ['Param_UC_','_Model_',num2str(Model_type),'.',fig_format];

%% BEM
evaluationType = 1;
[fit_Percent,fit_Distribution,~] = fun_parameterUncertainty(pars,strExperimentalData, sequenceParameters,Nrep,N_sims,evaluationType,modelReduction);
evaluationType = 2;
[~,~,fit_HT] = fun_parameterUncertainty(pars,strExperimentalData, sequenceParameters,Nrep,N_sims,evaluationType,modelReduction);
ftvals_best = [fit_Percent,fit_Distribution,fit_HT];
best_fit_obj = sum(ftvals_best);

Best_L3 = ftvals_best(3);  % Best ever Likelihood for HT data.
log_likelihood_dif_ToRejectValues =30;
Threshold = best_fit_obj + log_likelihood_dif_ToRejectValues;

%%
if justMakePlot ==0
    
    %% For loop that evaluates each parameter and salves the parameter if the objective function is within a 10% error.
    for ip=1:numberTotalEvaluations
        if ip==1
            delta = 0;
        else
            delta = 0.01; % Percentage of random changes to parameters
        end
        %% Deffing parameters.
        %%%%%%%%%%%%%%%%%%
        pars_new = pars.*(1+delta*randn(size(pars)));
        % this section makes zero the parameters for the cases with 3 states
        % and not-cross-over
        if Model_type ==1
            pars_new(8:12) =0;
        elseif Model_type ==2
            pars_new(8:11) =0;
        elseif Model_type ==3
            pars_new(12) =0;
        % This section make the specific association between parameters for each model.    
        elseif Model_type ==5
            pars_new(8) = pars_new(4);  % k_34 = k_12;
            pars_new(9) = pars_new(5);  % k_43 = k_21;
            pars_new(10) = pars_new(6); % k_24 = k_13;
            pars_new(11) = pars_new(7);  % k_42 = k_31;
            pars_new(12) =0;
        elseif Model_type ==6
            pars_new(8) = pars_new(4);  % k_34 = k_12;
            pars_new(9) = pars_new(5);  % k_43 = k_21;
            pars_new(10) = pars_new(6); % k_24 = k_13;
            pars_new(11) = pars_new(7);  % k_42 = k_31;
        elseif Model_type ==7
            pars_new(8) = pars_new(4);  % k_34 = k_12;
            pars_new(9) = pars_new(5);  % k_43 = k_21;
            pars_new(12) =0;
        elseif Model_type ==8
            pars_new(8) = pars_new(4);  % k_34 = k_12;
            pars_new(9) = pars_new(5);  % k_43 = k_21;
        elseif Model_type ==9
            pars_new(10) = pars_new(6); % k_24 = k_13;
            pars_new(11) = pars_new(7);  % k_42 = k_31;
            pars_new(12) =0;
        elseif Model_type ==10
            pars_new(10) = pars_new(6); % k_24 = k_13;
            pars_new(11) = pars_new(7);  % k_42 = k_31;
        elseif Model_type ==11
            pars_new(8) = pars_new(4);  % k_34 = k_12;
            pars_new(9) = pars_new(5);  % k_43 = k_21;
            pars_new(10) = pars_new(6); % k_24 = k_13;
            pars_new(12) =0;
        elseif Model_type ==12
            pars_new(8) = pars_new(4);  % k_34 = k_12;
            pars_new(9) = pars_new(5);  % k_43 = k_21;
            pars_new(10) = pars_new(6); % k_24 = k_13;
        elseif Model_type ==13
            pars_new(8) = pars_new(4);  % k_34 = k_12;
            pars_new(9) = pars_new(5);  % k_43 = k_21;
            pars_new(11) = pars_new(7);  % k_42 = k_31;
            pars_new(12) =0;
        elseif Model_type ==14
            pars_new(8) = pars_new(4);  % k_34 = k_12;
            pars_new(9) = pars_new(5);  % k_43 = k_21;
            pars_new(11) = pars_new(7);  % k_42 = k_31;
        end
        
        pars_new(pars_new<0)=0;
        pars_sv = [pars_sv;pars_new];
        evaluationType =1;
        [fit_Percent,fit_Distribution,~] = fun_parameterUncertainty(pars_new,strExperimentalData, sequenceParameters,Nrep,N_sims,evaluationType,modelReduction);
        
        %% Performing the Harringtonine Assays
        if ip==1 || (fit_Percent+fit_Distribution+Best_L3*0.95)<Threshold    %
            
            % Harringtonine assays
            evaluationType =2;
            [~,~,fit_HT] = fun_parameterUncertainty(pars_new,strExperimentalData, sequenceParameters,Nrep,N_sims,evaluationType,modelReduction);
            Ftvals = [Ftvals;fit_Percent, fit_Distribution, fit_HT];
            
            %% BEM
            % Updates best HT fit if it finds something better.
            if fit_HT<Best_L3
                Best_L3=fit_HT;
            end
            % Updates best overall fit and threshold if it finds something better.
            if (fit_Percent+fit_Distribution+fit_HT)<best_fit_obj
                best_fit_obj = (fit_Percent+fit_Distribution+Best_L3);
            end
            
            if  (fit_Percent+fit_Distribution+fit_HT)<Threshold % && ip>burninParameters
                pars = pars_new;
                inorout = [inorout,1];
            else
                inorout = [inorout,0];
            end
            
        else
            Ftvals = [Ftvals;fit_Percent, fit_Distribution, Inf];
            inorout = [inorout,0];
        end
    end
    %% reject
    inorout(1:burninParameters)= 0;
    %% matrixes that are saved during the analysis
    Parameter_AboveThreshold = pars_sv(inorout==1,:);
    total_tested_Parameters = pars_sv;                    % total tested parameters
    total_Objective_function_values = Ftvals(2:end,:);      % objective function values for all tested parameters
    total_aceptance_rejection_parameter_binary = inorout; % zeros is rejected one is accepted
    %% Saving the selected parameters to file.
    cd (folderName)
    S = dir('randomSearchData*.*');
    try
        numberOf_RS_Repetitions =  length (S);
    catch
        numberOf_RS_Repetitions = 0;
    end
    fname=[pwd,'/randomSearchData_',num2str(numberOf_RS_Repetitions+1),'.mat'];
    if ispc
        fname=[pwd,'\randomSearchData_',num2str(numberOf_RS_Repetitions+1),'.mat'];
    end
    save(fname,'Parameter_AboveThreshold','initPop','total_tested_Parameters','total_Objective_function_values','total_aceptance_rejection_parameter_binary','numberTotalEvaluations')
    cd ..
    cd (folderName)
    %% Load previous list of parameters that fullfil the selection criterion
    sel_param_list = zeros(1,numberOfFreeParameters);%[0,0,0,0,0];
    % load all repetitions to Plot data.
    S_load = dir('randomSearchData*.*');
    numberOf_RS_Repetitions =  length (S_load);
    for i = 1: numberOf_RS_Repetitions
        fileNames =  ['randomSearchData_',num2str(i),'.mat'];
        load (fileNames)
        pre_sel_param_list = Parameter_AboveThreshold;
        sel_param_list = [sel_param_list;pre_sel_param_list];
    end
    % removing rows with only zeros.
    sel_param_list = sel_param_list(any(sel_param_list,2),:);
    cd ..
    %% Plotting the covariance matrix
    [meanParameters_UC,std_Parameters_UC] = make_uncertainty_plot(sel_param_list(:,selPar),'nameplot',namePlot,'contour_plot',1,'ellipse',0,'parameter_names',paramNames,'kde',1,'true_parameters',pars_real(selPar),'figformat',fig_format);
    movefile(horzcat(namePlot),horzcat(folderName),'f');
    %% Saving the selected parameters to file.
    dlmwrite('mean_Parameters_UC.txt',meanParameters_UC','delimiter','\t')% ,'precision',3
    movefile ('mean_Parameters_UC.txt', folderName)
    dlmwrite('std_Parameters_UC.txt',std_Parameters_UC','delimiter','\t')% ,'precision',3
    movefile ('std_Parameters_UC.txt', folderName)
else
    cd (folderName)
    %% Load previous list of parameters that fullfil the selection criterion
    sel_param_list = zeros(1,numberOfFreeParameters);%[0,0,0,0,0];
    % load all repetitions to Plot data.
    S_load = dir('randomSearchData*.*');
    numberOf_RS_Repetitions =  length (S_load);
    for i = 1: numberOf_RS_Repetitions
        fileNames =  ['randomSearchData_',num2str(i),'.mat'];
        load (fileNames)
        pre_sel_param_list = Parameter_AboveThreshold;
        sel_param_list = [sel_param_list;pre_sel_param_list];
    end
    % removing rows with only zeros.
    sel_param_list = sel_param_list(any(sel_param_list,2),:);
    size(sel_param_list)
    cd ..
    %% Plotting the covariance matrix
    coef_int =80;
    [meanParameters_UC,std_Parameters_UC,CI_bestParameters_UC] = make_uncertainty_plot(sel_param_list(:,selPar),'nameplot',namePlot,'contour_plot',1,'ellipse',0,'parameter_names',paramNames,'kde',0,'true_parameters',pars_real(selPar),'figformat',fig_format,'coef_interval',coef_int);
    movefile(horzcat(namePlot),horzcat(folderName),'f');
    
    %% Saving the selected parameters to file.
    dlmwrite('mean_Parameters_UC.txt',meanParameters_UC','delimiter','\t')% ,'precision',3
    movefile ('mean_Parameters_UC.txt', folderName)
    
    dlmwrite('std_Parameters_UC.txt',std_Parameters_UC,'delimiter','\t')% ,'precision',3
    movefile ('std_Parameters_UC.txt', folderName)
    
    dlmwrite('CI_Parameters_UC.txt',CI_bestParameters_UC,'delimiter','\t')% ,'precision',3
    movefile ('CI_Parameters_UC.txt', folderName)

end
cd ..
end