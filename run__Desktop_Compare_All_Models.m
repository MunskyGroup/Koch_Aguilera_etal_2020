clear all; close all; clc;
justPlot = 1;             % 1 = plots the data.     REPeat  12 and 8
runAllModels =0;          % 1 = runs all simulations.
run_SpecificModel =[]; % only use a value if it is needed to run an specific model. If it is needed to run for all models leave empty.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This function is intended to perform a parameter uncertainty routine.
% INPUT:
%  - justPlot. An integer with values 0 or 1, where 1 indicates that the code only generates plots with the save data, and 0 performs simulations and make plots.
%  - runAllModels. Integer with values 0 or 1, where 1 indicates that the code evaluates all models. 0 indicates that the code only evaluates a single model.
%- run_SpecificModel. An integer between 1 and 4. Where the number indicates one of the following models:
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
% OUTPUT:
% The comparision for all models is saved in the folder Model_Comparison.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd Source_Code
mainName = 'Fit_OP_MT_';
folderName = 'Model_Comparison';
N_Models =14; % total number of evaluated models.
modelReduction =1;        % 1 = Uses the fast approximaiton. 0 = Uses the complete model.
N_runs_forStatistics =3;  % number of repetitions used to calculate SD.
Nrep = 10;                 % number of repetitions for HT, NaAS and DTT.
N_sims = 4000;            % number of independent SSAs.
if isempty(run_SpecificModel) ==0 && runAllModels ==1
    error('please select one model or all models')
end
labelOrder = [1,2,5,6,7,8,9,10,11,12,13,14,3,4];
labels_Model={'3S','3S_{C}','4S_{II}','4S_{IIC}','4S_{DI}','4S_{DIC}','4S_{ID}','4S_{IDC}', '4S_{Im1}','4S_{Im1C}','4S_{Im2}','4S_{Im2C}', '4S_{DD}','4S_{DDC}'};
labels_Model_notSorted={'3S','3S_{C}', '4S_{DD}','4S_{DDC}','4S_{II}','4S_{IIC}','4S_{DI}','4S_{DIC}','4S_{ID}','4S_{IDC}', '4S_{Im1}','4S_{Im1C}','4S_{Im2}','4S_{Im2C}'};
n_freeParameters = [7    8    11    12    7    8    9    10    9    10    8    9    8    9];
strExperimentalData = extractingExperimentalData; % loading experimental data
sequenceParameters = fun_sequenceParameters; % loading sequence parameters

selectedModel = 0;
if justPlot ==0
    if isempty(run_SpecificModel) == 1 && runAllModels == 1
        kk=1;
    elseif runAllModels ==0
        kk = N_Models;
    end
    while kk <= N_Models
        if runAllModels == 1
            selectedModel =0;
        end
        if isempty(run_SpecificModel) ==1
            j =0;
            while selectedModel ==0
                j =j+1;
                Files=dir(['Model',num2str(j),'_*']);
                if isempty(Files) ==1
                    selectedModel = j;
                end
            end
        else
            selectedModel = run_SpecificModel;
        end
        if  selectedModel ==1 &&  isempty(run_SpecificModel) ==1;  try; rmdir(folderName, 's'); catch;end; mkdir(folderName); end
        [param] = initialPopulation (selectedModel);
        
        if size(param,1)>1;    sum_fit = sum(fval(:,1:end-1),2);    [~,Ind]=sort(sum_fit); else;    Ind =1;end
        for k=1:N_runs_forStatistics
            % This section runs the SSA for the selected model
            %             subfolderName = ['Model',num2str(selectedModel),'_', num2str(k)]; mkdir(subfolderName);
            pre_Report_Matrix  = fun_plot_complete_sweep (param(Ind(1),:),strExperimentalData, sequenceParameters,0,1,Nrep,N_sims,modelReduction);
            name_model_rep = ['comp_',num2str(selectedModel),'_', num2str(k),'.mat'];
            save(name_model_rep, 'pre_Report_Matrix')
            movefile( name_model_rep,folderName,'f');
        end
        kk = kk+1;
    end
    
    %% test if all folder exist on the main folder. then move and generate complete matrix
    if selectedModel == N_Models
        cd (folderName)
        for k=1:N_runs_forStatistics
            for i =1:N_Models
                name_model_rep = ['comp_',num2str(i),'_', num2str(k),'.mat'];
                load (name_model_rep)
                Report_Matrix(:,i)  = pre_Report_Matrix;
            end
            vector_Report_Matrix (:,:,k) = Report_Matrix;
            vector_Report_Matrix_norm (:,:,k) = (Report_Matrix'-min(Report_Matrix)')';
        end
        save compareModels.mat vector_Report_Matrix vector_Report_Matrix_norm
        cd ..
    end
    
else
    cd (folderName)
    for k=1:N_runs_forStatistics
        for i =1:N_Models
            name_model_rep = ['comp_',num2str(i),'_', num2str(k),'.mat'];
            load (name_model_rep)
            Report_Matrix(:,i)  = pre_Report_Matrix;
        end
        vector_Report_Matrix (:,:,k) = Report_Matrix;
        vector_Report_Matrix_norm (:,:,k) = (Report_Matrix'-min(Report_Matrix)')';
    end
    save compareModels.mat vector_Report_Matrix vector_Report_Matrix_norm
    load compareModels.mat
    cd ..
    N_runs_forStatistics = size(vector_Report_Matrix,3);
    
end

if selectedModel == N_Models || justPlot ==1
    %% Calculating Likelihood functions
    mean_LL= median(vector_Report_Matrix,3);
    SD_LL = std(vector_Report_Matrix,[],3);
    sum_LL_total = sum(mean_LL);
    sum_LL_wo_Stress = sum(mean_LL(1:3,:));
    sum_LL = sum(mean_LL);
    
    mean_LL_normalized= mean_LL - min(mean_LL,[],2);
    sum_LL_norm_total = sum(mean_LL_normalized);
    sum_LL_norm_wo_Stress = sum(mean_LL_normalized(1:3,:));
    
    SD_LL_normalized = SD_LL;
    % mean values fot the Log-likelihood
    complete_Report_Matrix = [mean_LL;sum_LL_total];
    complete_Report_Matrix_normalized = [mean_LL_normalized;sum_LL_norm_total];
    complete_Report_Matrix_normalized_wo_Stress =[mean_LL_normalized(1:3,:);sum_LL_norm_wo_Stress];
    
    % calculating SD form the sums
    for k =1: N_runs_forStatistics
        sum_matrix_LL(k,:) = sum(vector_Report_Matrix(:,:,k));
    end
    sum_SD_LL = std(sum_matrix_LL);
    sum_SD_LL_nromalized = sum_SD_LL;
    SD_Report_Matrix = [SD_LL;sum_SD_LL];
    SD_Report_Matrix_normalized =[SD_LL_normalized;sum_SD_LL_nromalized];
    
    
    % calculating SD form the sums without stresses
    for k =1: N_runs_forStatistics
        sum_matrix_LL_wo_stress(k,:) = sum(vector_Report_Matrix(1:3,:,k));
    end
    sum_SD_LL_wo_stress = std(sum_matrix_LL_wo_stress)./sqrt(N_runs_forStatistics);
    SD_Report_Matrix_normalized_wo_stress = [SD_LL(1:3,:);sum_SD_LL_wo_stress];
    
    % fitValue(1)=LL_Percent;
    % fitValue(2)=fit_LL_Distribution;
    % fitValue(3)=LL_HT;
    %
    % fitValue=[fitValue,LL_NA_1]; % 4 - 7  % type I. k_c = [0,33, 66, 100].
    % fitValue=[fitValue,LL_NA_2]; % 8 - 11 % type II. k_c = [0,33, 66, 100]
    %
    % fitValue=[fitValue,LL_DTT_1]; % 16 - 19 % type I. k_c = [0,33, 66, 100]
    % fitValue=[fitValue,LL_DTT_2]; % 20 - 23 % type II. k_c = [0,33, 66, 100]
    
    % NaAS stress type I
    mean_prediction_Na_1 = complete_Report_Matrix(4:7,:);
    %    mean_prediction_Na_1 = complete_Report_Matrix_normalized(4:7,:);
    sd_prediction_Na_1 = SD_Report_Matrix_normalized(4:7,:)./sqrt(N_runs_forStatistics);
    
    % NaAS stress type II
    mean_prediction_Na_2 = complete_Report_Matrix(8:11,:);
    %    mean_prediction_Na_2 = complete_Report_Matrix_normalized(8:11,:);
    sd_prediction_Na_2 = SD_Report_Matrix_normalized(8:11,:)./sqrt(N_runs_forStatistics);
    
    % DTT stress type I
    mean_prediction_DTT_1 = complete_Report_Matrix(12:15,:);
    %mean_prediction_DTT_1 = complete_Report_Matrix_normalized(12:15,:);
    sd_prediction_DTT_1 = SD_Report_Matrix_normalized(16:19,:)./sqrt(N_runs_forStatistics);
    
    % DTT stress type III
    mean_prediction_DTT_2 = complete_Report_Matrix(16:19,:);
    % mean_prediction_DTT_2 = complete_Report_Matrix_normalized(16:19,:);
    sd_prediction_DTT_2 = SD_Report_Matrix_normalized(20:23,:)./sqrt(N_runs_forStatistics);
    
    namePlot = 'Comparison_Na';
    typeStress=1;
    plottingLikelihood (typeStress,complete_Report_Matrix_normalized_wo_Stress, SD_Report_Matrix_normalized_wo_stress,mean_prediction_Na_1,sd_prediction_Na_1,mean_prediction_Na_2,sd_prediction_Na_2,folderName,namePlot,labels_Model_notSorted,n_freeParameters)
    
    namePlot = 'Comparison_DTT';
    typeStress=2;
    plottingLikelihood (typeStress,complete_Report_Matrix_normalized_wo_Stress, SD_Report_Matrix_normalized_wo_stress,mean_prediction_DTT_1,sd_prediction_DTT_1,mean_prediction_DTT_2,sd_prediction_DTT_2,folderName,namePlot,labels_Model_notSorted,n_freeParameters)
    
    % selected model
    sel_model=13;
    mostComplexModel =4;
    
    %% Report for optimization Best Fit.
    L_op_sel_model= [sum_LL_norm_wo_Stress(sel_model),sum_SD_LL_wo_stress(sel_model)]
    L_op_complex_model= [sum_LL_norm_wo_Stress(mostComplexModel),sum_SD_LL_wo_stress(mostComplexModel)]
    
    %% Reproted values for NaAS
    %     Best_crossvalidation=1; % 100% inibition or k_c =0
    [~, Best_crossvalidation] =min(mean_prediction_Na_1(:,sel_model));
    L_NA_typeI =[mean_prediction_Na_1(Best_crossvalidation,sel_model), sd_prediction_Na_1(Best_crossvalidation,sel_model)]
    
    %    Best_crossvalidation=2; % 66% inibition or k_c =33
    [~, Best_crossvalidation] =min(mean_prediction_Na_2(:,sel_model));
    L_NA_typeII =[mean_prediction_Na_2(Best_crossvalidation,sel_model),sd_prediction_Na_2(Best_crossvalidation,sel_model)]
    
    % Reproted values for DTT
    % Best_crossvalidation=1; % 100% inibition or k_c =0
    [~, Best_crossvalidation] =min(mean_prediction_DTT_1(:,sel_model));
    L_DTT_typeI = [mean_prediction_DTT_1(Best_crossvalidation,sel_model),sd_prediction_DTT_1(Best_crossvalidation,sel_model)]
    %     Best_crossvalidation=2; % 66% inibition or k_c =33
    [~, Best_crossvalidation] =min(mean_prediction_DTT_2(:,sel_model));
    L_DTT_typeII = [mean_prediction_DTT_2(Best_crossvalidation,sel_model),sd_prediction_DTT_2(Best_crossvalidation,sel_model)]
    
end

cd ..

