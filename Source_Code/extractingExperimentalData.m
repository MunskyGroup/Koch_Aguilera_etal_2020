function strExperimentalData = extractingExperimentalData
cd Experimental_Data

%% Data representing the Fraction of spots
strExperimentalData.NumberForPercentage = 3748; % total number of spots

% Loading data for percentages
fileName = 'PercentPerCell_OriginalTag.xlsx';
[Data_Percentages,~,~]=xlsread(fileName);
percentage_Non_Translating = Data_Percentages(1,:);
percentage_Cap = Data_Percentages(2,:);
percentage_IRES = Data_Percentages(3,:);
percentage_Cap_IRES  = Data_Percentages(4,:);

% Loading number of spots per cell
fileName = 'NumbersPerCell_OriginalTag.xlsx';
[Data_Percentages_numberSpots,~,~]=xlsread(fileName);
no_spots_percentage_Non_Translating = Data_Percentages_numberSpots(1,:);
no_spots_percentage_Cap = Data_Percentages_numberSpots(2,:);
no_spots_percentage_IRES = Data_Percentages_numberSpots(3,:);
no_spots_percentage_Cap_IRES  = Data_Percentages_numberSpots(4,:);
total_no_spots_cell = sum(Data_Percentages_numberSpots);

% Normalizing mean values by the number of spots
norm_by_number_spots_mean1_percentage_Cap = mean(percentage_Cap./total_no_spots_cell);
norm_by_number_spots_mean1_percentage_IRES = mean(percentage_IRES./total_no_spots_cell);
norm_by_number_spots_mean1_percentage_Cap_IRES = mean(percentage_Cap_IRES./total_no_spots_cell);
norm_by_number_spots_mean1_percentage_Non_Translating =  mean(percentage_Non_Translating./total_no_spots_cell);
sum_mean1 = norm_by_number_spots_mean1_percentage_Cap+norm_by_number_spots_mean1_percentage_IRES+ norm_by_number_spots_mean1_percentage_Cap_IRES+ norm_by_number_spots_mean1_percentage_Non_Translating;

% Second Normalization by the total percentage
mean1_percentage_Cap = (norm_by_number_spots_mean1_percentage_Cap / sum_mean1)*100;
mean1_percentage_IRES = (norm_by_number_spots_mean1_percentage_IRES / sum_mean1 )*100;
mean1_percentage_Cap_IRES = (norm_by_number_spots_mean1_percentage_Cap_IRES / sum_mean1 )*100;
mean1_percentage_Non_Translating = (norm_by_number_spots_mean1_percentage_Non_Translating/ sum_mean1 )*100;

% Claculating SEM for the original data
sem_percentage_Cap = std(percentage_Cap)/sqrt(mean(no_spots_percentage_Cap));
sem_percentage_IRES = std(percentage_IRES)/sqrt(mean(no_spots_percentage_IRES));
sem_percentage_Cap_IRES = std(percentage_Cap_IRES)/sqrt(mean(no_spots_percentage_Cap_IRES));
sem_percentage_Non_Translating = std(percentage_Non_Translating)/sqrt(mean(no_spots_percentage_Non_Translating));

% Saving mean values and SEM for percentages.

strExperimentalData.Percentage_Cap = mean1_percentage_Cap;
strExperimentalData.Percentage_Ires = mean1_percentage_IRES;
strExperimentalData.Percentage_Cap_Ires = mean1_percentage_Cap_IRES;
strExperimentalData.Percentage_NT = mean1_percentage_Non_Translating;

strExperimentalData.Percentage_err_Cap = sem_percentage_Cap;
strExperimentalData.Percentage_err_Ires = sem_percentage_IRES;
strExperimentalData.Percentage_err_Cap_Ires = sem_percentage_Cap_IRES;
strExperimentalData.Percentage_err_NT = sem_percentage_Non_Translating;

%% Intensity Data
calibrated_No_Ribosomes = 14.652; % calibrated by Amanda.  Cap in Cap Only (the calibration)
fileName = 'Cap_Only_mRNA_Threshold.xls';
[Data_CAP,~,~]=xlsread(fileName);
Int_cap= Data_CAP(:,2);
% transforming intensities to UMP
mean_Intensity_Associated_to_Calibration = mean  (Int_cap);
% the scaling factor transforms intensities to UMP. It is used to convert
% units for both genes, beacuse the same tag (SM) is used. in the 'normal
% tag' and 'Switch Tag' constructs.
scalingFactor_FromIntensity_to_UMP = calibrated_No_Ribosomes/mean_Intensity_Associated_to_Calibration;
cap_ump= Int_cap .* scalingFactor_FromIntensity_to_UMP;
fileName = 'SwitchTag_IRESOnly_mRNA_Threshold.xls';
[Data_IRES,~,~]=xlsread(fileName);
Int_ires= Data_IRES(:,2);
ires_ump= Int_ires .* scalingFactor_FromIntensity_to_UMP;
strExperimentalData.Int_cap = cap_ump;
strExperimentalData.Int_ires = ires_ump;

%% Scatter Plot
fileName = 'Cap_IRES_mRNA_Threshold.xls';
[Data_Scatter,~,~]=xlsread(fileName);
Int_ST_CAP = Data_Scatter(:,2);
strExperimentalData.Int_ST_CAP = Int_ST_CAP.* scalingFactor_FromIntensity_to_UMP;
Int_ST_IRES = Data_Scatter(:,3);
strExperimentalData.Int_ST_IRES = Int_ST_IRES.* scalingFactor_FromIntensity_to_UMP;
% plot(strExperimentalData.Int_ST_CAP,strExperimentalData.Int_ST_IRES,'o','MarkerSize',10,'MarkerEdge','k','MarkerFaceColor','g')
% [rho,pval] = corr(strExperimentalData.Int_ST_IRES,strExperimentalData.Int_ST_CAP);

%% Harringtonine Data
% CAP Harringtonine Data. Notice that to make compatible simulations and
% experimental data the first time point is ignored.
fileName = 'CapIntPerCell_HT.xlsx';
[Data_HT_CAP,~,~]=xlsread(fileName);
strExperimentalData.timeExp_Harringtonine = Data_HT_CAP(1,2:36);
Data_HT_CAP=Data_HT_CAP(2:end,2:36);
original_Data_HT_CAP=Data_HT_CAP;

% % normalizing data with respect to the first 5 time points
for i=1:size(Data_HT_CAP,1)
    normalizer = mean (Data_HT_CAP(i,1:5));
    Data_HT_CAP(i,:) = Data_HT_CAP(i,:) ./ normalizer;
end
% % Calculating the error
sem_exp_HT_CAP = std(Data_HT_CAP)./sqrt(size(Data_HT_CAP,1));
% % Calculating the mean
mean_exp_HT_CAP = mean(Data_HT_CAP);
% % removing the non-responding fraction
mean_exp_HT_CAP = mean_exp_HT_CAP- mean(mean_exp_HT_CAP(end-6:end));
% % second normalization with respect to the first 5 time points
mean_exp_HT_CAP = mean_exp_HT_CAP./mean (mean_exp_HT_CAP(1:5));

% IRES Harringtonine Data. Notice that to make compatible simulations and
% experimental data the first time point is ignored.
fileName = 'IRESIntPerCell_HT.xlsx';
[Data_HT_IRES,~,~]=xlsread(fileName);
Data_HT_IRES=Data_HT_IRES(2:end,2:36);
original_Data_HT_IRES=Data_HT_IRES;
% % normalizing data with respect to the first 5 time points
for i=1:size(Data_HT_IRES,1)
    normalizer = mean (Data_HT_IRES(i,1:5));
    Data_HT_IRES(i,:) = Data_HT_IRES(i,:) ./ normalizer;
end
% % Calculating the error
sem_exp_HT_IRES = std(Data_HT_IRES)/sqrt(size(Data_HT_IRES,1));
% % Calculating the mean
mean_exp_HT_IRES = mean(Data_HT_IRES);
% % removing the non-responding fraction
mean_exp_HT_IRES = mean_exp_HT_IRES- mean(mean_exp_HT_IRES(end-6:end));
% % second normalization with respect to the first 5 time points
mean_exp_HT_IRES = mean_exp_HT_IRES./mean (mean_exp_HT_IRES(1:5));
% saving Harrintonine Data to Structure
strExperimentalData.trajectories_exp_CAP = Data_HT_CAP;
strExperimentalData.trajectories_exp_IRES = Data_HT_IRES;

strExperimentalData.mean_exp_HT_IRES = mean_exp_HT_IRES;
strExperimentalData.sem_exp_HT_IRES = sem_exp_HT_IRES;

strExperimentalData.mean_exp_HT_CAP = mean_exp_HT_CAP;
strExperimentalData.sem_exp_HT_CAP = sem_exp_HT_CAP;

strExperimentalData.repetitions_exp_HT_CAP = size(Data_HT_CAP,1);
strExperimentalData.repetitions_exp_HT_IRES = size(Data_HT_IRES,1);

%% NaAs Data
% NaAs Cap in cap-only
fileName = 'greenInYellowTotalInt_NaAs.xlsx';
[Data_NaAs_CAP,~,~]=xlsread(fileName);
strExperimentalData.mean_exp_NaAs_CAP = Data_NaAs_CAP(:,2)';
strExperimentalData.sem_exp_NaAs_CAP = Data_NaAs_CAP(:,3)';
% strExperimentalData.repetitions_exp_NaAs_CAP = Data_NaAs_CAP(:,4)';
strExperimentalData.timeExp_NaAs = Data_NaAs_CAP(:,1)';
% NaAs IRES in IRES-only
fileName = 'blueInPurpleTotalInt_NaAs.xlsx';
[Data_NaAs_IRES,~,~]=xlsread(fileName);
strExperimentalData.mean_exp_NaAs_IRES = Data_NaAs_IRES(:,2)';
strExperimentalData.sem_exp_NaAs_IRES = Data_NaAs_IRES(:,3)';
% strExperimentalData.repetitions_exp_NaAs_IRES = Data_NaAs_IRES(:,4)';
% NaAs CAP in CAP-IRES
fileName = 'greenInWhiteTotalInt_NaAs.xlsx';
[Data_NaAs_CAP_IRES,~,~]=xlsread(fileName);
strExperimentalData.mean_exp_NaAs_CAP_IRES = Data_NaAs_CAP_IRES(:,2)';
strExperimentalData.sem_exp_NaAs_CAP_IRES = Data_NaAs_CAP_IRES(:,3)';
% strExperimentalData.repetitions_exp_NaAs_CAP_IRES = Data_NaAs_CAP_IRES(:,4)';
% NaAs IRES in CAP-IRES
fileName = 'blueInWhiteTotalInt_NaAs.xlsx';
[Data_NaAs_IRES_CAPIRES,~,~]=xlsread(fileName);
strExperimentalData.mean_exp_NaAs_IRES_CAPIRES = Data_NaAs_IRES_CAPIRES(:,2)';
strExperimentalData.sem_exp_NaAs_IRES_CAPIRES = Data_NaAs_IRES_CAPIRES(:,3)';
% strExperimentalData.repetitions_exp_NaAs_IRES_CAPIRES = Data_NaAs_IRES_CAPIRES(:,4)';

%% DTT data
fileName = 'greenInYellowTotalInt_DTT.xlsx';
[Data_DTT_CAP,~,~]=xlsread(fileName);
strExperimentalData.mean_exp_DTT_CAP = Data_DTT_CAP(1:35,2)';
strExperimentalData.sem_exp_DTT_CAP = Data_DTT_CAP(1:35,3)';
strExperimentalData.timeExp_DTT = Data_DTT_CAP(1:35,1)';

fileName = 'blueInPurpleTotalInt_DTT.xlsx';
[Data_DTT_IRES,~,~]=xlsread(fileName);
strExperimentalData.mean_exp_DTT_IRES = Data_DTT_IRES(1:35,2)';
strExperimentalData.sem_exp_DTT_IRES = Data_DTT_IRES(1:35,3)';

fileName = 'greenInWhiteTotalInt_DTT.xlsx';
[Data_DTT_CAP_IRES,~,~]=xlsread(fileName);
strExperimentalData.mean_exp_DTT_CAP_IRES = Data_DTT_CAP_IRES(1:35,2)';
strExperimentalData.sem_exp_DTT_CAP_IRES = Data_DTT_CAP_IRES(1:35,3)';

fileName = 'blueInWhiteTotalInt_DTT.xlsx';
[Data_DTT_IRES_CAPIRES,~,~]=xlsread(fileName);
strExperimentalData.mean_exp_DTT_IRES_CAPIRES = Data_DTT_IRES_CAPIRES(1:35,2)';
strExperimentalData.sem_exp_DTT_IRES_CAPIRES = Data_DTT_IRES_CAPIRES(1:35,3)';

cd ..
end