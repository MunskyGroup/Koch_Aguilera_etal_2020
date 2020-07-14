clear all
clc

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

report_experimental_percentage = [mean1_percentage_Non_Translating,sem_percentage_Non_Translating;...
                        mean1_percentage_Cap,sem_percentage_Cap; ...
                        mean1_percentage_IRES, sem_percentage_IRES;...
                        mean1_percentage_Cap_IRES, sem_percentage_Cap_IRES ]


cd ..   


cd ('Selected_Model_13')
load data_Figure_5_B

report_simulation_percentage = [sim_NT,err_exp_NT;...
                        sim_1G,err_sim_1G; ...
                        sim_2G, err_sim_2G;...
                        sim_BG, err_sim_BG ]
                    
cd ..
                    



