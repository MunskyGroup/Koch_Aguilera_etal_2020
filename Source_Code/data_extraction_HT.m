clear all
clc

cd ('Selected_Model_13')
load data_Figure_5_D

report_CAP = [ds_expTime' , dataTrajectories_CAP',mean_sim_1G',err_sim_1G'];
 
load data_Figure_5_E

report_IRES = [ds_expTime' , dataTrajectories_IRES',mean_sim_2G',err_sim_2G'];

cd ..