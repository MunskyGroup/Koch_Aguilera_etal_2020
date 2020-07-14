clear all
close all
clc

%% Extracting Experimental Data

strExperimentalData = extractingExperimentalData; 
ds_expTime = strExperimentalData.timeExp_NaAs;

data_1G = strExperimentalData.mean_exp_NaAs_CAP;
data_2G = strExperimentalData.mean_exp_NaAs_IRES;
data_BG = strExperimentalData.mean_exp_NaAs_CAP_IRES;

err_data_1G = strExperimentalData.sem_exp_NaAs_CAP;
err_data_2G = strExperimentalData.sem_exp_NaAs_IRES;
err_data_BG = strExperimentalData.sem_exp_NaAs_CAP_IRES;

matrix_NaAs_experimental = [ds_expTime', data_1G', err_data_1G' , data_2G', err_data_2G', data_BG', err_data_BG'];

%% Extracting Experimental Data
ds_expTime_DTT = strExperimentalData.timeExp_DTT;

data_1G_DTT = strExperimentalData.mean_exp_DTT_CAP;
data_2G_DTT = strExperimentalData.mean_exp_DTT_IRES;
data_BG_DTT = strExperimentalData.mean_exp_DTT_CAP_IRES;

err_data_1G_DTT = strExperimentalData.sem_exp_DTT_CAP;
err_data_2G_DTT = strExperimentalData.sem_exp_DTT_IRES;
err_data_BG_DTT = strExperimentalData.sem_exp_DTT_CAP_IRES;

matrix_DTT_experimental = [ds_expTime_DTT', data_1G_DTT', err_data_1G_DTT' , data_2G_DTT', err_data_2G_DTT', data_BG_DTT', err_data_BG_DTT'];

%% Extracting Simulated Data for NaAs type II
cd ('Selected_Model_13')

load sim_data_NaAS_CAP_typeII.mat
% value_sweep
% sim_CAP_ONLY
% error_CAP_ONLY

load sim_data_NaAS_IRES_typeII.mat
% value_sweep
% sim_IRES_ONLY
% error_IRES_ONLY

load sim_data_NaAS_CAP_IRES_typeII.mat
% value_sweep
% sim_CI
% error_CI

matrix_NaAs_Simulations_typeII = [sim_CAP_ONLY',error_CAP_ONLY', sim_IRES_ONLY', error_IRES_ONLY',sim_CI',error_CI' ];

%% Extracting Simulated Data for DTT type I

load sim_data_DTT_CAP_typeI.mat
% value_sweep
% sim_CAP_ONLY
% error_CAP_ONLY

load sim_data_DTT_IRES_typeI.mat
% value_sweep
% sim_IRES_ONLY
% error_IRES_ONLY

load sim_data_DTT_CAP_IRES_typeI.mat
% value_sweep
% sim_CI
% error_CI

matrix_DTT_Simulations_typeI = [sim_CAP_ONLY',error_CAP_ONLY', sim_IRES_ONLY', error_IRES_ONLY',sim_CI',error_CI' ];

cd ..