clear all
close all
clc

cd ('Selected_Model_13')

load data_Figure_5_C_cap

probability_cap = mean_stairs_1G([1:2:end])';
sd_probability_cap = std_stairs_1G([1:2:end])';

report_cap = [probability_cap,sd_probability_cap]


load data_Figure_5_C_ires

probability_ires = mean_stairs_2G([1:2:end])';
sd_probability_ires = std_stairs_2G([1:2:end])';

report_ires = [probability_ires,sd_probability_ires]