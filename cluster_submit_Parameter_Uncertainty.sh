#!/bin/csh
#$ -cwd
#$ -N PU
#$ -o PU.txt
#$ -e PU.err
#$ -pe matlab 12

# #$ -q defaultfaculty.q@node*
# #$ -q munsky-gpu.q@gpu*
#$ -q munsky.q@node*
# #$ -q short.q@node*
#$ -V

# /usr/local/bin/matlab -r "run__Cluster_Parameter_Uncertainty (Model_type numberTotalEvaluations Nrep N_sims justMakePlot) "

# Model_type = integer between 1 and 14. Represents the specific model use. For the model description use the readme file.
# numberTotalEvaluations = integer larger than 0. This number is the total number of parameters to evaluate.
# N_sims = integer larger than 0. The number of independent trajectories or stochastic simulations for each parameter set. 
# Nrep = integer larger than 0. The number of times the complete N_sims is repeated. Recommend to use more than 3.
# justMakePlot  = [0 or 1]; if justMakePlot = 0 it runs the simulations. If justMakePlot = 1 it uses the stored results to create the plots.

/usr/local/bin/matlab -r "run__Cluster_Parameter_Uncertainty 13 10000 3 3000 0 "
