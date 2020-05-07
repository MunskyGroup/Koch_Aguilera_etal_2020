#!/bin/csh
#$ -cwd
#$ -N o14
#$ -o o14.txt
#$ -e o14.err
#$ -pe matlab 12

# #$ -q defaultfaculty.q@node*
# #$ -q munsky-gpu.q@gpu*
#$ -q munsky.q@node*
# #$ -q short.q@node*
#$ -V

# /usr/local/bin/matlab -r "run__Cluster_Compare_All_Models(run_SpecificModel,justPlot,runAllModels)"
# run_SpecificModel = integer between 0 and 14. If run_SpecificModel is empty or 0 it is ignored and can be used with runAllModels. If run_SpecificModel>0 represents the specific model use. For the model description use the readme file.
# justMakePlot  = [0 or 1]; if justMakePlot = 0 it runs the simulations. If justMakePlot = 1 it uses the stored results to create the plots.
# runAllModels = [0 or 1]; if justMakePlot = 1 it runs the simulations for all models.

/usr/local/bin/matlab -r "run__Cluster_Compare_All_Models 0 1 0 "
