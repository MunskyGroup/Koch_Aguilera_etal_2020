#!/bin/csh
#$ -cwd
#$ -N opt
#$ -o opt.txt
#$ -e opt.err
#$ -pe matlab 12

# #$ -q defaultfaculty.q@node*
# #$ -q munsky-gpu.q@gpu*
#$ -q munsky.q@node*
# #$ -q short.q@node*
#$ -V

# /usr/local/bin/matlab -r "run__Cluster_Optimization (Model_type,randInit, optimizationMethod, PopulationSize, generations, maxSimulationTime)"

# Model_type = integer between 1 and 14. Represents the specific model use. For the model description use the readme file.
# randInit = [0 or 1]. If randInit =1 it creates generates an random initial population.
# Optimization Metods 1 = Multi-onjective Optimization. 2 = GA. 3 = Pattern Search.
# PopulationSize = integer larger than 0. Populization size for the genetic algorithm.
# generations = integer larger than 0. The number of generations for the genetic algorithm.
# maxSimulationTime = integer larger than 0. This number is the total number of hours that the algorithm is allowed to run.

/usr/local/bin/matlab -r "run__Cluster_Optimization 13 1 2 30 20 10 "
