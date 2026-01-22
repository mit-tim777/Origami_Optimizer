#!/bin/bash


bash bash_scripts/restart.sh 100


cp input/*.pdb output.pdb

bash bash_scripts/one_optimization_iteration.sh 0
# bash bash_scripts/mutate.sh  # only works rn in input is already standardized

# bash bash_scripts/restart.sh 0


# bash bash_scripts/one_optimization_iteration.sh 1
# bash bash_scripts/mutate.sh

# bash bash_scripts/restart.sh 1


# bash bash_scripts/one_optimization_iteration.sh 2
