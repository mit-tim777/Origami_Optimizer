#!/bin/bash


bash bash_scripts/clear.sh

cp input/*.pdb output.pdb

bash bash_scripts/one_optimization_iteration.sh 0 "restrained" "nodebug"
bash bash_scripts/one_optimization_iteration.sh 1 "restrained" "nodebug"

cp input/*.pdb output.pdb

bash bash_scripts/one_optimization_iteration.sh 2 "unrestrained" "nodebug"
bash bash_scripts/one_optimization_iteration.sh 3 "unrestrained" "nodebug"



# bash bash_scripts/one_optimization_iteration.sh 1
# bash bash_scripts/mutate.sh

# bash bash_scripts/restart.sh 1


# bash bash_scripts/one_optimization_iteration.sh 2
