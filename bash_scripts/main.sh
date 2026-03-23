#!/bin/bash

set -e

bash bash_scripts/clear.sh


for i in {1..20}; do
    echo "-----------------------------------------------------------------------------------------------------------------------------------"
    cp template_pdb/*.pdb input/input.pdb

    tleap -f MD_simulation/scripts/standardise_pdb.in  >leap.out

    cp input/*.pdb output.pdb

    python Helix_separator/find_bound_double_strands.py

    python Mutate/assign_new_starting_sequence.py

    cp input/*.pdb output.pdb

    bash bash_scripts/one_optimization_iteration.sh 1 "restrained" "nodebug"
    bash bash_scripts/one_optimization_iteration.sh 2 "restrained" "nodebug"

    cp input/*.pdb output.pdb

    bash bash_scripts/one_optimization_iteration.sh 3 "unrestrained" "nodebug"
    bash bash_scripts/one_optimization_iteration.sh 4 "unrestrained" "nodebug"

    bash bash_scripts/restart_sequence.sh $i

done
