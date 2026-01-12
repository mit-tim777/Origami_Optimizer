#!/bin/bash


# conda activate AmberTools25

bash restart.sh

cp input/*.pdb output.pdb

bash one_optimization_iteration.sh 0

# for i in {1..1}; do
#     python Mutate/delete_residues.py
#     bash one_optimization_iteration.sh $i
# done

