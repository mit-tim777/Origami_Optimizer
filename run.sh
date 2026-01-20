#!/bin/bash


# conda activate AmberTools25

bash restart.sh

cp input/*.pdb output.pdb

bash one_optimization_iteration.sh 0
python Mutate/delete_residues.py

rm -rf MD_simulation/long_sim_data/*
mv MD_simulation/trajectorys MD_simulation/long_sim_data/
mv MD_simulation/info MD_simulation/long_sim_data/
mv MD_simulation/out MD_simulation/long_sim_data/
mv MD_simulation/restart_files MD_simulation/long_sim_data/
mv MD_simulation/input.prmtop MD_simulation/long_sim_data/
mkdir MD_simulation/trajectorys
mkdir MD_simulation/out
mkdir MD_simulation/info
mkdir MD_simulation/restart_files

bash one_optimization_iteration.sh 1

# for i in {1..1}; do
# done

