
rm -rf previous_iteration/iteration_$1/*
cp output.pdb previous_iteration/iteration_$1/

mv MD_simulation/trajectorys previous_iteration/iteration_$1/
mv MD_simulation/info previous_iteration/iteration_$1/
mv MD_simulation/out previous_iteration/iteration_$1/
mv MD_simulation/restart_files previous_iteration/iteration_$1/
mv display_energys.tcl previous_iteration/iteration_$1/
mv Mutate/mutation_information.txt previous_iteration/iteration_$1/
mv cpptraj_base_pairing.txt previous_iteration/iteration_$1/
mv MD_simulation/input.* previous_iteration/iteration_$1/
mv cpptraj.out previous_iteration/iteration_$1/
mv leap.out previous_iteration/iteration_$1/
mv leap.log previous_iteration/iteration_$1/
mv Offset_energy_calculator/MD_Results/ previous_iteration/iteration_$1/

mkdir MD_simulation/trajectorys
mkdir MD_simulation/out
mkdir MD_simulation/info
mkdir MD_simulation/restart_files
mkdir Offset_energy_calculator/MD_Results


