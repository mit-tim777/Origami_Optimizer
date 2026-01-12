mkdir backup/long_sim_data1
mv -f MD_simulation/long_sim_data/* backup/long_sim_data1

rm display_energys.tcl
rm output.pdb
rm cpptraj.out
rm leap.out
rm leap.log
rm Offset_energy_calculator/MD_Results/*
rm Mutate/mutation_information.txt
rm MD_simulation/input.*
rm MD_simulation/info/*
rm MD_simulation/out/*
rm MD_simulation/restart_files/*
rm MD_simulation/trajectorys/*
rm Helix_separator/cpptraj_base_pairing.txt
