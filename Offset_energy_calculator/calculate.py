import csv
from pathlib import Path
import re
import numpy as np
import os, os.path
import random
import copy

def print_helix_text_reprensentation(helix):
    
    for b in helix['strand_sequences'][0]:
        print(b, end = " ")
    print("")
    for b in helix['strand_sequences'][0]:
        print('|', end = " ")
    print("")
    for b in helix['strand_sequences'][1]:
        print(b, end = " ")
    print()
    print()


def extract_data(filename): # read out the averaged helical parameters for one helix

    with open(filename) as f:
        
        bp_params = []  # [ [shear,stretch,stagger,buckle,prop,open] , same for basepair 2 , ...]
        step_params = [] # [ [shift,slide,rise,tilt,roll,twist] , same for step 2 , ...]  ( bp1 step1 bp2 step2 ... ) 
        heli_params = [] # [ [x_disp,y_disp,hrise,incl,tip,htwist] , step 2 , ...]
        bp_params_sd = []  # standard deviations of the parameters across the trajectory
        step_params_sd = []
        heli_params_sd = []
        strands = []

        meta_data = f.readline().split()
        strands.append(meta_data[1])
        strands.append(meta_data[2])
        
        strand_res_inds = [[int(i.split('-')[j])-1 for i in meta_data[0].split(',')[0:-1]] for j in range(2)]

        next(f)
        for line in f:
            params_raw = [ float(i) for i in line.split()]
            if(params_raw == []):
                continue
            bp_params.append(params_raw[1:7])
            step_params.append(params_raw[7:13])
            heli_params.append(params_raw[13:19])
            bp_params_sd.append(params_raw[19:25])
            step_params_sd.append(params_raw[25:31])
            heli_params_sd.append(params_raw[31:37])
        
        step_params.pop(-1)
        heli_params.pop(-1)
        step_params_sd.pop(-1)
        heli_params_sd.pop(-1)
        
        helix = {
            'strand_sequences' : strands,
            'strand_res_inds' : strand_res_inds,
            'bp_params' : bp_params,
            'step_params' : step_params,
            'heli_params' : heli_params,
            'bp_params_sd' : bp_params_sd,
            'step_params_sd' : step_params_sd,
            'heli_params_sd' : heli_params_sd,
            'energys' : None,
            'displacements' : None,
            'stiffs' : None,
            'eq_params' : None
        }


    return helix

def safe_float(x):
    try:
        return float(x)
    except ValueError:
        return x 

def load_equalibrium_params(): # load the equalibrium parameters from the paper "Sequence-Dependent Shape and Stiffness of DNA and RNA Double Helices"
    with open(root_dir / "Offset_energy_calculator" / "hexamers_csv" / "DNA" / "coords_grooves_DNA_hexamers_table.csv") as f:
        data = csv.reader(f)
        equalibrium_step_params = { row[0] : [safe_float(i) for i in row[1:7]] for row in data }
    with open(root_dir / "Offset_energy_calculator" / "hexamers_csv" / "DNA" / "coords_grooves_DNA_hexamers_table.csv") as f:
        data = csv.reader(f)
        equalibrium_heli_params = { row[0] : [safe_float(i) for i in row[7:13]] for row in data }
    with open(root_dir / "Offset_energy_calculator" / "hexamers_csv" / "DNA" / "coords_grooves_DNA_heptamers_table.csv") as f:
        data = csv.reader(f)
        equalibrium_bp_params = { row[0] : [safe_float(i) for i in row[1:7]] for row in data }
    equalibrium_params = {
        'bp' : equalibrium_bp_params,
        'step' : equalibrium_step_params,
        'heli' : equalibrium_heli_params
    }
    return equalibrium_params

def load_stiffs():  # load the quadratic offset energy stiffness
    with open(root_dir / "Offset_energy_calculator" / "hexamers_csv" / "DNA" / "coords_stiffs_DNA_hexamers_table.csv") as f:
        data = csv.reader(f)
        step_stiffs = { row[0] : [safe_float(i) for i in row[1:7]] for row in data }
    with open(root_dir / "Offset_energy_calculator" / "hexamers_csv" / "DNA" / "coords_stiffs_DNA_hexamers_table.csv") as f:
        data = csv.reader(f)
        heli_stiffs = { row[0] : [safe_float(i) for i in row[7:13]] for row in data }
    with open(root_dir / "Offset_energy_calculator" / "hexamers_csv" / "DNA" / "coords_stiffs_DNA_heptamers_table.csv") as f:
        data = csv.reader(f)
        bp_stiffs = { row[0] : [safe_float(i) for i in row[1:7]] for row in data }
    stiffs = {
        'bp'   :   bp_stiffs,
        'step' : step_stiffs,
        'heli' : heli_stiffs
    }
    return stiffs

def get_all_possible_sequences(sequence):   # given any sequence of any length and '-' as placeholder, return all possible sequences by replacing '-' with A,T,C or G
    sequence = "".join(sequence)
    placeholder_inds = [ index for index, val in enumerate(sequence) if val == '-']
    possible_sequences = [sequence]
    for i in placeholder_inds:
        next_possible_sequences = []
        for seq in possible_sequences:
            for base in ['A','T','C','G']:
                seq2 = list(seq)
                seq2[i] = base
                next_possible_sequences.append(''.join(seq2))
        possible_sequences = next_possible_sequences.copy() 
    return possible_sequences

def get_equalibrium_params(param_type, sequence):   # for sequences with placeholders '-', return average equalibrium parameters ( nessasary since equalibrium parameters are always defined for hexamers or heptamers)
    possible_sequences = get_all_possible_sequences(sequence)
    avg_params = [0]*6

    for seq in possible_sequences:
        for i in range(6):
            avg_params[i] += equalibrium_params[param_type][seq][i]

    for i in range(6):
        avg_params[i] /= len(possible_sequences)

    return avg_params

def get_stiffs(param_type, sequence):   # for sequences with placeholders '-', return average stiffness parameters ( nessasary since stiffness parameters are always defined for hexamers or heptamers) 
    possible_sequences = get_all_possible_sequences(sequence)
    avg_params = [0]*6

    for seq in possible_sequences:
        for i in range(6):
            avg_params[i] += stiffs[param_type][seq][i]

    for i in range(6):
        avg_params[i] /= len(possible_sequences)

    return avg_params
  
def compare_steps_to_eql(helix):  # calculate the offset energy for all step parameters compared to equalibrium
    differences = [[0,0,0,0,0,0]] * len(helix['step_params'])
    stiffnesses = [[0,0,0,0,0,0]] * len(helix['step_params'])
    energys = [[0,0,0,0,0,0]] * len(helix['step_params'])
    equalibrium_params = [[0,0,0,0,0,0]] * len(helix['step_params'])

    seq = helix['strand_sequences'][0]
    for step_n in range(0, len(helix['step_params'])):
        diffs_of_step = []
        energys_of_step = []
        hex_seq = ''.join([ seq[i] if (i in range(len(seq))) else '-' for i in range(step_n-2,step_n+4)])
        step_eql_params = get_equalibrium_params('step', hex_seq)
        step_stiffs = get_stiffs('step', hex_seq)
        for j in range(6):
            diff = helix['step_params'][step_n][j] - step_eql_params[j]
            energy = 0.5 * step_stiffs[j] * pow(diff, 2)
            diffs_of_step.append(diff)
            energys_of_step.append(energy)
        energys[step_n] = energys_of_step
        differences[step_n] = diffs_of_step
        stiffnesses[step_n] = step_stiffs
        equalibrium_params[step_n] = step_eql_params
    return(energys, stiffnesses, equalibrium_params, differences)

def compare_heli_to_eql(helix):  # calculate the offset energy for all helical parameters compared to equalibrium
    differences = [[0,0,0,0,0,0]] * len(helix['heli_params'])
    stiffnesses = [[0,0,0,0,0,0]] * len(helix['heli_params'])
    energys = [[0,0,0,0,0,0]] * len(helix['heli_params'])
    equalibrium_params = [[0,0,0,0,0,0]] * len(helix['heli_params'])

    seq = helix['strand_sequences'][0]
    for step_n in range(0, len(helix['heli_params'])):
        diffs_of_step = []
        energy_of_step = []
        hex_seq = ''.join([ seq[i] if (i in range(len(seq))) else '-' for i in range(step_n-2,step_n+4)])
        step_eql_params = get_equalibrium_params('heli', hex_seq)
        step_stiffs = get_stiffs('heli', hex_seq)
        for j in range(6):
            diff = helix['heli_params'][step_n][j] - step_eql_params[j]
            energy = 0.5 * step_stiffs[j] * pow(diff, 2)
            diffs_of_step.append(diff)
            energy_of_step.append(energy)
        energys[step_n] = energy_of_step
        differences[step_n] = diffs_of_step
        stiffnesses[step_n] = step_stiffs
        equalibrium_params[step_n] = step_eql_params
    return(energys, stiffnesses, equalibrium_params, differences)

def compare_bp_to_eql(helix):   # calculate the offset energy for all base pair parameters compared to equalibrium
    differences = [[0,0,0,0,0,0]] * len(helix['bp_params'])
    stiffnesses = [[0,0,0,0,0,0]] * len(helix['bp_params'])
    energys = [[0,0,0,0,0,0]] * len(helix['bp_params'])
    equalibrium_params = [[0,0,0,0,0,0]] * len(helix['bp_params'])

    seq = helix['strand_sequences'][0]
    for bp_n in range(0, len(helix['bp_params'])):
        diffs_of_bp = []
        energys_of_bp = []
        hep_seq = ''.join([ seq[i] if (i in range(len(seq))) else '-' for i in range(bp_n-3,bp_n+4)])
        bp_eql_params = get_equalibrium_params('bp',hep_seq)
        bp_stiffs = get_stiffs('bp',hep_seq)
        for j in range(6):
            diff = helix['bp_params'][bp_n][j] - bp_eql_params[j]
            energy = 0.5 * bp_stiffs[j] * pow(diff, 2)
            diffs_of_bp.append(diff)
            energys_of_bp.append(energy)
        energys[bp_n] = energys_of_bp
        differences[bp_n] = diffs_of_bp
        stiffnesses[bp_n] = bp_stiffs
        equalibrium_params[bp_n] = bp_eql_params
    return(energys, stiffnesses, equalibrium_params, differences)


def calculate_displacement_energy(helix):  
    energys = {
        'bp'   : compare_bp_to_eql(helix)[0],
        'step' : compare_steps_to_eql(helix)[0], 
        'heli' : compare_heli_to_eql(helix)[0]
    }
    stiffs = {
        'bp'   : compare_bp_to_eql(helix)[1],
        'step' : compare_steps_to_eql(helix)[1],
        'heli' : compare_heli_to_eql(helix)[1]
    }
    equalibrium_params = {
        'bp'   : compare_bp_to_eql(helix)[2],
        'step' : compare_steps_to_eql(helix)[2], 
        'heli' : compare_heli_to_eql(helix)[2]
    }
    differences = {
        'bp'   : compare_bp_to_eql(helix)[3],
        'step' : compare_steps_to_eql(helix)[3], 
        'heli' : compare_heli_to_eql(helix)[3]
    }

    return (energys, stiffs, equalibrium_params, differences)

def calculate_displacement_energy_alternate_sequence(helix, sequence):   # assuming the helix had another sequence, then calculate the displacement energy
    complements = {
        'A' : 'T', 'T' : 'A', 'C' : 'G', 'G' : 'C'
    }
    complementary = "".join([complements[i] for i in sequence])
    helix_cpy = copy.deepcopy(helix)
    helix_cpy['strand_sequences'] = [ "".join(sequence) , complementary]

    energys = {
        'bp'   : compare_bp_to_eql(helix_cpy)[0],
        'step' : compare_steps_to_eql(helix_cpy)[0], 
        'heli' : compare_heli_to_eql(helix_cpy)[0]
    }
    return (energys)

def write_tcl_representation_script(helices):   # the display_energys.tcl file can be loadid with vmd to visualize the energy of the helix deformations

    COLORING_FACTOR = 0.3
    with open('display_energys.tcl', 'w') as f:

        f.write('mol new output.pdb \n')
        f.write('set molid [lindex [expr {[molinfo list]}] end] \n\n\n')
        
        for helix in helices:

            energy_sums = [sum(i) for i in helix['energys']['bp']]
            color_param_bp = [i*COLORING_FACTOR for i in energy_sums]

            energy_sums = [sum(helix['energys']['step'][i]) + sum(helix['energys']['heli'][i]) for i in range(len(helix['energys']['step']))]
            color_param_step = [i*COLORING_FACTOR for i in energy_sums]

            selection_string = '"backbone and not name OP1 OP2"'
            f.write('mol addrep $molid \n')
            f.write('set repindex [expr {[molinfo $molid get numreps] - 1}] \n')
            f.write('mol modselect $repindex $molid ' + selection_string + ' \n')
            f.write('mol modstyle $repindex $molid Licorice 0.8 12.0 \n')

            f.write('color scale method GWR \n')
            f.write('mol modcolor $repindex $molid user \n')

            for i in range(len(helix['strand_res_inds'][0])-1):
                selection_string = '"backbone and residue ' + str(helix['strand_res_inds'][0][i]) + " " + str(helix['strand_res_inds'][1][i]) + ' and not name OP1 OP2"'

                # f.write('mol addrep $molid \n')
                # f.write('set repindex [expr {[molinfo $molid get numreps] - 1}] \n')
                # f.write('mol modselect $repindex $molid ' + selection_string + ' \n')
                # f.write('mol modstyle $repindex $molid Licorice 0.8 12.0 \n')
                
                f.write('set sel [atomselect $molid ' + selection_string + '] \n')
                f.write('$sel set user '+ str((color_param_step[i])) +' \n')
                f.write('$sel delete \n')
                
                f.write(' \n')
            f.write(' \n')
            f.write(' \n')
            for i in range(len(helix['strand_res_inds'][0])):
                selection_string = '"not backbone and residue ' + str(helix['strand_res_inds'][0][i]) + " " + str(helix['strand_res_inds'][1][i]) + '"'
                
                f.write('mol addrep $molid \n')
                f.write('set repindex [expr {[molinfo $molid get numreps] - 1}] \n')
                f.write('mol modselect $repindex $molid ' + selection_string + ' \n')
                f.write('mol modstyle $repindex $molid VDW 0.5 10.0 \n')
                
                f.write('set sel [atomselect $molid ' + selection_string + '] \n')
                f.write('$sel set user '+ str((color_param_bp[i])) +' \n')
                f.write('$sel delete \n')

                f.write('color scale method GWR \n')
                f.write('mol modcolor $repindex $molid user \n')
                
                f.write(' \n')

            

            f.write(' \n')

def sum_all_offset_energys(offset_energys):  
    energy_sums_bp = [sum(i) for i in offset_energys['bp']]
    energy_sums_step = [(sum(offset_energys['step'][i]) + sum(offset_energys['heli'][i])) for i in range(len(offset_energys['step']))]
    offset_energys = energy_sums_bp[0] + energy_sums_step[0]
    return offset_energys

def find_energy_minimum_sequence(helices):  # for each helix and assuming the same average helix parameters as with old sequence find the sequence with minimal offset energy and write new seuence to mutation information
    complements = {
        'A' : 'T', 'T' : 'A', 'C' : 'G', 'G' : 'C'
    }
    with open('Mutate/mutation_information.txt', 'w') as f:
        for helix in helices:
            oldSeq = helix['strand_sequences'][0]
            possible_new_sequences = get_all_possible_sequences('-' * len(oldSeq))
            offset_energys_for_all_possible_sequences = []
            
            CONTENT_CHANGE_ALLOWED = 0
            old_AT_count = oldSeq.count('A') + oldSeq.count('T')
            for seq in possible_new_sequences[:]:
                new_AT_count = seq.count('A') + seq.count('T')
                if new_AT_count > old_AT_count+CONTENT_CHANGE_ALLOWED or new_AT_count < old_AT_count-CONTENT_CHANGE_ALLOWED:
                    possible_new_sequences.remove(seq)

            for seq in possible_new_sequences:
                offset_energys_for_all_possible_sequences.append(sum_all_offset_energys(calculate_displacement_energy_alternate_sequence(helix, seq)))

            new_seq = possible_new_sequences[offset_energys_for_all_possible_sequences.index(min(offset_energys_for_all_possible_sequences))] 
            old_Energy = calculate_displacement_energy(helix)[0]
            old_Energy = sum_all_offset_energys(old_Energy)
            new_Energy = offset_energys_for_all_possible_sequences[possible_new_sequences.index(new_seq)]
            
            for i, newRes in enumerate(new_seq):    # print to mutation information file in format <residue index> <new resname>
                if(newRes != oldSeq[i]):
                    f.write(str(helix['strand_res_inds'][0][i]+1) + " D" + newRes + '\n')
                    f.write(str(helix['strand_res_inds'][1][i]+1) + " D" + complements[newRes] + '\n')
            print(oldSeq + ' E='+str(old_Energy)+ " < old  |  new > "+ new_seq + ' E='+str(new_Energy)) 
            with open('energy_log.txt', 'a') as f2:
                f2.write(str(old_Energy) + "\n")

def get_helix_snippet(helix, start_bp, end_bp):
    if end_bp is None:
        end_bp = len(helix['strand_sequences'][0])
    helix_snippet = {
        'strand_sequences' : [ helix['strand_sequences'][0][start_bp:end_bp] , helix['strand_sequences'][1][start_bp:end_bp] ],
        'strand_res_inds' : [ helix['strand_res_inds'][0][start_bp:end_bp] , helix['strand_res_inds'][1][start_bp:end_bp] ],
        'bp_params' : helix['bp_params'][start_bp:end_bp],
        'step_params' : helix['step_params'][start_bp:end_bp-1],
        'heli_params' : helix['heli_params'][start_bp:end_bp-1],
        'energys' : None
    }
    return helix_snippet

def print_res_ind_of_Helix(helix):
    for res_ind in helix['strand_res_inds'][0]:
        print(res_ind, end=' ')
    for res_ind in helix['strand_res_inds'][1]:
        print(res_ind, end=' ')

def calculate_total_energy(helices):
    total_energy = 0
    for helix in helices:
        energy_sums_bp = [sum(i) for i in helix['energys']['bp']]
        energy_sums_step = [(sum(helix['energys']['step'][i]) + sum(helix['energys']['heli'][i])) for i in range(len(helix['energys']['step']))]
        total_energy += sum(energy_sums_bp) + sum(energy_sums_step)
    return total_energy


root_dir = Path(__file__).resolve().parents[1]
equalibrium_params = load_equalibrium_params()
stiffs = load_stiffs()

if __name__ == "__main__":
  
    # load in helices which were previously found by find_bound_double_strands.py
    helices = []
    md_results_dir = root_dir / "Offset_energy_calculator" / "MD_Results"
    for helix_file in md_results_dir.glob("MD_averaged_parameters_of_helix_*.dat"):
        helix = extract_data(helix_file)
        helix['energys'], helix['stiffs'], helix['eq_params'], helix['differences'] = calculate_displacement_energy(helix)
        helices.append(helix)

    print_helix_text_reprensentation(helices[0])
    
    total_energy = calculate_total_energy(helices)
    print("total Energy of Structure: " + str(total_energy))
    with open('energy_log.txt', 'a') as f2:
        f2.write(str(total_energy) + " ")

    helices = [get_helix_snippet(helices[0], 1, 16)]

    with open('energy_log.txt', 'a') as f2:
        f2.write(str(total_energy) + " ")
        
    helices = [get_helix_snippet(helices[0], 5, 12)]

    for helix in helices:
        helix['energys'], helix['stiffs'], helix['eq_params'], helix['differences'] = calculate_displacement_energy(helix)


    find_energy_minimum_sequence(helices)

    write_tcl_representation_script(helices)



