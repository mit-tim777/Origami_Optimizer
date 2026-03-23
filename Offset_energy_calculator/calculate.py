import csv
from pathlib import Path
import re
import numpy as np
import os, os.path
import random
import copy
from uncertainties import ufloat

root_dir = Path(__file__).resolve().parents[1] #/ "result_data" / "sequence_1" / "iteration_1"

csv_directory = Path(__file__).resolve().parents[1] / "Offset_energy_calculator" / "hexamers_csv" / "DNA"


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
            #heli_params.append(params_raw[13:19])
            bp_params_sd.append(params_raw[19:25])
            step_params_sd.append(params_raw[25:31])
            #heli_params_sd.append(params_raw[31:37])
        
        step_params.pop(-1)
        #heli_params.pop(-1)
        step_params_sd.pop(-1)
        #heli_params_sd.pop(-1)
        
        # Build ufloat containers upfront so that every downstream path uses the same uncertainty-aware values.
        bp_params_u = []
        for bp_row, bp_sd_row in zip(bp_params, bp_params_sd):
            bp_params_u.append([ufloat(v, sem_from_sd(sd)) for v, sd in zip(bp_row, bp_sd_row)])

        step_params_u = []
        for step_row, step_sd_row in zip(step_params, step_params_sd):
            step_params_u.append([ufloat(v, sem_from_sd(sd)) for v, sd in zip(step_row, step_sd_row)])

        helix = {
            'strand_sequences' : strands,
            'strand_res_inds' : strand_res_inds,
            'bp_params' : bp_params,
            'step_params' : step_params,
          #  'heli_params' : heli_params,
            'bp_params_sd' : bp_params_sd,
            'step_params_sd' : step_params_sd,
            'bp_params_u' : bp_params_u,
            'step_params_u' : step_params_u,
         #   'heli_params_sd' : heli_params_sd,
            'energys' : None,
            'stiffs' : None,
            'eq_params' : None
        }


    return helix

def safe_float(x):
    try:
        return float(x)
    except ValueError:
        return x 

MD_SAMPLE_SIZE = 1000  # number of frames used in the average; adjust if needed

def sem_from_sd(standard_dev, n=MD_SAMPLE_SIZE):
    if standard_dev is None:
        return 0.0
    return float(standard_dev) / np.sqrt(n)


def load_equalibrium_params(): # load the equalibrium parameters from the paper "Sequence-Dependent Shape and Stiffness of DNA and RNA Double Helices"
    with open(csv_directory / "coords_grooves_DNA_hexamers_table.csv") as f:
        data = csv.reader(f)
        equalibrium_step_params = { row[0] : [safe_float(i) for i in row[1:7]] for row in data }
    # with open(root_dir / "Offset_energy_calculator" / "hexamers_csv" / "DNA" / "coords_grooves_DNA_hexamers_table.csv") as f:
    #     data = csv.reader(f)
    #     equalibrium_heli_params = { row[0] : [safe_float(i) for i in row[7:13]] for row in data }
    with open(csv_directory / "coords_grooves_DNA_heptamers_table.csv") as f:
        data = csv.reader(f)
        equalibrium_bp_params = { row[0] : [safe_float(i) for i in row[1:7]] for row in data }
    equalibrium_params = {
        'bp' : equalibrium_bp_params,
        'step' : equalibrium_step_params,
        # 'heli' : equalibrium_heli_params
    }
    return equalibrium_params

def load_stiffs():  # load the quadratic offset energy stiffness
    with open(csv_directory / "coords_stiffs_DNA_hexamers_table.csv") as f:
        data = csv.reader(f)
        step_stiffs = { row[0] : [safe_float(i) for i in row[1:7]] for row in data }
    # with open(root_dir / "Offset_energy_calculator" / "hexamers_csv" / "DNA" / "coords_stiffs_DNA_hexamers_table.csv") as f:
    #     data = csv.reader(f)
    #     heli_stiffs = { row[0] : [safe_float(i) for i in row[7:13]] for row in data }
    with open(csv_directory / "coords_stiffs_DNA_heptamers_table.csv") as f:
        data = csv.reader(f)
        bp_stiffs = { row[0] : [safe_float(i) for i in row[1:7]] for row in data }
    stiffs = {
        'bp'   :   bp_stiffs,
        'step' : step_stiffs,
        # 'heli' : heli_stiffs
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
    differences = [[0]*6 for _ in range(len(helix['step_params']) )]
    stiffnesses = [[0]*6 for _ in range(len(helix['step_params']) )]
    energys = [[0]*6 for _ in range(len(helix['step_params']) )]
    equalibrium_params = [[0]*6 for _ in range(len(helix['step_params']) )]

    seq = helix['strand_sequences'][0]
    for step_n in range(0, len(helix['step_params'])):
        diffs_of_step = []
        energys_of_step = []
        hex_seq = ''.join([ seq[i] if (i in range(len(seq))) else '-' for i in range(step_n-2,step_n+4)])
        step_eql_params = get_equalibrium_params('step', hex_seq)
        step_stiffs = get_stiffs('step', hex_seq)
        for j in range(6):
            meas_u = helix['step_params_u'][step_n][j]
            eq_val = step_eql_params[j]
            diff = meas_u - eq_val
            energy = 0.5 * step_stiffs[j] * diff**2
            diffs_of_step.append(diff)
            energys_of_step.append(energy)
        energys[step_n] = energys_of_step
        differences[step_n] = diffs_of_step
        stiffnesses[step_n] = step_stiffs
        equalibrium_params[step_n] = step_eql_params
    return(energys, stiffnesses, equalibrium_params, differences)

def compare_bp_to_eql(helix):   # calculate the offset energy for all base pair parameters compared to equalibrium
    differences = [[0]*6 for _ in range(len(helix['bp_params']) )]
    stiffnesses = [[0]*6 for _ in range(len(helix['bp_params']) )]
    energys = [[0]*6 for _ in range(len(helix['bp_params']) )]
    equalibrium_params = [[0]*6 for _ in range(len(helix['bp_params']) )]

    seq = helix['strand_sequences'][0]
    for bp_n in range(0, len(helix['bp_params'])):
        diffs_of_bp = []
        energys_of_bp = []
        hep_seq = ''.join([ seq[i] if (i in range(len(seq))) else '-' for i in range(bp_n-3,bp_n+4)])
        bp_eql_params = get_equalibrium_params('bp',hep_seq)
        bp_stiffs = get_stiffs('bp',hep_seq)
        for j in range(6):
            meas_u = helix['bp_params_u'][bp_n][j]
            eq_val = bp_eql_params[j]
            diff = meas_u - eq_val
            energy = 0.5 * bp_stiffs[j] * diff**2
            diffs_of_bp.append(diff)
            energys_of_bp.append(energy)
        energys[bp_n] = energys_of_bp
        differences[bp_n] = diffs_of_bp
        stiffnesses[bp_n] = bp_stiffs
        equalibrium_params[bp_n] = bp_eql_params
    return(energys, stiffnesses, equalibrium_params, differences)


def calculate_displacement_energy(helix):  
    bp_res = compare_bp_to_eql(helix)
    step_res = compare_steps_to_eql(helix)

    energys = {
        'bp'   : bp_res[0],
        'step' : step_res[0], 
        # 'heli' : compare_heli_to_eql(helix)[0]
    }
    stiffs = {
        'bp'   : bp_res[1],
        'step' : step_res[1],
        # 'heli' : compare_heli_to_eql(helix)[1]
    }
    equalibrium_params = {
        'bp'   : bp_res[2],
        'step' : step_res[2], 
        # 'heli' : compare_heli_to_eql(helix)[2]
    }
    differences = {
        'bp'   : bp_res[3],
        'step' : step_res[3], 
        # 'heli' : compare_heli_to_eql(helix)[3]
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
        # 'heli' : compare_heli_to_eql(helix_cpy)[0]
    }
    return (energys)

def write_tcl_representation_script(helices):
    COLORING_FACTOR = 0.5
    with open('display_energys.tcl', 'w') as f:
        # Global Setup
        f.write('mol new output.pdb \n')
        f.write('set molid [lindex [expr {[molinfo list]}] end] \n')
        f.write('color scale method GWR \n')
        f.write('color scale min 0.3 \n')
        f.write('color scale midpoint 0.5 \n')
        f.write('color scale max 1 \n\n')

        for helix in helices:
            energy_sums_bp = [sum(i) for i in helix['energys']['bp']]
            energy_sums_bp = [float(e.nominal_value if hasattr(e, 'nominal_value') else e) for e in energy_sums_bp]
            color_param_bp = [i * COLORING_FACTOR for i in energy_sums_bp]

            energy_sums_step = [sum(helix['energys']['step'][i]) for i in range(len(helix['energys']['step']))]
            energy_sums_step = [float(e.nominal_value if hasattr(e, 'nominal_value') else e) for e in energy_sums_step]
            color_param_step = [i * COLORING_FACTOR for i in energy_sums_step]

            # 1. Create per-residue backbone representations
            for i in range(len(helix['strand_res_inds'][0]) - 1):
                selection_string = f'"backbone and residue {helix["strand_res_inds"][0][i]} {helix["strand_res_inds"][1][i]} and not name OP1 OP2"'
                
                f.write('mol addrep $molid \n')
                f.write('set repindex [expr {[molinfo $molid get numreps] - 1}] \n')
                f.write(f'mol modselect $repindex $molid {selection_string} \n')
                f.write(f'mol modstyle $repindex $molid Licorice 0.8 12.0 \n')
                f.write(f'mol modcolor $repindex $molid user \n')
                
                # --- LOCK THE COLOR RANGE ---
                f.write(f'mol colupdate $repindex $molid 0 \n')
                f.write(f'mol scaleminmax $molid $repindex 0.3 1.0 \n')
                
                f.write(f'set sel [atomselect $molid {selection_string}] \n')
                f.write(f'$sel set user {color_param_step[i]} \n')
                f.write(f'$sel delete \n\n')
            
            # 2. Create per-residue non-backbone representations
            for i in range(len(helix['strand_res_inds'][0])):
                selection_string = f'"not backbone and residue {helix["strand_res_inds"][0][i]} {helix["strand_res_inds"][1][i]} and not hydrogen"'
                
                f.write('mol addrep $molid \n')
                f.write('set repindex [expr {[molinfo $molid get numreps] - 1}] \n')
                f.write(f'mol modselect $repindex $molid {selection_string} \n')
                f.write(f'mol modstyle $repindex $molid VDW 0.5 10.0 \n')
                f.write(f'mol modcolor $repindex $molid user \n')

                # --- LOCK THE COLOR RANGE ---
                f.write(f'mol colupdate $repindex $molid 0 \n')
                f.write(f'mol scaleminmax $molid $repindex 0.3 1.0 \n')

                f.write(f'set sel [atomselect $molid {selection_string}] \n')
                f.write(f'$sel set user {color_param_bp[i]} \n')
                f.write(f'$sel delete \n\n')

            f.write('\n')

def sum_all_offset_energys(offset_energys):  
    # some callers may pass calculate_displacement_energy() tuple; take first item if needed.
    if not isinstance(offset_energys, dict) and len(offset_energys) > 0:
        offset_energys = offset_energys[0]

    energy_sums_bp = [sum(i) for i in offset_energys['bp']]
    energy_sums_step = [sum(offset_energys['step'][i]) for i in range(len(offset_energys['step']))] # + sum(offset_energys['heli'][i])
    return sum(energy_sums_bp) + sum(energy_sums_step)
     

def find_energy_minimum_sequence(helices, context_helices):  # for each helix and assuming the same average helix parameters as with old sequence find the sequence with minimal offset energy and write new seuence to mutation information
    complements = {
        'A' : 'T', 'T' : 'A', 'C' : 'G', 'G' : 'C'
    }
    with open('Mutate/mutation_information.txt', 'w') as f:
        for helix, context_helix in zip(helices, context_helices):
            oldSeq = helix['strand_sequences'][0]
            possible_new_sequences = get_all_possible_sequences('-' * len(oldSeq))
            offset_energys_for_all_possible_sequences = []
            
            CONTENT_CHANGE_ALLOWED = 0
            old_AT_count = oldSeq.count('A') + oldSeq.count('T')
            possible_new_sequences = [seq for seq in possible_new_sequences if abs(seq.count('A') + seq.count('T') - old_AT_count) <= CONTENT_CHANGE_ALLOWED]

            for seq in possible_new_sequences:
                seq_with_context = context_helix['strand_sequences'][0][:5] + seq + context_helix['strand_sequences'][0][-5:]
                offset_energys_for_all_possible_sequences.append(sum_all_offset_energys(calculate_displacement_energy_alternate_sequence(context_helix, seq_with_context)))

            best_index = min(range(len(offset_energys_for_all_possible_sequences)), key=lambda k: offset_energys_for_all_possible_sequences[k].nominal_value if hasattr(offset_energys_for_all_possible_sequences[k], 'nominal_value') else offset_energys_for_all_possible_sequences[k])
            new_seq = possible_new_sequences[best_index]

            old_Energy = sum_all_offset_energys(calculate_displacement_energy(context_helix))
            new_Energy = offset_energys_for_all_possible_sequences[best_index]
            
            # sequence_and_energy = list(zip(possible_new_sequences, offset_energys_for_all_possible_sequences))
            # sequence_and_energy.sort(key=lambda x: x[1])
            # print("All possible sequences and their energies (sorted):")
            # for seq, energy in sequence_and_energy:
            #     print(f"Sequence: {seq}, Energy: {energy}")

            for i, newRes in enumerate(new_seq):    # print to mutation information file in format <residue index> <new resname>
                if(newRes != oldSeq[i]):
                    f.write(str(helix['strand_res_inds'][0][i]+1) + " D" + newRes + '\n')
                    f.write(str(helix['strand_res_inds'][1][i]+1) + " D" + complements[newRes] + '\n')
            print(oldSeq + ' E='+str(old_Energy)+ " < old  |  new > "+ new_seq + ' E='+str(new_Energy)) 

            with open('energy_log.txt', 'a') as f2:
                f2.write(str(old_Energy) + " " + str(new_Energy) + "\n")

def get_helix_snippet(helix, start_bp, end_bp):
    if end_bp is None:
        end_bp = len(helix['strand_sequences'][0])
    if end_bp > 0:
        step_end_ind = end_bp-1
    else:
        step_end_ind = end_bp
    stiffs = {
            'bp' : helix['stiffs']['bp'][start_bp:end_bp],
            'step' : helix['stiffs']['step'][start_bp:step_end_ind],
        }
    equalibrium_params = {
            'bp' : helix['eq_params']['bp'][start_bp:end_bp],
            'step' : helix['eq_params']['step'][start_bp:step_end_ind],
        }
    helix_snippet = {
        'strand_sequences' : [ helix['strand_sequences'][0][start_bp:end_bp] , helix['strand_sequences'][1][start_bp:end_bp] ],
        'strand_res_inds' : [ helix['strand_res_inds'][0][start_bp:end_bp] , helix['strand_res_inds'][1][start_bp:end_bp] ],
        'bp_params' : helix['bp_params'][start_bp:end_bp],
        'step_params' : helix['step_params'][start_bp:step_end_ind],
        'bp_params_sd' : helix['bp_params_sd'][start_bp:end_bp],
        'step_params_sd' : helix['step_params_sd'][start_bp:step_end_ind],
        'bp_params_u' : helix['bp_params_u'][start_bp:end_bp],
        'step_params_u' : helix['step_params_u'][start_bp:step_end_ind],
        'stiffs' : stiffs, 
        'eq_params' : equalibrium_params,
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
        if helix['energys'] is None:
            helix['energys'], helix['stiffs'], helix['eq_params'], helix['differences'] = calculate_displacement_energy(helix)
        energy_sums_bp = [sum(i) for i in helix['energys']['bp']]
        energy_sums_step = [sum(helix['energys']['step'][i]) for i in range(len(helix['energys']['step']))] #+ sum(helix['energys']['heli'][i])
        helix_energy = sum(energy_sums_bp) + sum(energy_sums_step)
        total_energy += helix_energy
    return total_energy

def construct_stiffness_matrix(sequence):
    matrix_size = (6*len(sequence) + 6 * (len(sequence)-1))  # 6 parameters per base pair and 6 parameters per step, with one less step than base pairs
    K = np.zeros((matrix_size, matrix_size))
    K_layered = [[[] for _ in range(matrix_size)] for _ in range(matrix_size)]
    hexamers = [ sequence[i:i+6] for i in range(len(sequence)-5) ]

    with open(csv_directory / "K_intra_inter_DNA.csv") as Kf:
        data_raw = csv.reader(Kf)
        data = [row for row in data_raw]
        for n, hexamer in enumerate(hexamers): 
            # get all lines that start with hexamer and read the 66x66 block of stiffness values that follow
            K_block = []

            for row in data:
                if row[0] == hexamer:
                    K_block.append(list(map(float, row[1:67])))
            for i, row in enumerate(K_block):
                for j, value in enumerate(row):
                    K_layered[12*n+i][12*n+j].append(value)

    # Now average the values in K_layered and fill the final stiffness matrix K
    for i in range(matrix_size):
        for j in range(matrix_size):
            if K_layered[i][j]:  # if there are values to average
                K[i][j] = sum(K_layered[i][j]) / len(K_layered[i][j])
            else:
                K[i][j] = 0.0  # or some default value if no data is available
    
    P = get_projection_matrix(len(sequence))

    w_eq = np.zeros((matrix_size))  
    C = np.zeros((matrix_size))
    for step_n in range(len(sequence)-1):
        hex_seq = ''.join([ sequence[i] if (i in range(len(sequence))) else '-' for i in range(step_n-2,step_n+4)])
        w_eq[6+12*step_n : 12*(step_n+1)] = get_equalibrium_params('step', hex_seq)
    
    for bp_n in range(len(sequence)):
        hep_seq = ''.join([ sequence[i] if (i in range(len(sequence))) else '-' for i in range(bp_n-3,bp_n+4)])
        w_eq[12*bp_n : 6+12*bp_n] = get_equalibrium_params('bp',hep_seq)
    



    # for j, base in enumerate(sequence[:-1]):
    #     hex_seq = ''.join([ sequence[i] if (i in range(len(sequence))) else '-' for i in range(j-2,j+4)])
    #     hep_seq = ''.join([ sequence[i] if (i in range(len(sequence))) else '-' for i in range(j-3,j+4)])
    #     print( "position " + str(j) + " base " + base + " hexamer seq: " + hex_seq + " heptamer seq: " + hep_seq)
    #     for i, parmtype in enumerate(["shear","stretch","stagger","buckle","propeller","opening","shift","slide","rise","tilt","roll","twist"]):
    #         print(f"{parmtype}: {w_eq[i+12*j]}")

    # with open('projection_matrix.txt', 'w') as f:
    #     for y in range(len(P[0])):
    #         row = []
    #         for x in range(len(P)):
    #             row.append(str(int(P[x][y])))
    #         f.write(''.join(row) + '\n')
    return K

def get_projection_matrix(sequence_length):
    # This function should return the projection matrix P that maps the full parameter space to the sub
    P = np.zeros((12*sequence_length-12, 12*sequence_length-6))
    for i in range(12*sequence_length-18):
        P[i][i] = 1
    for i in range(sequence_length-2):
        for j in range(6):
            P[i*12+6+j][12*sequence_length-18+j] = -1
    for i in range(6):
        P[12*sequence_length-18+i][12*(sequence_length-1)+i] = 1
    
    return P


equalibrium_params = load_equalibrium_params()
stiffs = load_stiffs()

if __name__ == "__main__":
  
    # # load in helices which were previously found by find_bound_double_strands.py
    # helices = []
    # md_results_dir = root_dir / "Offset_energy_calculator" / "MD_Results"
    # for helix_file in md_results_dir.glob("MD_averaged_parameters_of_helix_*.dat"):
    #     helix = extract_data(helix_file)
    #     helix['energys'], helix['stiffs'], helix['eq_params'], helix['differences'] = calculate_displacement_energy(helix)
    #     helices.append(helix)


    # total_energy = calculate_total_energy(helices)

    # with open('energy_log.txt', 'a') as f2:
    #     f2.write(str(total_energy) + " ")

    # helices_temp = [get_helix_snippet(helices[0], 1, -1)]
    # total_energy = calculate_total_energy(helices_temp)
    
    # with open('energy_log.txt', 'a') as f2:     # emergy_log line structure: <total energy of whole helix> <total energy of helix without terminal base pairs> <total energy of central hexamer> <total energy of central hexamer with alternate sequence with minimal energy (unsimulated)>
    #     f2.write(str(total_energy) + " ")


    # helix_snip = [get_helix_snippet(helices[0], 5, -5)]

    # for helix in helix_snip:
    #     helix['energys'], helix['stiffs'], helix['eq_params'], helix['differences'] = calculate_displacement_energy(helix)


    # find_energy_minimum_sequence( helix_snip,helices )

    # write_tcl_representation_script(helices)

    construct_stiffness_matrix("ATCGGTAGCTGAGC")



