import csv
import numpy as np
import os, os.path
import random
import copy

def print_helix_data(helix):
    
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

    # print("base pair parameters:")
    # print("shear \tstretch stagger buckle propeller opening")

    # for bp_params in helix['bp_params']:
    #     for param in bp_params:
    #         print(param, end='\t')
    #     print()
    
    # print()

    # print("step paramenters:")
    # print("shift \t slide \trise \ttilt \troll \ttwist")
    # for step_params in helix['step_params']:
    #     for param in step_params:
    #         print(param, end='\t')
    #     print()
    # print()

    # print("helical parameters:")
    # print("x_disp	y_disp	hrise	incl	tip	htwist")
    # for heli_params in helix['heli_params']:
    #     for param in heli_params:
    #         print(param, end='\t')
    #     print()
    # print()


def extract_data(filename):

    with open(filename) as f:
        
        bp_params = []
        step_params = []
        heli_params = []
        strands = []
        
        meta_data = f.readline().split()
        strands.append(meta_data[1])
        strands.append(meta_data[2])
        
        strand_res_inds = [[int(i.split('-')[j])-1 for i in meta_data[0].split(',')[0:-1]] for j in range(2)]

        next(f)
        for line in f:
            params_raw = [ float(i) for i in line.split()]
            bp_params.append(params_raw[1:7])
            step_params.append(params_raw[7:13])
            heli_params.append(params_raw[13:21])
        
        step_params.pop(-1)
        heli_params.pop(-1)

        
        helix = {
            'strand_sequences' : strands,
            'strand_res_inds' : strand_res_inds,
            'step_params' : step_params,
            'bp_params' : bp_params,
            'heli_params' : heli_params,
            'energys' : None
        }


    return helix

def safe_float(x):
    try:
        return float(x)
    except ValueError:
        return x 

def load_equalibrium_params():
    with open(root+'hexamers_csv/DNA/coords_grooves_DNA_hexamers_table.csv') as f:
        data = csv.reader(f)
        equalibrium_step_params = { row[0] : [safe_float(i) for i in row[1:7]] for row in data }
    with open(root+'hexamers_csv/DNA/coords_grooves_DNA_hexamers_table.csv') as f:
        data = csv.reader(f)
        equalibrium_heli_params = { row[0] : [safe_float(i) for i in row[7:13]] for row in data }
    with open(root+'hexamers_csv/DNA/coords_grooves_DNA_heptamers_table.csv') as f:
        data = csv.reader(f)
        equalibrium_bp_params = { row[0] : [safe_float(i) for i in row[1:7]] for row in data }
    equalibrium_params = {
        'bp' : equalibrium_bp_params,
        'step' : equalibrium_step_params,
        'heli' : equalibrium_heli_params
    }
    return equalibrium_params

def load_stiffs():
    with open(root+'hexamers_csv/DNA/coords_stiffs_DNA_hexamers_table.csv') as f:
        data = csv.reader(f)
        step_stiffs = { row[0] : [safe_float(i) for i in row[1:7]] for row in data }
    with open(root+'hexamers_csv/DNA/coords_stiffs_DNA_hexamers_table.csv') as f:
        data = csv.reader(f)
        heli_stiffs = { row[0] : [safe_float(i) for i in row[7:13]] for row in data }
    with open(root+'hexamers_csv/DNA/coords_stiffs_DNA_heptamers_table.csv') as f:
        data = csv.reader(f)
        bp_stiffs = { row[0] : [safe_float(i) for i in row[1:7]] for row in data }
    stiffs = {
        'bp'   :   bp_stiffs,
        'step' : step_stiffs,
        'heli' : heli_stiffs
    }
    return stiffs

def get_all_possible_sequences(sequence):
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

def get_equalibrium_params(param_type, sequence):
    
    possible_sequences = get_all_possible_sequences(sequence)

    #calculate average paramenters for all possible sequences
    avg_params = [0]*6
    for seq in possible_sequences:
        for i in range(6):
            avg_params[i] += equalibrium_params[param_type][seq][i]

    for i in range(6):
        avg_params[i] /= len(possible_sequences)

    return avg_params

def get_stiffs(param_type, sequence):
    possible_sequences = get_all_possible_sequences(sequence)

        #calculate average paramenters for all possible sequences
    avg_params = [0]*6
    for seq in possible_sequences:
        for i in range(6):
            avg_params[i] += stiffs[param_type][seq][i]

    for i in range(6):
        avg_params[i] /= len(possible_sequences)

    return avg_params
  
def compare_steps_to_eql(helix):
    differences = [[0,0,0,0,0,0]] * len(helix['step_params'])

    seq = helix['strand_sequences'][0]
    for step_n in range(0, len(helix['step_params'])):
        diffs_of_step = []
        hex_seq = ''.join([ seq[i] if (i in range(len(seq))) else '-' for i in range(step_n-2,step_n+4)])
        step_eql_params = get_equalibrium_params('step', hex_seq)
        step_stiffs = get_stiffs('step', hex_seq)
        for j in range(6):
            diff = 0.5 * step_stiffs[j] * pow(helix['step_params'][step_n][j] - step_eql_params[j],2)
            diffs_of_step.append(diff)
        differences[step_n] = diffs_of_step
    return(differences)

def compare_heli_to_eql(helix):
    differences = [[0,0,0,0,0,0]] * len(helix['heli_params'])

    seq = helix['strand_sequences'][0]
    for step_n in range(0, len(helix['heli_params'])):
        diffs_of_step = []
        hex_seq = ''.join([ seq[i] if (i in range(len(seq))) else '-' for i in range(step_n-2,step_n+4)])
        step_eql_params = get_equalibrium_params('heli', hex_seq)
        step_stiffs = get_stiffs('heli', hex_seq)
        for j in range(6):
            diff = 0.5 * step_stiffs[j] * pow(helix['heli_params'][step_n][j] - step_eql_params[j],2)
            diffs_of_step.append(diff)
        differences[step_n] = diffs_of_step
    return(differences)

def compare_bp_to_eql(helix):
    differences = [[0,0,0,0,0,0]] * len(helix['bp_params'])

    seq = helix['strand_sequences'][0]
    for bp_n in range(0, len(helix['bp_params'])):
        diffs_of_bp = []
        hep_seq = ''.join([ seq[i] if (i in range(len(seq))) else '-' for i in range(bp_n-3,bp_n+4)])
        eql_params = get_equalibrium_params('bp',hep_seq)
        bp_stiffs = get_stiffs('bp',hep_seq)
        for j in range(6):
            diff = 0.5 * bp_stiffs[j] * pow(helix['bp_params'][bp_n][j] - eql_params[j], 2)
            diffs_of_bp.append(diff)
        differences[bp_n] = diffs_of_bp
    return(differences)

def calculate_displacement_energy(helix):
    coordinate_displacements = {
        'bp'   : compare_bp_to_eql(helix),
        'step' : compare_steps_to_eql(helix), 
        'heli' : compare_heli_to_eql(helix)
    }
    return coordinate_displacements

def calculate_displacement_energy2(helix, sequence):
    complements = {
        'A' : 'T', 'T' : 'A', 'C' : 'G', 'G' : 'C'
    }
    complementary = "".join([complements[i] for i in sequence])
    helix_cpy = copy.deepcopy(helix)
    helix_cpy['strand_sequences'] = [ "".join(sequence) , complementary]
    coordinate_displacements = {
        'bp'   : compare_bp_to_eql(helix_cpy),
        'step' : compare_steps_to_eql(helix_cpy), 
        'heli' : compare_heli_to_eql(helix_cpy)
    }
    return coordinate_displacements


def write_tcl_representation_script(helices):

    with open('display_energys.tcl', 'w') as f:
        #f.write('cd .. \n')
        #f.write('cd Helix_separator \n')
        f.write('mol new output.pdb \n')
        f.write('set molid [lindex [expr {[molinfo list]}] end] \n\n\n')
        
        for helix in helices:

            energy_sums = [np.sqrt(sum(i)) for i in helix['energys']['bp']]
            color_param_bp = [i/max(energy_sums) for i in energy_sums]

            energy_sums = [np.sqrt(sum(helix['energys']['step'][i]) + sum(helix['energys']['heli'][i])) for i in range(len(helix['energys']['step']))]
            color_param_step = [i/max(energy_sums) for i in energy_sums]

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

# for a MD averaged origami replacing the second last bp of each helix so the eqlb energy of the last bp is minimized
def write_mutation_information(helices):
    with open('Mutate/mutation_information.txt', 'w') as f:
        for helix in helices:

            possible_bps = [['A','T'],['T','A'],['C','G'],['G','C']]

            energy_sums_with_mutation = []
            bp_to_change = 1 # random.randint(1, len(helix['strand_sequences'][0])-2)
            
            placeholder_seq = list(helix['strand_sequences'][0])
            placeholder_seq[1] = '-'
            possible_mutated_sequences = get_all_possible_sequences(placeholder_seq)
            for seq in possible_mutated_sequences:
                offset_mutation = calculate_displacement_energy2(helix, seq)
                energy_sums_bp = [sum(i) for i in offset_mutation['bp']]
                energy_sums_step = [(sum(offset_mutation['step'][i]) + sum(offset_mutation['heli'][i])) for i in range(len(offset_mutation['step']))]
                energy_sums_with_mutation.append(energy_sums_bp[0] + energy_sums_step[0])
                # print( seq + " " + str(energy_sums_bp[0] + energy_sums_step[0]))
            #print(str(energy_sums_with_mutation) + "   " + helix['strand_sequences'][0][bp_to_change] + "-" + helix['strand_sequences'][1][bp_to_change])
            # print(helix['strand_sequences'][0])
            new_bp = possible_bps[energy_sums_with_mutation.index(min(energy_sums_with_mutation))] # [energy_sums_with_mutation.index(min(energy_sums_with_mutation))]
            
            old_bp = [helix['strand_sequences'][0][1],helix['strand_sequences'][1][1]]

            if new_bp == old_bp:
                return
            # for new_bp in possible_bps:
                # mutated_helix = copy.deepcopy(helix)
                # new_sequenceFow = list(copy.copy(mutated_helix['strand_sequences'][0]))
                # new_sequenceBac = list(copy.copy(mutated_helix['strand_sequences'][1]))
                # new_sequenceFow[bp_to_change] = new_bp[0]
                # new_sequenceBac[bp_to_change] = new_bp[1]
                # mutated_helix['strand_sequences'] = ["".join(new_sequenceFow) , "".join(new_sequenceBac)]
                # mutated_helix['energys'] = calculate_displacement_energy(mutated_helix)
                # energy_sums_bp = [sum(i) for i in mutated_helix['energys']['bp']]
                # energy_sums_step = [(sum(mutated_helix['energys']['step'][i]) + sum(mutated_helix['energys']['heli'][i])) for i in range(len(mutated_helix['energys']['step']))]
                # energy_sums_with_mutation.append(sum(energy_sums_bp) + sum(energy_sums_step))
            
            #print(str(energy_sums_with_mutation) + "   " + helix['strand_sequences'][0][bp_to_change] + "-" + helix['strand_sequences'][1][bp_to_change])
 
 
 
 
 

            # die bisherige sequenz ist immer am nächsten am equilibrium mit der bisherigen sequenz :/
 
 
 
 
 
 
            # energy_sums_bp[0] = 0
            # energy_sums_bp[-1] = 0
            # print(str(helix['strand_sequences'][0]) + "  :  Etot = " + str(sum(energy_sums_bp)))
            # bp_to_change = [helix['strand_res_inds'][i][energy_sums_bp.index(max(energy_sums_bp))]+1 for i in range(2)]
            # new_bp_names = possible_bps[random.randint(0,3)]
            
            f.write(str(helix['strand_res_inds'][0][bp_to_change]) + " " + new_bp[0] + '\n')
            f.write(str(helix['strand_res_inds'][1][bp_to_change]) + " " + new_bp[1] + '\n')

def find_energy_minimum_sequence(helices):
    complements = {
        'A' : 'T', 'T' : 'A', 'C' : 'G', 'G' : 'C'
    }
    with open('Mutate/mutation_information.txt', 'w') as f:
        for helix in helices:
            possible_bps = [['A','T'],['T','A'],['C','G'],['G','C']]

            energy_sums_with_mutation = []

            placeholder_seq = '-' * len(helix['strand_sequences'][0])
            possible_mutated_sequences = get_all_possible_sequences(placeholder_seq)
            for seq in possible_mutated_sequences:
                offset_mutation = calculate_displacement_energy2(helix, seq)
                energy_sums_bp = [sum(i) for i in offset_mutation['bp']]
                energy_sums_step = [(sum(offset_mutation['step'][i]) + sum(offset_mutation['heli'][i])) for i in range(len(offset_mutation['step']))]
                energy_sums_with_mutation.append(energy_sums_bp[0] + energy_sums_step[0])


            new_seq = possible_mutated_sequences[energy_sums_with_mutation.index(min(energy_sums_with_mutation))] # [energy_sums_with_mutation.index(min(energy_sums_with_mutation))]
            old_Energy = energy_sums_with_mutation[possible_mutated_sequences.index(helix['strand_sequences'][0])]
            new_Energy = energy_sums_with_mutation[possible_mutated_sequences.index(new_seq)]
            
            for i, newRes in enumerate(new_seq):
                if(newRes != helix['strand_sequences'][0][i]):
                    f.write(str(helix['strand_res_inds'][0][i]) + " D" + newRes + '\n')
                    f.write(str(helix['strand_res_inds'][1][i]) + " D" + complements[newRes] + '\n')
            print(helix['strand_sequences'][0] + ' E='+str(old_Energy)+ " < old  |  new > "+ new_seq + ' E='+str(new_Energy)) 

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


root = './Offset_energy_calculator/'
# root = './'
if __name__ == "__main__":
    equalibrium_params = load_equalibrium_params()
    stiffs = load_stiffs()
    
    
    helices = []
    for i in range(len([name for name in os.listdir(root+'MD_Results')])):
        helix = extract_data(root+'MD_Results/MD_averaged_parameters_of_helix_'+str(i)+'.dat')
        helix['energys'] = calculate_displacement_energy(helix)
        helices.append(helix)

    print("total Energy of Structure: " + str(calculate_total_energy(helices)))

    hexamers_next_to_transition = []
    hexamer_len = 6

    for helix in helices:       
        starting_hexamer = get_helix_snippet(helix, 0, hexamer_len)
        hexamers_next_to_transition.append(starting_hexamer)

        if(len(helix['strand_sequences'][0]) > hexamer_len):
            ending_hexamer = get_helix_snippet(helix, -hexamer_len, None)
            hexamers_next_to_transition.append(ending_hexamer)

    energys = []
    for helix in hexamers_next_to_transition:
        helix['energys'] = calculate_displacement_energy(helix)
        print_helix_data(helix)

# only for testing, outputs the residue indices so that you can paste it in vmd selection
    # for helix in hexamers_next_to_transition:
    #     print_res_ind_of_Helix(helix)

    find_energy_minimum_sequence(hexamers_next_to_transition)
    # write_mutation_information(helices)
    
    write_tcl_representation_script(helices)

