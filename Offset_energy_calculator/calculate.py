import csv
from pathlib import Path
import numpy as np
import copy

import scipy
from scipy import sparse
from scipy.optimize import NonlinearConstraint
# from uncertainties import ufloat

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
        
        step_params.pop(-1)

        helix = {
            'strand_sequences' : strands,
            'strand_res_inds' : strand_res_inds,
            'w' : build_w_vector(bp_params, step_params),
            'wKw' : None,
            'w_eq' : None,
            'K' : None
        }


    return helix

def safe_float(x):
    try:
        return float(x)
    except ValueError:
        return x 

def build_w_vector(bp_params, step_params):
    sequence_length = len(bp_params)
    if step_params and len(step_params) != sequence_length - 1:
        raise ValueError("Step parameters must contain one fewer entry than base-pair parameters.")

    matrix_size = 6 * sequence_length + 6 * max(0, sequence_length - 1)
    w = np.zeros(matrix_size)

    for bp_n in range(sequence_length):
        w[12 * bp_n : 6 + 12 * bp_n] = bp_params[bp_n]

    for step_n in range(max(0, sequence_length - 1)):
        w[6 + 12 * step_n : 12 * (step_n + 1)] = step_params[step_n]

    return w


def slice_w_vector(w_vector, start_bp, end_bp):
    if end_bp < 0:
        end_bp = ( len(w_vector) // 12 + 1 ) + end_bp

    snippet_length = end_bp - start_bp

    matrix_size = 6 * snippet_length + 6 * (snippet_length - 1)
    snippet_w = np.zeros(matrix_size)

    for bp_n in range(snippet_length):
        global_bp_n = start_bp + bp_n
        snippet_w[12 * bp_n : 6 + 12 * bp_n] = w_vector[12 * global_bp_n : 6 + 12 * global_bp_n]

    for step_n in range(max(0, snippet_length - 1)):
        global_step_n = start_bp + step_n
        snippet_w[6 + 12 * step_n : 12 * (step_n + 1)] = w_vector[6 + 12 * global_step_n : 12 * (global_step_n + 1)]

    return snippet_w

def build_equilibrium_w_vector(sequence):
    matrix_size = 6 * len(sequence) + 6 * max(0, len(sequence) - 1)
    w_eq = np.zeros(matrix_size)

    for bp_n in range(len(sequence)):
        hep_seq = ''.join([sequence[i] if (i in range(len(sequence))) else '-' for i in range(bp_n - 3, bp_n + 4)])
        w_eq[12 * bp_n : 6 + 12 * bp_n] = get_equalibrium_params('bp', hep_seq)

    for step_n in range(len(sequence) - 1):
        hex_seq = ''.join([sequence[i] if (i in range(len(sequence))) else '-' for i in range(step_n - 2, step_n + 4)])
        w_eq[6 + 12 * step_n : 12 * (step_n + 1)] = get_equalibrium_params('step', hex_seq)

    return w_eq

MD_SAMPLE_SIZE = 1000  # number of frames used in the average; adjust if needed

def sem_from_sd(standard_dev, n=MD_SAMPLE_SIZE):
    if standard_dev is None:
        return 0.0
    return float(standard_dev) / np.sqrt(n)

def load_equalibrium_params(): # load the equalibrium parameters from the paper "Sequence-Dependent Shape and Stiffness of DNA and RNA Double Helices"
    with open(csv_directory / "coords_grooves_DNA_hexamers_table.csv") as f:
        data = csv.reader(f)
        equalibrium_step_params = { row[0] : [safe_float(i) for i in row[1:7]] for row in data }
    with open(csv_directory / "coords_grooves_DNA_heptamers_table.csv") as f:
        data = csv.reader(f)
        equalibrium_bp_params = { row[0] : [safe_float(i) for i in row[1:7]] for row in data }
    equalibrium_params = {
        'bp' : equalibrium_bp_params,
        'step' : equalibrium_step_params,
    }
    return equalibrium_params

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
  
complements = {
    'A' : 'T', 'T' : 'A', 'C' : 'G', 'G' : 'C'
}
def find_energy_minimum_sequence(helices ):  # for each helix and assuming the same average helix parameters as with old sequence find the sequence with minimal offset energy and write new seuence to mutation information
    with open('Mutate/mutation_information.txt', 'w') as f:
        # for helix, context_helix in zip(helices, context_helices):
        for helix in helices:
            oldSeq = helix['strand_sequences'][0]
            possible_new_sequences = get_all_possible_sequences('-' * len(oldSeq))
            offset_energys_for_all_possible_sequences = []
            
            CONTENT_CHANGE_ALLOWED = 0
            old_AT_count = oldSeq.count('A') + oldSeq.count('T')
            possible_new_sequences = [seq for seq in possible_new_sequences if abs(seq.count('A') + seq.count('T') - old_AT_count) <= CONTENT_CHANGE_ALLOWED]

            for seq in possible_new_sequences:
                # seq_with_context = context_helix['strand_sequences'][0][:5] + seq + context_helix['strand_sequences'][0][-5:]
                offset_energys_for_all_possible_sequences.append(calculate_wKw(helix, sequence=seq))

            best_index = min(range(len(offset_energys_for_all_possible_sequences)), key=lambda k: offset_energys_for_all_possible_sequences[k].nominal_value if hasattr(offset_energys_for_all_possible_sequences[k], 'nominal_value') else offset_energys_for_all_possible_sequences[k])

            new_seq = possible_new_sequences[best_index]

            old_Energy = calculate_wKw(helix)
            new_Energy = offset_energys_for_all_possible_sequences[best_index]
            
            sequence_and_energy = list(zip(possible_new_sequences, offset_energys_for_all_possible_sequences))
            sequence_and_energy.sort(key=lambda x: x[1])
            # print("All possible sequences and their energies (sorted):")
            # for seq, energy in sequence_and_energy:
            #     print(f"Sequence: {seq}, Energy: {energy}")

            for i, newRes in enumerate(new_seq):    # print to mutation information file in format <residue index> <new resname>
                if(newRes != oldSeq[i]):
                    f.write(str(helix['strand_res_inds'][0][i]+1) + " D" + newRes + '\n')
                    f.write(str(helix['strand_res_inds'][1][i]+1) + " D" + complements[newRes] + '\n')
            print(oldSeq + ' E='+str(old_Energy)+ " < old  |  new > "+ new_seq + ' E='+str(new_Energy)) 

            # with open('energy_log.txt', 'a') as f2:
            #     f2.write(str(old_Energy) + " " + str(new_Energy) + "\n")

def get_helix_snippet(helix, start_bp, end_bp):

    snippet_sequence = helix['strand_sequences'][0][start_bp:end_bp]
    if(len(snippet_sequence) < 6):
        print("Warning: snippet sequence is less than 6 base pairs long")
        return None
    
    snippet_w = slice_w_vector(helix['w'], start_bp, end_bp)
    snippet_w_eq = slice_w_vector(helix['w_eq'], start_bp, end_bp)

    helix_snippet = {
        'strand_sequences' : [ snippet_sequence, helix['strand_sequences'][1][start_bp:end_bp] ],
        'strand_res_inds' : [ helix['strand_res_inds'][0][start_bp:end_bp] , helix['strand_res_inds'][1][start_bp:end_bp] ],
        'w' : snippet_w,
        'w_eq' : snippet_w_eq,
        'K' : None,
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
        total_energy += calculate_wKw(helix)
    return total_energy

def calculate_wKw(helix):  # calculate the quadratic energy w^T K w for a helix, where w is the vector of displacements from equalibrium and K is the stiffness matrix; if sequence is provided, use that sequence instead of the one in the helix to determine K and equalibrium parameters (useful for calculating energy of alternate sequences without needing to recalculate displacements)
    
    sequence = helix['strand_sequences'][0]
    matrix_size = (6*len(sequence) + 6 * (len(sequence)-1))  # 6 parameters per base pair and 6 parameters per step, with one less step than base pairs
    hexamers = [ sequence[i:i+6] for i in range(len(sequence)-5) ]

    block_size = 66
    block_rows = np.arange(block_size)
    rows = []
    cols = []
    values = []
    counts = []

    for block_index, hexamer in enumerate(hexamers):
        # Assemble overlapping 66x66 blocks directly as sparse COO entries.
        K_block = np.asarray(K_intra_inter_DNA[hexamer], dtype=float)
        if K_block.shape != (block_size, block_size):
            raise ValueError(
                f"Expected a 66x66 stiffness block for {hexamer}, "
                f"got {K_block.shape}"
            )

        row_indices, col_indices = np.meshgrid(
            12 * block_index + block_rows,
            12 * block_index + block_rows,
            indexing='ij',
        )
        rows.append(row_indices.ravel())
        cols.append(col_indices.ravel())
        values.append(K_block.ravel())
        counts.append(np.ones(block_size * block_size))

    if rows:
        row_indices = np.concatenate(rows)
        col_indices = np.concatenate(cols)
        K_sum = sparse.coo_matrix(
            (np.concatenate(values), (row_indices, col_indices)),
            shape=(matrix_size, matrix_size),
        ).tocsr()
        K_count = sparse.coo_matrix(
            (np.concatenate(counts), (row_indices, col_indices)),
            shape=(matrix_size, matrix_size),
        ).tocsr()
        K = K_sum.multiply(K_count.power(-1))
    else:
        K = sparse.csr_matrix((matrix_size, matrix_size))

    cutoff = 0.44
    K_dense = K.toarray()
    eigenvalues, eigenvectors = np.linalg.eigh(K_dense)
    corrected_eigenvalues = np.maximum(eigenvalues, cutoff)
    K_corrected = (eigenvectors * corrected_eigenvalues) @ eigenvectors.T

    # The eigendecomposition is dense; restore the half-bandwidth of the
    # overlapping 66-parameter blocks before converting back to sparse form.
    band_mask = np.abs(
        np.arange(matrix_size)[:, None] - np.arange(matrix_size)[None, :]
    ) < block_size
    K_corrected[~band_mask] = 0.0
    K = sparse.csr_matrix((K_corrected + K_corrected.T) / 2)

    helix['w_eq'] = build_equilibrium_w_vector(sequence)

    w_offset = helix['w'] - helix['w_eq']
    energy = float(w_offset @ K.dot(w_offset))
    helix['wKw'] = energy
    helix['K'] = K
    return energy

def get_parms_of_type(param_type, w):
    parnames = ['shift', 'slide', 'rise', 'tilt', 'roll', 'twist']
    for i in range(len(w)//12):
        yield w[6+12*i + parnames.index(param_type)]

def rodrigues_matrix(axis, angle):
    """Generates a 3x3 rotation matrix using Rodrigues' formula."""
    axis = np.array(axis)
    if np.linalg.norm(axis) < 1e-9:
        return np.eye(3)
    axis = axis / np.linalg.norm(axis)
    
    # Skew-symmetric matrix K
    K = np.array([
        [0, -axis[2], axis[1]],
        [axis[2], 0, -axis[0]],
        [-axis[1], axis[0], 0]
    ])
    
    return np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * np.dot(K, K)

def get_rotation_matrix(step_parms):

    tilt, roll, twist = np.deg2rad(step_parms[3:6])
    L = np.sqrt(tilt**2 + roll**2)
    o = np.arctan2(tilt, roll)
    
    # Axis of tilt in the T1 xy-plane (offset by o - w/2)
    theta_axis = o - twist/2
    u = np.array([np.sin(theta_axis), np.cos(theta_axis), 0.0])
    
    # Full Rotation R = R_tilt(L) * R_twist(w)
    R_tilt = rodrigues_matrix(u, L)
    R_twist = rodrigues_matrix([0, 0, 1], twist)
    R_full = np.dot(R_tilt, R_twist)
    
    return R_full

def get_half_rotation_matrix(step_parms):
    tilt, roll, twist = np.deg2rad(step_parms[3:6])
    L = np.sqrt(tilt**2 + roll**2)
    o = np.arctan2(tilt, roll)

    theta_axis = o - twist/2
    u = np.array([np.sin(theta_axis), np.cos(theta_axis), 0.0])
    
    R_half_tilt = rodrigues_matrix(u, L/2)
    R_half_twist = rodrigues_matrix([0, 0, 1], twist/2)
    R_half = np.dot(R_half_tilt, R_half_twist)

    return R_half
    
def get_step_matrix(step_parms):
    # Assemble 4x4 Matrix

    M = np.eye(4)
    M[:3, :3] = get_rotation_matrix(step_parms)  
    M[:3, 3] = np.dot(get_half_rotation_matrix(step_parms), np.array(step_parms[:3]))  # Apply half rotation to the local displacement
    return M

def mutate_bp(helix, bp_index, new_bp_type):

    helix = copy.deepcopy(helix)  # create a copy of the helix to avoid modifying the original

    sequence = list(helix['strand_sequences'][0])
    sequence[bp_index] = new_bp_type
    helix['strand_sequences'][0] = ''.join(sequence)
    helix['strand_sequences'][1] = ''.join([complements[base] for base in sequence])

    calculate_wKw(helix)  # recalculate wKw and K for the new sequence

    # Optimize step i-1, base pair i, and step i.
    B = np.arange(12 * bp_index - 6, 12 * bp_index + 12)
    first_step = np.arange(12 * bp_index - 6, 12 * bp_index)
    intra_bp = np.arange(12 * bp_index, 12 * bp_index + 6)
    second_step = np.arange(12 * bp_index + 6, 12 * bp_index + 12)

    K = helix['K']
    w = helix['w'] - helix['w_eq']

    def objective(local_w):
        trial_w = w.copy()
        trial_w[B] = local_w
        return float(trial_w @ K.dot(trial_w))

    original_step_parameters = [
        helix['w'][first_step],
        helix['w'][second_step],
    ]
    original_boundary = (
        get_step_matrix(original_step_parameters[0])
        @ get_step_matrix(original_step_parameters[1])
    )

    def boundary_constraint(local_w):
        trial_w = w.copy()
        trial_w[B] = local_w
        first_step_cords = trial_w[first_step] + helix['w_eq'][first_step]
        second_step_cords = trial_w[second_step] + helix['w_eq'][second_step]
        boundary = get_step_matrix(first_step_cords) @ get_step_matrix(second_step_cords)
        rotation_delta = boundary[:3, :3] @ original_boundary[:3, :3].T
        rotation_residual = 0.5 * np.array([
            rotation_delta[2, 1] - rotation_delta[1, 2],
            rotation_delta[0, 2] - rotation_delta[2, 0],
            rotation_delta[1, 0] - rotation_delta[0, 1],
        ])
        translation_residual = boundary[:3, 3] - original_boundary[:3, 3]
        return np.concatenate((translation_residual, rotation_residual))

    constraint = NonlinearConstraint(boundary_constraint, 0, 0)
    # step_change_bounds = [(-3,3), (-3,3), (-5,5), (-10,10), (-10,10), (-20,20)]
    # intra_change_bounds = [(-1,1), (-1,1), (-1,1), (-10,10), (-10,10), (-10,10)]
    # bounds = step_change_bounds + intra_change_bounds + step_change_bounds  # realistic bounds for the optimization variables
    result = scipy.optimize.minimize(
        objective,
        w[B].copy(),
        method='SLSQP',
        constraints=constraint,
        # bounds=bounds,
        options={'maxiter': 1000, 'ftol': 1e-10},
    )
    if not result.success:
        raise RuntimeError(f"Mutation optimization failed at base pair {bp_index}: {result.message}")

    helix['w'][B] = result.x + helix['w_eq'][B]
    calculate_wKw(helix)

    return helix

def write_rebuild_file(helix, filename):

    with open(filename + '.txt', 'w') as f:
        f.write("# Sequence: " + helix['strand_sequences'][0] + "\n")
        j = 0
        for i, val in enumerate(helix['w']):
            if j == 0:
                f.write(helix['strand_sequences'][0][i//12] + "-" + helix['strand_sequences'][1][i//12] + " ")
            j +=1
            if j > 6:
                f.write(f"{val:.6f} ")
                if (i + 1) % 12 == 0:
                    f.write("\n")
                    j=0


K_intra_inter_DNA = {}
with open(csv_directory / "K_intra_inter_DNA.csv") as Kf:
    data_raw = csv.reader(Kf)
    data = [list(row) for row in data_raw]
    for row in data:
        hexamer = row[0]
        values = list(map(float, row[1:67]))
        if hexamer not in K_intra_inter_DNA:
            K_intra_inter_DNA[hexamer] = []
        K_intra_inter_DNA[hexamer].append(values)
equalibrium_params = load_equalibrium_params()

param_types = ['shear', 'stretch', 'stagger', 'buckle', 'propeller', 'opening', 'shift', 'slide', 'rise', 'tilt', 'roll', 'twist']



if __name__ == "__main__":

    print("\n\n\nloaded\n\n\n")
    # load in helices which were previously found by find_bound_double_strands.py
    helices = []
    helix = None

    md_results_dir = root_dir / "Offset_energy_calculator" / "MD_Results"
    for helix_file in md_results_dir.glob("MD_averaged_parameters_of_helix_*.dat"):
        helix = extract_data(helix_file)
        calculate_wKw(helix)


    helix = get_helix_snippet(helix, 0, -1)

    print_helix_text_reprensentation(helix)
    print("wKw = " + str(calculate_wKw(helix)) + "\n\n\n")


    # write_rebuild_file(helix, "rebuild_helix_long_initial")

    old_helix = copy.deepcopy(helix)
    last_energy = calculate_wKw(helix)
    for iteration in range(3):
        print("Iteration " + str(iteration+1) + " \n")
        for i in range(3, len(helix['strand_sequences'][0])-3):
            helices = []
            helices.append(mutate_bp(helix, i, 'C'))
            helices.append(mutate_bp(helix, i, 'G'))
            helices.append(mutate_bp(helix, i, 'A'))
            helices.append(mutate_bp(helix, i, 'T'))
            energy_values = [calculate_wKw(h) for h in helices]
            helix = copy.deepcopy(helices[energy_values.index(min(energy_values))])
            print(str(round(calculate_wKw(helix), 2)) + " | " + str(helix['strand_sequences'][0]) + " | " + str(round((helix['w']-helix['w_eq']).min(), 2)) + " " + str(round((helix['w']-helix['w_eq']).max(), 2)))
            if( calculate_wKw(helix) > last_energy):
                print("fuck")
            last_energy = calculate_wKw(helix)

    with open('Mutate/mutation_information.txt', 'w') as f:
        for i, newRes in enumerate(helix['strand_sequences'][0]):    # print to mutation information file in format <residue index> <new resname>
            if(newRes != old_helix['strand_sequences'][0][i]):
                f.write(str(helix['strand_res_inds'][0][i]+1) + " D" + newRes + '\n')
                f.write(str(helix['strand_res_inds'][1][i]+1) + " D" + complements[newRes] + '\n')
        print(old_helix['strand_sequences'][0] + ' E='+str(old_helix['wKw'])+ " < old  |  new > "+ helix['strand_sequences'][0] + ' E='+str(helix['wKw'])) 

        print_helix_text_reprensentation(helix)


             











             
    # write_rebuild_file(helix, "rebuild_helix_optim_broken")   


    # helix = {
    #     'strand_sequences' : ['GTTCGCCGGCTTTCCCCGTCAAGCTCTAAA', None],
    #     'strand_res_inds' : None,
    #     'w' : build_equilibrium_w_vector('GTTCGCCGGCTTTCCCCGTCAAGCTCTAAA'),
    #     'wKw' : None,
    #     'w_eq' : None,
    #     'K' : None
    # }
    # helix['strand_sequences'][1] = ''.join([complements[base] for base in helix['strand_sequences'][0]])


    # helix['w'][10 + 12 * 11] += 15
    # helix['w'][10 + 12 * 12] += 15
    # helix['w'][10 + 12 * 13] += 15
    # helix['w'][10 + 12 * 14] += 15
    # helix['w'][10 + 12 * 15] -= 15
    # helix['w'][10 + 12 * 16] -= 15
    # helix['w'][10 + 12 * 17] -= 15
    # helix['w'][10 + 12 * 18] -= 15

    # write_rebuild_file(helix, "rebuild_helix_80")        

    # print()



# TTCCCTTCCTTTCTCGCCACGTTCGCCGGCTTTCCCCGTCAAGCTCTAAATCGGGGGCTCCCTTTAGGGTTCCGATTTAGTGCTTTACGGCACCTCGACCCCAAAAAACTTGATTTGGGTGATGGTTCACGTAGTGGGCCATCGCCCTGATAGACGGTTTTTCGCCCTTTGACGTTGGAGTCCACGTTCTTTAATAGTGGACTCTTGTTCCAAACTGGAACAACACTCAACCCTATCTCGGGCTATTCTTTTGATTTATAAGGGATTTTGCCGATTTCGGAACCACCATCAAACAGGATTTTCGCCTGCTGGGGCAAACCAGCGTGGACCGCTTGCTGCAACTCTCTCAGGGCCAGGCGGTGAAGGGCAATCAGCTGTTGCCCGTCTCACTGGTGAAAAGAAAAACCACCCTGGCGCCCAATACGCAAACCGCCTCTCCCCGCGCGTTGGCCGATTCATTAATGCAGCTGGCACGACAGGTTTCCCGACTGGAAAGCGGGCAGTGAGCGCAACGCAATTAATGTGAGTTAGCTCACTCATTAGGCACCCCAGGCTTTACACTTTATGCTTCCGGCTCGTATGTTGTGTGGAATTGTGAGCGGATAACAATTTCACACAGGAAACAGCTATGACCATGATTACGAATTCGAGCTCGGTACCCGGGGATCCTCTAGAGTCGACCTGCAGGCATGCAAGCTTGGCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGCGAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGGCGCTTTGCCTGGTTTCCGGCACCAGAAGCGGTGCCGGAAAGCTGGCTGGAGTGCGATCTTCCTGAGGCCGATACTGTCGTCGTCCCCTCAAACTGGCAGATGCACGGTTACGATGCGCCCATCTACACCAACGTGACCTATCCCATTACGGTCAATCCGCCGTTTGTTCCCACGGAGAATCCGACGGGTTGTTACTCGCTCACATTTAATGTTGATGAAAGCTGGCTACAGGAAGGCCAGACGCGAATTATTTTTGATGGCGTTCCTATTGGTTAAAAAATGAGCTGATTTAACAAAAATTTAATGCGAATTTTAACAAAATATTAACGTTTACAATTTAAATATTTGCTTATACAATCTTCCTGTTTTTGGGGCTTTTCTGATTATCAACCGGGGTACATATGATTGACATGCTAGTTTTACGATTACCGTTCATCGATTCTCTTGTTTGCTCCAGACTCTCAGGCAATGACCTGATAGCCTTTGTAGATCTCTCAAAAATAGCTACCCTCTCCGGCATTAATTTATCAGCTAGAACGGTTGAATATCATATTGATGGTGATTTGACTGTCTCCGGCCTTTCTCACCCTTTTGAATCTTTACCTACACATTACTCAGGCATTGCATTTAAAATATATGAGGGTTCTAAAAATTTTTATCCTTGCGTTGAAATAAAGGCTTCTCCCGCAAAAGTATTACAGGGTCATAATGTTTTTGGTACAACCGATTTAGCTTTATGCTCTGAGGCTTTATTGCTTAATTTTGCTAATTCTTTGCCTTGCCTGTATGATTTATTGGATGTTAATGCTACTACTATTAGTAGAATTGATGCCACCTTTTCAGCTCGCGCCCCAAATGAAAATATAGCTAAACAGGTTATTGACCATTTGCGAAATGTATCTAATGGTCAAACTAAATCTACTCGTTCGCAGAATTGGGAATCAACTGTTATATGGAATGAAACTTCCAGACACCGTACTTTAGTTGCATATTTAAAACATGTTGAGCTACAGCATTATATTCAGCAATTAAGCTCTAAGCCATCCGCAAAAATGACCTCTTATCAAAAGGAGCAATTAAAGGTACTCTCTAATCCTGACCTGTTGGAGTTTGCTTCCGGTCTGGTTCGCTTTGAAGCTCGAATTAAAACGCGATATTTGAAGTCTTTCGGGCTTCCTCTTAATCTTTTTGATGCAATCCGCTTTGCTTCTGACTATAATAGTCAGGGTAAAGACCTGATTTTTGATTTATGGTCATTCTCGTTTTCTGAACTGTTTAAAGCATTTGAGGGGGATTCAATGAATATTTATGACGATTCCGCAGTATTGGACGCTATCCAGTCTAAACATTTTACTATTACCCCCTCTGGCAAAACTTCTTTTGCAAAAGCCTCTCGCTATTTTGGTTTTTATCGTCGTCTGGTAAACGAGGGTTATGATAGTGTTGCTCTTACTATGCCTCGTAATTCCTTTTGGCGTTATGTATCTGCATTAGTTGAATGTGGTATTCCTAAATCTCAACTGATGAATCTTTCTACCTGTAATAATGTTGTTCCGTTAGTTCGTTTTATTAACGTAGATTTTTCTTCCCAACGTCCTGACTGGTATAATGAGCCAGTTCTTAAAATCGCATAAGGTAATTCACAATGATTAAAGTTGAAATTAAACCATCTCAAGCCCAATTTACTACTCGTTCTGGTGTTTCTCGTCAGGGCAAGCCTTATTCACTGAATGAGCAGCTTTGTTACGTTGATTTGGGTAATGAATATCCGGTTCTTGTCAAGATTACTCTTGATGAAGGTCAGCCAGCCTATGCGCCTGGTCTGTACACCGTTCATCTGTCCTCTTTCAAAGTTGGTCAGTTCGGTTCCCTTATGATTGACCGTCTGCGCCTCGTTCCGGCTAAGTAACATGGAGCAGGTCGCGGATTTCGACACAATTTATCAGGCGATGATACAAATCTCCGTTGTACTTTGTTTCGCGCTTGGTATAATCGCTGGGGGTCAAAGATGAGTGTTTTAGTGTATTCTTTTGCCTCTTTCGTTTTAGGTTGGTGCCTTCGTAGTGGCATTACGTATTTTACCCGTTTAATGGAAACTTCCTCATGAAAAAGTCTTTAGTCCTCAAAGCCTCTGTAGCCGTTGCTACCCTCGTTCCGATGCTGTCTTTCGCTGCTGAGGGTGACGATCCCGCAAAAGCGGCCTTTAACTCCCTGCAAGCCTCAGCGACCGAATATATCGGTTATGCGTGGGCGATGGTTGTTGTCATTGTCGGCGCAACTATCGGTATCAAGCTGTTTAAGAAATTCACCTCGAAAGCAAGCTGATAAACCGATACAATTAAAGGCTCCTTTTGGAGCCTTTTTTTTGGAGATTTTCAACGTGAAAAAATTATTATTCGCAATTCCTTTAGTTGTTCCTTTCTATTCTCACTCCGCTGAAACTGTTGAAAGTTGTTTAGCAAAATCCCATACAGAAAATTCATTTACTAACGTCTGGAAAGACGACAAAACTTTAGATCGTTACGCTAACTATGAGGGCTGTCTGTGGAATGCTACAGGCGTTGTAGTTTGTACTGGTGACGAAACTCAGTGTTACGGTACATGGGTTCCTATTGGGCTTGCTATCCCTGAAAATGAGGGTGGTGGCTCTGAGGGTGGCGGTTCTGAGGGTGGCGGTTCTGAGGGTGGCGGTACTAAACCTCCTGAGTACGGTGATACACCTATTCCGGGCTATACTTATATCAACCCTCTCGACGGCACTTATCCGCCTGGTACTGAGCAAAACCCCGCTAATCCTAATCCTTCTCTTGAGGAGTCTCAGCCTCTTAATACTTTCATGTTTCAGAATAATAGGTTCCGAAATAGGCAGGGGGCATTAACTGTTTATACGGGCACTGTTACTCAAGGCACTGACCCCGTTAAAACTTATTACCAGTACACTCCTGTATCATCAAAAGCCATGTATGACGCTTACTGGAACGGTAAATTCAGAGACTGCGCTTTCCATTCTGGCTTTAATGAGGATTTATTTGTTTGTGAATATCAAGGCCAATCGTCTGACCTGCCTCAACCTCCTGTCAATGCTGGCGGCGGCTCTGGTGGTGGTTCTGGTGGCGGCTCTGAGGGTGGTGGCTCTGAGGGTGGCGGTTCTGAGGGTGGCGGCTCTGAGGGAGGCGGTTCCGGTGGTGGCTCTGGTTCCGGTGATTTTGATTATGAAAAGATGGCAAACGCTAATAAGGGGGCTATGACCGAAAATGCCGATGAAAACGCGCTACAGTCTGACGCTAAAGGCAAACTTGATTCTGTCGCTACTGATTACGGTGCTGCTATCGATGGTTTCATTGGTGACGTTTCCGGCCTTGCTAATGGTAATGGTGCTACTGGTGATTTTGCTGGCTCTAATTCCCAAATGGCTCAAGTCGGTGACGGTGATAATTCACCTTTAATGAATAATTTCCGTCAATATTTACCTTCCCTCCCTCAATCGGTTGAATGTCGCCCTTTTGTCTTTGGCGCTGGTAAACCATATGAATTTTCTATTGATTGTGACAAAATAAACTTATTCCGTGGTGTCTTTGCGTTTCTTTTATATGTTGCCACCTTTATGTATGTATTTTCTACGTTTGCTAACATACTGCGTAATAAGGAGTCTTAATCATGCCAGTTCTTTTGGGTATTCCGTTATTATTGCGTTTCCTCGGTTTCCTTCTGGTAACTTTGTTCGGCTATCTGCTTACTTTTCTTAAAAAGGGCTTCGGTAAGATAGCTATTGCTATTTCATTGTTTCTTGCTCTTATTATTGGGCTTAACTCAATTCTTGTGGGTTATCTCTCTGATATTAGCGCTCAATTACCCTCTGACTTTGTTCAGGGTGTTCAGTTAATTCTCCCGTCTAATGCGCTTCCCTGTTTTTATGTTATTCTCTCTGTAAAGGCTGCTATTTTCATTTTTGACGTTAAACAAAAAATCGTTTCTTATTTGGATTGGGATAAATAATATGGCTGTTTATTTTGTAACTGGCAAATTAGGCTCTGGAAAGACGCTCGTTAGCGTTGGTAAGATTCAGGATAAAATTGTAGCTGGGTGCAAAATAGCAACTAATCTTGATTTAAGGCTTCAAAACCTCCCGCAAGTCGGGAGGTTCGCTAAAACGCCTCGCGTTCTTAGAATACCGGATAAGCCTTCTATATCTGATTTGCTTGCTATTGGGCGCGGTAATGATTCCTACGATGAAAATAAAAACGGCTTGCTTGTTCTCGATGAGTGCGGTACTTGGTTTAATACCCGTTCTTGGAATGATAAGGAAAGACAGCCGATTATTGATTGGTTTCTACATGCTCGTAAATTAGGATGGGATATTATTTTTCTTGTTCAGGACTTATCTATTGTTGATAAACAGGCGCGTTCTGCATTAGCTGAACATGTTGTTTATTGTCGTCGTCTGGACAGAATTACTTTACCTTTTGTCGGTACTTTATATTCTCTTATTACTGGCTCGAAAATGCCTCTGCCTAAATTACATGTTGGCGTTGTTAAATATGGCGATTCTCAATTAAGCCCTACTGTTGAGCGTTGGCTTTATACTGGTAAGAATTTGTATAACGCATATGATACTAAACAGGCTTTTTCTAGTAATTATGATTCCGGTGTTTATTCTTATTTAACGCCTTATTTATCACACGGTCGGTATTTCAAACCATTAAATTTAGGTCAGAAGATGAAATTAACTAAAATATATTTGAAAAAGTTTTCTCGCGTTCTTTGTCTTGCGATTGGATTTGCATCAGCATTTACATATAGTTATATAACCCAACCTAAGCCGGAGGTTAAAAAGGTAGTCTCTCAGACCTATGATTTTGATAAATTCACTATTGACTCTTCTCAGCGTCTTAATCTAAGCTATCGCTATGTTTTCAAGGATTCTAAGGGAAAATTAATTAATAGCGACGATTTACAGAAGCAAGGTTATTCACTCACATATATTGATTTATGTACTGTTTCCATTAAAAAAGGTAATTCAAATGAAATTGTTAAATGTAATTAATTTTGTTTTCTTGATGTTTGTTTCATCATCTTCTTTTGCTCAGGTAATTGAAATGAATAATTCGCCTCTGCGCGATTTTGTAACTTGGTATTCAAAGCAATCAGGCGAATCCGTTATTGTTTCTCCCGATGTAAAAGGTACTGTTACTGTATATTCATCTGACGTTAAACCTGAAAATCTACGCAATTTCTTTATTTCTGTTTTACGTGCAAATAATTTTGATATGGTAGGTTCTAACCCTTCCATTATTCAGAAGTATAATCCAAACAATCAGGATTATATTGATGAATTGCCATCATCTGATAATCAGGAATATGATGATAATTCCGCTCCTTCTGGTGGTTTCTTTGTTCCGCAAAATGATAATGTTACTCAAACTTTTAAAATTAATAACGTTCGGGCAAAGGATTTAATACGAGTTGTCGAATTGTTTGTAAAGTCTAATACTTCTAAATCCTCAAATGTATTATCTATTGACGGCTCTAATCTATTAGTTGTTAGTGCTCCTAAAGATATTTTAGATAACCTTCCTCAATTCCTTTCAACTGTTGATTTGCCAACTGACCAGATATTGATTGAGGGTTTGATATTTGAGGTTCAGCAAGGTGATGCTTTAGATTTTTCATTTGCTGCTGGCTCTCAGCGTGGCACTGTTGCAGGCGGTGTTAATACTGACCGCCTCACCTCTGTTTTATCTTCTGCTGGTGGTTCGTTCGGTATTTTTAATGGCGATGTTTTAGGGCTATCAGTTCGCGCATTAAAGACTAATAGCCATTCAAAAATATTGTCTGTGCCACGTATTCTTACGCTTTCAGGTCAGAAGGGTTCTATCTCTGTTGGCCAGAATGTCCCTTTTATTACTGGTCGTGTGACTGGTGAATCTGCCAATGTAAATAATCCATTTCAGACGATTGAGCGTCAAAATGTAGGTATTTCCATGAGCGTTTTTCCTGTTGCAATGGCTGGCGGTAATATTGTTCTGGATATTACCAGCAAGGCCGATAGTTTGAGTTCTTCTACTCAGGCAAGTGATGTTATTACTAATCAAAGAAGTATTGCTACAACGGTTAATTTGCGTGATGGACAGACTCTTTTACTCGGTGGCCTCACTGATTATAAAAACAC



















        


    # helices_temp = [get_helix_snippet(helices[4], 7, -6)]
    # total_energy = calculate_total_energy(helices_temp)

    # print_helix_text_reprensentation(helices_temp[0])
    # print(total_energy)


