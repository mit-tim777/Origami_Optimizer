import math
import sys

TOLERANCE_DISTANCE = 8
TOLERATED_UNMATCHES = 1
ALLOW_ONE_STRAND_HELIX = False

# trys to find all double strand helices in a pdb file and write their base pairing in a file 

def parse_helix_pdb(file):  # parse a pdb file into an array of atoms with their properties
    atoms = []
    
    chain_N = 0
    last_res_N = 0
    res_displacement = 0
    for line in file:
        if line.startswith(('ATOM')):
            chain_ID = line[21].strip()
            global_res_ID = int(line[22:26])
            last_res_N = int(line[22:26])
            if(chain_ID != ''):     # existence of chain_ID in pdb file depends on format
                chain_ID = ord(chain_ID) - ord('A')
                res_ID = int(line[22:26])
                if(chain_ID > 25):
                    chain_ID -= 6
            else:
                chain_ID = chain_N
                res_ID = int(line[22:26]) - res_displacement
            atom = {
                'atom_serial': int(line[6:11]) - chain_ID,
                
                'atom_name': line[12:16].strip(),
                'res_name': line[17:20].strip(),
                'res_ind': res_ID-1,
                'global_res_ind': global_res_ID-1,
                'chain_ind': chain_ID,
                'element': line[76:78].strip(),
                'x': float(line[30:38]),
                'y': float(line[38:46]),
                'z': float(line[46:54])
                
            }
            atoms.append(atom)
        if line.startswith(('TER')):
            chain_N += 1
            res_displacement = last_res_N
    
    if(len(atoms) <= 1):
        print('no atoms in pdb')
        exit()
    return atoms

def seperate_Atoms(atoms):  # sort the atoms into an array of chains that contain residues that contain atoms
    chains = [None] * (atoms[-1]['chain_ind']+1)
    for i in range(len(chains)):
        aa = [k for k in atoms if k['chain_ind'] == i]
        chains[i] = [None] * (aa[-1]['res_ind'] + 1)

    for i in range(len(chains)):
        for j in range(len(chains[i])):
            chains[i][j] = []

    for atom in atoms:
        chains[atom['chain_ind']][atom['res_ind']].append(atom)

    return chains

def pick_binding_atoms(chains): # pick for each residue only the atom that binds to a fitting counter residue

    dict = {
    'DC5': 'N3',
    'DC3': 'N3',
    'DC' : 'N3',
    'DT3': 'N3',
    'DT5': 'N3',
    'DT' : 'N3',
    'DA5': 'N1',
    'DA3': 'N1',
    'DA' : 'N1',
    'DG3': 'N1',
    'DG5': 'N1',
    'DG' : 'N1'
    }

    chains_red = [None] * len(chains)
    for chain_ind in range(len(chains)):
        chains_red[chain_ind] = [None] * len(chains[chain_ind])
    
    for chain_ind in range(len(chains)):
        for res_ind in range(len(chains[chain_ind])):
            chains_red[chain_ind][res_ind] = [ i for i in chains[chain_ind][res_ind] if i['atom_name'] == dict[i['res_name']] ][0]
    
    return chains_red


def could_be_binded(atom1, atom2):  # checks if two binding site atoms are close enough and of matching bases
    if(atom1 == None or atom2 == None): return False

    matching_pairs = {frozenset(("A", "T")), frozenset(("G", "C"))}
    if(not ALLOW_ONE_STRAND_HELIX and atom1['chain_ind'] == atom2['chain_ind'] ):
        return False
    
    close_enough = getDistance(atom1, atom2) < TOLERANCE_DISTANCE
    bases_match = frozenset( ( atom1['res_name'][1], atom2['res_name'][1] ) ) in matching_pairs

    return close_enough and bases_match

def getDistance(atom1, atom2):
    sum = 0
    for i in ['x','y','z']: sum += ((atom1[i]) - (atom2[i]))**2
    return math.sqrt(sum)

def find_closest_resis(chains, atom):   # find the three residues with smallest binding atom distance
    
    closest_dist = [float('inf')]   * 3
    closest_atoms =[None]           * 3
    
    chains_exept = []
    for i in chains:
        if( not atom['chain_ind'] == i[0]['chain_ind'] ): chains_exept.append(i)

    for chain in chains_exept:
        for atom_iter in chain:
            dist = getDistance( atom, atom_iter)
            if(dist < closest_dist[0]):
                closest_atoms[2] = closest_atoms[1]
                closest_atoms[1] = closest_atoms[0]
                closest_dist[2] = closest_dist[1]
                closest_dist[1] = closest_dist[0]
                closest_dist[0] = dist
                closest_atoms[0] = atom_iter
            elif(dist < closest_dist[1]):
                closest_atoms[2] = closest_atoms[1]
                closest_dist[2] = closest_dist[1]
                closest_dist[1] = dist
                closest_atoms[1] = atom_iter
            elif(dist < closest_dist[2]):
                closest_dist[2] = dist
                closest_atoms[2] = atom_iter
    return closest_atoms


def get_next(chains, atom, step):   # moving on the chain find the next residue with direction of step
    if(atom == None): return None
    if(atom['res_ind'] >= -step and atom['res_ind'] < -step+len(chains[atom['chain_ind']])):
        return chains[atom['chain_ind']][atom['res_ind']+step]
    else:
        return None

def get_last_possible_binding(chains, orig_atom1, orig_atom2):  # assuming the two input residues are bound, checks in all directions on the strands if helix conditions are met  

    directons = [[1,1],[-1,-1],[1,-1],[-1,1]]
    last_possible_of_direction = []

    trys = [None] *4

    for i, direc in enumerate(directons):
        atom1 = orig_atom1
        atom2 = orig_atom2
        try_len = 1
        while(True):
            trys[i] = [atom1, atom2, try_len]
            atom1 = get_next(chains, atom1, direc[0])
            atom2 = get_next(chains, atom2, direc[1])
            if(not could_be_binded(atom1, atom2)):
                break
            try_len += 1
    
    len_parrallel = trys[0][2] + trys[1][2]
    len_antiParal = trys[1][2] + trys[2][2]
    if  (len_parrallel > len_antiParal):
        helix = {
            'start_atoms' : (trys[0][0], trys[0][1]),
            'end_atoms'   : (trys[1][0], trys[1][1]),
            'length' : trys[0][2] + trys[1][2]
            }
        return helix
    elif(len_parrallel <= len_antiParal):
        helix = {
            'start_atoms' : (trys[2][0], trys[2][1]),
            'end_atoms'   : (trys[3][0], trys[3][1]),
            'length' : trys[2][2] + trys[3][2]
            }
        return helix

def get_helix_chains(chains, helix):   # helices are defined just by start/end atoms. this gives back the whole chains in between
    helix_chains = [[], []]

    chain1_start = helix['start_atoms'][0]['res_ind']
    chain1_end   = helix['end_atoms'][0]['res_ind']
    chain1_ind   = helix['start_atoms'][0]['chain_ind']

    step = 1 - 2*(chain1_start > chain1_end)
    for i in range(chain1_start, chain1_end+step, step):
        helix_chains[0].append(chains[chain1_ind][i])

    chain2_start = helix['start_atoms'][1]['res_ind']
    chain2_end   = helix['end_atoms'][1]['res_ind']
    chain2_ind   = helix['start_atoms'][1]['chain_ind']

    step = 1 - 2*(chain2_start > chain2_end)
    for i in range(chain2_start, chain2_end+step, step):
        helix_chains[1].append(chains[chain2_ind][i])
    
    if(chains[chain1_ind][0]['res_name'][2] == '5'):
        helix_chains[0],helix_chains[1] = helix_chains[1],helix_chains[0]

    return helix_chains

def atom_is_in_helix(atom, helix):
    chain1_ind   = helix['start_atoms'][0]['chain_ind']
    chain1_start = helix['start_atoms'][0]['res_ind']
    chain1_end   = helix['end_atoms'][0]['res_ind']
    chain2_start = helix['start_atoms'][1]['res_ind']
    chain2_end   = helix['end_atoms'][1]['res_ind']
    chain2_ind   = helix['start_atoms'][1]['chain_ind']
    if(atom['chain_ind'] == chain1_ind):
        if(chain1_start < chain1_end):
            return atom['res_ind'] >= chain1_start and atom['res_ind'] <= chain1_end
        else:
            return atom['res_ind'] <= chain1_start and atom['res_ind'] >= chain1_end
    elif(atom['chain_ind'] == chain2_ind):
        if(chain2_start < chain2_end):
            return atom['res_ind'] >= chain2_start and atom['res_ind'] <= chain2_end
        else:
            return atom['res_ind'] <= chain2_start and atom['res_ind'] >= chain2_end
    else:
        return False

def match_chains(chains):   # goes through all residues and constructs helices 
    helices = []

    for chain in chains:
        for res in chain:
            if(not any([atom_is_in_helix(res, helix) for helix in helices])):
                attempts = []
                binding_candidates = find_closest_resis(chains, res)    # sometimes in a helix the residues are closer so another residue than the residue that mathes the sequence. that why the closest 3 are tried for constructing the helix
                for candidate in binding_candidates:
                    if(could_be_binded(res, candidate)):
                        attempt = get_last_possible_binding(chains, res, candidate)
                        if(attempt != None): attempts.append(attempt)

                lengths = [i['length'] for i in attempts]
                if(len(lengths) != 0):
                    if(max(lengths) > 4):
                        helices.append(attempts[lengths.index(max(lengths))])
    return helices

def get_helix_sequences(chains, helix):
    helix_chains = get_helix_chains(chains, helix)
    sequences = [[],[]]
    one_letter_names = {
        'DC5': 'C',
        'DC3': 'C',
        'DC' : 'C',
        'DT3': 'T',
        'DT5': 'T',
        'DT' : 'T',
        'DA5': 'A',
        'DA3': 'A',
        'DA' : 'A',
        'DG3': 'G',
        'DG5': 'G',
        'DG' : 'G'
    }
    
    for bp1,bp2 in zip(helix_chains[0],helix_chains[1]):
        sequences[0].append(one_letter_names[bp1['res_name']])
        sequences[1].append(one_letter_names[bp2['res_name']])

    return sequences

if __name__ == "__main__":

    pdb_filename = str("output.pdb")
    origami = open(pdb_filename, 'r')

    chains_raw = seperate_Atoms(parse_helix_pdb(origami))
    chains = pick_binding_atoms(chains_raw)

    helices = match_chains(chains)

    with open("Helix_separator/cpptraj_base_pairing.txt", 'w') as f: # write base pairing file for cpptraj in the format '<residue id 1>-<residue id 2>,... sequence1 sequence2
        for helix in helices:
            helix_chains = get_helix_chains(chains, helix)

            for bp1,bp2 in zip(helix_chains[0],helix_chains[1]):
                f.write( str(bp1['global_res_ind']+1) +"-"+str(bp2['global_res_ind']+1)+ "," )


            sequences = get_helix_sequences(chains, helix)
            f.write( " " + "".join(sequences[0])+ " " + "".join(sequences[1]))

            f.write("\n")