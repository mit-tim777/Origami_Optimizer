#!/usr/bin/env python3
import sys

pdb_file = 'output.pdb'
atoms_to_keep = {  #                     backbone atoms beginning                       reference atoms                                                                     referece                backbone end
    'A': ["P","OP1","OP2","O5'","C5'","H5'","H5''","C4'","H4'","O4'","C1'","H1'",      "N9","C8","H8",   "C4",       "C3'","H3'","C2'","H2'","H2''","O3'"],
    'T': ["P","OP1","OP2","O5'","C5'","H5'","H5''","C4'","H4'","O4'","C1'","H1'",      "N1","C6","H6",   "C2",       "C3'","H3'","C2'","H2'","H2''","O3'"],
    'C': ["P","OP1","OP2","O5'","C5'","H5'","H5''","C4'","H4'","O4'","C1'","H1'",      "N1","C6","H6",   "C2",       "C3'","H3'","C2'","H2'","H2''","O3'"],
    'G': ["P","OP1","OP2","O5'","C5'","H5'","H5''","C4'","H4'","O4'","C1'","H1'",      "N9","C8","H8",   "C4",       "C3'","H3'","C2'","H2'","H2''","O3'"]
}
                # deleted atoms
# """ "N7","C5","C6","N6","H61","H62","N1","C2","H2","N3",       """ 
# """ "C5","C7","H71","H72","H73","C4","O4","N3","H3","O2",      """ 
# """ "N7","C5","C6","O6","N1","H1","C4","N2","H21","H22","N3",  """ 
# """ "C5","H5","N4","H41","H42","N3","C2","O2",                 """ 


change_names = {  #   hold list of atoms to change atoms from (oldname, newname),... when changing from change_names[old_residue][new_residue] 
    'A': { 'A' : [],                                                     'T' : [('N9','N1'),('C8','C6'),('H8','H6'),('C4','C2')],  'C' : [('N9','N1'),('C8','C6'),('H8','H6'),('C4','C2')], 'G' : []},



    'T': { 'A' : [('N1','N9'),('C6','C8'),('H6','H8'),('C2','C4')],      'T' : [],                                                 'C' : [],                                                'G' : [('N1','N9'),('C6','C8'),('H6','H8'),('C2','C4')]},



    'C': { 'A' : [('N1','N9'),('C6','C8'),('H6','H8'),('C2','C4')],      'T' : [],                                                 'C' : [],                                                'G' : [('N1','N9'),('C6','C8'),('H6','H8'),('C2','C4')]},



    'G': { 'A' : [],                                                     'T' : [('N9','N1'),('C8','C6'),('H8','H6'),('C4','C2')],  'C' : [('N9','N1'),('C8','C6'),('H8','H6'),('C4','C2')], 'G' : []}
}

def change_residue(res_index,new_name):

    out = []
    for line in open(pdb_file):
        
        if not line.startswith(("ATOM", "TER")):
            out.append(line)
            continue
        resid = int(line[22:26])
        resname = str(line[17:20]).strip()
        atom_name = str(line[12:16]).strip()

        if resid == res_index and len(resname) != 3:
            line = line[:17] + f"{'D'+new_name:>3}" + line[20:]
            
            change_list = [i[0] for i in change_names[resname[1]][new_name]]
            change_index = change_list.index(atom_name) if atom_name in change_list else -1
            if change_index >= 0:
                line = line[:12] + f"{change_names[resname[1]][new_name][change_index][1]:^4}" + line[16:]

            if atom_name in atoms_to_keep[resname[1]]:
                out.append(line)

        else:
            out.append(line)


    open("output.pdb", "w").writelines(out)

if __name__ == "__main__":
    with open('Mutate/mutation_information.txt', 'r') as f:
        for line in f.readlines():
            res_to_change = int(line.split()[0])
            res_new_name = line.split()[1]
            print("replacing residue " + str(res_to_change) + " with " + res_new_name )
            change_residue(res_to_change, res_new_name[1])


