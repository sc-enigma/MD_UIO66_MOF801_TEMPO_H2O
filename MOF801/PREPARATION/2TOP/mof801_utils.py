import numpy as np

from atom import Atom, calculate_angle

def remove_extra_oxygens(atoms):
    def getNeighbours(atom_idx):
        return [atoms[adj_idx] for adj_idx in atoms[atom_idx].adjacency]
    
    def getNeighbourElems(atom_idx):
        return np.sort([atom.name[0] for atom in getNeighbours(atom_idx)])
    
    shift_ids = {}
    removed_count = 0
    new_atoms = []
    for atom_idx in range(len(atoms)):
        neighbours = getNeighbourElems(atom_idx)
        if atoms[atom_idx].name[0] == 'O' and len(neighbours) == 1 and neighbours[0] == 'O':
            shift_ids[atom_idx] = -1
            removed_count += 1
        else:
            shift_ids[atom_idx] = atom_idx - removed_count
            new_atoms.append(atoms[atom_idx])
    print(np.size(atoms), np.size(new_atoms))
    atoms = new_atoms
    
    for atom_idx in range(len(atoms)):
        adjacency = atoms[atom_idx].adjacency
        for adj_idx in range(len(adjacency)):
            atoms[atom_idx].adjacency[adj_idx] = shift_ids[adjacency[adj_idx]]
            
        if -1 in atoms[atom_idx].adjacency:
            atoms[atom_idx].adjacency.remove(-1)
    return atoms
    
def define_mof801_atom_types(atoms):
    def getNeighbours(atom_idx):
        return [atoms[adj_idx] for adj_idx in atoms[atom_idx].adjacency]
    
    def getNeighbourElems(atom_idx):
        return np.sort([atom.name[0] for atom in getNeighbours(atom_idx)])
    
    def getElem(atom_idx):
        return atoms[atom_idx].name[0]
        
    # define Zr: itcFF_Zr
    for atom_idx in range(len(atoms)):
        if atoms[atom_idx].name[0] == 'Z':
            atoms[atom_idx].atom_type = 'itcFF_Zr'
    
    # define H1: itcFF_H1
    for atom_idx in range(len(atoms)):
        if atoms[atom_idx].mol2name == 'H':
            atoms[atom_idx].atom_type = 'itcFF_H1'
        
    # define H7: itcFF_H7
    for atom_idx in range(len(atoms)):
        if atoms[atom_idx].mol2name == 'H7':
            atoms[atom_idx].atom_type = 'itcFF_H7'
    
    # define O1: itcFF_O1
    for atom_idx in range(len(atoms)):
        if atoms[atom_idx].name[0] == 'O':
            neighbours = getNeighbourElems(atom_idx)
            if len(neighbours) == 2 and neighbours[0] == 'C' and neighbours[1] == 'Z':
                atoms[atom_idx].atom_type = 'itcFF_O1'
    
    # define O2: itcFF_O2
    for atom_idx in range(len(atoms)):
        if atoms[atom_idx].name[0] == 'O':
            neighbours = getNeighbourElems(atom_idx)
            if len(neighbours) == 3 and neighbours[0] == 'Z' and neighbours[1] == 'Z' and neighbours[2] == 'Z':
                atoms[atom_idx].atom_type = 'itcFF_O2'
            if len(neighbours) == 4 and neighbours[0] == 'H' and neighbours[1] == 'Z' and neighbours[2] == 'Z' and neighbours[3] == 'Z':
                atoms[atom_idx].atom_type = 'itcFF_O2'
                    
    # define C1: itcFF_C1
    for atom_idx in range(len(atoms)):
        if atoms[atom_idx].name[0] == 'C':
            neighbours = getNeighbourElems(atom_idx)
            if len(neighbours) == 3 and neighbours[0] == 'C' and neighbours[1] == 'O' and neighbours[2] == 'O':
                atoms[atom_idx].atom_type = 'itcFF_C1'
                
    # define C2: itcFF_C2
    for atom_idx in range(len(atoms)):
        if atoms[atom_idx].name[0] == 'C':
            neighbours = getNeighbourElems(atom_idx)
            if len(neighbours) == 3 and neighbours[0] == 'C' and neighbours[1] == 'C' and neighbours[2] == 'H':
                atoms[atom_idx].atom_type = 'itcFF_C2'
    return atoms

def check_mof801_atom_names(atoms):
    def getNeighbours(atom_idx):
        return [atoms[adj_idx] for adj_idx in atoms[atom_idx].adjacency]
    
    def getNeighbourElems(atom_idx):
        return np.sort([atom.name[0] for atom in getNeighbours(atom_idx)])
    
    def getElem(atom_idx):
        return atoms[atom_idx].name[0]
    
    for atom_idx in range(len(atoms)):
        if atoms[atom_idx].atom_type == 'unk_type':
            print('ERROR', atom_idx, atoms[atom_idx].name, atoms[atom_idx].mol2name)            
            print(getElem(atom_idx), getNeighbourElems(atom_idx))

def define_mof801_atom_names(atoms):  
    for atom_idx in range(len(atoms)):
        atoms[atom_idx].name = atoms[atom_idx].atom_type.replace('itcFF_', '')
        atoms[atom_idx].resid_idx = 1
        atoms[atom_idx].resid = 'MOF'
        atoms[atom_idx].atom_idx = atom_idx
    return atoms