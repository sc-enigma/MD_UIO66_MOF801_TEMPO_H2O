import numpy as np
import matplotlib.pyplot as plt

class Atom:
    def __init__(self, name, r):
        self.name = name        
        self.r = np.array(r)
        
        # topology properties
        self.resid_idx = 1
        self.resid = 'UNK'
        self.atom_idx = 0
        self.atom_type = 'unk_type'
        self.mol2name = ''
        
        # internal properties
        self.r_internal = np.zeros(3)
        self.adjacency = []

def mol2_to_atoms(data):
    atoms = []
    
    data_atoms = data['ATOM']
    for atom_idx in range(len(data_atoms)):
        split = data_atoms[atom_idx].split()
        if len(split) < 5:
            continue
        name = split[1]
        r = [float(c) for c in split[2:5]]
        atoms.append(Atom(name, r))
        atoms[-1].atom_idx = atom_idx
        atoms[-1].mol2name = split[5]
    
    if 'BOND' in data.keys():
        data_bonds = data['BOND']
        for bond_idx in range(len(data_bonds)):
            split = data_bonds[bond_idx].split()
            connected = [int(split[1]) - 1, int(split[2]) - 1]
            atoms[connected[0]].adjacency.append(connected[1])
            atoms[connected[1]].adjacency.append(connected[0])
    return atoms

def gro_to_atoms(data):
    atoms = []
    for atom_idx in range(len(data)):
        split = data[atom_idx].split()
        name = split[1]
        r = [float(c) for c in split[3:6]]
        atoms.append(Atom(name, r))
        atoms[-1].resid = data[atom_idx][5:8]
        atoms[-1].atom_idx = atom_idx
    return atoms

def copy_atom(atom):
    a = Atom(atom.name, np.array(atom.r))
    a.resid_idx = atom.resid_idx
    a.resid = atom.resid
    a.atom_idx = atom.atom_idx
    a.atom_type = atom.atom_type
    a.mol2name = atom.mol2name
    a.r_internal = np.array(a.r_internal)
    a.adjacency = atom.adjacency.copy()
    return a

def count_atoms(atoms):
    count = {}
    for atom_idx in range(len(atoms)):
        atom_type = atoms[atom_idx].atom_type
        if atom_type in count.keys():
            count[atom_type] += 1
        else:
            count[atom_type] = 1
    for atom_type in count.keys():
        print(f'{atom_type}    {count[atom_type]}')
        
def dump_atoms(atoms, dim_1=0, dim_2=1, dumpBonds=False):
    points_single = []
    point_pairs = []
    for atom_idx in range(len(atoms)):
        points_single.append(atoms[atom_idx].r)
        for adj_idx in atoms[atom_idx].adjacency:
            point_pairs.append([atoms[atom_idx].r, atoms[adj_idx].r])
        
    points_single = np.transpose(np.array(points_single))
    s = np.ones(len(points_single[0])) * 0.5
    plt.scatter(points_single[dim_1], points_single[dim_2], s=s)
    
    if dumpBonds:
        for pair_idx in range(len(point_pairs)):
            pair = point_pairs[pair_idx]
            plt.plot([pair[0][dim_1], pair[1][dim_1]], [pair[0][dim_2], pair[1][dim_2]])        
    plt.show()
    
def select_atoms(atoms, atom_type):
    selected_ids = []
    for atom_idx in range(len(atoms)):
        if atoms[atom_idx].atom_type == atom_type:
            selected_ids.append(atom_idx)
    return selected_ids

def remove_atoms(atoms, ids_to_remove):
    if len(ids_to_remove) == 0:
        return atoms
    shift_ids = {}
    removed_count = 0
    for atom_idx in range(len(atoms)):
        if atom_idx in ids_to_remove:
            shift_ids[atom_idx] = -1
            removed_count += 1
        else:
            shift_ids[atom_idx] = atom_idx - removed_count
            
    atoms = [copy_atom(atoms[atom_idx]) for atom_idx in range(len(atoms)) if not (atom_idx in ids_to_remove)]

    for atom_idx in range(len(atoms)):
        adjacency = atoms[atom_idx].adjacency
        for adj_idx in range(len(adjacency)):
            atoms[atom_idx].adjacency[adj_idx] = shift_ids[adjacency[adj_idx]]
            
        while -1 in atoms[atom_idx].adjacency:
            atoms[atom_idx].adjacency.remove(-1)
    return atoms
    
def remove_non_bonded_atoms(atoms):
    ids_to_remove = []
    for atom_idx in range(len(atoms)):
        if len(atoms[atom_idx].adjacency) == 0:
            ids_to_remove.append(atom_idx)
    return remove_atoms(atoms, ids_to_remove)

def shift_atoms(atoms, vec):
    for atom_idx in range(len(atoms)):
        atoms[atom_idx].r += vec
    return atoms

def calculate_distance(atom1, atom2):
    v1 = atom1.r - atom2.r
    return np.linalg.norm(v1) / 10.

def calculate_angle(atom1, atom2, atom3):
    v1 = atom1.r - atom2.r
    v2 = atom3.r - atom2.r
    cosine_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    cosine_angle = np.clip(cosine_angle, -1.0, 1.0)
    return np.arccos(cosine_angle) / np.pi * 180.0

def calculate_dihedral(atom1, atom2, atom3, atom4):
    v1 = atom2.r - atom1.r
    v2 = atom3.r - atom2.r
    v3 = atom4.r - atom3.r
    n1 = np.cross(v1, v2)
    n2 = np.cross(v2, v3)
    cosine_angle = np.dot(n1, n2) / (np.linalg.norm(n1) * np.linalg.norm(n2))
    cosine_angle = np.clip(cosine_angle, -1.0, 1.0)
    return np.arccos(cosine_angle)

def set_residue(atoms, res):
    for atom_idx in range(len(atoms)):
        atoms[atom_idx].resid = res
    return atoms

def add_atoms(atoms1, atoms2):
    atoms = atoms1 + atoms2
    for atom_idx in range(len(atoms1), len(atoms)):
        for adjacency_idx in range(len(atoms[atom_idx].adjacency)):
            atoms[atom_idx].adjacency[adjacency_idx] += len(atoms1)
    return atoms

    