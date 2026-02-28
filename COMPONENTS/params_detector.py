import numpy as np
from atom import Atom, calculate_distance, calculate_angle, calculate_dihedral

def detect_bonds(atoms):
    unique_bonds = set()
    length = {}
    for atom_idx in range(len(atoms)):
        for adj_idx in atoms[atom_idx].adjacency:
            bond = [atoms[atom_idx].atom_type, atoms[adj_idx].atom_type]
            if bond[0] > bond[-1]:
                bond.reverse()
            unique_bonds.add((bond[0], bond[1]))
            length[(bond[0], bond[1])] = calculate_distance(atoms[atom_idx], atoms[adj_idx])
    print('BONDS')
    for rec in unique_bonds:
        print(rec[0], rec[1], length[rec])
    
def detect_angles(atoms):
    unique_angles = set()
    values = {}
    for atom_idx in range(len(atoms)):
        for adj_idx1 in atoms[atom_idx].adjacency:
            for adj_idx2 in atoms[adj_idx1].adjacency:
                if len(set([atom_idx, adj_idx1, adj_idx2])) != 3:
                    continue
                angle = [atoms[atom_idx].atom_type, atoms[adj_idx1].atom_type, atoms[adj_idx2].atom_type]
                if angle[0] > angle[-1]:
                    angle.reverse()
                unique_angles.add((angle[0], angle[1], angle[2]))
                values[(angle[0], angle[1], angle[2])] = calculate_angle(atoms[atom_idx], atoms[adj_idx1], atoms[adj_idx2])
    print('ANGLES')
    for rec in unique_angles:
        print(rec[0], rec[1], rec[2], round(values[rec] / np.pi * 180, 4))

def detect_dihedrals(atoms):
    unique_dihedrals = set()
    values = {}
    for atom_idx in range(len(atoms)):
        for adj_idx1 in atoms[atom_idx].adjacency:
            for adj_idx2 in atoms[adj_idx1].adjacency:
                for adj_idx3 in atoms[adj_idx2].adjacency:
                    if len(set([atom_idx, adj_idx1, adj_idx2, adj_idx3])) != 4:
                       continue
                    dihedral = [atoms[atom_idx].atom_type, atoms[adj_idx1].atom_type, atoms[adj_idx2].atom_type, atoms[adj_idx3].atom_type]
                    if dihedral[0] > dihedral[-1]:
                        dihedral.reverse()
                    unique_dihedrals.add((dihedral[0], dihedral[1], dihedral[2], dihedral[3]))
                    values[(dihedral[0], dihedral[1], dihedral[2], dihedral[3])] = calculate_dihedral(atoms[atom_idx], atoms[adj_idx1], atoms[adj_idx2], atoms[adj_idx3])
    print('DIHEDRALS')
    for rec in unique_dihedrals:
        print(rec[0], rec[1], rec[2], rec[3], values[rec] / np.pi * 180)