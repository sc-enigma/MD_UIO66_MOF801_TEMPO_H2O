import sys
import numpy as np
import pickle

sys.path.append('../ITP')
from atom import Atom

def add_unique(vecs, new_vec):
    for vec in vecs:
        if np.linalg.norm(new_vec - vec) < 1e-6:
            return vecs
    vecs.append(new_vec)
    return vecs
        
def add_hydrogens(atoms):
    x_shift1 = np.array([-8.9174, -0.4238, 8.9174])
    x_shift2 = np.array([8.9174, -0.4238, 8.9174])
    y_shift1 = np.array([0.397002, -8.9174, 8.9174])
    y_shift2 = np.array([0.397002, 8.9174, 8.9174])
    z_shift1 = np.array([0, 0, -17.8348])
    
    oxygens_in_one_cluster = np.array([
        [25.387800, 28.387699, 27.906099],
        [27.906099, 25.387800, 28.387699],
        [28.116600, 25.116699, 25.598301],
        [25.598301, 28.116600, 25.116699]
        ])
    
    shifts = []
    for idx_x_shift1 in range(-2, 2):
        for idx_x_shift2 in range(-2, 2):
            for idx_y_shift1 in range(-2, 2):
                for idx_y_shift2 in range(-2, 2):
                    for idx_z_shift1 in range(-2, 3):
                        new_shift = x_shift1 * idx_x_shift1 + x_shift2 * idx_x_shift2 + y_shift1 * idx_y_shift1 + y_shift2 * idx_y_shift2 + z_shift1 * idx_z_shift1 
                        shifts = add_unique(shifts, new_shift)
    
    ids_selected = []
    for idx in range(len(atoms)):
        if atoms[idx].mol2name == 'O3':
            for o_atom in oxygens_in_one_cluster:
                for shift in shifts:
                    if np.linalg.norm(o_atom + shift - atoms[idx].r) < 1.0:
                        ids_selected.append(idx)
    ids_selected = list(set(ids_selected))
    ids_selected = sorted(ids_selected, reverse=True)
    
    for idx in ids_selected:
        bestShiftVector = np.zeros(3)
        bestDistance = sys.float_info.max
        for thetta in np.linspace(0.0, np.pi, 10):
            for phi in np.linspace(0.0, 2 * np.pi, 20):
                shiftVector = np.array([1.09 * np.sin(thetta) * np.cos(phi), 1.09 * np.sin(thetta) * np.sin(phi), 1.09 * np.cos(thetta)])
                distance = 0.0
                for atom in atoms:
                    if atom.mol2name == 'Zr':
                        distance += 1 / np.linalg.norm(atoms[idx].r + shiftVector - atom.r)
                if distance < bestDistance:
                    bestShiftVector = shiftVector
                    bestDistance = distance
                    
        atoms.append(Atom('H7', atoms[idx].r + bestShiftVector))
        atoms[-1].atom_idx = len(atoms) - 1
        atoms[-1].mol2name = 'H7'
        atoms[-1].adjacency.append(idx)
        atoms[idx].adjacency.append(len(atoms) - 1)
        
    return atoms
