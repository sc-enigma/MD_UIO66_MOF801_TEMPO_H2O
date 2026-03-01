import sys
import numpy as np
import pickle

sys.path.append('../../../COMPONENTS/')

from atom import Atom

def find_linkers(atoms, forbidden_names):
    max_linker_size = 10
    ids_linkers = []
    used = np.full(len(atoms), False)
    
    for i in range(len(atoms)):
        if used[i] or atoms[i].name in forbidden_names:
            continue
        ids_linker = [i]
        for iter in range(max_linker_size):
            for j in ids_linker:
                if used[j]:
                    continue
                for k in atoms[j].adjacency:
                    if not used[k] and not atoms[k].name in forbidden_names and not k in ids_linker:
                        ids_linker.append(k)
                used[j] = True
        ids_linkers.append(ids_linker)
    return ids_linkers

def find_connected_atoms(atoms, target_names, linker_ids):
    all_ids_connected = []
    for i in linker_ids:
        ids_connected = []
        for j in atoms[i].adjacency:
            if not j in linker_ids and atoms[j].name in target_names:
                ids_connected.append(j)
        if len(ids_connected) > 0:
            all_ids_connected.append(ids_connected)
    return all_ids_connected

def select_best_r(atoms, center, radius, ids_ignored):
    best_r = center + np.array([0, 0, radius])
    best_dist = 0.0
    for thetta in np.linspace(0, np.pi, 4):
        for phi in np.linspace(0, 2 * np.pi, 7):
            r = center + np.array([radius * np.sin(thetta) * np.cos(phi), radius * np.sin(thetta) * np.sin(phi), radius * np.cos(thetta)])
            dist = 100
            for i in range(len(atoms)):
                if not i in ids_ignored:
                    dist = min(dist, np.linalg.norm(atoms[i].r - r))
            if dist > best_dist:
                best_r = r
                best_dist = dist
    return best_r
                    
def add_h(atoms, idx_site, ids_ignored):
    r_h = select_best_r(atoms, atoms[idx_site].r, 0.960, ids_ignored)
    atoms[idx_site].adjacency.append(len(atoms))
    atoms.append(Atom('H7', r_h))
    atoms[-1].resid_idx = atoms[-2].resid_idx
    atoms[-1].resid = atoms[-2].resid
    atoms[-1].atom_idx = atoms[-2].atom_idx + 1
    atoms[-1].atom_type = 'repl_H'
    atoms[-1].adjacency.append(idx_site)
    return atoms

def add_oh(atoms, idx_site, ids_ignored):
    r_o = select_best_r(atoms, atoms[idx_site].r, 2.098, ids_ignored)
    atoms[idx_site].adjacency.append(len(atoms))
    atoms.append(Atom('O7', r_o))
    atoms[-1].resid_idx = atoms[-2].resid_idx
    atoms[-1].resid = atoms[-2].resid
    atoms[-1].atom_idx = atoms[-2].atom_idx + 1
    atoms[-1].atom_type = 'repl_O'
    atoms[-1].adjacency.append(idx_site)
    
    r_h = select_best_r(atoms, atoms[-1].r, 0.960, ids_ignored)
    atoms[-1].adjacency.append(len(atoms))
    atoms.append(Atom('H7', r_h))
    atoms[-1].resid_idx = atoms[-2].resid_idx
    atoms[-1].resid = atoms[-2].resid
    atoms[-1].atom_idx = atoms[-2].atom_idx + 1
    atoms[-1].atom_type = 'repl_H'
    atoms[-1].adjacency.append(len(atoms) - 2)
    return atoms