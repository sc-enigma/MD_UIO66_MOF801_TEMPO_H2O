import sys
import numpy as np
import pickle

sys.path.append('../../COMPONENTS')
from atom import Atom, select_atoms, remove_atoms, mol2_to_atoms, count_atoms, shift_atoms, dump_atoms, count_atoms, shift_atoms
from read_utils import read_mol2_file
from write_utils import write_gro_file, write_mol2_file
from write_utils import write_atoms, write_bonds, write_angles, write_dihedrals, compose_itp_files
from pbc import calculate_lattice_vectors, select_in_box
from add_hydrogens import add_hydrogens

# STEP 1. Read data
# Set lower bound in Mercury calculate packing = 0.0
a     = 17.8348
b     = 17.8348
c     = 17.8348
alpha = 90.00 / 180 * np.pi
beta  = 90.00 / 180 * np.pi
gamma = 90.00 / 180 * np.pi
atoms = mol2_to_atoms(read_mol2_file('__source.mol2'))

# STEP 2. Select atoms
ids_selected = []
# For 2x2x2
for idx in range(len(atoms)):
    if atoms[idx].r[0] > 17 and atoms[idx].r[0] < 37 and atoms[idx].r[1] > 17 and atoms[idx].r[1] < 37 and atoms[idx].r[2] > 17 and atoms[idx].r[2] < 37:
        ids_selected.append(idx)
    if atoms[idx].mol2name == 'Zr' and atoms[idx].r[0] > 12 and atoms[idx].r[0] < 42 and atoms[idx].r[1] > 12 and atoms[idx].r[1] < 42 and atoms[idx].r[2] > 12 and atoms[idx].r[2] < 42:
        ids_selected.append(idx)
    if atoms[idx].mol2name == 'O3' and atoms[idx].r[0] > 12 and atoms[idx].r[0] < 42 and atoms[idx].r[1] > 12 and atoms[idx].r[1] < 42 and atoms[idx].r[2] > 12 and atoms[idx].r[2] < 42:
        ids_selected.append(idx)
# For 4x4x4
'''
for idx in range(len(atoms)):
    if atoms[idx].r[0] > 8.0826 and atoms[idx].r[0] < 45.9174 and atoms[idx].r[1] > 8.0826 and atoms[idx].r[1] < 45.9174 and atoms[idx].r[2] > 8.0826 and atoms[idx].r[2] < 45.9174:
        ids_selected.append(idx)
    if atoms[idx].mol2name == 'Zr' and atoms[idx].r[0] > 3.0826 and atoms[idx].r[0] < 50.9174 and atoms[idx].r[1] > 3.0826 and atoms[idx].r[1] < 50.9174 and atoms[idx].r[2] > 3.0826 and atoms[idx].r[2] < 50.9174:
        ids_selected.append(idx)
    if atoms[idx].mol2name == 'O3' and atoms[idx].r[0] > 3.0826 and atoms[idx].r[0] < 50.9174 and atoms[idx].r[1] > 3.0826 and atoms[idx].r[1] < 50.9174 and atoms[idx].r[2] > 3.0826 and atoms[idx].r[2] < 50.9174:
        ids_selected.append(idx)
'''

atoms = remove_atoms(atoms, list(set(list(np.arange(len(atoms)))) ^ set(ids_selected)))

# STEP 3. Add hydrogens
atoms = add_hydrogens(atoms)

# For pacmof2
ids_selected = []
for idx in range(len(atoms)):
    if atoms[idx].r[0] > 17 and atoms[idx].r[0] < 36 and atoms[idx].r[1] > 17 and atoms[idx].r[1] < 36 and atoms[idx].r[2] > 17 and atoms[idx].r[2] < 36:
        ids_selected.append(idx)
atoms = remove_atoms(atoms, list(set(list(np.arange(len(atoms)))) ^ set(ids_selected)))
atoms = shift_atoms(atoms, np.array([17.8348, 17.8348, 17.8348]))

# STEP 4. Write data
# write_mol2_file(atoms, 'mof801_222.mol2', a, b, c, alpha, beta, gamma)
# write_mol2_file(atoms, 'mof801_444.mol2', a, b, c, alpha, beta, gamma)
write_mol2_file(atoms, 'mof801_222_pacmof.mol2', a, b, c, alpha, beta, gamma)
