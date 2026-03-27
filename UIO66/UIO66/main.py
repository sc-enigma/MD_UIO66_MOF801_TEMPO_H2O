import sys
import numpy as np
import pickle
import random
import itertools

sys.path.append('../../COMPONENTS/')
from atom import Atom, mol2_to_atoms, count_atoms, remove_non_bonded_atoms, remove_atoms, shift_atoms, remove_duplicates
from linker_utils import find_linkers, find_connected_atoms, add_h, add_oh
from read_utils import read_mol2_file
from write_utils import write_gro_file, write_mol2_file
from write_utils import write_atoms, write_bonds, write_angles, write_dihedrals, write_restraints, compose_itp_files

from uio66_pbc import apply_pbc, calculate_lattice_vectors
from uio66_utils import define_uio66_atom_types, define_uio66_atom_names, remove_Zr_Zr_bond, remove_periodic_bonds, select_atoms
from uio66_params import get_uio66_params

a     = 20.7465
b     = 20.7465
c     = 20.7465
alpha = 90.00 / 180 * np.pi
beta  = 90.00 / 180 * np.pi
gamma = 90.00 / 180 * np.pi
bounds_a = [0.0, 2.0]
bounds_b = [0.0, 2.0]
bounds_c = [0.0, 2.0]

# STEP 1. Preprocess data
### CASHING START ###
'''
atoms = mol2_to_atoms(read_mol2_file('__uio66_source.mol2'))
atoms = remove_duplicates(atoms)
atoms = remove_Zr_Zr_bond(atoms)
atoms = remove_non_bonded_atoms(atoms)
atoms = remove_periodic_bonds(atoms)
atoms = define_uio66_atom_types(atoms)
atoms = define_uio66_atom_names(atoms)
write_mol2_file(atoms, '__tmp/uio66_source_processed1.mol2', 38, 38, 38, alpha, beta, gamma, True)
'''
### CASHING END ###

### SELECTION START ###
'''
x < 44 and x > 7 and y < 44 and y > 7 and z < 44 and z > 7 and type Zr
or
x < 44 and x > 7 and y < 44 and y > 7 and z < 44 and z > 7 and  type O3
or
x > 9.5 and x < 43 and y > 9.5 and y < 43 and z > 9.5 and z < 43
-> __tmp/uio66_source_processed2.mol2
'''
### SELECTION END ###

### CASHING START ###
'''
atoms = mol2_to_atoms(read_mol2_file('__tmp/uio66_source_processed2.mol2'))
atoms = shift_atoms(atoms, np.array([-6.89, -6.89, -6.89]))
atoms = define_uio66_atom_types(atoms)
atoms = define_uio66_atom_names(atoms)
with open('__tmp/atoms_uio66.pickle', 'wb') as handle:
    pickle.dump(atoms, handle, protocol=pickle.HIGHEST_PROTOCOL)
'''
### CASHING END ###
with open('__tmp/atoms_uio66.pickle', 'rb') as handle:
    atoms = pickle.load(handle)

# STEP 2. Remove some linkers
linkers = find_linkers(atoms, ['Zr', 'O1', 'O3'], 12)
ids_linkers_to_remove = random.sample(range(0, len(linkers)), int(0.25 * len(linkers)))

# a. Small pore - center (20.772, 20.772, 21.6001)
# May be removed [0, 2, 6, 73, 76, 102]
'''for idx in [0, 2, 6, 73, 76, 102]:
    if idx in ids_linkers_to_remove:
        ids_linkers_to_remove.remove(idx)
for idx in [73, 76, 102]:
    if not idx in ids_linkers_to_remove:
        ids_linkers_to_remove.append(idx)
ids_atoms_to_remove = list(itertools.chain.from_iterable([linkers[i] for i in ids_linkers_to_remove]))
'''
# b. Big pore - center (12.3761, 12.3761, 12.3761)
'''
# May be removed [0, 6, 73, 60]
# Must remain [1, 9, 16, 17, 18, 19, 61, 74]
for idx in [1, 9, 16, 17, 18, 19, 61, 74, 60]:
    if idx in ids_linkers_to_remove:
        ids_linkers_to_remove.remove(idx)
for idx in [0, 6, 73]:
    if not idx in ids_linkers_to_remove:
        ids_linkers_to_remove.append(idx)
ids_atoms_to_remove = list(itertools.chain.from_iterable([linkers[i] for i in ids_linkers_to_remove]))
'''

# Add h / oh groups
for i in ids_linkers_to_remove:
    connected_oxygens = find_connected_atoms(atoms, ['O1'], linkers[i])
    atoms = add_h(atoms, connected_oxygens[0][0], ids_atoms_to_remove)
    atoms = add_h(atoms, connected_oxygens[1][0], ids_atoms_to_remove)
    atoms = add_oh(atoms, connected_oxygens[0][1], ids_atoms_to_remove)
    atoms = add_oh(atoms, connected_oxygens[1][1], ids_atoms_to_remove)
atoms = remove_atoms(atoms, ids_atoms_to_remove)

# STEP 3. Write .gro and .mol2 files
write_gro_file(atoms, 'uio66.gro', 38, 38, 38, alpha, beta, gamma)
write_mol2_file(atoms, 'uio66.mol2', 38, 38, 38, alpha, beta, gamma, True)

# STEP 4. Write .itp files
mass, charge, bond_params, angle_params, dihedral_params = get_uio66_params()
write_atoms(atoms, charge, mass, 'UIO', 'atoms.itp')
write_bonds(atoms, bond_params, 'bonds.itp')
write_angles(atoms, angle_params, 'angles.itp')
write_dihedrals(atoms, dihedral_params, 'dihedrals.itp')
write_restraints(atoms, ['flexFF_Zr'], 'restraints.itp', 5000)
compose_itp_files(['moleculetype.itp', 'atoms.itp', 'bonds.itp', 'angles.itp', 'dihedrals.itp', 'restraints.itp'], 'uio66.itp')
