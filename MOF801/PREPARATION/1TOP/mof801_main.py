import sys
import numpy as np
import pickle
import random
import itertools

sys.path.append('../../../COMPONENTS')
from atom import Atom, select_atoms, remove_atoms, mol2_to_atoms, count_atoms, shift_atoms, dump_atoms, count_atoms, shift_atoms
from linker_utils import find_linkers, find_connected_atoms, add_h, add_oh
from read_utils import read_mol2_file
from write_utils import write_gro_file, write_mol2_file
from write_utils import write_atoms, write_bonds, write_angles, write_dihedrals, write_restraints, compose_itp_files
from params_detector import detect_bonds, detect_angles, detect_dihedrals

from mof801_utils import remove_extra_oxygens, define_mof801_atom_types, define_mof801_atom_names, check_mof801_atom_names
from mof801_params import get_mof801_params

a     = 17.8348
b     = 17.8348
c     = 17.8348
alpha = 90.00 / 180 * np.pi
beta  = 90.00 / 180 * np.pi
gamma = 90.00 / 180 * np.pi
bounds_a = [0.0, 2.5]
bounds_b = [0.0, 2.5]
bounds_c = [0.0, 2.5]

# STEP 1. Preprocess data
### CASHING START ###
'''
atoms = mol2_to_atoms(read_mol2_file('__source_mof801_444.mol2'))
atoms = shift_atoms(atoms, np.array([-4.4587, -4.4587, -4.4587]))
atoms = define_mof801_atom_types(atoms)
atoms = define_mof801_atom_names(atoms)
check_mof801_atom_names(atoms)
write_mol2_file(atoms, '__tmp/mof801_source_processed1.mol2', 44.5870, 44.5870, 44.5870, alpha, beta, gamma, True)
'''
### CASHING END ###

### SELECTION START ###
'''
x < 32.5 and y < 32.5 and z < 32.5
or
x < 36 and y < 36 and z < 36 and type O2
or
x < 36 and y < 36 and z < 36 and type H7
or
x < 36 and y < 36 and z < 36 and type Zr
-> __tmp/mof801_source_processed2.mol2
'''
### SELECTION END ###

### CASHING START ###
'''
atoms = mol2_to_atoms(read_mol2_file('__tmp/mof801_source_processed2.mol2'))
atoms = define_mof801_atom_types(atoms)
atoms = define_mof801_atom_names(atoms)
atoms = shift_atoms(atoms, np.array([-1,-1,-1]))
with open('__tmp/atoms_mof801.pickle', 'wb') as handle:
    pickle.dump(atoms, handle, protocol=pickle.HIGHEST_PROTOCOL)
'''
### CASHING END ###
with open('__tmp/atoms_mof801.pickle', 'rb') as handle:
    atoms = pickle.load(handle)
    
# STEP 2. Remove some linkers
linkers = find_linkers(atoms, ['Zr', 'O1', 'O2'], 6)
ids_linkers_to_remove = np.unique([random.randint(0, len(linkers)-1) for _ in range(int(len(linkers) * 0.25))])
ids_atoms_to_remove = list(itertools.chain.from_iterable([linkers[i] for i in ids_linkers_to_remove]))
for i in ids_linkers_to_remove:
    connected_oxygens = find_connected_atoms(atoms, ['O1'], linkers[i])
    atoms = add_h(atoms, connected_oxygens[0][0], ids_atoms_to_remove)
    atoms = add_h(atoms, connected_oxygens[1][0], ids_atoms_to_remove)
    atoms = add_oh(atoms, connected_oxygens[0][1], ids_atoms_to_remove)
    atoms = add_oh(atoms, connected_oxygens[1][1], ids_atoms_to_remove)
atoms = remove_atoms(atoms, ids_atoms_to_remove)

# STEP 3. Write .gro and .mol2 files
write_gro_file(atoms, 'mof801.gro', 33.67, 33.67, 33.67, alpha, beta, gamma)
write_mol2_file(atoms, 'mof801.mol2', 33.67, 33.67, 33.67, alpha, beta, gamma)

# STEP 4. Write .itp files
check_mof801_atom_names(atoms)
mass, charge, bond_params, angle_params, dihedral_params = get_mof801_params()
write_atoms(atoms, charge, mass, 'MOF', 'atoms.itp')
write_bonds(atoms, bond_params, 'bonds.itp')
write_angles(atoms, angle_params, 'angles.itp')
write_dihedrals(atoms, dihedral_params, 'dihedrals.itp')
write_restraints(atoms, ['itcFF_Zr'], 'restraints.itp', 5000)
compose_itp_files(['moleculetype.itp', 'atoms.itp', 'bonds.itp', 'angles.itp', 'dihedrals.itp', 'restraints.itp'], 'mof801.itp')

