import sys
import numpy as np
import pickle

sys.path.append('../../../COMPONENTS')
from atom import Atom, select_atoms, remove_atoms, mol2_to_atoms, count_atoms, shift_atoms, dump_atoms, count_atoms, shift_atoms
from read_utils import read_mol2_file
from write_utils import write_gro_file, write_mol2_file
from write_utils import write_atoms, write_bonds, write_angles, write_dihedrals, write_restraints, compose_itp_files
from params_detector import detect_bonds, detect_angles, detect_dihedrals

from mof801_utils import remove_extra_oxygens, define_mof801_atom_types, define_mof801_atom_names, check_mof801_atom_names
from mof801_params import get_mof801_params

# STEP 1. Read data
# Set lower bound in Mercury calculate packing = 0.0
a     = 17.8348
b     = 17.8348
c     = 17.8348
alpha = 90.00 / 180 * np.pi
beta  = 90.00 / 180 * np.pi
gamma = 90.00 / 180 * np.pi
bounds_a, bounds_b, bounds_c = [0.0, 2.5], [0.0, 2.5], [0.0, 2.5]
atoms = mol2_to_atoms(read_mol2_file('__source_mof801_444.mol2'))
atoms = shift_atoms(atoms, np.array([-4.4587, -4.4587, -4.4587]))

# STEP 2. Define atom types and names
atoms = define_mof801_atom_types(atoms)
atoms = define_mof801_atom_names(atoms)
check_mof801_atom_names(atoms)

# STEP 3. Write .gro and .mol2 files
# with open('__tmp/atoms_zif7_lp.pickle', 'wb') as handle:
#     pickle.dump(atoms, handle, protocol=pickle.HIGHEST_PROTOCOL)
write_gro_file(atoms, 'mof801.gro', a, b, c, alpha, beta, gamma, bounds_a, bounds_b, bounds_c)
write_mol2_file(atoms, 'mof801.mol2', a, b, c, alpha, beta, gamma)

# STEP 4. Write .itp files
check_mof801_atom_names(atoms)
mass, charge, bond_params, angle_params, dihedral_params = get_mof801_params()
write_atoms(atoms, charge, mass, 'MOF', 'atoms.itp')
write_bonds(atoms, bond_params, 'bonds.itp')
write_angles(atoms, angle_params, 'angles.itp')
write_dihedrals(atoms, dihedral_params, 'dihedrals.itp')
write_restraints(atoms, ['itcFF_Zr'], 'restraints.itp', 5000)
compose_itp_files(['moleculetype.itp', 'atoms.itp', 'bonds.itp', 'angles.itp', 'dihedrals.itp', 'restraints.itp'], 'mof801.itp')

