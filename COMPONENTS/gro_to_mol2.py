import sys
import numpy as np
import pickle

sys.path.append('../COMPONENTS')
from atom import Atom, select_atoms, remove_atoms, gro_to_atoms, mol2_to_atoms, count_atoms, shift_atoms, dump_atoms, count_atoms, shift_atoms
from read_utils import read_mol2_file, read_gro_file
from write_utils import write_gro_file, write_mol2_file
from write_utils import write_atoms, write_bonds, write_angles, write_dihedrals, write_restraints, compose_itp_files
from params_detector import detect_bonds, detect_angles, detect_dihedrals

a     = 17.8348
b     = 17.8348
c     = 17.8348
alpha = 90.00 / 180 * np.pi
beta  = 90.00 / 180 * np.pi
gamma = 90.00 / 180 * np.pi
bounds_a, bounds_b, bounds_c = [0.0, 2.5], [0.0, 2.5], [0.0, 2.5]

atoms_from_mol2 = mol2_to_atoms(read_mol2_file('../TOP/MOF801/mof801.mol2'))
atoms_from_gro = gro_to_atoms(read_gro_file('../TOP/MOF801/mof801.gro'))

for atom_idx in range(len(atoms_from_mol2)):
    atoms_from_mol2[atom_idx].r = atoms_from_gro[atom_idx].r * 10.0

write_mol2_file(atoms_from_mol2, 'abc.mol2', a, b, c, alpha, beta, gamma)