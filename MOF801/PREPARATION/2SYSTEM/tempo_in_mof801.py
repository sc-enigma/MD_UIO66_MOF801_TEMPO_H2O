import sys
import numpy as np
import pickle

sys.path.append('../../../COMPONENTS')
from atom import Atom, select_atoms, remove_atoms, mol2_to_atoms, count_atoms, shift_atoms, dump_atoms, count_atoms, shift_atoms, set_residue, add_atoms
from read_utils import read_mol2_file
from write_utils import write_gro_file, write_mol2_file
from write_utils import write_atoms, write_bonds, write_angles, write_dihedrals, write_restraints, compose_itp_files

a     = 17.8348
b     = 17.8348
c     = 17.8348
alpha = 90.00 / 180 * np.pi
beta  = 90.00 / 180 * np.pi
gamma = 90.00 / 180 * np.pi
bounds_a, bounds_b, bounds_c = [0.0, 2.5], [0.0, 2.5], [0.0, 2.5]

atoms_mof801 = mol2_to_atoms(read_mol2_file('../2TOP/mof801.mol2'))
atoms_tempo = mol2_to_atoms(read_mol2_file('tempo.mol2'))
atoms_mof801 = set_residue(atoms_mof801, "MOF")
atoms_tempo = set_residue(atoms_tempo, "TMP")

atoms_tempo = shift_atoms(atoms_tempo, 0.5 * (atoms_mof801[1770].r + atoms_mof801[2272].r) - atoms_tempo[0].r)

atoms = add_atoms(atoms_mof801, atoms_tempo)

write_mol2_file(atoms, 'mof801_tempo.mol2', a, b, c, alpha, beta, gamma)
write_gro_file(atoms, 'mof801_tempo.gro', a, b, c, alpha, beta, gamma)

