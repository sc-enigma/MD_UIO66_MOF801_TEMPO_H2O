import sys
import numpy as np
import pickle

sys.path.append('../../../COMPONENTS/')
from atom import Atom, set_residue, mol2_to_atoms, shift_atoms, add_atoms
from read_utils import read_mol2_file
from write_utils import write_gro_file, write_mol2_file

a     = 33.6700
b     = 33.6700
c     = 33.6700
alpha = 90.00 / 180 * np.pi
beta  = 90.00 / 180 * np.pi
gamma = 90.00 / 180 * np.pi

atoms_mof801 = mol2_to_atoms(read_mol2_file('mof801/3_D/mof801.mol2'))
atoms_tempo = mol2_to_atoms(read_mol2_file('tempo/tempo.mol2'))

atoms_mof801 = set_residue(atoms_mof801, "MOF")
atoms_tempo = set_residue(atoms_tempo, "TMP")
atoms_tempo = shift_atoms(atoms_tempo, np.array([12.3761, 12.3761, 12.3761]) - atoms_tempo[0].r)

atoms = add_atoms(atoms_mof801, atoms_tempo)

write_mol2_file(atoms, 'mof801+tempo/3_D/mof801_tempo.mol2', a, b, c, alpha, beta, gamma)
write_gro_file(atoms, 'mof801+tempo/3_D/mof801_tempo.gro', a, b, c, alpha, beta, gamma)
