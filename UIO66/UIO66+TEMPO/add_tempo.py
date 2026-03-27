import sys
import numpy as np
import pickle

sys.path.append('../../COMPONENTS/')
from atom import Atom, set_residue, mol2_to_atoms, shift_atoms, add_atoms
from read_utils import read_mol2_file
from write_utils import write_gro_file, write_mol2_file

a     = 38.0
b     = 38.0
c     = 38.0
alpha = 90.00 / 180 * np.pi
beta  = 90.00 / 180 * np.pi
gamma = 90.00 / 180 * np.pi

atoms_mof801 = mol2_to_atoms(read_mol2_file('uio66.mol2'))
atoms_tempo = mol2_to_atoms(read_mol2_file('tempo/tempo.mol2'))

atoms_mof801 = set_residue(atoms_mof801, "UIO")
atoms_tempo = set_residue(atoms_tempo, "TMP")
# a. Small pore - center (20.872, 20.872, 20.0001)
# b. Big pore - center (12.3761, 12.3761, 12.3761)
atoms_tempo = shift_atoms(atoms_tempo, np.array([12.3761, 12.3761, 12.3761]) - atoms_tempo[0].r)

atoms = add_atoms(atoms_mof801, atoms_tempo)

write_mol2_file(atoms, 'uio66_tempo.mol2', a, b, c, alpha, beta, gamma)
write_gro_file(atoms, 'uio66_tempo.gro', a, b, c, alpha, beta, gamma)
