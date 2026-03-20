# https://mdtraj.readthedocs.io/en/stable/load_functions.html
import mdtraj as md
from plot_data import plot_data

def unpack_load_params(load_params):
    chunk = 1000
    stride = None
    skip = 0
    atom_indices = []
    if 'chunk' in load_params:
        chunk = load_params['chunk']
    if 'stride' in load_params:
        stride = load_params['stride']
    if 'skip' in load_params:
        skip = load_params['skip']
    if 'atom_indices' in load_params:
        atom_indices = load_params['atom_indices']
    return chunk, stride, skip, atom_indices

def process(task, trajectory:str, task_params=None, load_params={}):
    chunk, stride, skip, atom_indices = unpack_load_params(load_params)
    data = plot_data()
    chunk = md.load(f'{trajectory}/traj_comp.xtc', top=f'{trajectory}/prod.gro', stride=stride, atom_indices=atom_indices)
    data.add(task(chunk, params))
    return data

def iterprocess(task, trajectory:str, task_params=None, load_params={}):
    chunk, stride, skip, atom_indices = unpack_load_params(load_params)
    data = plot_data()
    for chunk in md.iterload(f'{trajectory}/traj_comp.xtc', top=f'{trajectory}/prod.gro', chunk=chunk, stride=stride, skip=skip, atom_indices=atom_indices):
        data.add(task(chunk, task_params))
    return data
