# https://mdtraj.readthedocs.io/en/stable/load_functions.html
import os
import mdtraj as md
import numpy as np
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
    traj, top = f'{trajectory}/traj_comp.xtc', f'{trajectory}/prod.gro'    
    chunk, stride, skip, atom_indices = unpack_load_params(load_params)
    traj_name = trajectory.split('/')[-1]
    data = plot_data()
    if os.path.exists(traj) and os.path.exists(top):
        print(f'{traj_name} 0 ps')
        coords = md.load(f'{trajectory}/traj_comp.xtc', top=f'{trajectory}/prod.gro', stride=stride, atom_indices=atom_indices)
        x, y = task(coords, params)
        data.add(x, y)
    return data

def iterprocess(task, trajectory:str, task_params=None, load_params={}):
    traj, top = f'{trajectory}/traj_comp.xtc', f'{trajectory}/prod.gro'
    chunk, stride, skip, atom_indices = unpack_load_params(load_params)
    traj_name = trajectory.split('/')[-1]
    data = plot_data()
    if os.path.exists(traj) and os.path.exists(top):
        idx_chunk = 0
        for coords in md.iterload(f'{trajectory}/traj_comp.xtc', top=f'{trajectory}/prod.gro', chunk=chunk, stride=stride, skip=skip, atom_indices=atom_indices):
            print(f'{traj_name} {idx_chunk * chunk / 10} ps')
            if np.shape(coords.xyz)[0] < chunk:
                break
            x, y = task(coords, task_params)
            data.add(x, y)
            idx_chunk += 1
    return data
