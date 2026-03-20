from process import process, iterprocess
from plot_manager import plot_manager
from rotacf import rotacf
from utils import make_name

data_path = '/home/sc_enigma/Projects/MD_UIO66_MOF801_TEMPO_H2O/PROD/DATA/'
prod_path = '/home/sc_enigma/Projects/MD_UIO66_MOF801_TEMPO_H2O/PROD/'
trajectories = [\
    prod_path + 'UIO66+TEMPO',
    prod_path + 'UIO66+TEMPO+H2O_0.1',
    prod_path + 'UIO66+TEMPO+H2O_0.2',
    prod_path + 'UIO66+TEMPO+H2O_0.3',
    prod_path + 'UIO66+TEMPO+H2O_0.4',
    prod_path + 'UIO66+TEMPO+H2O_0.5',
    prod_path + 'UIO66+TEMPO+H2O_0.6',
    prod_path + 'UIO66+TEMPO+H2O_1.0'
    ]

# rotacf
man_rotacf = plot_manager(data_path+'rotacf', 'rotacf')
for trajectory in trajectories:
    data = iterprocess(rotacf, trajectory, [], {'atom_indices':[2020,2030]})
    man_rotacf.update(make_name(trajectory), data)
    break

man_rotacf.dump()
man_rotacf.save()