import os
import time as tm
from datetime import datetime

prod_dir = '/home/sc_enigma/Projects/MD_UIO66_MOF801_TEMPO_H2O/PROD/'

def run_md(path):
    if not os.path.exists(path):
        return 
    os.chdir(path)
    os.system('gmx_mpi mdrun -s *.tpr -v -o prod -c prod')
    while True: 
        tm.sleep(1)
        if 'prod.gro' in os.listdir('.'):
            break
    tm.sleep(10)

def write_log(path):
    now = datetime.now()
    current_time = now.strftime("%m/%d/%Y, %H:%M:%S")
    with open(prod_dir + 'log.txt', 'a') as file:
        file.write(f'{current_time} gmx_mpi mdrun {path}\n')

path_list = [\
    prod_dir + 'UIO66+TEMPO/0/',
    prod_dir + 'MOF801+TEMPO/0/'
    ]

for path in path_list:
    write_log(path)
    run_md(path)
