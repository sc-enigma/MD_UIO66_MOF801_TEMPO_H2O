import numpy as np
from scipy.interpolate import splrep, BSpline
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def read_xvg(file_name):
    file = open(file_name, 'r')
    lines = [line.replace('\n', '') for line in file]
    lines = [line for line in lines if len(line) != 0 \
             if line[0] != '@' and line[0] != '#' and line[0] != '&']
    file.close()

    x = [float(line.split()[0]) for line in lines]
    y = [float(line.split()[1]) for line in lines]
    return x, y

def smooth(x, y):
    tck = splrep(x, y, s=0)
    return x, BSpline(*tck)(x)

def find_tau(x, y):
    def exp(x,a,b,c):
        return a * np.exp(-b*x) + c
    popt, pcov = curve_fit(exp, x, y)
    a, b, c = popt
    return 1/b # ps
    
def plot_tau(x, y):
    plt.cla()
    def exp(x,a,b,c):
        return a * np.exp(-b*x) + c
    plt.plot(x, y)
    plt.plot(x, [a * np.exp(-b*v) + c for v in x ])
    plt.show()

trajectories = [
    'UIO66+TEMPO+H2O_B3_0.0.xvg',
    'UIO66+TEMPO+H2O_B3_0.1.xvg',
    'UIO66+TEMPO+H2O_B3_0.2.xvg',
    'UIO66+TEMPO+H2O_B3_0.3.xvg',
    'UIO66+TEMPO+H2O_B3_0.4.xvg',
    'UIO66+TEMPO+H2O_B3_0.5.xvg',
    'UIO66+TEMPO+H2O_B3_0.6.xvg',
    ]

for trajectory in trajectories:
    x, y = read_xvg(trajectory)
    x, y = smooth(x, y)
    x, y = x[:1000], y[:1000]
    print(trajectory, find_tau(x,y))
    
    