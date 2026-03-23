# 10 frames = 1 ps, 10000 frames = 1 ns
import numpy as np
from scipy.interpolate import splrep, BSpline
from numba import njit

import matplotlib.pyplot as plt

@njit
def rotacf(x, y, dt_max):
    j_max = np.shape(y)[0]
    for j in range(j_max):
        y[j] = y[j] / np.linalg.norm(y[j])
    out = np.zeros(dt_max)
    for dt in range(0, dt_max):
        i_max = np.shape(y)[0] - dt*10 - 1
        for i in range(i_max):
            out[dt] += y[i].dot(y[i+dt*10]) / i_max
    out = out / out[0]
    return out

def task_rotacf(coords, params):
    v1 = coords.xyz[:,0,:] - coords.xyz[:,1,:]
    v2 = coords.xyz[:,1,:] - coords.xyz[:,2,:]
    n = np.cross(v1, v2)
    x = np.arange(np.shape(n)[0])
    s = 0
    tck = [splrep(x, n[:,0], s=s), splrep(x, n[:,1], s=s), splrep(x, n[:,2], s=s)]
    n_smoothed = np.array([BSpline(*tck[0])(x), BSpline(*tck[1])(x), BSpline(*tck[2])(x)]).T
    t_max = min(int(np.shape(n)[0] / 25), 5000) # 5000 ns
    t = np.arange(0, t_max)
    r = rotacf(x, n_smoothed, t_max)
    return t, r