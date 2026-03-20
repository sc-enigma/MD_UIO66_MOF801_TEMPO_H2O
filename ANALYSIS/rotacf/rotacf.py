import math
import numpy as np
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

def truncate_by_y_val(x, y, y_val):
    i = 0
    while i + 1 < len(y):
        i += 1
        if y[i] < y_val:
            break
    return x[:i], y[:i]

def exp_fit(x, y, min_y_val, max_y_val):
    x_for_opt, y_for_opt = [], []
    search_threshold = max_y_val
    for i in range(len(y)):
        if y[i] < search_threshold:
            x_for_opt.append(x[i])
            y_for_opt.append(y[i])
            search_threshold -= 0.1
        if search_threshold <= min_y_val:
            break
    plt.scatter(x_for_opt, y_for_opt, color='k')
    def exp_func(x, a, b):
        return a * np.exp(-b * x)
    params, covariance = curve_fit(exp_func, x_for_opt, y_for_opt)
    return params
    # b = math.log(y_for_opt[0] / y_for_opt[-1]) / (x_for_opt[-1] - x_for_opt[0])
    # a = y_for_opt[0] / np.exp(-b * x_for_opt[0])
    # return [a, b]

x, y = read_xvg('test.xvg')
x, y = truncate_by_y_val(x, y, 0.1)
params = exp_fit(x, y, 0.1, 0.8)
y_exp = [params[0] * np.exp(-params[1] * val) for val in x]

t_corr = 1 / params[1] / 1000
print(f't_corr = {t_corr} ns')

plt.plot(x, y, color='k')
plt.plot(x, y_exp, color='r')
plt.show()
