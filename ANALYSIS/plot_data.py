import numpy as np

class plot_data:
    def __init__(self, i_x=[], i_y=[]):
        self.x = np.array(i_x)
        self.y = np.array(i_y)
        self.n = int(len(self.y) > 0)
    def append(self, i_y):
        self.y = np.append(self.y, np.array(i_y))
        self.n = int(len(self.y) > 0)
    def add(self, i_y):
        if self.n == 0:
            self.y = np.array(i_y)
            self.n = 1
        else:
            assert(len(self.y) == len(i_y))
            self.y = (self.y*self.n + np.array(i_y)) / (self.n+1)
            self.n += 1

