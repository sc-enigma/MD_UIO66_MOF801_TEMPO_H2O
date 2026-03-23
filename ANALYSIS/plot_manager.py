import numpy as np
import pickle

from plot_data import plot_data

class plot_manager:
    def __init__(self, i_path:str, i_name:str):
        self._storage = {}
        self._path = i_path
        self._name = i_name
    def load(self):
        with open(f'{self._path}/{self._name}.pickle', 'rb') as handle:
            proxy = pickle.load(handle)
            self._storage = proxy._storage
    def save(self):
        with open(f'{self._path}/{self._name}.pickle', 'wb') as handle:
            pickle.dump(self, handle, protocol=pickle.HIGHEST_PROTOCOL)
    def update(self, i_title:str, i_data:plot_data):
        self._storage[i_title] = i_data
    def dump(self):
        for key in self._storage.keys():
            print(f'{key} x:{np.shape(self._storage[key].x)} y:{np.shape(self._storage[key].y)}')
    def get(self, i_title:str):
        for i_title in self._storage.keys():
            return self._storage[i_title]
        return None