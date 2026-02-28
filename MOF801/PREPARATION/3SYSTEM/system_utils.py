import numpy as np

def optimize(atoms, other_atoms, delta, steps = 5):
    r = np.array([atom.r for atom in atoms])
    other_r = np.array([atom.r for atom in other_atoms])
    best_shift = np.zeros(3)
    best_dist = 0.0
    for dx in np.linspace(-delta, delta, steps):
        for dy in np.linspace(-delta, delta, steps):
            for dz in np.linspace(-delta, delta, steps):
                shift = np.array([dx, dy, dz])
                r += shift
                dist = 1e10
                for ri in r:
                    for rj in other_r:
                        dist = min(dist, np.linalg.norm(ri - rj))
                r -= shift
                if dist > best_dist:
                    best_dist = dist
                    best_shift = shift
    return best_shift
    