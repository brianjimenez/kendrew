import numpy as np


# Some common colors
white = np.array([1.0, 1.0, 1.0, 1.0])
black = np.array([0.0, 0.0, 0.0, 1.0])
carbon = np.array([0.17, 0.17, 0.18, 1.0])
red = np.array([0.95, 0.03, 0.01, 1.0])
blue = np.array([0.01, 0.03, 0.95, 1.0])
sky = np.array([0.233, 0.686, 1.0, 1.0])
yellow = np.array([0.95, 0.95, 0.0, 1.0])
green = np.array([0.0, 0.53, 0.0, 1.0])
pink = np.array([0.53, 0.12, 0.36, 1.0])
darkRed = np.array([0.59, 0.13, 0.0, 1.0])
violet = np.array([0.46, 0.0, 1.0, 1.0])
darkviolet = np.array([0.39, 0.0, 0.73, 1.0])
cyan = np.array([0.0, 1.0, 1.0, 1.0])
orange = np.array([1.0, 0.59, 0.0, 1.0])
peach = np.array([1.0, 0.66, 0.46, 1.0])
darkGreen = np.array([0.0, 0.46, 0.0, 1.0])
gray = np.array([0.59, 0.59, 0.59, 1.0])
darkorange = np.array([0.86, 0.46, 0.0, 1.0])
golden = np.array([0.75, 0.625, 0.38, 1.0])

# Atom colors
atom_colors = {'H':white, 'C':carbon, 'N':sky, 'O':red, 'S':yellow, 'AU':golden, 'ZN':violet}

# Atom VdW radii
vdw_radius = {'C':1.9080, 'N':1.8240, 'O':1.6612, 'S':2.0000, 'AU':1.37, 'H':1.0, 'ZN':1.39}


def centroid(arr):
    """Calculates the centroid of a given array"""
    length = arr.shape[0]
    sum_x = np.sum(arr[:, 0])
    sum_y = np.sum(arr[:, 1])
    sum_z = np.sum(arr[:, 2])
    return sum_x/length, sum_y/length, sum_z/length
