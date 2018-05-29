import numpy as np


# Some common colors
_white = np.array([1.0, 1.0, 1.0, 1.0])
_black = np.array([0.0, 0.0, 0.0, 1.0])
_carbon = np.array([0.17, 0.17, 0.18, 1.0])
_red = np.array([0.95, 0.03, 0.01, 1.0])
_blue = np.array([0.01, 0.03, 0.95, 1.0])
_sky = np.array([0.233, 0.686, 1.0, 1.0])
_yellow = np.array([0.95, 0.95, 0.0, 1.0])
_green = np.array([0.0, 0.53, 0.0, 1.0])
_pink = np.array([0.53, 0.12, 0.36, 1.0])
_dark_red = np.array([0.59, 0.13, 0.0, 1.0])
_violet = np.array([0.46, 0.0, 1.0, 1.0])
_dark_violet = np.array([0.39, 0.0, 0.73, 1.0])
_cyan = np.array([0.0, 1.0, 1.0, 1.0])
_orange = np.array([1.0, 0.59, 0.0, 1.0])
_peach = np.array([1.0, 0.66, 0.46, 1.0])
_dark_green = np.array([0.0, 0.46, 0.0, 1.0])
_gray = np.array([0.59, 0.59, 0.59, 1.0])
_dark_orange = np.array([0.86, 0.46, 0.0, 1.0])
_golden = np.array([0.75, 0.625, 0.38, 1.0])
_brown = np.array([0.4, 0.0, 0.0, 1.0])
_turquoise = np.array([0.0, 0.78, 0.84, 1.0])
_salmon = np.array([1.0, 0.75, 0.51, 1.0])
_purple = np.array([0.35, 0.0, 0.59, 1.0])

chain_colors = [_red, _blue, _sky, _yellow, _green, _pink, _dark_red, _violet,
            _peach, _orange, _dark_violet, _cyan, _salmon, _golden, _dark_orange]

# Atom colors
_atom_colors = {'H':_white, 'C':_carbon, 'N':_sky, 'O':_red, 'S':_yellow, 'AU':_golden,
                'ZN':_violet, 'F':_green, 'CL':_green, 'BR':_brown, 'I':_dark_violet,
                'HE':_turquoise, 'NE':_turquoise, 'AR':_turquoise, 'XE':_turquoise,'KR':_turquoise,
                'P':_orange, 'B':_salmon, 'LI':_purple, 'NA':_purple, 'K':_purple, 'RB':_purple, 'CS':_purple,
                'BE':_dark_green, 'MG':_dark_green, 'SR':_dark_green, 'BA':_dark_green, 'RA':_dark_green,
                'TI':_gray, 'FE':_dark_orange}
_default_atom_color = _peach

# Atom VdW radii
_vdw_radius = {'C':1.9080, 'N':1.8240, 'O':1.6612, 'S':2.0000, 'AU':1.37,
               'H':1.0, 'ZN':1.39, 'F':1.47,'CL':1.75, 'BR':1.85, 'I':1.98,
               'HE':1.40, 'NE':1.54, 'AR':1.88, 'XE':2.16, 'KR':2.02, 'P':1.8,
               'S':1.8, 'B':1.92, 'LI':1.82, 'NA':2.27,'K':2.75,'RB':3.03,
               'CS':3.43, 'BE':1.53, 'MG':1.73, 'Ca':2.31, 'SR':2.49, 'BA':2.68,
               'RA':2.83}
_default_vdw_radius = 1.5


def centroid(coordinates):
    """Calculates the centroid of a given array"""
    length = coordinates.shape[0]
    sum_x = np.sum(coordinates[:, 0])
    sum_y = np.sum(coordinates[:, 1])
    sum_z = np.sum(coordinates[:, 2])
    return sum_x/length, sum_y/length, sum_z/length


def get_atom_color(atom_type):
    """Returns the color corresponding to this atom type"""
    try:
        return _atom_colors[atom_type]
    except:
        return _default_atom_color


def get_vdw_radius(atom_tyoe):
    """Returns the VdWaals radius of this atom type"""
    try:
        return _vdw_radius[atom_type]
    except:
        return _default_vdw_radius


def get_color_by_chain(chain_id):
    """Get a color for this chain_id"""
    return chain_colors[ord(chain_id) % 15]


def get_color_by_model(model_id):
    """Get a color for this model_id"""
    return chain_colors[model_id % 15]
