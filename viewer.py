import numpy as np
from Bio.PDB.PDBParser import PDBParser

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

from vispy import gloo
from vispy import app
from vispy.util.transforms import perspective, translate, rotate
from vispy.io import load_data_file


class Viewer(object):
    def __init__(self, pdb_file_name):
        self._parser = PDBParser(QUIET=True, PERMISSIVE=True)
        self.structure = self._parser.get_structure('structure', pdb_file_name)

    def __str__(self):
        return str(self.structure)

    def show(self):
        return NotImplementedError()


class Viewer2D(Viewer):

    visualization_modes = {'ca':'_ca', 'cpk':'_cpk'}

    def __init__(self, pdb_file_name, mode='ca'):
        super(Viewer2D, self).__init__(pdb_file_name)
        if mode not in Viewer2D.visualization_modes.keys():
            raise Exception('Not recognized visualization mode %s' % mode)
        self.mode = mode

    def show(self):
        fig = plt.figure()
        getattr(self, Viewer2D.visualization_modes[self.mode])(fig)
        plt.show()

    def _ca(self, fig):
        """ Draws CA atoms linked by lines """
        ax = fig.add_subplot(111, projection='3d')
        for chain in self.structure.get_chains():
            ca_atoms = [atom for atom in chain.get_atoms() if atom.get_name() == 'CA']
            ca_coordinates = np.array([atom.coord for atom in ca_atoms])
            xs = ca_coordinates[...,0]
            ys = ca_coordinates[...,1]
            zs = ca_coordinates[...,2]
            chain_color = np.random.rand(3,1)
            ax.scatter(xs, ys, zs, c=chain_color, marker='o')
            for i, atom in enumerate(ca_atoms[:-1]):
                a1 = ca_atoms[i].coord
                a2 = ca_atoms[i+1].coord
                ax.plot([a1[0], a2[0]], [a1[1], a2[1]], [a1[2], a2[2]], c=chain_color)

    cpk = {'C':'black', 'N':'blue', 'O':'red'}
    def _cpk(self, fig):
        """ Draws atoms using CPK notation """
        ax = fig.add_subplot(111, projection='3d')
        for atom_type, cpk_color in Viewer2D.cpk.iteritems():
            atoms = [atom for atom in self.structure.get_atoms() if atom.get_id() == atom_type]
            coordinates = np.array([atom.coord for atom in atoms])
            xs = coordinates[...,0]
            ys = coordinates[...,1]
            zs = coordinates[...,2]
            ax.scatter(xs, ys, zs, c=cpk_color, marker='o')

