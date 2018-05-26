import numpy as np
import argparse
from glumpy import app, gloo, gl
from glumpy.transforms import Position, Trackball
from glumpy.graphics.filter import Filter
from Bio.PDB.PDBParser import PDBParser
import kendrew.atomic as atomic
from kendrew.shaders import dgoodsell, bonds


_radii_ratio = 10.0


# Get command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("pdb_file", help="PDB file name")
parser.add_argument('--ca', help='Display only CA atoms', action='store_true')
parser.add_argument('--bonds', help='Draw atomic bonds', action='store_true')
args = parser.parse_args()

# Read atomic information
pdb_file_name = args.pdb_file
parser = PDBParser(QUIET=True, PERMISSIVE=True)
structure = parser.get_structure('structure', pdb_file_name)
if args.ca:
    atoms = [atom for atom in structure.get_atoms() if atom.id == 'CA']
else:
    atoms = [atom for atom in structure.get_atoms()]
coordinates = np.array([atom.coord for atom in atoms])
center = atomic.centroid(coordinates)
coordinates -= center
coordinates *= 0.05

window = app.Window(width=1600, height=1600, color=(1,1,1,1))

protein = gloo.Program(dgoodsell.vertex, dgoodsell.fragment)
protein['light_position'] = 0., 0., 2.
protein["transform"] = Trackball(Position())

a_data = []
colors = [ atomic.get_atom_color(atom.element) for atom in atoms ]
radius = [ atomic.get_vdw_radius(atom.element) for atom in atoms ]

a_data = np.zeros(len(colors), [("position", np.float32, 3),
                                ("color",    np.float32, 4),
                                ("radius",    np.float32, 1)])

for i in range(len(colors)):
    a_data[i]['position'] = [coordinates[i][0], coordinates[i][1], coordinates[i][2]]
    a_data[i]['color'] = colors[i]
    a_data[i]['radius'] = radius[i] / _radii_ratio

protein.bind(a_data.view(gloo.VertexBuffer))
protein['color'] *= .25
protein['color'] += .75

bonds = gloo.Program(bonds.vertex, bonds.fragment)
bonds["transform"] = protein["transform"]
bonds.bind(a_data.view(gloo.VertexBuffer))
bonds['color'] = protein['color']

window.attach(protein["transform"])

@window.event
def on_draw(dt):
    window.clear()
    protein.draw(gl.GL_POINTS)
    if args.bonds:
        bonds.draw(gl.GL_LINE_STRIP)

@window.event
def on_init():
    gl.glEnable(gl.GL_DEPTH_TEST)
    gl.glEnable(gl.GL_LINE_SMOOTH)
    gl.glLineWidth(50.0)

@window.event
def on_key_press(key, modifiers):
    global protein
    if key == 65:
        protein['color'] += .01
    if key == 83:
        protein['color'] -= .01


app.run()
