import sys
import numpy as np
from glumpy import app, gloo, gl
from glumpy.transforms import Position, Trackball
from glumpy.graphics.filter import Filter
from Bio.PDB.PDBParser import PDBParser
from atomic import *


only_CA = False
draw_bonds = False
radii_ratio = 10.0


vertex = """
uniform vec3 light_position;

attribute vec3 position;
attribute vec3 color;
attribute float radius;

varying float v_size;
varying vec3 v_color;
varying float v_radius;
varying vec4 v_eye_position;
varying vec3 v_light_direction;

void main (void)
{
    v_color = color;
    v_radius = radius;
    v_eye_position = <transform.trackball_view> *
                     <transform.trackball_model> *
                     vec4(position,1.0);
    v_light_direction = normalize(light_position);
    gl_Position = <transform(position)>;
    vec4 p = <transform.trackball_projection> *
             vec4(radius, radius, v_eye_position.z, v_eye_position.w);
    v_size = 512.0 * p.x / p.w;
    gl_PointSize = v_size;
}
"""

fragment = """
#include "antialias/outline.glsl"


varying float v_size;
varying vec3 v_color;
varying float v_radius;
varying vec4 v_eye_position;
varying vec3 v_light_direction;

void main()
{
    vec2 P = gl_PointCoord.xy - vec2(0.5,0.5);
    float point_size = v_size  + 5.0;
    float distance = length(P*point_size) - v_size/2;

    vec2 texcoord = gl_PointCoord* 2.0 - vec2(1.0);
    float x = texcoord.x;
    float y = texcoord.y;
    float d = 1.0 - x*x - y*y;

    if (d <= 0.0) discard;

    float z = sqrt(d);
    vec4 pos = v_eye_position;
    pos.z += v_radius*z;
    vec3 pos2 = pos.xyz;

    pos = <transform.trackball_projection> * pos;
    gl_FragDepth = 0.5*(pos.z / pos.w)+0.5;
    vec3 normal = vec3(x,y,z);
    float diffuse = clamp(dot(normal, v_light_direction), 0.0, 1.0);

    vec4 color = vec4((0.5 + 0.5*diffuse)*v_color, 1.0);
    gl_FragColor = outline(distance, 1.0, 1.0, vec4(0,0,0,1), color);
    // gl_FragColor = color;
}
"""


line_shader = """
#include "math/constants.glsl"

attribute vec3 position;
attribute vec3 color;

varying float v_size;
varying vec3 v_color;
varying vec4 v_eye_position;

void main (void) {
    v_color = color;
    v_eye_position = <transform.trackball_view> *
                     <transform.trackball_model> *
                     vec4(position,1.0);
    gl_Position = <transform(position)>;
    vec4 p = <transform.trackball_projection> *
             vec4(1., 1., v_eye_position.z, v_eye_position.w);
    v_size = 512.0 * p.x / p.w;
    gl_PointSize = v_size;
}
"""

line_fragment = """
varying vec3 v_color;

void main()
{
    gl_FragColor = vec4(v_color, 1);
}
"""


window = app.Window(width=1600, height=1600, color=(1,1,1,1))

protein = gloo.Program(vertex, fragment)
protein['light_position'] = 0., 0., 2.
protein["transform"] = Trackball(Position())

pdb_file_name = sys.argv[1]
parser = PDBParser(QUIET=True, PERMISSIVE=True)
structure = parser.get_structure('structure', pdb_file_name)
if only_CA:
    atoms = [atom for atom in structure.get_atoms() if atom.id == 'CA']
else:
    atoms = [atom for atom in structure.get_atoms()]
coordinates = np.array([atom.coord for atom in atoms])
center = centroid(coordinates)
coordinates -= center
coordinates *= 0.05

a_data = []
colors = [ atom_colors[atom.element] if atom_colors.has_key(atom.element) else white for atom in atoms ]
radius = [ vdw_radius[atom.element] if vdw_radius.has_key(atom.element) else 1.0 for atom in atoms]

a_data = np.zeros(len(colors), [("position", np.float32, 3),
                                ("color",    np.float32, 4),
                                ("radius",    np.float32, 1)])

for i in range(len(colors)):
    a_data[i]['position'] = [coordinates[i][0], coordinates[i][1], coordinates[i][2]]
    a_data[i]['color'] = colors[i]
    a_data[i]['radius'] = radius[i] / radii_ratio

#a_data = np.load('./protein.npy')

protein.bind(a_data.view(gloo.VertexBuffer))
protein['color'] *= .25
protein['color'] += .75


bonds = gloo.Program(line_shader, line_fragment)
bonds["transform"] = protein["transform"]
bonds.bind(a_data.view(gloo.VertexBuffer))
bonds['color'] = protein['color']


@window.event
def on_draw(dt):
    window.clear()
    protein.draw(gl.GL_POINTS)
    if draw_bonds:
        bonds.draw(gl.GL_LINE_STRIP)

@window.event
def on_init():
    gl.glEnable(gl.GL_DEPTH_TEST)
    gl.glEnable(gl.GL_LINE_SMOOTH)
    gl.glLineWidth(50.0)


window.attach(protein["transform"])
app.run()
