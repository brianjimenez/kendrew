# -*- coding: utf-8 -*-
# vispy: gallery 30
# -----------------------------------------------------------------------------
# 2014, Aurore Deschildre, Gael Goret, Cyrille Rossant, Nicolas P. Rougier.
# Distributed under the terms of the new BSD License.
# -----------------------------------------------------------------------------
import numpy as np
import sys
from vispy import gloo
from vispy import app
from vispy.util.transforms import perspective, translate, rotate
from vispy.io import load_data_file
from Bio.PDB.PDBParser import PDBParser

"""
White = Color(1.0, 1.0, 1.0)
Black = Color(0.0, 0.0, 0.0)
Carbon = Color(0.17, 0.17, 0.18)
Red = Color(0.95, 0.03, 0.01)
Blue = Color(0.01, 0.03, 0.95)
Sky = Color(0.233, 0.686, 1.0)
Yellow = Color(1.0, 1.0, 0.0)
Green = Color(0.0, 0.53, 0.0)
Pink = Color(0.53, 0.12, 0.36)
DarkRed = Color(0.59, 0.13, 0.0)
Violet = Color(0.46, 0.0, 1.0)
DarkViolet = Color(0.39, 0.0, 0.73)
Cyan = Color(0.0, 1.0, 1.0)
Orange = Color(1.0, 0.59, 0.0)
Peach = Color(1.0, 0.66, 0.46)
DarkGreen = Color(0.0, 0.46, 0.0)
Gray = Color(0.59, 0.59, 0.59)
DarkOrange = Color(0.86, 0.46, 0.0)
"""

atom_colors = {'H':np.array([1.0, 1.0, 1.0]), 'C':np.array([0.17, 0.17, 0.18]), 'N':np.array([0.233, 0.686, 1.0]),
               'O':np.array([0.95, 0.03, 0.01]), 'S':np.array([0.95, 0.95, 0.00]), 'AU':np.array([0.75, 0.625, 0.38])}

vdw_radius = {'C':1.9080, 'N':1.8240, 'O':1.6612, 'S':2.0000, 'AU':1.37, 'H':1.0}

vertex = """
#version 120

uniform mat4 u_model;
uniform mat4 u_view;
uniform mat4 u_projection;
uniform vec3 u_light_position;
uniform vec3 u_light_spec_position;

attribute vec3  a_position;
attribute vec3  a_color;
attribute float a_radius;

varying vec3  v_color;
varying vec4  v_eye_position;
varying float v_radius;
varying vec3  v_light_direction;

void main (void) {
    v_radius = a_radius;
    v_color = a_color;

    v_eye_position = u_view * u_model * vec4(a_position,1.0);
    v_light_direction = normalize(u_light_position);
    float dist = length(v_eye_position.xyz);

    gl_Position = u_projection * v_eye_position;

    // stackoverflow.com/questions/8608844/...
    //  ... resizing-point-sprites-based-on-distance-from-the-camera
    vec4  proj_corner = u_projection * vec4(a_radius, a_radius, v_eye_position.z, v_eye_position.w);  // # noqa
    gl_PointSize = 512.0 * proj_corner.x / proj_corner.w;
}
"""

fragment = """
#version 120

uniform mat4 u_model;
uniform mat4 u_view;
uniform mat4 u_projection;
uniform vec3 u_light_position;
uniform vec3 u_light_spec_position;

varying vec3  v_color;
varying vec4  v_eye_position;
varying float v_radius;
varying vec3  v_light_direction;
void main()
{
    // r^2 = (x - x0)^2 + (y - y0)^2 + (z - z0)^2
    vec2 texcoord = gl_PointCoord* 2.0 - vec2(1.0);
    float x = texcoord.x;
    float y = texcoord.y;
    float d = 1.0 - x*x - y*y;
    if (d <= 0.0)
        discard;

    float z = sqrt(d);
    vec4 pos = v_eye_position;
    pos.z += v_radius*z;
    vec3 pos2 = pos.xyz;
    pos = u_projection * pos;
    // gl_FragDepth = 0.5*(pos.z / pos.w)+0.5;
    vec3 normal = vec3(x,y,z);
    float diffuse = clamp(dot(normal, v_light_direction), 0.0, 1.0);

    // Specular lighting.
    vec3 M = pos2.xyz;
    vec3 O = v_eye_position.xyz;
    vec3 L = u_light_spec_position;
    vec3 K = normalize(normalize(L - M) + normalize(O - M));
    // WARNING: abs() is necessary, otherwise weird bugs may appear with some
    // GPU drivers...
    float specular = clamp(pow(abs(dot(normal, K)), 40.), 0.0, 1.0);
    vec3 v_light = vec3(1., 1., 1.);
    gl_FragColor.rgb = (.15*v_color + .55*diffuse * v_color
                        + .35*specular * v_light);
}
"""

def centroid(arr):
    length = arr.shape[0]
    sum_x = np.sum(arr[:, 0])
    sum_y = np.sum(arr[:, 1])
    sum_z = np.sum(arr[:, 2])
    return sum_x/length, sum_y/length, sum_z/length


def generate_colors_for_models(num_models):
    color_models = [np.random.rand(3) for i in range(num_models)]
    return color_models


class Canvas(app.Canvas):

    def __init__(self, pdb_file_name, mode='cpk'):
        app.Canvas.__init__(self, title='Molecular viewer',
                            keys='interactive', size=(1200, 800))
        self.ps = self.pixel_scale

        self.parser = PDBParser(QUIET=True, PERMISSIVE=True)
        self.structure = self.parser.get_structure('structure', pdb_file_name)
        self.num_models = len(self.structure)
        self.model_colors = generate_colors_for_models(self.num_models)

        self.atoms = [atom for atom in self.structure.get_atoms()] #if atom.get_name() == 'CA']
        self.coordinates = np.array([atom.coord for atom in self.atoms])
        self.center = centroid(self.coordinates)
        self.coordinates -= self.center
        self.translate = 120
        self.program = gloo.Program(vertex, fragment)
        self.view = translate((0, 0, -self.translate))
        self.model = np.eye(4, dtype=np.float32)
        self.projection = np.eye(4, dtype=np.float32)
        self.mode = mode


        self.apply_zoom()

        self.load_molecule()
        self.load_data()

        self.theta = 0
        self.phi = 0

        gloo.set_state(depth_test=True, clear_color='black')
        self.timer = app.Timer('auto', connect=self.on_timer, start=True)

        self.show()

    def load_molecule(self):
        self._nAtoms = len(self.atoms)

        # The array that will store the color and alpha scale for all the atoms
        #self.atomsColours = np.array([1.0,2.0,3.0]*self._nAtoms).reshape(self._nAtoms,3)
        self.atomsColours = []
        if self.mode == 'model':
            for atom in self.atoms:
                residue = atom.get_parent()
                chain = residue.get_parent()
                model = chain.get_parent()
                self.atomsColours.append(self.model_colors[model.id])
        else:
            for atom in self.atoms:
                self.atomsColours.append(atom_colors[atom.element])
        self.atomsColours = np.array(self.atomsColours)

        # The array that will store the scale for all the atoms.
        self.atomsScales = []
        for atom in self.atoms:
            self.atomsScales.append(vdw_radius[atom.element])
        self.atomsScales = np.array(self.atomsScales)  #np.random.uniform(low=12.0, high=12.0, size=(self._nAtoms,))

    def load_data(self):
        n = self._nAtoms

        data = np.zeros(n, [('a_position', np.float32, 3),
                            ('a_color', np.float32, 3),
                            ('a_radius', np.float32, 1)])

        data['a_position'] = self.coordinates
        data['a_color'] = self.atomsColours
        data['a_radius'] = self.atomsScales*self.ps

        self.program.bind(gloo.VertexBuffer(data))

        self.program['u_model'] = self.model
        self.program['u_view'] = self.view
        self.program['u_light_position'] = 0., 0., 2.
        self.program['u_light_spec_position'] = -5., 5., -5.

    def on_key_press(self, event):
        if event.text == ' ':
            if self.timer.running:
                self.timer.stop()
            else:
                self.timer.start()
        if event.text == 'a':
            self.theta += 5.0
            self.model = np.dot(rotate(self.theta, (0, 0, 1)),
                            rotate(self.phi, (0, 1, 0)))
            self.program['u_model'] = self.model
            self.update()
        if event.text == 's':
            self.theta -= 5.0
            self.model = np.dot(rotate(self.theta, (0, 0, 1)),
                            rotate(self.phi, (0, 1, 0)))
            self.program['u_model'] = self.model
            self.update()

    def on_timer(self, event):
        self.theta += .25
        self.phi += .25
        self.model = np.dot(rotate(self.theta, (0, 0, 1)),
                            rotate(self.phi, (0, 1, 0)))
        self.program['u_model'] = self.model
        self.update()

    def on_resize(self, event):
        width, height = event.size

    def apply_zoom(self):
        width, height = self.physical_size
        gloo.set_viewport(0, 0, width, height)
        self.projection = perspective(95.0, width / float(height), 1.0, 400.0)
        self.program['u_projection'] = self.projection

    def on_mouse_wheel(self, event):
        self.translate -= event.delta[1]
        self.translate = max(-1, self.translate)
        self.view = translate((0, 0, -self.translate))

        self.program['u_view'] = self.view
        self.update()

    def on_draw(self, event):
        gloo.clear()
        self.program.draw('points')


if __name__ == '__main__':
    pdb_file_name = sys.argv[1]
    try:
        mode = sys.argv[2]
    except:
        mode = 'cpk'
    mvc = Canvas(pdb_file_name, mode)
    app.run()
