import numpy as np
import sys
from vispy import gloo
from vispy import app
from vispy.util.transforms import perspective, translate, rotate
from Bio.PDB.PDBParser import PDBParser


W, H = 1200, 800


VERT_SHADER = """
uniform mat4 u_model;
uniform mat4 u_view;
uniform mat4 u_projection;
attribute vec3 a_position;

void main (void) {
    gl_Position = u_projection * u_view * u_model * vec4(a_position,1.0);
}
"""

FRAG_SHADER = """

void main()
{
    gl_FragColor = vec4(0,0,0,1);
}
"""


def centroid(arr):
    length = arr.shape[0]
    sum_x = np.sum(arr[:, 0])
    sum_y = np.sum(arr[:, 1])
    sum_z = np.sum(arr[:, 2])
    return sum_x/length, sum_y/length, sum_z/length


class Canvas(app.Canvas):

    def __init__(self, pdb_file_name):
        app.Canvas.__init__(self, keys='interactive', size=(W, H))

        self.program = gloo.Program(VERT_SHADER, FRAG_SHADER)

        self.parser = PDBParser(QUIET=True, PERMISSIVE=True)
        self.structure = self.parser.get_structure('structure', pdb_file_name)
        self.atoms = [atom for atom in self.structure.get_atoms() if atom.id == 'CA']
        self.coordinates = np.array([atom.coord for atom in self.atoms])
        self.center = centroid(self.coordinates)
        self.coordinates -= self.center

        self.program['a_position'] = gloo.VertexBuffer(self.coordinates)

        self.translate = 50
        self.view = translate((0, 0, -self.translate), dtype=np.float32)
        self.model = np.eye(4, dtype=np.float32)

        gloo.set_viewport(0, 0, self.physical_size[0], self.physical_size[1])
        self.projection = perspective(45.0, self.size[0] /
                                      float(self.size[1]), 1.0, 1000.0)
        self.program['u_projection'] = self.projection

        self.program['u_model'] = self.model
        self.program['u_view'] = self.view

        self.theta = 0
        self.phi = 0

        self.context.set_clear_color('white')
        self.context.set_state('translucent')

        self.timer = app.Timer('auto', connect=self.on_timer)

        self.show()

    def on_key_press(self, event):
        if event.text == ' ':
            if self.timer.running:
                self.timer.stop()
            else:
                self.timer.start()
        if event.text == 'a':
            self.theta += 5.0
            self.rotate_molecule()
        if event.text == 's':
            self.theta -= 5.0
            self.rotate_molecule()
        if event.text == 'q':
            self.phi += 5.0
            self.rotate_molecule()
        if event.text == 'w':
            self.phi -= 5.0
            self.rotate_molecule()

    def rotate_molecule(self):
        self.model = np.dot(rotate(self.theta, (0, 0, 1)),
                            rotate(self.phi, (0, 1, 0)))
        self.program['u_model'] = self.model
        self.update()

    def on_timer(self, event):
        self.phi += .25
        self.rotate_molecule()

    def on_resize(self, event):
        gloo.set_viewport(0, 0, event.physical_size[0], event.physical_size[1])
        self.projection = perspective(45.0, event.size[0] /
                                      float(event.size[1]), 1.0, 1000.0)
        self.program['u_projection'] = self.projection

    def on_mouse_wheel(self, event):
        self.translate += event.delta[1]
        self.translate = max(2, self.translate)
        self.view = translate((0, 0, -self.translate))
        self.program['u_view'] = self.view
        self.update()

    def on_draw(self, event):
        self.context.clear()
        self.program.draw('line_strip')


if __name__ == '__main__':
    pdb_file_name = sys.argv[1]
    c = Canvas(pdb_file_name)
    app.run()
