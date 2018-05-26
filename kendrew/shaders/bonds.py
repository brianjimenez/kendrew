
vertex = """
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

fragment = """
varying vec3 v_color;

void main()
{
    gl_FragColor = vec4(v_color, 1);
}
"""
