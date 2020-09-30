// diffuse.vs
//
// set up interpolants for diffuse lighting

uniform vec3 lightPos[2];
varying vec3 N, L0, L1;

void main(void)
{
    // vertex MVP transform
    gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;

    // eye-space normal
    N = gl_NormalMatrix * gl_Normal;

    // eye-space light vector
    vec4 V = gl_ModelViewMatrix * gl_Vertex;    
    L0 = lightPos[0] - V.xyz;
    L1 = lightPos[1] - V.xyz;

    // Copy the primary color
    gl_FrontColor = gl_Color;
}
