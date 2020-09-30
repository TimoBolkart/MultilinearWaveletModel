// TemplateDisparityProject.vs
//
// Projects mesh into template disparity (left) image frame
//
// Alan Brunton, NRC-IIT VIT July 2009



uniform float baseline;
uniform float focal;
uniform float pixelSize;

varying vec3 viewCoords;


void main(void)
{
    // vertex MVP transform
    gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;

	vec4 vecView;
	vecView = gl_ModelViewMatrix * gl_Vertex;
	
	viewCoords = vecView.xyz;
	
    // Copy the primary color
    gl_FrontColor = gl_Color;
}
