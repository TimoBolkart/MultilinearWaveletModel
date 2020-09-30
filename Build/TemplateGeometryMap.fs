// TemplateDisparityMap.fs
//
// Compute template disparity map from 3D mesh
//
// Alan Burnton, NRC-IIT VIT July 2009


uniform float baseline;
uniform float focal;
uniform float pixelSize;

varying vec3 viewCoords;


void main(void)
{
	gl_FragData[0] = vec4(viewCoords, 1.0);
//	gl_FragData[0] = vec4(1.0);
}
