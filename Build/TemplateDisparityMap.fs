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
	float disparity = focal * baseline / abs(viewCoords.z);
//	disparity /= pixelSize;
	
	gl_FragData[0] = disparity;
}
