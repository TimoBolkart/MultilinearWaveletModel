// TemplateDisparityProject.vs
//
// Projects mesh into template disparity (left) image frame
//
// Alan Brunton, NRC-IIT VIT July 2009



uniform mat4	matRigidAlign;

varying vec3	viewCoords;


void main(void)
{
	vec4 vecTemplatePos = vec4(0.0, 0.0, 0.0, 1.0);
	vec4 vecView, vecAligned;

	//apply ModelView transform to vertex coordinates, pass to fragment shader	
	vecView = gl_ModelViewMatrix * gl_Vertex;
//	viewCoords = vecView.xyz;
//	viewCoords = (matRigidAlign * gl_Vertex).xyz;
//	viewCoords = gl_Vertex.xyz;
	vecAligned = matRigidAlign * gl_Vertex;
	viewCoords = vecAligned.xyz;
	
    //2D orthographic projection
	vecTemplatePos.xy = gl_MultiTexCoord0.xy;
//	vecTemplatePos.z = abs(vecView.z);
//	vecTemplatePos.z = 0.0;
//	vecTemplatePos.z = vecView.z + 1.0;
//	vecTemplatePos.z = vecView.z - 1.0;
	gl_Position = gl_ProjectionMatrix * vecTemplatePos;

    // Copy the primary color
    gl_FrontColor = gl_Color;
}
