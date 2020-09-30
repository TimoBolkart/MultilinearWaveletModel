// DisplayImage.fs
//
// display an image

#version 120

uniform sampler2D texDisparity;

uniform float scale, bias;

uniform float baseline;
uniform float focal;
uniform float pixelSize;

uniform mat4x3 matProjLeft;
uniform mat4x3 matProjRight;

uniform vec4 gridTemplate;

void main(void)
{
#if 0
	const float eps_disp = 1e-8;
    const float rhoD = 0.7;
	
	const vec4 light_pos = vec4(0.0, 1.0, -1.0, 1.0);
	
	float rcpFocal = 1.0 / focal;
	float rcpPixelSize = 1.0 / pixelSize;
	float dispCenter, dispRight, dispUp, dispLeft, dispDown;
	float intensity;
	
	vec2 gradDisp, gradZ;
	vec2 center = gridTemplate.zw * 0.5;
	
	vec3 N, L;
	vec4 v = vec4(0.0, 0.0, 0.0, 1.0);
	vec4 surfaceColor = vec4(0.8, 0.8, 0.8, 1.0);
	
	dispCenter = texture2D(texDisparity, gl_TexCoord[0].xy).x;
	dispRight = texture2D(texDisparity, gl_TexCoord[0].xy + vec2(gridTemplate.x, 0.0)).x;
	dispUp = texture2D(texDisparity, gl_TexCoord[0].xy + vec2(0.0, -gridTemplate.y)).x;
	dispLeft = texture2D(texDisparity, gl_TexCoord[0].xy + vec2(-gridTemplate.x, 0.0)).x;
	dispDown = texture2D(texDisparity, gl_TexCoord[0].xy + vec2(0.0, gridTemplate.y)).x;
	
	gradDisp.x = 0.5 * (dispRight - dispLeft);
	gradDisp.y = 0.5 * (dispUp - dispDown);
	
	gradZ = -focal * baseline / (dispCenter*dispCenter) * rcpPixelSize * gradDisp;
	
	v.xy = gl_FragCoord.xy;
	v.y = gridTemplate.w - v.y;
	v.z = focal * baseline / max(dispCenter, eps_disp);
	v.w = 1.0;
	
	v.xy -= center;
	v.xy *= rcpFocal;
	v.xy *= v.z;
	
	N = normalize(vec3(-gradZ.x, -gradZ.y, -1.0));
	L = normalize(light_pos.xyz - v.xyz);
	
    intensity = max(0.0, dot(N, L));
    
	surfaceColor = vec4(dispCenter) * scale + bias;
	gl_FragColor = surfaceColor * rhoD * intensity;
#else
	gl_FragColor = texture2D(texDisparity, gl_TexCoord[0].xy) * scale + bias;
#endif
}
