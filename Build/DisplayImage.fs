// DisplayImage.fs
//
// display an image

uniform sampler2D texImage;

varying vec2 imagePos;

void main(void)
{
	gl_FragColor = abs(texture2D(texImage, gl_TexCoord[0].xy));
//	gl_FragColor = texture2D(texImage, gl_TexCoord[0].xy) * 0.5 + 0.5;
}
