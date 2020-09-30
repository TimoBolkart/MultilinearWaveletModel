///////////////////////////////////////////////////////////////////////////////
//
//	StereoFace.cpp
//
//	Main source module for the StereoFace program
//
//	Alan Brunton, October 2008
//
///////////////////////////////////////////////////////////////////////////////


#include "StereoFace.h"

// Load shader from disk into a null-terminated string
GLchar *LoadShaderText(const char *fileName)
{
  GLchar *shaderText = NULL;
  GLint shaderLength = 0;
  FILE *fp;

  fp = fopen(fileName, "r");
  if (fp != NULL) {
    while (fgetc(fp) != EOF) {
      shaderLength++;
    }
    rewind(fp);
    shaderText = (GLchar *)malloc(shaderLength+1);
    if (shaderText != NULL) {
      fread(shaderText, 1, shaderLength, fp);
    }
    shaderText[shaderLength] = '\0';
    fclose(fp);
  }

  return shaderText;
}

///////////////////////////////////////////////////////////////////////////////
//OpenGL shader and framebuffer functions
///////////////////////////////////////////////////////////////////////////////

void fbInit(GLuint & nfb)
{
	if (nfb == 0)
	{
		glGenFramebuffersEXT(1, &nfb);
		_ASSERT(nfb != 0);
	}
}

void fbTerm(GLuint & nfb)
{
	if (nfb != 0)
	{
		fbUnbind();
		glDeleteFramebuffersEXT(1, &nfb);
		nfb = 0;
	}
}

void fbBind(GLuint & nfb)
{
	_ASSERT(nfb != 0);
	GLint fb;

	glGetIntegerv(GL_FRAMEBUFFER_BINDING_EXT, &fb);
//	if (fb != nfb)
	if (fb == 0)
		glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, nfb);
}

void fbUnbind()
{
	GLint fb;

	glGetIntegerv( GL_FRAMEBUFFER_BINDING_EXT, &fb);
	if (fb != 0)
		glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
}

int fbMaxAttachments()
{
	int maxAttach = 0;

	glGetIntegerv( GL_MAX_COLOR_ATTACHMENTS_EXT, &maxAttach );

	return maxAttach;
}

void fbAttachTexture(GLuint & nfb, GLenum attachment, GLenum texTarget, GLuint texId, int mipLevel, int zSlice)
{
	_ASSERT(nfb != 0);

	switch (texTarget)
	{
	case GL_TEXTURE_1D:
		glFramebufferTexture1DEXT( GL_FRAMEBUFFER_EXT, attachment, GL_TEXTURE_1D, texId, mipLevel );
		break;
	case GL_TEXTURE_3D:
		glFramebufferTexture3DEXT( GL_FRAMEBUFFER_EXT, attachment, GL_TEXTURE_3D, texId, mipLevel, zSlice );
		break;
	default:
		glFramebufferTexture2DEXT( GL_FRAMEBUFFER_EXT, attachment, texTarget, texId, mipLevel );
		break;
	}
}

void fbDetach(GLenum attachment, GLuint & nfb)
{
	_ASSERT(nfb != 0);

	fbAttachTexture(nfb, attachment, GL_TEXTURE_2D, 0);
}

void fbDetachAll(GLuint & nfb)
{
	_ASSERT(nfb != 0);

	int i, nAttach = fbMaxAttachments();
	for (i = 0; i < nAttach; i++)
		fbDetach(GL_COLOR_ATTACHMENT0_EXT + i, nfb);

	fbDetach(GL_DEPTH_ATTACHMENT_EXT, nfb);
}

GLuint loadVertexShader(const char* szFilename)
{
	GLuint shader = 0;
	const GLchar* aszShaderPtr[1];

	GLchar* szShaderString = LoadShaderText(szFilename);
	if (szShaderString == NULL)
	{
		fprintf(stderr, "Cannot load shader %s\n", szFilename);
		return shader;
	}

	shader = glCreateShader(GL_VERTEX_SHADER);
	aszShaderPtr[0] = szShaderString;
	glShaderSource(shader, 1, aszShaderPtr, NULL);

	free(szShaderString);
	return shader;
}

GLuint loadFragmentShader(const char* szFilename)
{
	GLuint shader = 0;
	const GLchar* aszShaderPtr[1];

	GLchar* szShaderString = LoadShaderText(szFilename);
	if (szShaderString == NULL)
	{
		fprintf(stderr, "Cannot load shader %s\n", szFilename);
		return shader;
	}

	shader = glCreateShader(GL_FRAGMENT_SHADER);
	aszShaderPtr[0] = szShaderString;
	glShaderSource(shader, 1, aszShaderPtr, NULL);

	free(szShaderString);
	return shader;
}

void programError(const char* szMessage, GLuint program)
{
	GLchar infoLog[MAX_INFO_LOG_SIZE];
	glGetProgramInfoLog(program, MAX_INFO_LOG_SIZE, NULL, infoLog);
	fprintf(stderr, "%s:\n", szMessage);
	fprintf(stderr, "%s\n", infoLog);
}

void shaderError(const char* szMessage, GLuint shader)
{
	GLchar infoLog[MAX_INFO_LOG_SIZE];
	glGetShaderInfoLog(shader, MAX_INFO_LOG_SIZE, NULL, infoLog);
	fprintf(stderr, "%s:\n", szMessage);
	fprintf(stderr, "%s\n", infoLog);
}

GLenum checkGLError(const char* szLastCommand)
{
	GLenum err = glGetError();

	switch (err)
	{
	case GL_INVALID_ENUM:
		fprintf(stderr, "GL ERROR: ");
		if (szLastCommand != NULL)
			fprintf(stderr, "%s; ", szLastCommand);
		fprintf(stderr, "Invalid enum\n");
		break;
	case GL_INVALID_VALUE:
		fprintf(stderr, "GL ERROR: ");
		if (szLastCommand != NULL)
			fprintf(stderr, "%s; ", szLastCommand);
		fprintf(stderr, "Invalid value\n");
		break;
	case GL_INVALID_OPERATION:
		fprintf(stderr, "GL ERROR: ");
		if (szLastCommand != NULL)
			fprintf(stderr, "%s; ", szLastCommand);
		fprintf(stderr, "Invalid operation\n");
		break;
	case GL_STACK_OVERFLOW:
		fprintf(stderr, "GL ERROR: ");
		if (szLastCommand != NULL)
			fprintf(stderr, "%s; ", szLastCommand);
		fprintf(stderr, "Stack overflow\n");
		break;
	case GL_STACK_UNDERFLOW:
		fprintf(stderr, "GL ERROR: ");
		if (szLastCommand != NULL)
			fprintf(stderr, "%s; ", szLastCommand);
		fprintf(stderr, "Stack underflow\n");
		break;
	case GL_OUT_OF_MEMORY:
		fprintf(stderr, "GL ERROR: ");
		if (szLastCommand != NULL)
			fprintf(stderr, "%s; ", szLastCommand);
		fprintf(stderr, "Out of memory\n");
		break;
	case GL_TABLE_TOO_LARGE:
		fprintf(stderr, "GL ERROR: ");
		if (szLastCommand != NULL)
			fprintf(stderr, "%s; ", szLastCommand);
		fprintf(stderr, "Table too large\n");
		break;
	case GL_NO_ERROR:
		break;
	}

	fflush(stderr);

	return err;
}


///////////////////////////////////////////////////////////////////////////////
//miscelaneous math functions
///////////////////////////////////////////////////////////////////////////////

//use SVD to make a rotation matrix orthogonal
void svdOrthogonalize(double* pmatR)
{
	char jobu = 'A';
	char jobvt = 'A';
	clapack::integer m = 3, n = 3;
	double* pmatA = new double[m * n];
	double* pS = new double[min(m, n)];
	double* pmatU = new double[m * m];
	double* pmatVtrans = new double[n * n];
	clapack::integer nWork = 10 * m * n;
	double* pWork = new double[nWork];
	clapack::integer nInfo;

	//column-major for clapack
	pmatA[0] = pmatR[0];
	pmatA[1] = pmatR[3];
	pmatA[2] = pmatR[6];
	pmatA[3] = pmatR[1];
	pmatA[4] = pmatR[4];
	pmatA[5] = pmatR[7];
	pmatA[6] = pmatR[2];
	pmatA[7] = pmatR[5];
	pmatA[8] = pmatR[8];

	clapack::dgesvd_(&jobu, &jobvt, &m, &n, pmatA, &m, pS, pmatU, &m, pmatVtrans, &n, pWork, &nWork, &nInfo);

	//row major for my code
	//R=UIV^t=UV^t
	mul_3x3_3x3(pmatA, pmatU, true, pmatVtrans, true);
	pmatR[0] = pmatA[0];
	pmatR[1] = pmatA[1];
	pmatR[2] = pmatA[2];
	pmatR[3] = pmatA[3];
	pmatR[4] = pmatA[4];
	pmatR[5] = pmatA[5];
	pmatR[6] = pmatA[6];
	pmatR[7] = pmatA[7];
	pmatR[8] = pmatA[8];
//	::memcpy(pmatR, pmatA, 9 * sizeof(double));

	delete [] pmatA;
	delete [] pS;
	delete [] pmatU;
	delete [] pmatVtrans;
	delete [] pWork;
}

///////////////////////////////////////////////////////////////////////////////
//helper functions
///////////////////////////////////////////////////////////////////////////////

//c = a x b
void cross_prod_3(double a[3], double b[3], double c[3])
{
	c[0] = a[1] * b[2] - a[2] * b[1];
	c[1] = a[2] * b[0] - a[0] * b[2];
	c[2] = a[0] * b[2] - a[1] * b[0];
}

//in-place normalization
void normalize(double a[3])
{
	double x2 = a[0]*a[0], y2 = a[1]*a[1], z2 = a[2]*a[2];
	double rcpmag = 1.0 / sqrt(x2 + y2 + z2);
	a[0] *= rcpmag;
	a[1] *= rcpmag;
	a[2] *= rcpmag;
}

void mul_3x3_3x3(double* pOut, double* pA, bool bTransA, double* pB, bool bTransB)
{
	if (bTransA && bTransB)
	{
		pOut[0] = pA[0] * pB[0] + pA[3] * pB[1] + pA[6] * pB[2];
		pOut[1] = pA[0] * pB[3] + pA[3] * pB[4] + pA[6] * pB[5];
		pOut[2] = pA[0] * pB[6] + pA[3] * pB[7] + pA[6] * pB[8];
		pOut[3] = pA[1] * pB[0] + pA[4] * pB[1] + pA[7] * pB[2];
		pOut[4] = pA[1] * pB[3] + pA[4] * pB[4] + pA[7] * pB[5];
		pOut[5] = pA[1] * pB[6] + pA[4] * pB[7] + pA[7] * pB[8];
		pOut[6] = pA[2] * pB[0] + pA[5] * pB[1] + pA[8] * pB[2];
		pOut[7] = pA[2] * pB[3] + pA[5] * pB[4] + pA[8] * pB[5];
		pOut[8] = pA[2] * pB[6] + pA[5] * pB[7] + pA[8] * pB[8];
	}
	else if (bTransA)
	{
		pOut[0] = pA[0] * pB[0] + pA[3] * pB[3] + pA[6] * pB[6];
		pOut[1] = pA[0] * pB[1] + pA[3] * pB[4] + pA[6] * pB[7];
		pOut[2] = pA[0] * pB[2] + pA[3] * pB[5] + pA[6] * pB[8];
		pOut[3] = pA[1] * pB[0] + pA[4] * pB[3] + pA[7] * pB[6];
		pOut[4] = pA[1] * pB[1] + pA[4] * pB[4] + pA[7] * pB[7];
		pOut[5] = pA[1] * pB[2] + pA[4] * pB[5] + pA[7] * pB[8];
		pOut[6] = pA[2] * pB[0] + pA[5] * pB[3] + pA[8] * pB[6];
		pOut[7] = pA[2] * pB[1] + pA[5] * pB[4] + pA[8] * pB[7];
		pOut[8] = pA[2] * pB[2] + pA[5] * pB[5] + pA[8] * pB[8];
	}
	else if (bTransB)
	{
		pOut[0] = pA[0] * pB[0] + pA[1] * pB[1] + pA[2] * pB[2];
		pOut[1] = pA[0] * pB[3] + pA[1] * pB[4] + pA[2] * pB[5];
		pOut[2] = pA[0] * pB[6] + pA[1] * pB[7] + pA[2] * pB[8];
		pOut[3] = pA[3] * pB[0] + pA[4] * pB[1] + pA[5] * pB[2];
		pOut[4] = pA[3] * pB[3] + pA[4] * pB[4] + pA[5] * pB[5];
		pOut[5] = pA[3] * pB[6] + pA[4] * pB[7] + pA[5] * pB[8];
		pOut[6] = pA[6] * pB[0] + pA[7] * pB[1] + pA[8] * pB[2];
		pOut[7] = pA[6] * pB[3] + pA[7] * pB[4] + pA[8] * pB[5];
		pOut[8] = pA[6] * pB[6] + pA[7] * pB[7] + pA[8] * pB[8];
	}
	else
	{
		pOut[0] = pA[0] * pB[0] + pA[1] * pB[3] + pA[2] * pB[6];
		pOut[1] = pA[0] * pB[1] + pA[1] * pB[4] + pA[2] * pB[7];
		pOut[2] = pA[0] * pB[2] + pA[1] * pB[5] + pA[2] * pB[8];
		pOut[3] = pA[3] * pB[0] + pA[4] * pB[3] + pA[5] * pB[6];
		pOut[4] = pA[3] * pB[1] + pA[4] * pB[4] + pA[5] * pB[7];
		pOut[5] = pA[3] * pB[2] + pA[4] * pB[5] + pA[5] * pB[8];
		pOut[6] = pA[6] * pB[0] + pA[7] * pB[3] + pA[8] * pB[6];
		pOut[7] = pA[6] * pB[1] + pA[7] * pB[4] + pA[8] * pB[7];
		pOut[8] = pA[6] * pB[2] + pA[7] * pB[5] + pA[8] * pB[8];
	}
}

void mul_4x4_4x4(double* pOut, double* pA, bool bTransA, double* pB, bool bTransB)
{
	if (bTransA && bTransB)
	{
		pOut[0]		= pA[0] * pB[0] + pA[4] * pB[1] + pA[8] * pB[2] + pA[12] * pB[3];
		pOut[1]		= pA[0] * pB[4] + pA[4] * pB[5] + pA[8] * pB[6] + pA[12] * pB[7];
		pOut[2]		= pA[0] * pB[8] + pA[4] * pB[9] + pA[8] * pB[10] + pA[12] * pB[11];
		pOut[3]		= pA[0] * pB[12] + pA[4] * pB[13] + pA[8] * pB[14] + pA[12] * pB[15];
		pOut[4]		= pA[1] * pB[0] + pA[5] * pB[1] + pA[9] * pB[2] + pA[13] * pB[3];
		pOut[5]		= pA[1] * pB[4] + pA[5] * pB[5] + pA[9] * pB[6] + pA[13] * pB[7];
		pOut[6]		= pA[1] * pB[8] + pA[5] * pB[9] + pA[9] * pB[10] + pA[13] * pB[11];
		pOut[7]		= pA[1] * pB[12] + pA[5] * pB[13] + pA[9] * pB[14] + pA[13] * pB[15];
		pOut[8]		= pA[2] * pB[0] + pA[6] * pB[1] + pA[10] * pB[2] + pA[14] * pB[3];
		pOut[9]		= pA[2] * pB[4] + pA[6] * pB[5] + pA[10] * pB[6] + pA[14] * pB[7];
		pOut[10]	= pA[2] * pB[8] + pA[6] * pB[9] + pA[10] * pB[10] + pA[14] * pB[11];
		pOut[11]	= pA[2] * pB[12] + pA[6] * pB[13] + pA[10] * pB[14] + pA[14] * pB[15];
		pOut[12]	= pA[3] * pB[0] + pA[7] * pB[1] + pA[11] * pB[2] + pA[15] * pB[3];
		pOut[13]	= pA[3] * pB[4] + pA[7] * pB[5] + pA[11] * pB[6] + pA[15] * pB[7];
		pOut[14]	= pA[3] * pB[8] + pA[7] * pB[9] + pA[11] * pB[10] + pA[15] * pB[11];
		pOut[15]	= pA[3] * pB[12] + pA[7] * pB[13] + pA[11] * pB[14] + pA[15] * pB[15];
	}
	else if (bTransA)
	{
		pOut[0]		= pA[0] * pB[0] + pA[4] * pB[4] + pA[8] * pB[8] + pA[12] * pB[12];
		pOut[1]		= pA[0] * pB[1] + pA[4] * pB[5] + pA[8] * pB[9] + pA[12] * pB[13];
		pOut[2]		= pA[0] * pB[2] + pA[4] * pB[6] + pA[8] * pB[10] + pA[12] * pB[14];
		pOut[3]		= pA[0] * pB[3] + pA[4] * pB[7] + pA[8] * pB[11] + pA[12] * pB[15];
		pOut[4]		= pA[1] * pB[0] + pA[5] * pB[4] + pA[9] * pB[8] + pA[13] * pB[12];
		pOut[5]		= pA[1] * pB[1] + pA[5] * pB[5] + pA[9] * pB[9] + pA[13] * pB[13];
		pOut[6]		= pA[1] * pB[2] + pA[5] * pB[6] + pA[9] * pB[10] + pA[13] * pB[14];
		pOut[7]		= pA[1] * pB[3] + pA[5] * pB[7] + pA[9] * pB[11] + pA[13] * pB[15];
		pOut[8]		= pA[2] * pB[0] + pA[6] * pB[4] + pA[10] * pB[8] + pA[14] * pB[12];
		pOut[9]		= pA[2] * pB[1] + pA[6] * pB[5] + pA[10] * pB[9] + pA[14] * pB[13];
		pOut[10]	= pA[2] * pB[2] + pA[6] * pB[6] + pA[10] * pB[10] + pA[14] * pB[14];
		pOut[11]	= pA[2] * pB[3] + pA[6] * pB[7] + pA[10] * pB[11] + pA[14] * pB[15];
		pOut[12]	= pA[3] * pB[0] + pA[7] * pB[4] + pA[11] * pB[8] + pA[15] * pB[12];
		pOut[13]	= pA[3] * pB[1] + pA[7] * pB[5] + pA[11] * pB[9] + pA[15] * pB[13];
		pOut[14]	= pA[3] * pB[2] + pA[7] * pB[6] + pA[11] * pB[10] + pA[15] * pB[14];
		pOut[15]	= pA[3] * pB[3] + pA[7] * pB[7] + pA[11] * pB[11] + pA[15] * pB[15];
	}
	else if (bTransB)
	{
		pOut[0]		= pA[0] * pB[0] + pA[1] * pB[1] + pA[2] * pB[2] + pA[3] * pB[3];
		pOut[1]		= pA[0] * pB[4] + pA[1] * pB[5] + pA[2] * pB[6] + pA[3] * pB[7];
		pOut[2]		= pA[0] * pB[8] + pA[1] * pB[9] + pA[2] * pB[10] + pA[3] * pB[11];
		pOut[3]		= pA[0] * pB[12] + pA[1] * pB[13] + pA[2] * pB[14] + pA[3] * pB[15];
		pOut[4]		= pA[4] * pB[0] + pA[5] * pB[1] + pA[6] * pB[2] + pA[7] * pB[3];
		pOut[5]		= pA[4] * pB[4] + pA[5] * pB[5] + pA[6] * pB[6] + pA[7] * pB[7];
		pOut[6]		= pA[4] * pB[8] + pA[5] * pB[9] + pA[6] * pB[10] + pA[7] * pB[11];
		pOut[7]		= pA[4] * pB[12] + pA[5] * pB[13] + pA[6] * pB[14] + pA[7] * pB[15];
		pOut[8]		= pA[8] * pB[0] + pA[9] * pB[1] + pA[10] * pB[2] + pA[11] * pB[3];
		pOut[9]		= pA[8] * pB[4] + pA[9] * pB[5] + pA[10] * pB[6] + pA[11] * pB[7];
		pOut[10]	= pA[8] * pB[8] + pA[9] * pB[9] + pA[10] * pB[10] + pA[11] * pB[11];
		pOut[11]	= pA[8] * pB[12] + pA[9] * pB[13] + pA[10] * pB[14] + pA[11] * pB[15];
		pOut[12]	= pA[12] * pB[0] + pA[13] * pB[1] + pA[14] * pB[2] + pA[15] * pB[3];
		pOut[13]	= pA[12] * pB[4] + pA[13] * pB[5] + pA[14] * pB[6] + pA[15] * pB[7];
		pOut[14]	= pA[12] * pB[8] + pA[13] * pB[9] + pA[14] * pB[10] + pA[15] * pB[11];
		pOut[15]	= pA[12] * pB[12] + pA[13] * pB[13] + pA[14] * pB[14] + pA[15] * pB[15];
	}
	else
	{
		pOut[0]		= pA[0] * pB[0] + pA[1] * pB[4] + pA[2] * pB[8] + pA[3] * pB[12];
		pOut[1]		= pA[0] * pB[1] + pA[1] * pB[5] + pA[2] * pB[9] + pA[3] * pB[13];
		pOut[2]		= pA[0] * pB[2] + pA[1] * pB[6] + pA[2] * pB[10] + pA[3] * pB[14];
		pOut[3]		= pA[0] * pB[3] + pA[1] * pB[7] + pA[2] * pB[11] + pA[3] * pB[15];
		pOut[4]		= pA[4] * pB[0] + pA[5] * pB[4] + pA[6] * pB[8] + pA[7] * pB[12];
		pOut[5]		= pA[4] * pB[1] + pA[5] * pB[5] + pA[6] * pB[9] + pA[7] * pB[13];
		pOut[6]		= pA[4] * pB[2] + pA[5] * pB[6] + pA[6] * pB[10] + pA[7] * pB[14];
		pOut[7]		= pA[4] * pB[3] + pA[5] * pB[7] + pA[6] * pB[11] + pA[7] * pB[15];
		pOut[8]		= pA[8] * pB[0] + pA[9] * pB[4] + pA[10] * pB[8] + pA[11] * pB[12];
		pOut[9]		= pA[8] * pB[1] + pA[9] * pB[5] + pA[10] * pB[9] + pA[11] * pB[13];
		pOut[10]	= pA[8] * pB[2] + pA[9] * pB[6] + pA[10] * pB[10] + pA[11] * pB[14];
		pOut[11]	= pA[8] * pB[3] + pA[9] * pB[7] + pA[10] * pB[11] + pA[11] * pB[15];
		pOut[12]	= pA[12] * pB[0] + pA[13] * pB[4] + pA[14] * pB[8] + pA[15] * pB[12];
		pOut[13]	= pA[12] * pB[1] + pA[13] * pB[5] + pA[14] * pB[9] + pA[15] * pB[13];
		pOut[14]	= pA[12] * pB[2] + pA[13] * pB[6] + pA[14] * pB[10] + pA[15] * pB[14];
		pOut[15]	= pA[12] * pB[3] + pA[13] * pB[7] + pA[14] * pB[11] + pA[15] * pB[15];
	}
}


///////////////////////////////////////////////////////////////////////////////
//functions for reading custom point cloud data formats
///////////////////////////////////////////////////////////////////////////////

void readPointCloudFile(std::string strFile, std::vector<SStereoPatch>& patches, float dataScale, bool& bNormals)
{
	int i, nPoints;

	abutil::C3Vectorf vpoint;
	unsigned char color[4];

	SStereoPatch patch;

	std::string strExt;
	std::ifstream file;

	size_t extpos = strFile.find_last_of('.') + 1;
	if (extpos < strFile.length())
	{
		strExt = strFile.substr(extpos);
	}

	bNormals = false;

	if (strcmp(strExt.c_str(), "ptcl") == 0)
	{
		file.open(strFile.c_str(), std::ios_base::binary | std::ios_base::in);
		if (!file.good())
		{
			std::cout << "readPointCloudFile(...): unable to open point cloud file " << strFile.c_str() << std::endl;
			return;
		}

		file.read((char*)&nPoints, sizeof(int));
		if (!file.good())
		{
			std::cout << "readPointCloudFile(...): error reading number of points from file " << strFile.c_str() << std::endl;
			return;
		}

		if (nPoints <= 0)
		{
			std::cout << "readPointCloudFile(...): invalid number of points " << nPoints << " in file " << strFile.c_str() << std::endl;
			return;
		}

		for (i = 0; i < nPoints && !file.eof(); i++)
		{
			file.read((char*)&vpoint, sizeof(abutil::C3Vectorf));
			if (!file.good())
			{
				std::cout << "readPointCloudFile(...): error reading position of point " << i << std::endl;
				break;
			}

			file.read((char*)color, 4 * sizeof(unsigned char));
			if (!file.good())
			{
				std::cout << "readPointCloudFile(...): error reading color of point " << i << std::endl;
				break;
			}

			patch.px = vpoint.x;
			patch.py = vpoint.y;
			patch.pz = vpoint.z;

			patch.r = (float)color[0] / 255.f;
			patch.g = (float)color[1] / 255.f;
			patch.b = (float)color[2] / 255.f;

			patch.nx = 0.f;
			patch.ny = 0.f;
			if (patch.pz > 0.f)
				patch.nz = -1.f;
			else
				patch.nz = 1.f;
			
			patch.tx = 1.f;
			patch.ty = 0.f;
			patch.tz = 0.f;

			patch.bx = 0.f;
			if (patch.pz > 0.f)
				patch.by = -1.f;
			else
				patch.by = 1.f;
			patch.bz = 0.f;

			patch.u = 0.f;
			patch.v = 0.f;

			patch.bLandmark = false;

			patches.push_back(patch);
		}

		file.close();
	}
	else
	{
		std::cout << "readPointCloudFile(...): unrecognized point cloud filename extension " << strExt << std::endl;
	}
}



