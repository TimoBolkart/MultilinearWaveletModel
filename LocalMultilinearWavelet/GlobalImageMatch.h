///////////////////////////////////////////////////////////////////////////////
//
//	GlobalImageMatch.h
//
//	Header file for the CGlobalImageMatch class
//
//	Alan Brunton
//
///////////////////////////////////////////////////////////////////////////////


#ifndef __GLOBALIMAGEMATCH_H__
#define __GLOBALIMAGEMATCH_H__


#include "StereoFace.h"


class CGlobalImageMatch
{
private:
	//surface member variables
	int										m_nVertices;
	float*									m_pVertices;
	float*									m_pNormals;
	int										m_nIndices;
	GLuint*									m_pIndices;
	GLenum									m_ePrimType;

	//surface member variables using buffer objects
	bool									m_bUseBufferObjects;
	bool									m_bDeleteBufferObjects;
	GLuint									m_nbufVertex;
	GLuint									m_nbufNormal;
	GLuint									m_nbufIndex;

public:

	CGlobalImageMatch()
	{
		setSurface(0, 0, GL_NONE, NULL, NULL);
		m_bUseBufferObjects = false;
		m_bDeleteBufferObjects = false;
		m_nbufVertex = 0;
		m_nbufNormal = 0;
		m_nbufIndex = 0;
	}

	~CGlobalImageMatch()
	{
		clearSurface();
	}

	void setSurface(int nVertices, int nIndices, GLenum eType, float* pVertices, GLuint* pIndices, bool bUseBufferObjects = false);
	void setSurface(int nVertices, int nIndices, GLenum eType, float* pVertices, float* pNormals, GLuint* pIndices, bool bUseBufferObjects = false);

	GLuint getVertexBuffer()				{ return m_nbufVertex; }
	GLuint getNormalBuffer()				{ return m_nbufNormal; }

	void clearSurface();

	///////////////////////////////////////////////////////////////////////////
	//public static members
	///////////////////////////////////////////////////////////////////////////

	static const GLuint						gim_indexMaskValue	= (GLuint)(-1);
};


#endif //__GLOBALIMAGEMATCH_H__

