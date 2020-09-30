///////////////////////////////////////////////////////////////////////////////
//
//	GlobalImageMatch.cpp
//
//	Source file for the CGlobalImageMatch class
//
//	Alan Brunton
//
///////////////////////////////////////////////////////////////////////////////


#include "GlobalImageMatch.h"


///////////////////////////////////////////////////////////////////////////////
//member functions
///////////////////////////////////////////////////////////////////////////////

void CGlobalImageMatch::setSurface(int nVertices, int nIndices, GLenum eType, float *pVertices, GLuint *pIndices, bool bUseBufferObjects)
{
	int i;

//	_ASSERT(pVertices != NULL);
//	_ASSERT(pIndices != NULL);

	m_nVertices = nVertices;
	m_nIndices = nIndices;
	m_ePrimType = eType;
	m_pVertices = pVertices;
	m_pIndices = pIndices;
	m_bUseBufferObjects = bUseBufferObjects;

	if (m_bUseBufferObjects)
	{
		glEnableClientState(GL_VERTEX_ARRAY);

		glGenBuffers(1, &m_nbufVertex);
		glBindBuffer(GL_ARRAY_BUFFER, m_nbufVertex);
		glBufferData(GL_ARRAY_BUFFER, 3 * m_nVertices * sizeof(float), m_pVertices, GL_DYNAMIC_DRAW);

		glBindBuffer(GL_ARRAY_BUFFER, 0);

		glGenBuffers(1, &m_nbufIndex);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_nbufIndex);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, m_nIndices * sizeof(GLuint), NULL, GL_STATIC_DRAW);

		nIndices = 0;
		for (i = 0; i < m_nIndices; i++)
		{
			if (m_pIndices[i] != gim_indexMaskValue)
			{
				glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, nIndices * sizeof(GLuint), sizeof(GLuint), m_pIndices + i);
				nIndices++;
			}
		}

		m_nIndices = nIndices;

		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

		m_bDeleteBufferObjects = true;
	}

}

void CGlobalImageMatch::setSurface(int nVertices, int nIndices, GLenum eType, float* pVertices, float* pNormals, GLuint* pIndices, bool bUseBufferObjects)
{
	int i;

	_ASSERT(pVertices != NULL);
	_ASSERT(pIndices != NULL);
	_ASSERT(pNormals != NULL);

	m_nVertices = nVertices;
	m_nIndices = nIndices;
	m_ePrimType = eType;
	m_pVertices = pVertices;
	m_pNormals = pNormals;
	m_pIndices = pIndices;
	m_bUseBufferObjects = bUseBufferObjects;

	if (m_bUseBufferObjects)
	{
		glEnableClientState(GL_VERTEX_ARRAY);
		glEnableClientState(GL_NORMAL_ARRAY);

		glGenBuffers(1, &m_nbufVertex);
		glBindBuffer(GL_ARRAY_BUFFER, m_nbufVertex);
		glBufferData(GL_ARRAY_BUFFER, 3 * m_nVertices * sizeof(float), m_pVertices, GL_DYNAMIC_DRAW);//GL_STREAM_COPY);//

		glBindBuffer(GL_ARRAY_BUFFER, 0);

		glGenBuffers(1, &m_nbufNormal);
		glBindBuffer(GL_ARRAY_BUFFER, m_nbufNormal);
		glBufferData(GL_ARRAY_BUFFER, 3 * m_nVertices * sizeof(float), m_pNormals, GL_DYNAMIC_DRAW);

		glBindBuffer(GL_ARRAY_BUFFER, 0);

		glGenBuffers(1, &m_nbufIndex);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_nbufIndex);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, m_nIndices * sizeof(GLuint), NULL, GL_STATIC_DRAW);

		nIndices = 0;
		for (i = 0; i < m_nIndices; i++)
		{
			if (m_pIndices[i] != gim_indexMaskValue)
			{
				glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, nIndices * sizeof(GLuint), sizeof(GLuint), m_pIndices + i);
				nIndices++;
			}
		}

		m_nIndices = nIndices;

		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

		m_bDeleteBufferObjects = true;
	}

}

void CGlobalImageMatch::clearSurface()
{
	m_pVertices = NULL;
	m_pNormals = NULL;
	m_pIndices = NULL;
	m_nVertices = 0;
	m_nIndices = 0;
	m_bUseBufferObjects = false;
	if (m_bDeleteBufferObjects)
	{
		if (m_nbufVertex != 0)
		{
			glDeleteBuffers(1, &m_nbufVertex);
			m_nbufVertex = 0;
		}
		if (m_nbufNormal != 0)
		{
			glDeleteBuffers(1, &m_nbufNormal);
			m_nbufNormal = 0;
		}
		if (m_nbufIndex != 0)
		{
			glDeleteBuffers(1, &m_nbufIndex);
			m_nbufIndex = 0;
		}
		m_bDeleteBufferObjects = false;
	}
}