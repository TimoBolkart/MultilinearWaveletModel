///////////////////////////////////////////////////////////////////////////////
//
//	NearestNeighborAssistant.cpp
//
//	Source file for the CNearestNeighborAssistant class
//
//	Alan Brunton, September 2010
//
///////////////////////////////////////////////////////////////////////////////


#include "NearestNeighborAssistant.h"

///////////////////////////////////////////////////////////////////////////////
//member functions
///////////////////////////////////////////////////////////////////////////////

void CNearestNeighborAssistant::clearInternalMemory()
{
	int i;

	if (m_pkdtree != NULL)
	{
		delete m_pkdtree;
		m_pkdtree = NULL;
	}
	if (m_pRefPoints != NULL)
	{
		for (i = 0; i < m_nRefPoints; i++)
		{
			if (m_pRefPoints[i] != NULL)
			{
				delete [] m_pRefPoints[i];
			}
		}
		delete [] m_pRefPoints;
		m_pRefPoints = NULL;
	}
	if (m_pQueryPoints != NULL)
	{
		delete [] m_pQueryPoints;
		m_pQueryPoints = NULL;
	}
}

void CNearestNeighborAssistant::setReferencePoints(int nPoints, NNAReal *pPoints)
{
	_ASSERT(pPoints != NULL);

	int i, j, iread, iwrite;

	//setting new reference points has the effect of essentially resetting the class
	//setQueryPoints must be called again
	clearInternalMemory();

	m_nRefPoints = nPoints;

	m_pRefPoints = new ANNpoint[m_nRefPoints];
	_ASSERT(m_pRefPoints != NULL);

	for (i = 0; i < m_nRefPoints; i++)
	{
		m_pRefPoints[i] = new ANNcoord[m_dim];
		_ASSERT(m_pRefPoints[i] != NULL);

		for (j = 0; j < m_dim; j++)
		{
			m_pRefPoints[i][j] = pPoints[i * m_dim + j];
		}
	}

	m_pkdtree = new ANNkd_tree(m_pRefPoints, m_nRefPoints, m_dim);
	_ASSERT(m_pkdtree != NULL);
}

void CNearestNeighborAssistant::setQueryPoints(int nPoints, NNAReal *pPoints, NNAReal* pDist, NNAIndex* pIndices, int nNNPerQuery)
{
	_ASSERT(pDist != NULL);
	_ASSERT(pIndices != NULL);
	_ASSERT(pPoints != NULL);

	int i, j, iread, iwrite;

	m_pDistanceOutput = pDist;
	m_pIndexOutput = pIndices;

	m_nNNPerQuery = nNNPerQuery;

	if (m_pQueryPoints != NULL && nPoints != m_nQueryPoints)
	{
		delete [] m_pQueryPoints;
		m_pQueryPoints = NULL;
	}

	if (m_pQueryPoints == NULL)
	{
		m_nQueryPoints = nPoints;
		m_pQueryPoints = new ANNcoord[m_dim * m_nQueryPoints];
		_ASSERT(m_pQueryPoints != NULL);
	}

	for (i = 0; i < m_nQueryPoints; i++)
	{
		for (j = 0; j < m_dim; j++)
		{
			m_pQueryPoints[i * m_dim + j] = pPoints[i * m_dim + j];
		}
	}
}

void CNearestNeighborAssistant::compute()
{
	_ASSERT(m_pRefPoints != NULL);
	_ASSERT(m_pQueryPoints != NULL);
	_ASSERT(m_pDistanceOutput != NULL);
	_ASSERT(m_pIndexOutput != NULL);

	ANNdist* pDistances = new ANNdist[m_nQueryPoints * m_nNNPerQuery];
	_ASSERT(pDistances != NULL);

	int i, k;
	for (i = 0; i < m_nQueryPoints; i++)
	{
		m_pkdtree->annkSearch(m_pQueryPoints + i * m_dim, m_nNNPerQuery, m_pIndexOutput + i * m_nNNPerQuery, pDistances + i * m_nNNPerQuery);

		for (k = 0; k < m_nNNPerQuery; k++)
			m_pDistanceOutput[i * m_nNNPerQuery + k] = sqrt(pDistances[i * m_nNNPerQuery + k]);
	}

	delete [] pDistances;
}

