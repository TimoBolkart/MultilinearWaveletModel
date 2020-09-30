///////////////////////////////////////////////////////////////////////////////
//
//	NearestNeighborAssistant.h
//
//	Header file for the CNearestNeighborAssistant class
//
//	Alan Brunton, September 2010
//
///////////////////////////////////////////////////////////////////////////////


#ifndef __NEARESTNEIGHBORASSISTANT_H__
#define __NEARESTNEIGHBORASSISTANT_H__


//#include "StereoFace.h"
#include "ANN/ANN.h"


typedef float NNAReal;
typedef int NNAIndex;


class CNearestNeighborAssistant
{
private:

	int										m_dim;
	int										m_nRefPoints;
	int										m_nQueryPoints;

	ANNpointArray							m_pRefPoints;
	ANNpoint								m_pQueryPoints;

	NNAReal*								m_pDistanceOutput;
	NNAIndex*								m_pIndexOutput;

	ANNkd_tree*								m_pkdtree;

	int										m_nNNPerQuery;
	
	void clearInternalMemory();


public:

	CNearestNeighborAssistant(int dimension = 3)
	{
		m_dim = dimension;
		m_nRefPoints = 0;
		m_nQueryPoints = 0;
		m_pRefPoints = NULL;
		m_pQueryPoints = NULL;
		m_pDistanceOutput = NULL;
		m_pIndexOutput = NULL;
		m_pkdtree = NULL;
		m_nNNPerQuery = 1;
	}

	~CNearestNeighborAssistant()
	{
		clearInternalMemory();
	}
	
	void setReferencePoints(int nPoints, NNAReal* pPoints);
	
	void setQueryPoints(int nPoints, NNAReal* pPoints, NNAReal* pDist, NNAIndex* pIndices, int nNNPerQuery = 1);
	
	void compute();
};


#endif //__NEARESTNEIGHBORASSISTANT_H__


