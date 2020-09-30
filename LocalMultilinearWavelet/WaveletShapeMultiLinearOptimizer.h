///////////////////////////////////////////////////////////////////////////////
//
//	WaveletShapeMultiLinearOptimizer.h
//
//	Header file for the CWaveletShapeMultiLinearOptimizer class
//	optimizes a local multi-linear model on wavelet coefficients to fit a shape to point cloud or other data
//
//	Alan Brunton, January 2013
//
///////////////////////////////////////////////////////////////////////////////


#ifndef __WAVELETSHAPEMULTILINEAROPTIMIZER_H__
#define __WAVELETSHAPEMULTILINEAROPTIMIZER_H__


#include "WaveletShapeFitter.h"
#include "MultilinearModelHandler.h"
#include "NearestNeighborAssistant.h"


class CWaveletShapeMultiLinearOptimizer : public CWaveletShapeFitter
{
protected:

	size_t											m_nBaseWidth, m_nBaseHeight, m_nLevels;
	size_t											m_nFullWidth, m_nFullHeight;

	std::vector<CBSplineGridWavelet<C3Vectorf>*>	m_activeWaveletModel;
	size_t											m_nModes;

	bool											m_bFitLandmarksOnly;

	MultilinearModelHandler*						m_pMultilinearModels;
	std::vector<double>*							m_pCurrentWeights;

	std::vector<double>*							m_pTransformCoeffs;
	std::vector<size_t>*							m_pTransformIndices;

	abutil::C4x4Matrixf								m_matModelToData;
	abutil::C4x4Matrixf								m_matDataToModel;
	abutil::C3Vectorf*								m_pTransformedRecon;
	CNearestNeighborAssistant*						m_pNNA;
	NNAReal*										m_pNNDists;
	NNAIndex*										m_pNNIdx;
	NNAReal*										m_pObservedPoints;
	size_t											m_nObservedPoints;
	abutil::C3Vectorf*								m_pObservedNormals;

	//abutil::C3Vectorf*								m_pPrediction;
	//float											m_predWeight;
	std::vector<abutil::C3Vectorf*>					m_predictions;
	std::vector<float>								m_predWeights;

	float											m_nnthresh;
	float											m_allowDeviation;
	float											m_landmarkRelativeWeight;
	float											m_inPlaneThresh;
	float											m_smoothRelativeWeight;

	size_t											m_nInnerIters;
	size_t											m_nRefineInnerIters;

	std::vector<bool>								m_fixMode;

	size_t											m_nLandmarks;
	abutil::C3Vectorf*								m_pModelLandmarks;
	abutil::C3Vectorf*								m_pDataLandmarks;
	size_t*											m_pLandmarkVerts;
	std::vector<double>								m_landmarkWeights;

	std::vector<std::vector<size_t>>		m_gridNeighbors;
	std::vector<bool>							m_invalidSmoothVertices;
	double										m_neighborSmoothWeight;
	std::vector<std::vector<size_t>>		m_precomputedSmoothIndices;
	std::vector<std::vector<double>>		m_precomputedSmoothWeights;

	CProfile										m_profCoeffOpt;
	CProfile										m_profObsUpdate;
	CProfile										m_profLevelOpt;
	CProfile										m_profLoadModel;
	CProfile										m_profRefine;
	
	void initActiveWaveletModel();
	void clear();

	void reconstructIntermediateModel() {}
	void reconstructIntermediateModel(int iCoeff) {}
	
	void updateObservation(std::vector<double>& obsData, std::vector<double>& obsMask);
	
	void updateObservations(size_t& numObs, std::vector<double>& obsWeights, std::vector<double>& obsData, std::vector<double>& obsMask);

	void computeGridNeighbors(std::vector<std::vector<size_t>>& neighbors);

public:
	
	CWaveletShapeMultiLinearOptimizer(): CWaveletShapeFitter()
	{
		m_nBaseWidth = m_nBaseHeight = m_nLevels = 0;
		m_nFullWidth = m_nFullHeight = 0;
		m_nModes = 0;
		m_bFitLandmarksOnly = false;

		m_pMultilinearModels = NULL;
		m_pCurrentWeights = NULL;

		m_pTransformCoeffs = NULL;
		m_pTransformIndices = NULL;

		m_matModelToData.setIdentity();
		m_matDataToModel.setIdentity();
		m_pTransformedRecon = NULL;
		m_pNNA = NULL;
		m_pNNDists = NULL;
		m_pNNIdx = NULL;
		m_pObservedPoints = NULL;
		m_nObservedPoints = 0;
		m_pObservedNormals = NULL;

		m_nnthresh = 1.f;
		m_allowDeviation = 1.f;
		m_landmarkRelativeWeight = 1.f;
		m_inPlaneThresh = 1.f;
		m_smoothRelativeWeight = 1.f;

		m_nInnerIters = 1;
		m_nRefineInnerIters = 1;

		m_nLandmarks = 0;
		m_pModelLandmarks = NULL;
		m_pDataLandmarks = NULL;
		m_pLandmarkVerts = NULL;

		m_neighborSmoothWeight = 0.0;

		m_profCoeffOpt.setName("optimize individual coefficient");
		m_profLevelOpt.setName("optimize full level");
		m_profObsUpdate.setName("update observation");
		m_profLoadModel.setName("load model");
		m_profRefine.setName("refine surface");
	}

	void setNumModes(size_t nModes)
	{
		clear();
		m_nModes = nModes;
		m_activeWaveletModel.resize(m_nModes, NULL);
	}
	
	void loadModels(std::string strModelFile);

	void loadModelsText(std::string strModelFile);

	void addObservation(CWaveletShapeObservation* pObs)
	{
		_ASSERT(pObs != NULL);
		pObs->m_nVertices = m_activeWaveletModel[0]->getNumCoefficients();
		pObs->m_nWidth = m_activeWaveletModel[0]->getFullResWidth();
		pObs->m_nHeight = m_activeWaveletModel[0]->getFullResHeight();
		pObs->m_pGeometry = m_pReconVerts;
		pObs->m_pGeometryMask = m_pVertexMask;
		pObs->init();
		m_observations.push_back(pObs);
	}

	float getAllowedDeviation()	const		{ return m_allowDeviation; }
	
	void setAllowedDeviation(float dev)		{ m_allowDeviation = fabs(dev); }

	float getLandmarkRelativeWeight() const	{ return m_landmarkRelativeWeight; }
	void setLandmarkRelativeWeight(float w)	{ m_landmarkRelativeWeight = fabs(w); }

	float getSmoothRelativeWeight() const	{ return m_smoothRelativeWeight; }
	
	void setSmoothRelativeWeight(float w)	{ m_smoothRelativeWeight = fabs(w); }

	size_t getNumInnerIters() const			{ return m_nInnerIters; }
	void setNumInnerIters(size_t iters)		{ m_nInnerIters = iters; }

	size_t getNumRefineInnerIters() const	{ return m_nRefineInnerIters; }
	
	void setNumRefineInnerIters(size_t n)	{ m_nRefineInnerIters = n; }

	bool isModeFixed(size_t mode)
	{
		if (mode == 2)	return m_fixMode[0];
		if (mode == 3)	return m_fixMode[1];
		return false;
	}
	void fixMode(size_t mode, bool bfix)
	{
		if (mode == 2)	m_fixMode[0] = bfix;
		if (mode == 3)	m_fixMode[1] = bfix;
	}

	void setNNA(CNearestNeighborAssistant* pNNA, NNAReal* pObsPoints, size_t nObsPoints, double* pmatModelToData, double* pmatDataToModel, float distThresh)
	{
		m_pNNA = pNNA;
		m_pObservedPoints = pObsPoints;
		m_nObservedPoints = nObsPoints;
		m_nnthresh = distThresh;
		m_pObservedNormals = NULL;

		m_matModelToData._11 = pmatModelToData[0];	m_matModelToData._12 = pmatModelToData[1];	m_matModelToData._13 = pmatModelToData[2];	m_matModelToData._14 = pmatModelToData[3];
		m_matModelToData._21 = pmatModelToData[4];	m_matModelToData._22 = pmatModelToData[5];	m_matModelToData._23 = pmatModelToData[6];	m_matModelToData._24 = pmatModelToData[7];
		m_matModelToData._31 = pmatModelToData[8];	m_matModelToData._32 = pmatModelToData[9];	m_matModelToData._33 = pmatModelToData[10];	m_matModelToData._34 = pmatModelToData[11];
		m_matModelToData._41 = pmatModelToData[12];	m_matModelToData._42 = pmatModelToData[13];	m_matModelToData._43 = pmatModelToData[14];	m_matModelToData._44 = pmatModelToData[15];

		m_matDataToModel._11 = pmatDataToModel[0];	m_matDataToModel._12 = pmatDataToModel[1];	m_matDataToModel._13 = pmatDataToModel[2];	m_matDataToModel._14 = pmatDataToModel[3];
		m_matDataToModel._21 = pmatDataToModel[4];	m_matDataToModel._22 = pmatDataToModel[5];	m_matDataToModel._23 = pmatDataToModel[6];	m_matDataToModel._24 = pmatDataToModel[7];
		m_matDataToModel._31 = pmatDataToModel[8];	m_matDataToModel._32 = pmatDataToModel[9];	m_matDataToModel._33 = pmatDataToModel[10];	m_matDataToModel._34 = pmatDataToModel[11];
		m_matDataToModel._41 = pmatDataToModel[12];	m_matDataToModel._42 = pmatDataToModel[13];	m_matDataToModel._43 = pmatDataToModel[14];	m_matDataToModel._44 = pmatDataToModel[15];
	}

	
	void setNNA(CNearestNeighborAssistant* pNNA, NNAReal* pObsPoints, NNAReal* pObsNormals, size_t nObsPoints, double* pmatModelToData, double* pmatDataToModel, float distThresh, float inPlaneThresh)
	{
		m_pNNA = pNNA;
		m_pObservedPoints = pObsPoints;
		m_nObservedPoints = nObsPoints;
		m_nnthresh = distThresh;
		m_inPlaneThresh = inPlaneThresh;

		m_matModelToData._11 = pmatModelToData[0];	m_matModelToData._12 = pmatModelToData[1];	m_matModelToData._13 = pmatModelToData[2];	m_matModelToData._14 = pmatModelToData[3];
		m_matModelToData._21 = pmatModelToData[4];	m_matModelToData._22 = pmatModelToData[5];	m_matModelToData._23 = pmatModelToData[6];	m_matModelToData._24 = pmatModelToData[7];
		m_matModelToData._31 = pmatModelToData[8];	m_matModelToData._32 = pmatModelToData[9];	m_matModelToData._33 = pmatModelToData[10];	m_matModelToData._34 = pmatModelToData[11];
		m_matModelToData._41 = pmatModelToData[12];	m_matModelToData._42 = pmatModelToData[13];	m_matModelToData._43 = pmatModelToData[14];	m_matModelToData._44 = pmatModelToData[15];

		m_matDataToModel._11 = pmatDataToModel[0];	m_matDataToModel._12 = pmatDataToModel[1];	m_matDataToModel._13 = pmatDataToModel[2];	m_matDataToModel._14 = pmatDataToModel[3];
		m_matDataToModel._21 = pmatDataToModel[4];	m_matDataToModel._22 = pmatDataToModel[5];	m_matDataToModel._23 = pmatDataToModel[6];	m_matDataToModel._24 = pmatDataToModel[7];
		m_matDataToModel._31 = pmatDataToModel[8];	m_matDataToModel._32 = pmatDataToModel[9];	m_matDataToModel._33 = pmatDataToModel[10];	m_matDataToModel._34 = pmatDataToModel[11];
		m_matDataToModel._41 = pmatDataToModel[12];	m_matDataToModel._42 = pmatDataToModel[13];	m_matDataToModel._43 = pmatDataToModel[14];	m_matDataToModel._44 = pmatDataToModel[15];

		m_pObservedNormals = new abutil::C3Vectorf[m_nObservedPoints];
		for (size_t i = 0; i < m_nObservedPoints; i++)
		{
			const size_t j = i * 3;
			m_pObservedNormals[i].set(pObsNormals + j);
		}
	}

	
	void setVertexMask(float* pMask)
	{
		if (pMask != NULL && m_pVertexMask != NULL && m_pIntermediateModel != NULL)
			memcpy(m_pVertexMask, pMask, m_pIntermediateModel->getNumCoefficients() * sizeof(float));
	}
	
	void setLandmarks(const size_t numLandmarks, const std::vector<size_t>& modelLandmarks, const std::vector<double>& dataLandmarks)
	{
		size_t i;

		if (m_pIntermediateModel == NULL || m_pReconVerts == NULL)
			return;

		if(numLandmarks != modelLandmarks.size() || 3*numLandmarks != dataLandmarks.size())
		{
			return;
		}

		delete [] m_pDataLandmarks;
		delete [] m_pModelLandmarks;
		m_nLandmarks = numLandmarks;

		m_pDataLandmarks = new abutil::C3Vectorf[m_nLandmarks];
		m_pModelLandmarks = new abutil::C3Vectorf[m_nLandmarks];
		m_pLandmarkVerts = new size_t[m_nLandmarks];
		m_landmarkWeights.resize(m_nLandmarks);

		for (i = 0; i < numLandmarks; i++)
		{
			m_pModelLandmarks[i] = m_pReconVerts[modelLandmarks[i]];
			m_pDataLandmarks[i].set(dataLandmarks[i * 3 + 0], dataLandmarks[i * 3 + 1], dataLandmarks[i * 3 + 2]);
			m_pLandmarkVerts[i] = modelLandmarks[i];
			m_landmarkWeights[i] = 1.0;
		}
	}

	void setLandmarks(const size_t numLandmarks, const size_t* pModelLandmarks, const double* pDataLandmarks, const double* pLandmarkWeights)
	{
		size_t i;

		if (m_pIntermediateModel == NULL || m_pReconVerts == NULL)
			return;

		delete [] m_pDataLandmarks;
		delete [] m_pModelLandmarks;
		m_nLandmarks = numLandmarks;

		m_pDataLandmarks = new abutil::C3Vectorf[m_nLandmarks];
		m_pModelLandmarks = new abutil::C3Vectorf[m_nLandmarks];
		m_pLandmarkVerts = new size_t[m_nLandmarks];
		m_landmarkWeights.resize(m_nLandmarks);

		for (i = 0; i < numLandmarks; i++)
		{
			m_pModelLandmarks[i] = m_pReconVerts[pModelLandmarks[i]];
			m_pDataLandmarks[i].set(pDataLandmarks[i * 3 + 0], pDataLandmarks[i * 3 + 1], pDataLandmarks[i * 3 + 2]);
			m_pLandmarkVerts[i] = pModelLandmarks[i];
			m_landmarkWeights[i] = pLandmarkWeights[i];
		}
	}

	//void setPrediction(abutil::C3Vectorf* pPrediction = NULL, const float weight = 1.f)
	//{
	//	m_pPrediction = pPrediction;
	//	m_predWeight = weight;
	//}
	void setPredictions(const std::vector<abutil::C3Vectorf*>& predictions = std::vector<abutil::C3Vectorf*>(NULL), const std::vector<float>& weights = std::vector<float>(0))
	{
		m_predictions.clear();
		m_predWeights.clear();
		const size_t numPreds = predictions.size();

		m_predictions.resize(numPreds);
		m_predWeights.resize(numPreds);
		for (size_t i = 0; i < numPreds; i++)
		{
			m_predictions[i] = predictions[i];
			m_predWeights[i] = weights[i];
		}
	}
	
	void getInternalReconstructionDataCoords(abutil::C3Vectorf* pTransformedVerts)
	{
		if (pTransformedVerts != NULL && m_pReconVerts != NULL && m_pIntermediateModel != NULL)
		{
			const size_t nverts = m_pIntermediateModel->getNumCoefficients();
			for (size_t i = 0; i < nverts; i++)
			{
				pTransformedVerts[i] = m_matModelToData * m_pReconVerts[i];
			}
		}
	}

	bool isFitToLandmarksOnly()						{ return m_bFitLandmarksOnly; }
	
	void setFitToLandmarksOnly(bool bLandmarksOnly)	{ m_bFitLandmarksOnly = bLandmarksOnly; }
	
	void setNeighborSmoothWeight(double smoothWeight) { m_neighborSmoothWeight = smoothWeight; }

	void optimizeActiveWaveletModel();
	
	void refineSurface();
	
	void refineSurface(abutil::C3Vectorf* pRefinedVerts);

	void highlightCoefficient(const abutil::C3Vectorf& baseColor, const abutil::C3Vectorf& highlightColor, const size_t icoeff, abutil::C3Vectorf* pResult);
	
	void computeGridNeighbors()
	{
		computeGridNeighbors(m_gridNeighbors);
	}

	void precomputeSmoothingWeights();

	void reportProfiling(FILE* pfOutput);
};


#endif //__WAVELETSHAPEMULTILINEAROPTIMIZER_H__


