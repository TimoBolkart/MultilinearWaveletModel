///////////////////////////////////////////////////////////////////////////////
//
//	WaveletShapeFitter.h
//
//	Base class for set of classes that fit a wavelet-based statistical model
//	to point cloud data, or possibly other data
//
//	Alan Brunton, January 2013
//
///////////////////////////////////////////////////////////////////////////////


#ifndef __WAVELETSHAPEFITTER_H__
#define __WAVELETSHAPEFITTER_H__


#include "StereoFace.h"
#include "BSGWFactory.h"
#include "Profile.h"

class CWaveletShapePrior
{
public:

	CBSplineGridWavelet<C3Vectorf>*			m_pMean;
	CBSplineGridWavelet<C3Vectorf>*			m_pStdDev;
	float*									m_pmatReductTranspose;

	float									m_weight;

	//Parameter of how many standard deviations we sample from the mean
	float									m_sampleRegion;

	CWaveletShapePrior()
	{
		m_pMean = NULL;
		m_pStdDev = NULL;
		m_pmatReductTranspose = NULL;
		m_weight = 1.0;
		m_sampleRegion = 1.0;
	}

	void inverseTransformCoefficient(int iCoeff, const C3Vectorf& vb, C3Vectorf& vc)
	{
		_ASSERT(m_pMean != NULL);

		vc = m_pMean->getCoefficient(iCoeff);

		if (m_pmatReductTranspose != NULL)
		{
			int imat = iCoeff * 9;
			vc.x += m_pmatReductTranspose[imat] * vb.x + m_pmatReductTranspose[imat + 1] * vb.y + m_pmatReductTranspose[imat + 2] * vb.z;
			vc.y += m_pmatReductTranspose[imat + 3] * vb.x + m_pmatReductTranspose[imat + 4] * vb.y + m_pmatReductTranspose[imat + 5] * vb.z;
			vc.z += m_pmatReductTranspose[imat + 6] * vb.x + m_pmatReductTranspose[imat + 7] * vb.y + m_pmatReductTranspose[imat + 8] * vb.z;
		}
		else
			vc += vb;
	}

	void inverseTransformSamples(int iCoeff, int iComponent, int nSamples, float* pSamples)
	{
		_ASSERT(m_pMean != NULL);
		_ASSERT(m_pmatReductTranspose == NULL);

		int i;
		float mu;
		C3Vectorf vmean = m_pMean->getCoefficient(iCoeff);

		switch (iComponent)
		{
		case 0:
			mu = vmean.x;
			break;
		case 1:
			mu = vmean.y;
			break;
		case 2:
			mu = vmean.z;
			break;
		}

		for (i = 0; i < nSamples; i++)
			pSamples[i] += mu;
	}

	void forwardTransformCoefficient(int iCoeff, const C3Vectorf& vc, C3Vectorf& vb)
	{
		_ASSERT(m_pMean != NULL);

		vb = vc - m_pMean->getCoefficient(iCoeff);

		if (m_pmatReductTranspose != NULL)
		{
			int imat = iCoeff * 9;
			vb =	C3Vectorf(m_pmatReductTranspose[imat], m_pmatReductTranspose[imat + 1], m_pmatReductTranspose[imat + 2]) * vb.x +
					C3Vectorf(m_pmatReductTranspose[imat + 3], m_pmatReductTranspose[imat + 4], m_pmatReductTranspose[imat + 5]) * vb.y +
					C3Vectorf(m_pmatReductTranspose[imat + 6], m_pmatReductTranspose[imat + 7], m_pmatReductTranspose[imat + 8]) * vb.z;
		}
	}

	void forwardTransformSamples(int iCoeff, int iComponent, int nSamples, float* pSamples)
	{
		_ASSERT(m_pMean != NULL);
		_ASSERT(m_pmatReductTranspose == NULL);

		int i;
		float mu;
		C3Vectorf vmean = m_pMean->getCoefficient(iCoeff);

		switch (iComponent)
		{
		case 0:
			mu = vmean.x;
			break;
		case 1:
			mu = vmean.y;
			break;
		case 2:
			mu = vmean.z;
			break;
		}

		for (i = 0; i < nSamples; i++)
			pSamples[i] -= mu;
	}

	void evaluateSamples(int iCoeff, int iComponent, int nSamples, float* pSamples, float* pCost, float* pTotalPriorCost = NULL)
	{
		_ASSERT(m_pStdDev != NULL);
		_ASSERT(iCoeff < m_pStdDev->getNumCoefficients());
		_ASSERT(iComponent >= 0 && iComponent < 3);

		int i;
		float sigma, rcpTwoSigma2;
		C3Vectorf vstd = m_pStdDev->getCoefficient(iCoeff);

		switch (iComponent)
		{
		case 0:
			sigma = vstd.x;
			break;
		case 1:
			sigma = vstd.y;
			break;
		case 2:
			sigma = vstd.z;
			break;
		}

		rcpTwoSigma2 = 1.f / (2.f * sigma*sigma);

		if (pTotalPriorCost != NULL)
		{
			for (i = 0; i < nSamples; i++)
			{
				pCost[i] = pSamples[i]*pSamples[i] * rcpTwoSigma2 * m_weight;
				pTotalPriorCost[i] += pCost[i];
			}
		}
		else
		{
			for (i = 0; i < nSamples; i++)
				pCost[i] = pSamples[i]*pSamples[i] * rcpTwoSigma2 * m_weight;
		}
	}
};


class CWaveletShapeObservation
{
public:

	int										m_nVertices;
	int										m_nWidth, m_nHeight;//in case the geometry is organized in a grid
	C3Vectorf*								m_pGeometry;
	float*									m_pGeometryMask;

	float									m_weight;

	CWaveletShapeObservation()
	{
		m_nVertices = m_nWidth = m_nHeight = 0;
		m_pGeometry = NULL;
		m_pGeometryMask = NULL;
		m_weight = 1.f;
	}

	virtual void init() = 0;
	virtual float computeDataCost() = 0;
	virtual void updateObserveration() = 0;
	virtual void reportProfiling(FILE* pfOutput) = 0;
};


class CWaveletShapeFitter
{
protected:

	std::vector<CWaveletShapePrior*>		m_priors;
	std::vector<CWaveletShapeObservation*>	m_observations;
	std::vector<CWaveletShapeObservation*>	m_refinements;

	CBSplineGridWavelet<C3Vectorf>*			m_pIntermediateModel;
	C3Vectorf*								m_pReconVerts;
	float*									m_pVertexMask;

	int										m_iStartCoefficient;
	int										m_level;


	virtual void initActiveWaveletModel() = 0;
	virtual void clear() = 0;

	virtual void reconstructIntermediateModel() = 0;
	virtual void reconstructIntermediateModel(int iCoeff) = 0;

	inline void reconstructGeometry()
	{
		m_pIntermediateModel->reconstructCopy(m_pIntermediateModel->getNumLevels() - 1, m_pReconVerts, m_pIntermediateModel->getFullResWidth(), 1);
	}


public:

	CWaveletShapeFitter()
	{
		m_pIntermediateModel = NULL;
		m_pReconVerts = NULL;
		m_pVertexMask = NULL;
		m_iStartCoefficient = 0;
		m_level = 0;
	}
	virtual ~CWaveletShapeFitter()
	{
	}

	virtual void resetFitting(const int level = 0)
	{
		if (level >= 0 && m_pIntermediateModel != NULL)
		{
			m_level = level;
			m_iStartCoefficient = level == 0 ? 0 : m_pIntermediateModel->getLevelHeight(level - 1) * m_pIntermediateModel->getLevelWidth(level - 1);
		}
	}

	virtual void addPrior(CWaveletShapePrior* pPrior) 
	{
		_ASSERT(pPrior != NULL);
		m_priors.push_back(pPrior);
		if (m_priors.size() == 1)
		{
			//this is the primary prior
			initActiveWaveletModel();
		}
	}
	virtual void addObservation(CWaveletShapeObservation* pObs) = 0;

	virtual void copyIntermediateWaveletModel(CBSplineGridWavelet<C3Vectorf>* pCopy)
	{
		_ASSERT(pCopy != NULL);
		_ASSERT(m_pIntermediateModel != NULL);
#if WSS_GPU
		((CCUDABSplineGridWavelet3f*)m_pIntermediateModel)->download();
#endif //WSS_GPU
		pCopy->copy(m_pIntermediateModel);
	}
	
	virtual void getInternalReconstruction(abutil::C3Vectorf* pReconVerts)
	{
		memcpy(pReconVerts, m_pReconVerts, m_pIntermediateModel->getNumCoefficients() * sizeof(abutil::C3Vectorf));
	}

	virtual void optimizeActiveWaveletModel() = 0;

	virtual void reportProfiling(FILE* pfOutput) = 0;
};


#endif //__WAVELETSHAPEFITTER_H__


