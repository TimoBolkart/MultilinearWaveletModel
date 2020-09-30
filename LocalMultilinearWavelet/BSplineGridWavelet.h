///////////////////////////////////////////////////////////////////////////////
//
//	BSplineGridWavelet.h
//
//	Header file defining the CBSplineGridWavelet class
//	implements generalized b-spline subdivision wavelet [Bertram et al. 2004]
//	on a regular grid
//
//	Alan Brunton, June 2009
//
///////////////////////////////////////////////////////////////////////////////


#ifndef __BSPLINEGRIDWAVELET_H__
#define __BSPLINEGRIDWAVELET_H__


#include <memory.h>


enum EBSGWType
{
	eBSGWFloat = 0,
	eBSGWDouble,
	eBSGWFloat3,
	eBSGWNumTypes
};


template<typename T>
class CBSplineGridWavelet
{
protected:

	int										m_nBaseWidth, m_nBaseHeight;
	int										m_iLevel, m_nLevels;
	int										m_nFullWidth, m_nFullHeight;

	T*										m_pData;

	EBSGWType								m_eType;
	bool									m_useCubicBasis;


public:

	CBSplineGridWavelet(int nBaseWidth, int nBaseHeight, int nLevels, EBSGWType eType):	m_nBaseWidth(nBaseWidth), m_nBaseHeight(nBaseHeight), 
																						m_nLevels(nLevels), m_eType(eType) 
	{
		m_pData = NULL;
		m_iLevel = 0;
		m_useCubicBasis = false;
		init();
	}
	~CBSplineGridWavelet()
	{
		clear();
	}

	void init();
	void clear();

	int getBaseResWidth()					{ return m_nBaseWidth; }
	int getBaseResHeight()					{ return m_nBaseHeight; }
	
	int getFullResWidth()					{ return m_nFullWidth; }
	
	int getFullResHeight()					{ return m_nFullHeight; }

	inline int getLevelWidth(int level)		{ return (1 << level) * m_nBaseWidth - (1 << level) + 1; }
	inline int getLevelHeight(int level)	{ return (1 << level) * m_nBaseHeight - (1 << level) + 1; }

	int getWidth()							{ return getLevelWidth(m_iLevel); }
	int getHeight()							{ return getLevelHeight(m_iLevel); }

	int getStep(int level)					{ return 1 << (m_nLevels - 1 - level); }

	int getNumLevels()						{ return m_nLevels; }
	int getCurrentLevel()					{ return m_iLevel; }
	void setLevel(int level)				{ m_iLevel = level; }

	void getLevelCoordinates(int iCoeff, int level, int* px, int* py)
	{
		_ASSERT(px != NULL);
		_ASSERT(py != NULL);

		int s = getStep(level);
		int y = iCoeff / m_nFullWidth;
		int x = iCoeff % m_nFullWidth;
		x /= s;
		y /= s;

		*px = x;
		*py = y;
	}

	
	virtual void setFullResolutionData(T* pFullResData, int nChannels = 1)
	{
		_ASSERT(m_pData != NULL);
		_ASSERT(pFullResData != NULL);

		if (nChannels == 1)
			memcpy(m_pData, pFullResData, m_nFullWidth * m_nFullHeight * sizeof(T));
		else
		{
			int i, j, nElements = m_nFullWidth * m_nFullHeight;
			for (j = i = 0; i < nElements; i++, j += nChannels)
				m_pData[i] = pFullResData[j];
		}

		m_iLevel = m_nLevels - 1;
	}

	void getCoefficientStats(T* pCoeffMean, T* pCoeffVar, T* pCoeffMin, T* pCoeffMax)
	{
		_ASSERT(pCoeffMin != NULL);
		_ASSERT(pCoeffMax != NULL);

		int i, nData = m_nFullHeight * m_nFullWidth;
		T val;
		T valMin = (T) 1e300, valMax = (T) -1e300;
		T valMean = (T) 0.0, valVar = 0.0;

		for (i = 0; i < nData; i++)
		{
			val = m_pData[i]
			valMean += val;
			if (val < valMin)
				valMin = val;
			if (val > valMax)
				valMax = val;
		}

		valMean /= (T) nData;

		for (i = 0; i < nData; i++)
		{
			val = m_pData[i];
			val -= valMean;
			valVar += val*val;
		}

		valVar /= (T) (nData - 1);

		*pCoeffMean = valMean;
		*pCoeffVar = valVar;
		*pCoeffMin = valMin;
		*pCoeffMax = valMax;
	}

	int getNumCoefficients()				{ return m_nFullHeight * m_nFullWidth; }
	virtual T getCoefficient(int iCoeff)			{ return m_pData[iCoeff]; }
	
	virtual void setCoefficient(int iCoeff, T val)	{ m_pData[iCoeff] = val; }
	T* getCoefficientStore()				{ return m_pData; }

	int getSerialIndex(int iSerial)
	{
		int i, j, r, c;
		int nCurLevel, nPrevLevel = 0;

		for (j = 0; j < m_nLevels; j++)
		{
			nCurLevel = getLevelWidth(j) * getLevelHeight(j);
			if (iSerial < nCurLevel)
				break;
			nPrevLevel = nCurLevel;
		}

		r = 0;
		if (nPrevLevel > 0)
		{
			c = 1;
			for (i = nPrevLevel; i < iSerial; i++)
			{
				if (r % 2 == 1)
					c++;
				else
					c += 2;
				if (c >= getLevelWidth(j))
				{
					r++;
					if (r % 2 == 1)
						c = 0;
					else
						c = 1;
				}
			}
		}
		else
		{
			c = 0;
			for (i = nPrevLevel; i < iSerial; i++)
			{
				c++;
				if (c >= getLevelWidth(j))
				{
					r++;
					c = 0;
				}
			}
		}

		return (r * m_nFullWidth + c) * getStep(j);
	}

	virtual void copy(CBSplineGridWavelet<T>* pbsgw)
	{
		_ASSERT(m_pData != NULL);
		_ASSERT(pbsgw != NULL);
		_ASSERT(m_nFullHeight == pbsgw->m_nFullHeight);
		_ASSERT(m_nFullWidth == pbsgw->m_nFullWidth);
		pbsgw->getRawCopy(m_pData);
		m_iLevel = pbsgw->m_iLevel;
	}

	void scale(float factor)
	{
		_ASSERT(m_pData != NULL);
		int i, nCoeff = m_nFullHeight * m_nFullWidth;
		for (i = 0; i < nCoeff; i++)
			m_pData[i] *= factor;
	}
	void scale(double factor)
	{
		_ASSERT(m_pData != NULL);
		int i, nCoeff = m_nFullHeight * m_nFullWidth;
		for (i = 0; i < nCoeff; i++)
			m_pData[i] *= factor;
	}

	void add(CBSplineGridWavelet<T>* pbsgw)
	{
		_ASSERT(m_pData != NULL);
		_ASSERT(pbsgw != NULL);
		_ASSERT(m_nFullHeight == pbsgw->m_nFullHeight);
		_ASSERT(m_nFullWidth == pbsgw->m_nFullWidth);
		int i, nCoeff = m_nFullHeight * m_nFullWidth;
		for (i = 0; i < nCoeff; i++)
			m_pData[i] += pbsgw->m_pData[i];
	}

	void sub(CBSplineGridWavelet<T>* pbsgw)
	{
		_ASSERT(m_pData != NULL);
		_ASSERT(pbsgw != NULL);
		_ASSERT(m_nFullHeight == pbsgw->m_nFullHeight);
		_ASSERT(m_nFullWidth == pbsgw->m_nFullWidth);
		int i, nCoeff = m_nFullHeight * m_nFullWidth;
		for (i = 0; i < nCoeff; i++)
			m_pData[i] -= pbsgw->m_pData[i];
	}

	void mul(CBSplineGridWavelet<T>* pbsgw)
	{
		_ASSERT(m_pData != NULL);
		_ASSERT(pbsgw != NULL);
		_ASSERT(m_nFullHeight == pbsgw->m_nFullHeight);
		_ASSERT(m_nFullWidth == pbsgw->m_nFullWidth);
		int i, nCoeff = m_nFullHeight * m_nFullWidth;
		for (i = 0; i < nCoeff; i++)
			m_pData[i] *= pbsgw->m_pData[i];
	}

	void sqrt()
	{
		_ASSERT(m_pData != NULL);
		int i, nCoeff = m_nFullHeight * m_nFullWidth;
		for (i = 0; i < nCoeff; i++)
			m_pData[i] = (T)::sqrt(m_pData[i]);
	}

	virtual void decompose(int level);
	virtual void reconstruct(int level);

	virtual void reconstructCopy(int level, T* pDest, int nPitch, int nStep);
	void reconstructCopyLocal(int iCoeff, int level, T* pDest, int nPitch, int nStep);

	void save(char* szFilename);

	void getRawCopy(T* pDest)
	{
		memcpy(pDest, m_pData, getFullResWidth() * getFullResHeight() * sizeof(T));
	}


	///////////////////////////////////////////////////////////////////////////
	//static lifting methods
	///////////////////////////////////////////////////////////////////////////

	static inline void s_lift_v(const double a, const double b, T* pCoeff, const int i, const int j, const int nPitch, const int nStep)
	{
		int iv, iv1;
		const double a2 = a*a, b2 = b*b, ab = a * b;
		T faceCentroid;
		T edgeCentroid;

		int nPitch1 = nPitch * nStep;

		//update v -> v'
		const int i1 = 2 * i;
		const int j1 = 2 * j;
		iv = iv1 = i1 * nPitch1 + j1 * nStep;
		pCoeff[iv1] *= b2;

		iv1 -= nPitch1;
		edgeCentroid = pCoeff[iv1];
		iv1 += nStep;
		faceCentroid = pCoeff[iv1];
		iv1 += nPitch1;
		edgeCentroid += pCoeff[iv1];
		iv1 += nPitch1;
		faceCentroid += pCoeff[iv1];
		iv1 -= nStep;
		edgeCentroid += pCoeff[iv1];
		iv1 -= nStep;
		faceCentroid += pCoeff[iv1];
		iv1 -= nPitch1;
		edgeCentroid += pCoeff[iv1];
		iv1 -= nPitch1;
		faceCentroid += pCoeff[iv1];

		edgeCentroid *= 0.25;
		faceCentroid *= 0.25;

//		pCoeff[iv] += 4.0 * a2 * faceCentroid + 4.0 * ab * edgeCentroid;
		pCoeff[iv] += faceCentroid * 4.0 * a2 + edgeCentroid * 4.0 * ab;
	}

	static inline void s_lift_e(const double a, const double b, T* pCoeff, const int i, const int j, const int nPitch, const int nStep)
	{
		int iv1, ie1;
		const double a2 = a*a, b2 = b*b, ab = a * b;
		T faceEdgeCentroid = 0.0;

		int nPitch1 = nPitch * nStep;

		//update e -> e'
		int i1 = 2 * i;
		int j1 = 2 * j + 1;
		ie1 = i1 * nPitch1 + j1 * nStep;
		pCoeff[ie1] *= b;
		iv1 = ie1 - nPitch1;
		faceEdgeCentroid = pCoeff[iv1];
		iv1 = ie1 + nPitch1;
		faceEdgeCentroid += pCoeff[iv1];

		faceEdgeCentroid *= 0.5;
//		pCoeff[ie1] += 2.0 * a * faceEdgeCentroid;
		pCoeff[ie1] += faceEdgeCentroid * 2.0 * a;

		//update e -> e'
		i1 = 2 * i + 1;
		j1 = 2 * j;
		ie1 = i1 * nPitch1 + j1 * nStep;
		pCoeff[ie1] *= b;
		iv1 = ie1 - nStep;
		faceEdgeCentroid = pCoeff[iv1];
		iv1 = ie1 + nStep;
		faceEdgeCentroid += pCoeff[iv1];

		faceEdgeCentroid *= 0.5;
//		pCoeff[ie1] += 2.0 * a * faceEdgeCentroid;
		pCoeff[ie1] += faceEdgeCentroid * 2.0 * a;
	}

	static inline void s_lift_corner_v(const double a, const double b, T* pCoeff, const int i, const int j, const int nPitch, const int nStep)
	{
		int iv, iv1;
		const double a2 = a*a, b2 = b*b, ab = a * b;
		T faceCentroid = 0.0;
		T edgeCentroid = 0.0;

		int nPitch1 = nPitch * nStep;

		//update v -> v'
		const int i1 = 2 * i;
		const int j1 = 2 * j;
		iv = iv1 = i1 * nPitch1 + j1 * nStep;
		pCoeff[iv1] *= b2;

		if (i != 0)
		{
			iv1 = iv - nPitch1;
			edgeCentroid += pCoeff[iv1];
			if (j == 0)
			{
				iv1 += nStep;
				faceCentroid += pCoeff[iv1];
			}
			else
			{
				iv1 -= nStep;
				faceCentroid += pCoeff[iv1];
			}
			iv1 += nPitch1;
			edgeCentroid += pCoeff[iv1];
		}
		else
		{
			iv1 = iv + nPitch1;
			edgeCentroid += pCoeff[iv1];
			if (j == 0)
			{
				iv1 += nStep;
				faceCentroid += pCoeff[iv1];
			}
			else
			{
				iv1 -= nStep;
				faceCentroid += pCoeff[iv1];
			}
			iv1 -= nPitch1;
			edgeCentroid += pCoeff[iv1];
		}

		//2 edge neighbors, 1 face neighbor
		edgeCentroid *= 0.5;

//		pCoeff[iv] += 4.0 * a2 * faceCentroid + 4.0 * ab * edgeCentroid;
		pCoeff[iv] += faceCentroid * 4.0 * a2 + edgeCentroid * 4.0 * ab;
	}

	static inline void s_lift_corner_e(const double a, const double b, T* pCoeff, const int i, const int j, const int nPitch, const int nStep)
	{
		int iv1, ie1;
		const double a2 = a*a, b2 = b*b, ab = a * b;
		T faceEdgeCentroid = 0.0;

		int nPitch1 = nPitch * nStep;

		//update e -> e'
		if (j == 0)
		{
			const int i1 = 2 * i;
			const int j1 = 2 * j + 1;
			ie1 = i1 * nPitch1 + j1 * nStep;
			pCoeff[ie1] *= b;
			if (i == 0)
			{
				iv1 = ie1 + nPitch1;
				faceEdgeCentroid = pCoeff[iv1];
			}
			else
			{
				iv1 = ie1 - nPitch1;
				faceEdgeCentroid = pCoeff[iv1];
			}

//			pCoeff[ie1] += 2.0 * a * faceEdgeCentroid;
			pCoeff[ie1] += faceEdgeCentroid * 2.0 * a;
		}

		//update e -> e'
		if (i == 0)
		{
			const int i1 = 2 * i + 1;
			const int j1 = 2 * j;
			ie1 = i1 * nPitch1 + j1 * nStep;
			pCoeff[ie1] *= b;
			if (j == 0)
			{
				iv1 = ie1 + nStep;
				faceEdgeCentroid = pCoeff[iv1];
			}
			else
			{
				iv1 = ie1 - nStep;
				faceEdgeCentroid = pCoeff[iv1];
			}

//			pCoeff[ie1] += 2.0 * a * faceEdgeCentroid;
			pCoeff[ie1] += faceEdgeCentroid * 2.0 * a;
		}
	}

	static inline void s_lift_vbound_v(const double a, const double b, T* pCoeff, const int i, const int j, const int nPitch, const int nStep)
	{
		int iv, iv1;
		const double a2 = a*a, b2 = b*b, ab = a * b;
		T faceCentroid = 0.0;
		T edgeCentroid = 0.0;

		int nPitch1 = nPitch * nStep;

		//update v -> v'
		const int i1 = 2 * i;
		const int j1 = 2 * j;
		iv = iv1 = i1 * nPitch1 + j1 * nStep;
		pCoeff[iv1] *= b2;

		if (i == 0)
		{
			iv1 += nStep;
			edgeCentroid = pCoeff[iv1];
			iv1 += nPitch1;
			faceCentroid += pCoeff[iv1];
			iv1 -= nStep;
			edgeCentroid += pCoeff[iv1];
			iv1 -= nStep;
			faceCentroid += pCoeff[iv1];
			iv1 -= nPitch1;
			edgeCentroid += pCoeff[iv1];
		}
		else
		{
			iv1 -= nPitch1;
			edgeCentroid = pCoeff[iv1];
			iv1 += nStep;
			faceCentroid = pCoeff[iv1];
			iv1 += nPitch1;
			edgeCentroid += pCoeff[iv1];
			iv1 -= nStep;
			iv1 -= nStep;
			edgeCentroid += pCoeff[iv1];
			iv1 -= nPitch1;
			faceCentroid += pCoeff[iv1];
		}

//		edgeCentroid *= 0.25;
		edgeCentroid /= 3.0;
//		faceCentroid *= 0.25;
		faceCentroid *= 0.5;

//		pCoeff[iv] += 4.0 * a2 * faceCentroid + 4.0 * ab * edgeCentroid;
		pCoeff[iv] += faceCentroid * 4.0 * a2 + edgeCentroid * 4.0 * ab;
	}

	static inline void s_lift_vbound_e(const double a, const double b, T* pCoeff, const int i, const int j, const int nPitch, const int nStep)
	{
		int i1, j1;
		int iv1, ie1;
		const double a2 = a*a, b2 = b*b, ab = a * b;
		T faceEdgeCentroid = 0.0;

		bool bzero = (i == 0);

		const int nPitch1 = nPitch * nStep;

		//update e -> e'
		i1 = 2 * i;
		j1 = 2 * j + 1;
		ie1 = i1 * nPitch1 + j1 * nStep;
		pCoeff[ie1] *= b;
		if (bzero)
		{
			iv1 = ie1 + nPitch1;
			faceEdgeCentroid += pCoeff[iv1];
		}
		else
		{
			iv1 = ie1 - nPitch1;
			faceEdgeCentroid += pCoeff[iv1];
		}

//		faceEdgeCentroid *= 0.5;
//		pCoeff[ie1] += 2.0 * a * faceEdgeCentroid;
		pCoeff[ie1] += faceEdgeCentroid * 2.0 * a;

		if (bzero)
		{
			//update e -> e'
			i1 = 2 * i + 1;
			j1 = 2 * j;
			ie1 = i1 * nPitch1 + j1 * nStep;
			pCoeff[ie1] *= b;
			iv1 = ie1 - nStep;
			faceEdgeCentroid = pCoeff[iv1];
			iv1 = ie1 + nStep;
			faceEdgeCentroid += pCoeff[iv1];

			faceEdgeCentroid *= 0.5;
//			pCoeff[ie1] += 2.0 * a * faceEdgeCentroid;
			pCoeff[ie1] += faceEdgeCentroid * 2.0 * a;
		}
	}

	static inline void s_lift_hbound_v(const double a, const double b, T* pCoeff, const int i, const int j, const int nPitch, const int nStep)
	{
		int iv, iv1;
		const double a2 = a*a, b2 = b*b, ab = a * b;
		T faceCentroid = 0.0;
		T edgeCentroid = 0.0;

		const int nPitch1 = nPitch * nStep;

		//update v -> v'
		const int i1 = 2 * i;
		const int j1 = 2 * j;
		iv = iv1 = i1 * nPitch1 + j1 * nStep;
		pCoeff[iv1] *= b2;

		if (j == 0)
		{
			iv1 -= nPitch1;
			edgeCentroid = pCoeff[iv1];
			iv1 += nStep;
			faceCentroid = pCoeff[iv1];
			iv1 += nPitch1;
			edgeCentroid += pCoeff[iv1];
			iv1 += nPitch1;
			faceCentroid += pCoeff[iv1];
			iv1 -= nStep;
			edgeCentroid += pCoeff[iv1];
		}
		else
		{
			iv1 -= nPitch1;
			edgeCentroid = pCoeff[iv1];
			iv1 += nPitch1;
			iv1 += nPitch1;
			edgeCentroid += pCoeff[iv1];
			iv1 -= nStep;
			faceCentroid += pCoeff[iv1];
			iv1 -= nPitch1;
			edgeCentroid += pCoeff[iv1];
			iv1 -= nPitch1;
			faceCentroid += pCoeff[iv1];
		}

//		edgeCentroid *= 0.25;
//		faceCentroid *= 0.25;
		edgeCentroid /= 3.0;
		faceCentroid *= 0.5;

//		pCoeff[iv] += 4.0 * a2 * faceCentroid + 4.0 * ab * edgeCentroid;
		pCoeff[iv] += faceCentroid * 4.0 * a2 + edgeCentroid * 4.0 * ab;
	}

	static inline void s_lift_hbound_e(const double a, const double b, T* pCoeff, const int i, const int j, const int nPitch, const int nStep)
	{
		int i1, j1;
		int iv1, ie1;
		const double a2 = a*a, b2 = b*b, ab = a * b;
		T faceEdgeCentroid = 0.0;

		bool bzero = (j == 0);

		const int nPitch1 = nPitch * nStep;

		//update e -> e'
		if (bzero)
		{
			i1 = 2 * i;
			j1 = 2 * j + 1;
			ie1 = i1 * nPitch1 + j1 * nStep;
			pCoeff[ie1] *= b;
			iv1 = ie1 - nPitch1;
			faceEdgeCentroid = pCoeff[iv1];
			iv1 = ie1 + nPitch1;
			faceEdgeCentroid += pCoeff[iv1];

			faceEdgeCentroid *= 0.5;
//			pCoeff[ie1] += 2.0 * a * faceEdgeCentroid;
			pCoeff[ie1] += faceEdgeCentroid * 2.0 * a;
		}

		//update e -> e'
		i1 = 2 * i + 1;
		j1 = 2 * j;
		ie1 = i1 * nPitch1 + j1 * nStep;
		pCoeff[ie1] *= b;
		if (bzero)
		{
			iv1 = ie1 + nStep;
			faceEdgeCentroid = pCoeff[iv1];
		}
		else
		{
			iv1 = ie1 - nStep;
			faceEdgeCentroid = pCoeff[iv1];
		}

//		faceEdgeCentroid *= 0.5;
//		pCoeff[ie1] += 2.0 * a * faceEdgeCentroid;
		pCoeff[ie1] += faceEdgeCentroid * 2.0 * a;
	}

	static inline void w_lift_f(const double a, const double b, T* pCoeff, const int i, const int j, const int nPitch, const int nStep)
	{
		int iv1, iface1;
		const double a2 = a*a, b2 = b*b, ab = a * b;
		T faceCentroid = 0.0;
		T faceEdgeCentroid = 0.0;

		const int nPitch1 = nPitch * nStep;

		//update f -> f'
		const int i1 = 2 * i + 1;
		const int j1 = 2 * j + 1;
		iface1 = i1 * nPitch1 + j1 * nStep;
		pCoeff[iface1] *= b2;
		iv1 = iface1 - nPitch1 - nStep;
		faceCentroid = pCoeff[iv1];
		iv1 += nStep;
		faceEdgeCentroid = pCoeff[iv1];
		iv1 += nStep;
		faceCentroid += pCoeff[iv1];
		iv1 += nPitch1;
		faceEdgeCentroid += pCoeff[iv1];
		iv1 += nPitch1;
		faceCentroid += pCoeff[iv1];
		iv1 -= nStep;
		faceEdgeCentroid += pCoeff[iv1];
		iv1 -= nStep;
		faceCentroid += pCoeff[iv1];
		iv1 -= nPitch1;
		faceEdgeCentroid += pCoeff[iv1];

		faceCentroid *= 0.25;
		faceEdgeCentroid *= 0.25;
//		pCoeff[iface1] += 4.0 * a2 * faceCentroid + 4.0 * ab * faceEdgeCentroid;
		pCoeff[iface1] += faceCentroid * 4.0 * a2 + faceEdgeCentroid * 4.0 * ab;
	}

	static inline void w_lift_e(const double a, const double b, T* pCoeff, const int i, const int j, const int nPitch, const int nStep)
	{
		int i1, j1;
		int iv1, ie1;
		const double a2 = a*a, b2 = b*b, ab = a * b;
		T edgeCentroid = 0.0;

		const int nPitch1 = nPitch * nStep;

		//update e -> e'
		i1 = 2 * i;
		j1 = 2 * j + 1;
		ie1 = i1 * nPitch1 + j1 * nStep;
		pCoeff[ie1] *= b;
		iv1 = ie1 - nStep;
		edgeCentroid = pCoeff[iv1];
		iv1 = ie1 + nStep;
		edgeCentroid += pCoeff[iv1];

		edgeCentroid *= 0.5;
//		pCoeff[ie1] += 2.0 * a * edgeCentroid;
		pCoeff[ie1] += edgeCentroid * 2.0 * a;

		//update e -> e'
		i1 = 2 * i + 1;
		j1 = 2 * j;
		ie1 = i1 * nPitch1 + j1 * nStep;
		pCoeff[ie1] *= b;
		iv1 = ie1 - nPitch1;
		edgeCentroid = pCoeff[iv1];
		iv1 = ie1 + nPitch1;
		edgeCentroid += pCoeff[iv1];

		edgeCentroid *= 0.5;
//		pCoeff[ie1] += 2.0 * a * edgeCentroid;
		pCoeff[ie1] += edgeCentroid * 2.0 * a;
	}

	static inline void w_lift_corner_e(const double a, const double b, T* pCoeff, const int i, const int j, const int nPitch, const int nStep)
	{
		int i1, j1;
		int iv1, ie1;
		const double a2 = a*a, b2 = b*b, ab = a * b;
		T edgeCentroid = 0.0;

		const int nPitch1 = nPitch * nStep;

		if (j == 0)
		{
			//update e -> e'
			i1 = 2 * i;
			j1 = 2 * j + 1;
			ie1 = i1 * nPitch1 + j1 * nStep;
			pCoeff[ie1] *= b;
			iv1 = ie1 - nStep;
			edgeCentroid = pCoeff[iv1];
			iv1 = ie1 + nStep;
			edgeCentroid += pCoeff[iv1];

			edgeCentroid *= 0.5;
//			pCoeff[ie1] += 2.0 * a * edgeCentroid;
			pCoeff[ie1] += edgeCentroid * 2.0 * a;
		}

		if (i == 0)
		{
			//update e -> e'
			i1 = 2 * i + 1;
			j1 = 2 * j;
			ie1 = i1 * nPitch1 + j1 * nStep;
			pCoeff[ie1] *= b;
			iv1 = ie1 - nPitch1;
			edgeCentroid = pCoeff[iv1];
			iv1 = ie1 + nPitch1;
			edgeCentroid += pCoeff[iv1];

			edgeCentroid *= 0.5;
//			pCoeff[ie1] += 2.0 * a * edgeCentroid;
			pCoeff[ie1] += edgeCentroid * 2.0 * a;
		}
	}

	static inline void w_lift_vbound_e(const double a, const double b, T* pCoeff, const int i, const int j, const int nPitch, const int nStep)
	{
		int i1, j1;
		int iv1, ie1;
		const double a2 = a*a, b2 = b*b, ab = a * b;
		T edgeCentroid = 0.0;

		const int nPitch1 = nPitch * nStep;

		//update e -> e'
		i1 = 2 * i;
		j1 = 2 * j + 1;
		ie1 = i1 * nPitch1 + j1 * nStep;
		pCoeff[ie1] *= b;
		iv1 = ie1 - nStep;
		edgeCentroid = pCoeff[iv1];
		iv1 = ie1 + nStep;
		edgeCentroid += pCoeff[iv1];

		edgeCentroid *= 0.5;
//		pCoeff[ie1] += 2.0 * a * edgeCentroid;
		pCoeff[ie1] += edgeCentroid * 2.0 * a;

		if (i == 0)
		{
			//update e -> e'
			i1 = 2 * i + 1;
			j1 = 2 * j;
			ie1 = i1 * nPitch1 + j1 * nStep;
			pCoeff[ie1] *= b;
			iv1 = ie1 - nPitch1;
			edgeCentroid = pCoeff[iv1];
			iv1 = ie1 + nPitch1;
			edgeCentroid += pCoeff[iv1];

			edgeCentroid *= 0.5;
//			pCoeff[ie1] += 2.0 * a * edgeCentroid;
			pCoeff[ie1] += edgeCentroid * 2.0 * a;
		}
	}

	static inline void w_lift_hbound_e(const double a, const double b, T* pCoeff, const int i, const int j, const int nPitch, const int nStep)
	{
		int i1, j1;
		int iv1, ie1;
		const double a2 = a*a, b2 = b*b, ab = a * b;
		T edgeCentroid = 0.0;

		const int nPitch1 = nPitch * nStep;

		if (j == 0)
		{
			//update e -> e'
			i1 = 2 * i;
			j1 = 2 * j + 1;
			ie1 = i1 * nPitch1 + j1 * nStep;
			pCoeff[ie1] *= b;
			iv1 = ie1 - nStep;
			edgeCentroid = pCoeff[iv1];
			iv1 = ie1 + nStep;
			edgeCentroid += pCoeff[iv1];

			edgeCentroid *= 0.5;
//			pCoeff[ie1] += 2.0 * a * edgeCentroid;
			pCoeff[ie1] += edgeCentroid * 2.0 * a;
		}

		//update e -> e'
		i1 = 2 * i + 1;
		j1 = 2 * j;
		ie1 = i1 * nPitch1 + j1 * nStep;
		pCoeff[ie1] *= b;
		iv1 = ie1 - nPitch1;
		edgeCentroid = pCoeff[iv1];
		iv1 = ie1 + nPitch1;
		edgeCentroid += pCoeff[iv1];

		edgeCentroid *= 0.5;
//		pCoeff[ie1] += 2.0 * a * edgeCentroid;
		pCoeff[ie1] += edgeCentroid * 2.0 * a;
	}


	///////////////////////////////////////////////////////////////////////////
	//static full passes over grid
	///////////////////////////////////////////////////////////////////////////

	static void s_lift_vertex_pass(double a, double b, T* pCoeff, int nWidth, int nHeight, int nPitch1, int nStep)
	{
		int i, j;

		int nLastRow = nHeight - 1, nRowStop = nWidth - 1;
		int nWidth1 = nWidth * 2 - 1;
		int nHeight1 = nHeight * 2 - 1;

		i = 0;
		j = 0;

		s_lift_corner_v(a, b, pCoeff, i, j, nPitch1, nStep);

		for (j = 1; j < nRowStop; j++)
			s_lift_vbound_v(a, b, pCoeff, i, j, nPitch1, nStep);

		s_lift_corner_v(a, b, pCoeff, i, j, nPitch1, nStep);

		for (i = 1; i < nLastRow; i++)
		{
			j = 0;

			s_lift_hbound_v(a, b, pCoeff, i, j, nPitch1, nStep);

			for (j = 1; j < nRowStop; j++)
				s_lift_v(a, b, pCoeff, i, j, nPitch1, nStep);

			s_lift_hbound_v(a, b, pCoeff, i, j, nPitch1, nStep);
		}

		j = 0;
		s_lift_corner_v(a, b, pCoeff, i, j, nPitch1, nStep);

		for (j = 1; j < nRowStop; j++)
			s_lift_vbound_v(a, b, pCoeff, i, j, nPitch1, nStep);

		s_lift_corner_v(a, b, pCoeff, i, j, nPitch1, nStep);
	}

	static void s_lift_edge_pass(double a, double b, T* pCoeff, int nWidth, int nHeight, int nPitch1, int nStep)
	{
		int i, j;

		int nWidth1 = nWidth * 2 - 1;
		int nHeight1 = nHeight * 2 - 1;
		int nLastRow = nHeight - 1, nRowStop = nWidth - 1;

		i = 0;
		j = 0;
		s_lift_corner_e(a, b, pCoeff, i, j, nPitch1, nStep);

		for (j = 1; j < nRowStop; j++)
			s_lift_vbound_e(a, b, pCoeff, i, j, nPitch1, nStep);

		s_lift_corner_e(a, b, pCoeff, i, j, nPitch1, nStep);

		for (i = 1; i < nLastRow; i++)
		{
			j = 0;
			s_lift_hbound_e(a, b, pCoeff, i, j, nPitch1, nStep);

			for (j = 1; j < nRowStop; j++)
				s_lift_e(a, b, pCoeff, i, j, nPitch1, nStep);

			s_lift_hbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
		}

		j = 0;
		s_lift_corner_e(a, b, pCoeff, i, j, nPitch1, nStep);

		for (j = 1; j < nRowStop; j++)
			s_lift_vbound_e(a, b, pCoeff, i, j, nPitch1, nStep);

		s_lift_corner_e(a, b, pCoeff, i, j, nPitch1, nStep);
	}

	static void w_lift_face_pass(double a, double b, T* pCoeff, int nWidth, int nHeight, int nPitch1, int nStep)
	{
		int i, j;

		int nWidth1 = 2 * nWidth - 1;
		int nHeight1 = 2 * nHeight - 1;
		int nLastRow = nHeight - 1, nRowStop = nWidth - 1;

		for (i = 0; i < nLastRow; i++)
		{
			for (j = 0; j < nRowStop; j++)
				w_lift_f(a, b, pCoeff, i, j, nPitch1, nStep);
		}
	}

	static void w_lift_edge_pass(double a, double b, T* pCoeff, int nWidth, int nHeight, int nPitch1, int nStep)
	{
		int i, j;

		int nWidth1 = 2 * nWidth - 1;
		int nHeight1 = 2 * nHeight - 1;
		int nLastRow = nHeight - 1, nRowStop = nWidth - 1;

		i = 0;
		j = 0;
		w_lift_corner_e(a, b, pCoeff, i, j, nPitch1, nStep);

		for (j = 1; j < nRowStop; j++)
			w_lift_vbound_e(a, b, pCoeff, i, j, nPitch1, nStep);

		w_lift_corner_e(a, b, pCoeff, i, j, nPitch1, nStep);

		for (i = 1; i < nLastRow; i++)
		{
			j = 0;
			w_lift_hbound_e(a, b, pCoeff, i, j, nPitch1, nStep);

			for (j = 1; j < nRowStop; j++)
				w_lift_e(a, b, pCoeff, i, j, nPitch1, nStep);

			w_lift_hbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
		}

		j = 0;
		w_lift_corner_e(a, b, pCoeff, i, j, nPitch1, nStep);
		
		for (j = 1; j < nRowStop; j++)
			w_lift_vbound_e(a, b, pCoeff, i, j, nPitch1, nStep);

		w_lift_corner_e(a, b, pCoeff, i, j, nPitch1, nStep);
	}

	//passes over sub-grids
	static void s_lift_vertex_subpass(double a, double b, T* pCoeff, int nWidth, int nHeight, int nPitch1, int nStep, int x0, int y0, int x1, int y1)
	{
		int i, j, k, i0, j0, i1, j1, iv0, iv1, ie1, iface1;

		int nLastRow = nHeight - 1, nRowStop = nWidth - 1;
		int nWidth1 = nWidth * 2 - 1;
		int nHeight1 = nHeight * 2 - 1;

		x0 = max(x0, 0);
		y0 = max(y0, 0);
		x1 = min(x1, nRowStop);
		y1 = min(y1, nLastRow);

		i = 0;
		j = 0;

		if (y0 == 0)
		{
			if (x0 == 0)
			{
				s_lift_corner_v(a, b, pCoeff, i, j, nPitch1, nStep);
				if (x1 == nRowStop)
				{
					for (j = 1; j < nRowStop; j++)
						s_lift_vbound_v(a, b, pCoeff, i, j, nPitch1, nStep);
					s_lift_corner_v(a, b, pCoeff, i, j, nPitch1, nStep);
					if (y1 == nLastRow)
					{
						for (i = 1; i < nLastRow; i++)
						{
							j = 0;
							s_lift_hbound_v(a, b, pCoeff, i, j, nPitch1, nStep);
							for (j = 1; j < nRowStop; j++)
								s_lift_v(a, b, pCoeff, i, j, nPitch1, nStep);
							s_lift_hbound_v(a, b, pCoeff, i, j, nPitch1, nStep);
						}
						j = 0;
						s_lift_corner_v(a, b, pCoeff, i, j, nPitch1, nStep);
						for (j = 1; j < nRowStop; j++)
							s_lift_vbound_v(a, b, pCoeff, i, j, nPitch1, nStep);
					}
					else
					{
						for (i = 1; i <= y1; i++)
						{
							j = 0;
							s_lift_hbound_v(a, b, pCoeff, i, j, nPitch1, nStep);
							for (j = 1; j < nRowStop; j++)
								s_lift_v(a, b, pCoeff, i, j, nPitch1, nStep);
							s_lift_hbound_v(a, b, pCoeff, i, j, nPitch1, nStep);
						}
					}
				}
				else
				{
					for (j = 1; j <= x1; j++)
						s_lift_vbound_v(a, b, pCoeff, i, j, nPitch1, nStep);
					if (y1 == nLastRow)
					{
						for (i = 1; i < nLastRow; i++)
						{
							j = 0;
							s_lift_hbound_v(a, b, pCoeff, i, j, nPitch1, nStep);
							for (j = 1; j <= x1; j++)
								s_lift_v(a, b, pCoeff, i, j, nPitch1, nStep);
						}
						j = 0;
						s_lift_corner_v(a, b, pCoeff, i, j, nPitch1, nStep);
						for (j = 1; j <= x1; j++)
							s_lift_v(a, b, pCoeff, i, j, nPitch1, nStep);
					}
					else
					{
						for (i = 1; i <= y1; i++)
						{
							j = 0;
							s_lift_hbound_v(a, b, pCoeff, i, j, nPitch1, nStep);
							for (j = 1; j <= x1; j++)
								s_lift_v(a, b, pCoeff, i, j, nPitch1, nStep);
						}
					}
				}
			}
			else
			{
				if (x1 == nRowStop)
				{
					for (j = x0; j < nRowStop; j++)
						s_lift_vbound_v(a, b, pCoeff, i, j, nPitch1, nStep);
					s_lift_corner_v(a, b, pCoeff, i, j, nPitch1, nStep);
					if (y1 == nLastRow)
					{
						for (i = 1; i < nLastRow; i++)
						{
							for (j = x0; j < nRowStop; j++)
								s_lift_v(a, b, pCoeff, i, j, nPitch1, nStep);
							s_lift_hbound_v(a, b, pCoeff, i, j, nPitch1, nStep);
						}
						for (j = x0; j < nRowStop; j++)
							s_lift_vbound_v(a, b, pCoeff, i, j, nPitch1, nStep);
						s_lift_corner_v(a, b, pCoeff, i, j, nPitch1, nStep);
					}
					else
					{
						for (i = 1; i <= y1; i++)
						{
							for (j = x0; j < nRowStop; j++)
								s_lift_v(a, b, pCoeff, i, j, nPitch1, nStep);
							s_lift_hbound_v(a, b, pCoeff, i, j, nPitch1, nStep);
						}
					}
				}
				else
				{
					for (j = x0; j <= x1; j++)
						s_lift_vbound_v(a, b, pCoeff, i, j, nPitch1, nStep);
					if (y1 == nLastRow)
					{
						for (i = 1; i < nLastRow; i++)
						{
							for (j = x0; j <= x1; j++)
								s_lift_v(a, b, pCoeff, i, j, nPitch1, nStep);
						}
						for (j = x0; j <= x1; j++)
							s_lift_vbound_v(a, b, pCoeff, i, j, nPitch1, nStep);
					}
					else
					{
						for (i = 1; i <= y1; i++)
						{
							for (j = x0; j <= x1; j++)
								s_lift_v(a, b, pCoeff, i, j, nPitch1, nStep);
						}
					}
				}
			}
		}
		else
		{
			if (x0 == 0)
			{
				if (x1 == nRowStop)
				{
					if (y1 == nLastRow)
					{
						for (i = y0; i < nLastRow; i++)
						{
							j = 0;
							s_lift_hbound_v(a, b, pCoeff, i, j, nPitch1, nStep);
							for (j = 1; j < nRowStop; j++)
								s_lift_v(a, b, pCoeff, i, j, nPitch1, nStep);
							s_lift_hbound_v(a, b, pCoeff, i, j, nPitch1, nStep);
						}
						j = 0;
						s_lift_corner_v(a, b, pCoeff, i, j, nPitch1, nStep);
						for (j = 1; j < nRowStop; j++)
							s_lift_vbound_v(a, b, pCoeff, i, j, nPitch1, nStep);
						s_lift_corner_v(a, b, pCoeff, i, j, nPitch1, nStep);
					}
					else
					{
						for (i = y0; i <= y1; i++)
						{
							j = 0;
							s_lift_hbound_v(a, b, pCoeff, i, j, nPitch1, nStep);
							for (j = 1; j < nRowStop; j++)
								s_lift_v(a, b, pCoeff, i, j, nPitch1, nStep);
							s_lift_hbound_v(a, b, pCoeff, i, j, nPitch1, nStep);
						}
					}
				}
				else
				{
					if (y1 == nLastRow)
					{
						for (i = y0; i < nLastRow; i++)
						{
							j = 0;
							s_lift_hbound_v(a, b, pCoeff, i, j, nPitch1, nStep);
							for (j = 1; j <= x1; j++)
								s_lift_v(a, b, pCoeff, i, j, nPitch1, nStep);
						}
						j = 0;
						s_lift_corner_v(a, b, pCoeff, i, j, nPitch1, nStep);
						for (j = 1; j <= x1; j++)
							s_lift_vbound_v(a, b, pCoeff, i, j, nPitch1, nStep);
					}
					else
					{
						for (i = y0; i <= y1; i++)
						{
							j = 0;
							s_lift_hbound_v(a, b, pCoeff, i, j, nPitch1, nStep);
							for (j = 1; j <= x1; j++)
								s_lift_v(a, b, pCoeff, i, j, nPitch1, nStep);
						}
					}
				}
			}
			else
			{
				if (x1 == nRowStop)
				{
					if (y1 == nLastRow)
					{
						for (i = y0; i < nLastRow; i++)
						{
							for (j = x0; j < nRowStop; j++)
								s_lift_v(a, b, pCoeff, i, j, nPitch1, nStep);
							s_lift_hbound_v(a, b, pCoeff, i, j, nPitch1, nStep);
						}
						for (j = x0; j < nRowStop; j++)
							s_lift_vbound_v(a, b, pCoeff, i, j, nPitch1, nStep);
						s_lift_corner_v(a, b, pCoeff, i, j, nPitch1, nStep);
					}
					else
					{
						for (i = y0; i <= y1; i++)
						{
							for (j = x0; j < nRowStop; j++)
								s_lift_v(a, b, pCoeff, i, j, nPitch1, nStep);
							s_lift_hbound_v(a, b, pCoeff, i, j, nPitch1, nStep);
						}
					}
				}
				else
				{
					if (y1 == nLastRow)
					{
						for (i = y0; i < nLastRow; i++)
						{
							for (j = x0; j <= x1; j++)
								s_lift_v(a, b, pCoeff, i, j, nPitch1, nStep);
						}
						for (j = x0; j <= x1; j++)
							s_lift_vbound_v(a, b, pCoeff, i, j, nPitch1, nStep);
					}
					else
					{
						for (i = y0; i <= y1; i++)
						{
							for (j = x0; j <= x1; j++)
								s_lift_v(a, b, pCoeff, i, j, nPitch1, nStep);
						}
					}
				}
			}
		}
	}

	static void s_lift_edge_subpass(double a, double b, T* pCoeff, int nWidth, int nHeight, int nPitch1, int nStep, int x0, int y0, int x1, int y1)
	{
		int i, j, k, i0, j0, i1, j1, iv0, iv1, ie1, iface1;

		int nWidth1 = nWidth * 2 - 1;
		int nHeight1 = nHeight * 2 - 1;
		int nLastRow = nHeight - 1, nRowStop = nWidth - 1;

		x0 = max(x0, 0);
		y0 = max(y0, 0);
		x1 = min(x1, nRowStop);
		y1 = min(y1, nLastRow);

		i = 0;
		j = 0;
		if (y0 == 0)
		{
			if (x0 == 0)
			{
				s_lift_corner_e(a, b, pCoeff, i, j, nPitch1, nStep);
				if (x1 == nRowStop)
				{
					for (j = 1; j < nRowStop; j++)
						s_lift_vbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
					s_lift_corner_e(a, b, pCoeff, i, j, nPitch1, nStep);
					if (y1 == nLastRow)
					{
						for (i = 1; i < nLastRow; i++)
						{
							j = 0;
							s_lift_hbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
							for (j = 1; j < nRowStop; j++)
								s_lift_e(a, b, pCoeff, i, j, nPitch1, nStep);
							s_lift_hbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
						}
						j = 0;
						s_lift_corner_e(a, b, pCoeff, i, j, nPitch1, nStep);
						for (j = 1; j < nRowStop; j++)
							s_lift_vbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
						s_lift_corner_e(a, b, pCoeff, i, j, nPitch1, nStep);
					}
					else
					{
						for (i = 1; i <= y1; i++)
						{
							j = 0;
							s_lift_hbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
							for (j = 1; j < nRowStop; j++)
								s_lift_e(a, b, pCoeff, i, j, nPitch1, nStep);
							s_lift_hbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
						}
					}
				}
				else
				{
					for (j = 1; j <= x1; j++)
						s_lift_vbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
					if (y1 == nLastRow)
					{
						for (i = 1; i < nLastRow; i++)
						{
							j = 0;
							s_lift_hbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
							for (j = 1; j <= x1; j++)
								s_lift_e(a, b, pCoeff, i, j, nPitch1, nStep);
						}
						j = 0;
						s_lift_corner_e(a, b, pCoeff, i, j, nPitch1, nStep);
						for (j = 1; j <= x1; j++)
							s_lift_vbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
					}
					else
					{
						for (i = 1; i <= y1; i++)
						{
							j = 0;
							s_lift_hbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
							for (j = 1; j <= x1; j++)
								s_lift_e(a, b, pCoeff, i, j, nPitch1, nStep);
						}
					}
				}
			}
			else
			{
				if (x1 == nRowStop)
				{
					for (j = x0; j < nRowStop; j++)
						s_lift_vbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
					s_lift_corner_e(a, b, pCoeff, i, j, nPitch1, nStep);
					if (y1 == nLastRow)
					{
						for (i = 1; i < nLastRow; i++)
						{
							for (j = x0; j < nRowStop; j++)
								s_lift_e(a, b, pCoeff, i, j, nPitch1, nStep);
							s_lift_hbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
						}
						for (j = x0; j < nRowStop; j++)
							s_lift_vbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
						s_lift_corner_e(a, b, pCoeff, i, j, nPitch1, nStep);
					}
					else
					{
						for (i = 1; i <= y1; i++)
						{
							for (j = x0; j < nRowStop; j++)
								s_lift_e(a, b, pCoeff, i, j, nPitch1, nStep);
							s_lift_hbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
						}
					}
				}
				else
				{
					for (j = x0; j <= x1; j++)
						s_lift_vbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
					if (y1 == nLastRow)
					{
						for (i = 1; i < nLastRow; i++)
						{
							for (j = x0; j <= x1; j++)
								s_lift_e(a, b, pCoeff, i, j, nPitch1, nStep);
						}
						for (j = x0; j <= x1; j++)
							s_lift_vbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
					}
					else
					{
						for (i = 1; i <= y1; i++)
						{
							for (j = x0; j <= x1; j++)
								s_lift_e(a, b, pCoeff, i, j, nPitch1, nStep);
						}
					}
				}
			}
		}
		else
		{
			if (x0 == 0)
			{
				if (x1 == nRowStop)
				{
					if (y1 == nLastRow)
					{
						for (i = y0; i < nLastRow; i++)
						{
							j = 0;
							s_lift_hbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
							for (j = 1; j < nRowStop; j++)
								s_lift_e(a, b, pCoeff, i, j, nPitch1, nStep);
							s_lift_hbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
						}
						j = 0;
						s_lift_corner_e(a, b, pCoeff, i, j, nPitch1, nStep);
						for (j = 1; j < nRowStop; j++)
							s_lift_vbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
						s_lift_corner_e(a, b, pCoeff, i, j, nPitch1, nStep);
					}
					else
					{
						for (i = y0; i <= y1; i++)
						{
							j = 0;
							s_lift_hbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
							for (j = 1; j < nRowStop; j++)
								s_lift_e(a, b, pCoeff, i, j, nPitch1, nStep);
							s_lift_hbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
						}
					}
				}
				else
				{
					if (y1 == nLastRow)
					{
						for (i = y0; i < nLastRow; i++)
						{
							j = 0;
							s_lift_hbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
							for (j = 1; j <= x1; j++)
								s_lift_e(a, b, pCoeff, i, j, nPitch1, nStep);
						}
						j = 0;
						s_lift_corner_e(a, b, pCoeff, i, j, nPitch1, nStep);
						for (j = 1; j <= x1; j++)
							s_lift_vbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
					}
					else
					{
						for (i = y0; i <= y1; i++)
						{
							j = 0;
							s_lift_hbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
							for (j = 1; j <= x1; j++)
								s_lift_e(a, b, pCoeff, i, j, nPitch1, nStep);
						}
					}
				}
			}
			else
			{
				if (x1 == nRowStop)
				{
					if (y1 == nLastRow)
					{
						for (i = y0; i < nLastRow; i++)
						{
							for (j = x0; j < nRowStop; j++)
								s_lift_e(a, b, pCoeff, i, j, nPitch1, nStep);
							s_lift_hbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
						}
						for (j = x0; j < nRowStop; j++)
							s_lift_vbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
						s_lift_corner_e(a, b, pCoeff, i, j, nPitch1, nStep);
					}
					else
					{
						for (i = y0; i <= y1; i++)
						{
							for (j = x0; j < nRowStop; j++)
								s_lift_e(a, b, pCoeff, i, j, nPitch1, nStep);
							s_lift_hbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
						}
					}
				}
				else
				{
					if (y1 == nLastRow)
					{
						for (i = y0; i < nLastRow; i++)
						{
							for (j = x0; j <= x1; j++)
								s_lift_e(a, b, pCoeff, i, j, nPitch1, nStep);
						}
						for (j = x0; j <= x1; j++)
							s_lift_vbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
					}
					else
					{
						for (i = y0; i <= y1; i++)
						{
							for (j = x0; j <= x1; j++)
								s_lift_e(a, b, pCoeff, i, j, nPitch1, nStep);
						}
					}
				}
			}
		}
	}

	static void w_lift_face_subpass(double a, double b, T* pCoeff, int nWidth, int nHeight, int nPitch1, int nStep, int x0, int y0, int x1, int y1)
	{
		int i, j, k;

		int nWidth1 = 2 * nWidth - 1;
		int nHeight1 = 2 * nHeight - 1;
		int nLastRow = nHeight - 1, nRowStop = nWidth - 1;

		x0 = max(x0, 0);
		y0 = max(y0, 0);
		x1 = min(x1, nRowStop);
		y1 = min(y1, nLastRow);

		for (i = y0; i < y1; i++)
		{
			for (j = x0; j < x1; j++)
				w_lift_f(a, b, pCoeff, i, j, nPitch1, nStep);
		}
	}

	static void w_lift_edge_subpass(double a, double b, T* pCoeff, int nWidth, int nHeight, int nPitch1, int nStep, int x0, int y0, int x1, int y1)
	{
		int i, j, k;

		int nWidth1 = 2 * nWidth - 1;
		int nHeight1 = 2 * nHeight - 1;
		int nLastRow = nHeight - 1, nRowStop = nWidth - 1;

		x0 = max(x0, 0);
		y0 = max(y0, 0);
		x1 = min(x1, nRowStop);
		y1 = min(y1, nLastRow);

		i = 0;
		j = 0;
		if (y0 == 0)
		{
			if (x0 == 0)
			{
				w_lift_corner_e(a, b, pCoeff, i, j, nPitch1, nStep);
				if (x1 == nRowStop)
				{
					for (j = 1; j < nRowStop; j++)
						w_lift_vbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
					w_lift_corner_e(a, b, pCoeff, i, j, nPitch1, nStep);
					if (y1 == nLastRow)
					{
						for (i = 1; i < nLastRow; i++)
						{
							j = 0;
							w_lift_hbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
							for (j = 1; j < nRowStop; j++)
								w_lift_e(a, b, pCoeff, i, j, nPitch1, nStep);
							w_lift_hbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
						}
						j = 0;
						w_lift_corner_e(a, b, pCoeff, i, j, nPitch1, nStep);
						for (j = 1; j < nRowStop; j++)
							w_lift_vbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
						w_lift_corner_e(a, b, pCoeff, i, j, nPitch1, nStep);
					}
					else
					{
						for (i = 1; i <= y1; i++)
						{
							j = 0;
							w_lift_hbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
							for (j = 1; j < nRowStop; j++)
								w_lift_e(a, b, pCoeff, i, j, nPitch1, nStep);
							w_lift_hbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
						}
					}
				}
				else
				{
					for (j = 1; j <= x1; j++)
						w_lift_vbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
					if (y1 == nLastRow)
					{
						for (i = 1; i < nLastRow; i++)
						{
							j = 0;
							w_lift_hbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
							for (j = 1; j <= x1; j++)
								w_lift_e(a, b, pCoeff, i, j, nPitch1, nStep);
						}
						j = 0;
						w_lift_corner_e(a, b, pCoeff, i, j, nPitch1, nStep);
						for (j = 1; j <= x1; j++)
							w_lift_vbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
					}
					else
					{
						for (i = 1; i <= y1; i++)
						{
							j = 0;
							w_lift_hbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
							for (j = 1; j <= x1; j++)
								w_lift_e(a, b, pCoeff, i, j, nPitch1, nStep);
						}
					}
				}
			}
			else
			{
				if (x1 == nRowStop)
				{
					for (j = x0; j < nRowStop; j++)
						w_lift_vbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
					w_lift_corner_e(a, b, pCoeff, i, j, nPitch1, nStep);
					if (y1 == nLastRow)
					{
						for (i = 1; i < nLastRow; i++)
						{
							for (j = x0; j < nRowStop; j++)
								w_lift_e(a, b, pCoeff, i, j, nPitch1, nStep);
							w_lift_hbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
						}
						for (j = x0; j < nRowStop; j++)
							w_lift_vbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
						w_lift_corner_e(a, b, pCoeff, i, j, nPitch1, nStep);
					}
					else
					{
						for (i = 1; i <= y1; i++)
						{
							for (j = x0; j < nRowStop; j++)
								w_lift_e(a, b, pCoeff, i, j, nPitch1, nStep);
							w_lift_hbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
						}
					}
				}
				else
				{
					for (j = x0; j <= x1; j++)
						w_lift_vbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
					if (y1 == nLastRow)
					{
						for (i = 1; i < nLastRow; i++)
						{
							for (j = x0; j <= x1; j++)
								w_lift_e(a, b, pCoeff, i, j, nPitch1, nStep);
						}
						for (j = x0; j <= x1; j++)
							w_lift_vbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
					}
					else
					{
						for (i = 1; i <= y1; i++)
						{
							for (j = x0; j <= x1; j++)
								w_lift_e(a, b, pCoeff, i, j, nPitch1, nStep);
						}
					}
				}
			}
		}
		else
		{
			if (x0 == 0)
			{
				if (x1 == nRowStop)
				{
					if (y1 == nLastRow)
					{
						for (i = y0; i < nLastRow; i++)
						{
							j = 0;
							w_lift_hbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
							for (j = 1; j < nRowStop; j++)
								w_lift_e(a, b, pCoeff, i, j, nPitch1, nStep);
							w_lift_hbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
						}
						j = 0;
						w_lift_corner_e(a, b, pCoeff, i, j, nPitch1, nStep);
						for (j = 1; j < nRowStop; j++)
							w_lift_vbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
						w_lift_corner_e(a, b, pCoeff, i, j, nPitch1, nStep);
					}
					else
					{
						for (i = y0; i <= y1; i++)
						{
							j = 0;
							w_lift_hbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
							for (j = 1; j < nRowStop; j++)
								w_lift_e(a, b, pCoeff, i, j, nPitch1, nStep);
							w_lift_hbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
						}
					}
				}
				else
				{
					if (y1 == nLastRow)
					{
						for (i = y0; i < nLastRow; i++)
						{
							j = 0;
							w_lift_hbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
							for (j = 1; j <= x1; j++)
								w_lift_e(a, b, pCoeff, i, j, nPitch1, nStep);
						}
						j = 0;
						w_lift_corner_e(a, b, pCoeff, i, j, nPitch1, nStep);
						for (j = 1; j <= x1; j++)
							w_lift_vbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
					}
					else
					{
						for (i = y0; i <= y1; i++)
						{
							j = 0;
							w_lift_hbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
							for (j = 1; j <= x1; j++)
								w_lift_e(a, b, pCoeff, i, j, nPitch1, nStep);
						}
					}
				}
			}
			else
			{
				if (x1 == nRowStop)
				{
					if (y1 == nLastRow)
					{
						for (i = y0; i < nLastRow; i++)
						{
							for (j = x0; j < nRowStop; j++)
								w_lift_e(a, b, pCoeff, i, j, nPitch1, nStep);
							w_lift_hbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
						}
						for (j = x0; j < nRowStop; j++)
							w_lift_vbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
						w_lift_corner_e(a, b, pCoeff, i, j, nPitch1, nStep);
					}
					else
					{
						for (i = y0; i <= y1; i++)
						{
							for (j = x0; j < nRowStop; j++)
								w_lift_e(a, b, pCoeff, i, j, nPitch1, nStep);
							w_lift_hbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
						}
					}
				}
				else
				{
					if (y1 == nLastRow)
					{
						for (i = y0; i < nLastRow; i++)
						{
							for (j = x0; j <= x1; j++)
								w_lift_e(a, b, pCoeff, i, j, nPitch1, nStep);
						}
						for (j = x0; j <= x1; j++)
							w_lift_vbound_e(a, b, pCoeff, i, j, nPitch1, nStep);
					}
					else
					{
						for (i = y0; i <= y1; i++)
						{
							for (j = x0; j <= x1; j++)
								w_lift_e(a, b, pCoeff, i, j, nPitch1, nStep);
						}
					}
				}
			}
		}
	}
};


template<typename T>
void CBSplineGridWavelet<T>::init()
{
	int j = m_nLevels - 1;

	clear();

	m_nFullWidth = getLevelWidth(j);
	m_nFullHeight = getLevelHeight(j);

	m_pData = new T[m_nFullWidth * m_nFullHeight];
	_ASSERT(m_pData != NULL);

	memset(m_pData, 0, m_nFullWidth * m_nFullHeight * sizeof(T));
}

template<typename T>
void CBSplineGridWavelet<T>::clear()
{
	if (m_pData != NULL)
	{
		delete [] m_pData;
		m_pData = NULL;
	}
}

template<typename T>
void CBSplineGridWavelet<T>::decompose(int level)
{
	if (level >= m_nLevels)
	{
		decompose(m_nLevels - 1);
		return;
	}
	else if (level < 0)
	{
		decompose(0);
		return;
	}
	else if (m_iLevel > level + 1)
		decompose(level + 1);

	if (m_iLevel == level + 1)
	{
		int nWidth = getLevelWidth(level);
		int nHeight = getLevelHeight(level);
		int nPitch1 = getFullResWidth();
		int nStep = getStep(m_iLevel);

		//linear basis function
		if(!m_useCubicBasis)
		{
			//w-lift faces
			w_lift_face_pass(-0.5, 1.0, m_pData, nWidth, nHeight, nPitch1, nStep);

			//w-lift edges
			w_lift_edge_pass(-0.5, 1.0, m_pData, nWidth, nHeight, nPitch1, nStep);

			//s-lift vertices
			s_lift_vertex_pass(0.25, 1.0, m_pData, nWidth, nHeight, nPitch1, nStep);

			//s-lift edges
			s_lift_edge_pass(0.25, 1.0, m_pData, nWidth, nHeight, nPitch1, nStep);
		}
		//cubic basis function
		else
		{
			//s-lift vertices
			s_lift_vertex_pass(-0.25, 1.0, m_pData, nWidth, nHeight, nPitch1, nStep);

			//s-lift edges
			s_lift_edge_pass(-0.25, 1.0, m_pData, nWidth, nHeight, nPitch1, nStep);

			//w-lift faces
			w_lift_face_pass(-1.0, 1.0, m_pData, nWidth, nHeight, nPitch1, nStep);

			//w-lift edges
			w_lift_edge_pass(-1.0, 1.0, m_pData, nWidth, nHeight, nPitch1, nStep);

			//s-lift vertices
			s_lift_vertex_pass(0.375, 2.0, m_pData, nWidth, nHeight, nPitch1, nStep);

			//s-lift edges
			s_lift_edge_pass(0.375, 2.0, m_pData, nWidth, nHeight, nPitch1, nStep);
		}

		m_iLevel = level;
	}
}

template<typename T>
void CBSplineGridWavelet<T>::reconstruct(int level)
{
	if (level >= m_nLevels)
	{
		reconstruct(m_nLevels - 1);
		return;
	}
	else if (level < 0)
	{
		reconstruct(0);
		return;
	}
	else if (m_iLevel < level - 1)
		reconstruct(level - 1);

	if (m_iLevel == level - 1)
	{
		//linear basis function
		if(!m_useCubicBasis)
		{
			//s-lift vertices
			s_lift_vertex_pass(-0.25, 1.0, m_pData, getWidth(), getHeight(), getFullResWidth(), getStep(level));

			//s-lift edges
			s_lift_edge_pass(-0.25, 1.0, m_pData, getWidth(), getHeight(), getFullResWidth(), getStep(level));

			//w-lift faces
			w_lift_face_pass(0.5, 1.0, m_pData, getWidth(), getHeight(), getFullResWidth(), getStep(level));

			//w-lift edges
			w_lift_edge_pass(0.5, 1.0, m_pData, getWidth(), getHeight(), getFullResWidth(), getStep(level));
		}
		//cubic basis function
		else
		{
			//s-lift vertices
			s_lift_vertex_pass(-0.1875, 0.5, m_pData, getWidth(), getHeight(), getFullResWidth(), getStep(level));

			//s-lift edges
			s_lift_edge_pass(-0.1875, 0.5, m_pData, getWidth(), getHeight(), getFullResWidth(), getStep(level));

			//w-lift faces
			w_lift_face_pass(1.0, 1.0, m_pData, getWidth(), getHeight(), getFullResWidth(), getStep(level));

			//w-lift edges
			w_lift_edge_pass(1.0, 1.0, m_pData, getWidth(), getHeight(), getFullResWidth(), getStep(level));

			//s-lift vertices
			s_lift_vertex_pass(0.25, 1.0, m_pData, getWidth(), getHeight(), getFullResWidth(), getStep(level));

			//s-lift edges
			s_lift_edge_pass(0.25, 1.0, m_pData, getWidth(), getHeight(), getFullResWidth(), getStep(level));
		}

		m_iLevel = level;
	}
}

template<typename T>
void CBSplineGridWavelet<T>::reconstructCopy(int level, T *pDest, int nPitch, int nStep)
{
	int i, j, iSrc, iDest;
	int nWidth, nHeight;
	int nSrcStep, nSrcPitchStep, nPitchStep;

	if (level >= m_nLevels)
	{
		reconstructCopy(m_nLevels - 1, pDest, nPitch, nStep);
		return;
	}
	else if (level < 0)
	{
		reconstructCopy(0, pDest, nPitch, nStep);
		return;
	}
#if 0
	else if (m_iLevel < level - 1)
		reconstructCopy(level - 1, pDest, nPitch, 2 * nStep);
	else if (m_iLevel == level - 1)
	{
		nWidth = getLevelWidth(level);
		nHeight = getLevelHeight(level);
		nSrcStep = getStep(level);
		nPitchStep = nStep * nPitch;
		nSrcPitchStep = nSrcStep * getFullResWidth();
		for (i = 0; i < nHeight; i++)
		{
			iSrc = i * nSrcPitchStep;
			iDest = i * nPitchStep;
			for (j = 0; j < nWidth; j++)
			{
				pDest[iDest] = m_pData[iSrc];
				iSrc += nSrcStep;
				iDest += nStep;
			}
		}
	}
	else if (m_iLevel == level)
	{
		nWidth = getLevelWidth(level);
		nHeight = getLevelHeight(level);
		nSrcStep = getStep(level);
		nPitchStep = nStep * nPitch;
		nSrcPitchStep = nSrcStep * getFullResWidth();
		for (i = 0; i < nHeight; i++)
		{
			iSrc = i * nSrcPitchStep;
			iDest = i * nPitchStep;
			for (j = 0; j < nWidth; j++)
			{
				pDest[iDest] = m_pData[iSrc];
				iSrc += nSrcStep;
				iDest += nStep;
			}
		}
		return;
	}
	else
		return;

	//s-lift vertices
	s_lift_vertex_pass(-0.25, 1.0, pDest, getLevelWidth(level - 1), getLevelHeight(level - 1), nPitch, nStep);

	//s-lift edges
	s_lift_edge_pass(-0.25, 1.0, pDest, getLevelWidth(level - 1), getLevelHeight(level - 1), nPitch, nStep);

	//w-lift faces
	w_lift_face_pass(0.5, 1.0, pDest, getLevelWidth(level - 1), getLevelHeight(level - 1), nPitch, nStep);

	//w-lift edges
	w_lift_edge_pass(0.5, 1.0, pDest, getLevelWidth(level - 1), getLevelHeight(level - 1), nPitch, nStep);
#else
	nWidth = getLevelWidth(level);
	nHeight = getLevelHeight(level);
	nSrcStep = getStep(level);
	nPitchStep = nStep * nPitch;
	nSrcPitchStep = nSrcStep * getFullResWidth();
	for (i = 0; i < nHeight; i++)
	{
		iSrc = i * nSrcPitchStep;
		iDest = i * nPitchStep;
		for (j = 0; j < nWidth; j++)
		{
			pDest[iDest] = m_pData[iSrc];
			iSrc += nSrcStep;
			iDest += nStep;
		}
	}

	j = m_iLevel;
	while (j < level)
	{
		i = nStep * (1 << (level - j - 1));
		nWidth = getLevelWidth(j);
		nHeight = getLevelHeight(j);
		//linear basis function
		if(!m_useCubicBasis)
		{
			//s-lift vertices
			s_lift_vertex_pass(-0.25, 1.0, pDest, nWidth, nHeight, nPitch, i);

			//s-lift edges
			s_lift_edge_pass(-0.25, 1.0, pDest, nWidth, nHeight, nPitch, i);

			//w-lift faces
			w_lift_face_pass(0.5, 1.0, pDest, nWidth, nHeight, nPitch, i);

			//w-lift edges
			w_lift_edge_pass(0.5, 1.0, pDest, nWidth, nHeight, nPitch, i);
		}
		//cubic basis function
		else
		{
			//s-lift vertices
			s_lift_vertex_pass(-0.1875, 0.5, pDest, nWidth, nHeight, nPitch, i);

			//s-lift edges
			s_lift_edge_pass(-0.1875, 0.5, pDest, nWidth, nHeight, nPitch, i);

			//w-lift faces
			w_lift_face_pass(1.0, 1.0, pDest, nWidth, nHeight, nPitch, i);

			//w-lift edges
			w_lift_edge_pass(1.0, 1.0, pDest, nWidth, nHeight, nPitch, i);

			//s-lift vertices
			s_lift_vertex_pass(0.25, 1.0, pDest, nWidth, nHeight, nPitch, i);

			//s-lift edges
			s_lift_edge_pass(0.25, 1.0, pDest, nWidth, nHeight, nPitch, i);
		}

		j++;
	}
#endif
}

template<typename T>
void CBSplineGridWavelet<T>::reconstructCopyLocal(int iCoeff, int level, T* pDest, int nPitch, int nStep)
{
	int i, j, iSrc, iDest;
	int x, y, xstart, xstop, ystart, ystop;
	int xstart1, xstop1, ystart1, ystop1;
	int nWidth, nHeight;
	int nSrcStep, nSrcPitchStep, nPitchStep;

	if (level >= m_nLevels)
	{
		reconstructCopy(m_nLevels - 1, pDest, nPitch, nStep);
		return;
	}
	else if (level < 0)
	{
		reconstructCopy(0, pDest, nPitch, nStep);
		return;
	}

	getLevelCoordinates(iCoeff, m_iLevel + 1, &x, &y);
	xstart1 = x - 1;
	xstop1 = x + 1;
	ystart1 = y - 1;
	ystop1 = y + 1;

	for (j = m_iLevel; j < level; j++)
	{
		xstart1 = xstart1 * 2 - 1;
		xstop1 = xstop1 * 2 + 1;
		ystart1 = ystart1 * 2 - 1;
		ystop1 = ystop1 * 2 + 1;
	}

	nWidth = getLevelWidth(level);
	nHeight = getLevelHeight(level);
	nSrcStep = getStep(level);
	nPitchStep = nStep * nPitch;
	nSrcPitchStep = nSrcStep * getFullResWidth();
	for (i = ystart; i <= ystop && i < nHeight; i++)
	{
		iSrc = i * nSrcPitchStep;
		iDest = i * nPitchStep;
		for (j = xstart; j <= xstop && j < nWidth; j++)
		{
			pDest[iDest] = m_pData[iSrc];
			iSrc += nSrcStep;
			iDest += nStep;
		}
	}

	xstart1 = x - 1;
	xstop1 = x + 1;
	ystart1 = y - 1;
	ystop1 = y + 1;
	j = m_iLevel;
	while (j < level)
	{
		i = nStep * (1 << (level - j - 1));
		nWidth = getLevelWidth(j);
		nHeight = getLevelHeight(j);

		xstart = xstart1 / 2;
		xstop = xstart1 / 2;
		ystart = ystart1 / 2;
		ystop = ystop1 / 2;

		//s-lift vertices
		s_lift_vertex_subpass(-0.25, 1.0, pDest, nWidth, nHeight, nPitch, i, xstart, ystart, xstop, ystop);

		//s-lift edges
		s_lift_edge_subpass(-0.25, 1.0, pDest, nWidth, nHeight, nPitch, i, xstart, ystart, xstop, ystop);

		//w-lift faces
		w_lift_face_subpass(0.5, 1.0, pDest, nWidth, nHeight, nPitch, i, xstart, ystart, xstop, ystop);

		//w-lift edges
		w_lift_edge_subpass(0.5, 1.0, pDest, nWidth, nHeight, nPitch, i, xstart, ystart, xstop, ystop);

		xstart1 = xstart1 * 2 - 1;
		xstop1 = xstop1 * 2 + 1;
		ystart1 = ystart1 * 2 - 1;
		ystop1 = ystop1 * 2 + 1;
		j++;
	}
}

template<typename T>
void CBSplineGridWavelet<T>::save(char *szFilename)
{
	FILE* pf = fopen(szFilename, "wb");
	_ASSERT(pf != NULL);

//	fprintf(pf, "BSGW %08i %08i %08i %08i\n", m_nBaseWidth, m_nBaseHeight, m_nLevels, m_eType);
//	fprintf(pf, "BSGW %8i %8i %8i %8i\n", m_nBaseWidth, m_nBaseHeight, m_nLevels, m_eType);
//	printf("BSGW %08i %08i %08i %08i\n", m_nBaseWidth, m_nBaseHeight, m_nLevels, m_eType);
//	printf("BSGW %8i %8i %8i %8i\n", m_nBaseWidth, m_nBaseHeight, m_nLevels, m_eType);
	fwrite(&m_nBaseWidth, sizeof(int), 1, pf);
	fwrite(&m_nBaseHeight, sizeof(int), 1, pf);
	fwrite(&m_nLevels, sizeof(int), 1, pf);
	fwrite(&m_eType, sizeof(EBSGWType), 1, pf);
	fwrite(m_pData, sizeof(T), m_nFullHeight * m_nFullWidth, pf);

	fclose(pf);
}


#endif //__BSPLINEGRIDWAVELET_H__


