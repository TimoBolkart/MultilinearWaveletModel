///////////////////////////////////////////////////////////////////////////////
//
//	MMProjectionCostFunction_Local.cpp
//
//	source file for the class MMProjectionCostFunction_Local
//	subclass of vnl_cost_function
//
//	Alan Brunton, September 2013
//
///////////////////////////////////////////////////////////////////////////////


#include "MMProjectionCostFunction_Local.h"
#include "MultilinearModel.h"


//MMProjectionCostFunction_Local::MMProjectionCostFunction_Local(Tensor* pMM, const std::vector<double>& mean, const std::vector<double>& reconOffset, const std::vector<double>& target, const std::vector<double>& targetMask, const std::vector<double>& transCoeffs, const std::vector<size_t>& transIndices): 
//	vnl_cost_function(pMM->getModeDimension(2) + pMM->getModeDimension(3)), m_pMultilinearModel(pMM)
//{
//	m_mean = mean;
//	m_pReconOffset = &reconOffset;
//	m_nTargets = 1;
//	m_nTargetDim = targetMask.size();
//	m_pTarget = &target;
//	m_pTargetMask = &targetMask;
//	m_transformCoeffs = transCoeffs;
//	m_transformIndices = transIndices;
//	m_modeFixed.resize(2, false);
//	m_fixedModeValues.resize(pMM->getModeDimension(2) + pMM->getModeDimension(3), 0.0);
//}
//
//MMProjectionCostFunction_Local::MMProjectionCostFunction_Local(Tensor* pMM, const std::vector<double>& mean, const std::vector<double>& reconOffset, const size_t numTargets, const size_t targetDim, const std::vector<double>& target, const std::vector<double>& targetMask, const std::vector<double>& transCoeffs, const std::vector<size_t>& transIndices): 
//	vnl_cost_function(pMM->getModeDimension(2) + pMM->getModeDimension(3)), m_pMultilinearModel(pMM)
//{
//	m_mean = mean;
//	m_pReconOffset = &reconOffset;
//	m_nTargets = numTargets;
//	m_nTargetDim = targetDim;
//	m_pTarget = &target;
//	m_pTargetMask = &targetMask;
//	m_transformCoeffs = transCoeffs;
//	m_transformIndices = transIndices;
//	m_modeFixed.resize(2, false);
//	m_fixedModeValues.resize(pMM->getModeDimension(2) + pMM->getModeDimension(3), 0.0);
//}

MMProjectionCostFunction_Local::MMProjectionCostFunction_Local(Tensor* pMM, const std::vector<double>& mean, const std::vector<double>& reconOffset, const size_t numTargets, const size_t targetDim, const std::vector<double>& target, const std::vector<double>& targetMask
																					, const std::vector<double>& transCoeffs, const std::vector<size_t>& transIndices, const double smoothNeighborWeight, const std::vector<std::vector<size_t>>& precomputedSmoothIndices, const std::vector<std::vector<double>>& precomputedSmoothWeights)
: vnl_cost_function(pMM->getModeDimension(2) + pMM->getModeDimension(3)), m_pMultilinearModel(pMM)
{
	m_mean = mean;
	m_pReconOffset = &reconOffset;
	m_nTargets = numTargets;
	m_nTargetDim = targetDim;
	m_pTarget = &target;
	m_pTargetMask = &targetMask;
	m_transformCoeffs = transCoeffs;
	m_transformIndices = transIndices;
	m_modeFixed.resize(2, false);
	m_fixedModeValues.resize(pMM->getModeDimension(2) + pMM->getModeDimension(3), 0.0);
	//m_pNeighborIndices = &neighborIndices;
	//m_pInvalidSmoothNeighbors = &invalidSmoothNeighbors;
	m_smoothNeighborWeight = smoothNeighborWeight;
	m_precomputedSmoothIndices = &precomputedSmoothIndices;
	m_precomputedSmoothWeights = &precomputedSmoothWeights;

	m_tmpRecon = *m_pReconOffset;
}

MMProjectionCostFunction_Local::~MMProjectionCostFunction_Local()
{
}

void MMProjectionCostFunction_Local::fixMode(size_t mode, const std::vector<double>& values)
{
	const size_t m2(m_pMultilinearModel->getModeDimension(2));
	const size_t m3(m_pMultilinearModel->getModeDimension(3));

	if (mode == 2)
	{
		m_modeFixed[0] = true;
		set_number_of_unknowns(m3);
		for (size_t i = 0; i < m2; i++)
			m_fixedModeValues[i] = values[i];
	}
	else if (mode == 3)
	{
		m_modeFixed[1] = true;
		set_number_of_unknowns(m2);
		for (size_t i = 0; i < m3; i++)
			m_fixedModeValues[i + m2] = values[i];
	}
}

void MMProjectionCostFunction_Local::compute(const vnl_vector<double>& x, double *f, vnl_vector<double>* g)
{
	const size_t d1(m_pMultilinearModel->getModeDimension(1));
	const size_t m2(m_pMultilinearModel->getModeDimension(2));
	const size_t m3(m_pMultilinearModel->getModeDimension(3));
	const size_t numWeights = m2 + m3;

	const std::vector<double>& target = *m_pTarget;
	const std::vector<double>& targetMask = *m_pTargetMask;
	const std::vector<double>& reconOffset = *m_pReconOffset;

	const size_t numTargs = m_nTargets;
	const size_t targStep = m_nTargetDim * d1;
	const size_t maskStep = m_nTargetDim;

	std::vector<double> w2, w3;
	const bool bM2Fixed = isModeFixed(2);
	const bool bM3Fixed = isModeFixed(3);
	const size_t offsetM3X = bM2Fixed ? 0 : m2;

	double energyValue(0.0);

	if (bM2Fixed)
	{
		for (size_t i = 0; i < m2; ++i)
			w2.push_back(m_fixedModeValues[i]);
	}
	else
	{
		for (size_t i = 0; i < m2; ++i)
			w2.push_back(x[i]);
	}

	if (bM3Fixed)
	{
		for (size_t i = 0; i < m3; ++i)
			w3.push_back(m_fixedModeValues[i + m2]);
	}
	else
	{
		for (size_t i = 0; i < m3; i++)
			w3.push_back(x[i + offsetM3X]);
	}

	Tensor M3;
	m_pMultilinearModel->modeMultiply(w2, "T", m2, 1, 2, M3);

	Tensor coeff;
	M3.modeMultiply(w3, "T", m3, 1, 3, coeff);

	const size_t numIndices = m_transformIndices.size();
	vnl_vector<double> tmp(numIndices * d1, 0.0);

	for (size_t i = 0; i < numIndices; ++i)
	{
		const size_t idx = m_transformIndices[i];
		const double tc = m_transformCoeffs[i];
		for (size_t j = 0; j < d1; j++)
		{
			const size_t m = idx * d1 + j;
			m_tmpRecon[m] = tc * (coeff.getElement(j, 0, 0, 0) + m_mean[j]) + reconOffset[m];
		}
	}

	//Recon is the reconstruction in x,y,z domain
	//Compute the smoothing energy and the gradient in x,y,z domain
	double neighborSmoothingEnergy(0.0);
	std::vector<double> neighborSmoothingGradient;
	const bool bSmooth = m_smoothNeighborWeight > 0.0 ? computeBiLaplacianEnergy(neighborSmoothingEnergy, neighborSmoothingGradient) : false;

	//Transform gradient to just contain all gradient values of vertices, that are influenced by current coefficient.
	std::vector<double> transformedSmoothingGradient;
	if(bSmooth)
	{
		transformedSmoothingGradient.resize(3*numIndices);

		for(size_t i = 0; i < numIndices; ++i)
		{
			const size_t idx = m_transformIndices[i];
			
			for (size_t j = 0; j < 3; j++)
			{
				const size_t k = i * 3 + j;
				const size_t m = idx * 3 + j;
				transformedSmoothingGradient[k] = neighborSmoothingGradient[m];
			}
		}
	}

	for (size_t t = 0; t < numTargs; t++)
	{
		const size_t targOff = t * targStep;
		const size_t maskOff = t * maskStep;
		for (size_t i = 0; i < numIndices; i++)
		{
			const size_t idx = m_transformIndices[i];
			const double mask = targetMask[idx + maskOff];
			//const double mask = targetMask[idx];
			for (size_t j = 0; j < d1; j++)
			{
				const size_t k = i * d1 + j;
				const size_t m1 = idx * d1 + j;
				const size_t m2 = idx * d1 + j + targOff;

				const double tmpDiff = m_tmpRecon[m1] - target[m2];
				energyValue += tmpDiff*tmpDiff * mask;
				tmp[k] += tmpDiff * mask;
			}
		}
	}

	*f = energyValue;

	if(bSmooth)
	{
		*f += m_smoothNeighborWeight*neighborSmoothingEnergy;
	}

	Tensor M2;
	m_pMultilinearModel->modeMultiply(w3, "T", m3, 1, 3, M2);

	if (!bM2Fixed)
	{
		if(bSmooth)
		{
			for(size_t i = 0; i < m2; ++i)
			{
				double out = 0.0;
				for(size_t j = 0; j < numIndices; ++j)
				{
					const size_t idx = m_transformIndices[j];
					const double tc = m_transformCoeffs[j];
					for (size_t k = 0; k < d1; k++)
					{
						const size_t m = j * d1 + k;
						const double mlc = M2.getElement(k, i, 0, 0);
						out += (tmp[m] + m_smoothNeighborWeight*transformedSmoothingGradient[m]) * tc * mlc;
					}
				}

				(*g)[i] = 2*out;
			}
		}
		else
		{
			for(size_t i = 0; i < m2; ++i)
			{
				double out = 0.0;
				for(size_t j = 0; j < numIndices; ++j)
				{
					const size_t idx = m_transformIndices[j];
					const double tc = m_transformCoeffs[j];
					for (size_t k = 0; k < d1; k++)
					{
						const size_t m = j * d1 + k;
						out += tmp[m] * tc * M2.getElement(k, i, 0, 0);
					}
				}

				(*g)[i] = 2*out;
			}
		}
	}

	if (!bM3Fixed)
	{
		if(bSmooth)
		{
			for(size_t i = 0; i < m3; ++i)
			{
				double out = 0.0;
				for(size_t j = 0; j < numIndices; ++j)
				{
					const size_t idx = m_transformIndices[j];
					const double tc = m_transformCoeffs[j];
					for (size_t k = 0; k < d1; k++)
					{
						const size_t m = j * d1 + k;
						const double mlc = M3.getElement(k, 0, i, 0);
						out += (tmp[m] + m_smoothNeighborWeight*transformedSmoothingGradient[m]) * tc * mlc;
					}
				}

				(*g)[i + offsetM3X] = 2*out;
			}
		}
		else
		{
			for(size_t i = 0; i < m3; ++i)
			{
				double out = 0.0;
				for(size_t j = 0; j < numIndices; ++j)
				{
					const size_t idx = m_transformIndices[j];
					const double tc = m_transformCoeffs[j];
					for (size_t k = 0; k < d1; k++)
					{
						const size_t m = j * d1 + k;
						out += tmp[m] * tc * M3.getElement(k, 0, i, 0);
					}
				}

				(*g)[i + offsetM3X] = 2*out;
			}
		}
	}
}

bool MMProjectionCostFunction_Local::computeBiLaplacianEnergy(double& neighborSmoothingEnergy, std::vector<double>& neighborSmoothingGradient)
{
	neighborSmoothingEnergy = 0.0;
	neighborSmoothingGradient.resize(3*m_nTargetDim, 0.0);

	const std::vector<std::vector<size_t>>& precomputedSmoothIndices = *m_precomputedSmoothIndices;
	const std::vector<std::vector<double>>& precomputedSmoothWeights = *m_precomputedSmoothWeights;
	if(precomputedSmoothIndices.empty() || precomputedSmoothIndices.size() != precomputedSmoothWeights.size())
	{
		return false;
	}
	
	const size_t numIndices = m_transformIndices.size();
	for (size_t i = 0; i < numIndices; ++i)
	{
		//Iterate over the indices of the vertices, that are influenced by currently processed coefficient
		const size_t pIndex = m_transformIndices[i];

		const std::vector<size_t>& currSmoothIndices = precomputedSmoothIndices[pIndex];
		if(currSmoothIndices.empty())
		{
			continue;
		}

		const std::vector<double>& currSmoothWeights = precomputedSmoothWeights[pIndex];

		double tmpx = 0.0, tmpy = 0.0, tmpz = 0.0;

		const size_t numCurrSmoothIndices = currSmoothIndices.size();
		for(size_t j = 0; j < numCurrSmoothIndices; ++j)
		{
			const size_t currSmoothStartIndex = 3*currSmoothIndices[j];
			const double currSmoothWeight = currSmoothWeights[j];

			tmpx += currSmoothWeight*m_tmpRecon[currSmoothStartIndex];
			tmpy += currSmoothWeight*m_tmpRecon[currSmoothStartIndex+1];
			tmpz += currSmoothWeight*m_tmpRecon[currSmoothStartIndex+2];
		}

		//neighborSmoothingEnergy += (tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2]);
		neighborSmoothingEnergy += (tmpx*tmpx + tmpy*tmpy + tmpz*tmpz);

		for(size_t j = 0; j < numCurrSmoothIndices; ++j)
		{
			//const size_t currSmoothIndex = currSmoothIndices[j];
			const size_t currSmoothStartIndex = 3*currSmoothIndices[j];
			const double currSmoothWeight = currSmoothWeights[j];

			neighborSmoothingGradient[currSmoothStartIndex] += /*2**/currSmoothWeight*tmpx;
			neighborSmoothingGradient[currSmoothStartIndex+1] += /*2**/currSmoothWeight*tmpy;
			neighborSmoothingGradient[currSmoothStartIndex+2] += /*2**/currSmoothWeight*tmpz;
		}
	}

	return true;
}
