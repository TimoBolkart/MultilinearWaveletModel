///////////////////////////////////////////////////////////////////////////////
//
//	MMProjectionCostFunction_Local.h
//
//	Header file for the MMProjectionCostFunction_Local class
//	subclass of vnl_cost_function
//	used for optimizing a local multilinear model wrt some input data
//
//	Alan Brunton, September 2013
//
///////////////////////////////////////////////////////////////////////////////


#ifndef __MMPROJECTIONCOSTFUNCTIONLOCAL_H__
#define __MMPROJECTIONCOSTFUNCTIONLOCAL_H__


#include <vnl/vnl_vector.h>
#include <vnl/vnl_cost_function.h>

#include <vector>


class Tensor;


class MMProjectionCostFunction_Local : public vnl_cost_function
{
private:

	Tensor*									m_pMultilinearModel;
	size_t									m_nTargets;
	size_t									m_nTargetDim;
	const std::vector<double>*			m_pTarget;
	const std::vector<double>*			m_pTargetMask;
	const std::vector<double>*			m_pReconOffset;
	std::vector<double>					m_mean;

	std::vector<double>					m_tmpRecon;

	double												m_smoothNeighborWeight;
	const std::vector<std::vector<size_t>>*	m_precomputedSmoothIndices;
	const std::vector<std::vector<double>>*	m_precomputedSmoothWeights;

	//for fixing one or more modes
	std::vector<bool>						m_modeFixed;
	std::vector<double>					m_fixedModeValues;

	//members for local-specific stuff
	std::vector<double>						m_transformCoeffs;
	std::vector<size_t>						m_transformIndices;

public:

	//MMProjectionCostFunction_Local(Tensor* pMM, const std::vector<double>& mean, const std::vector<double>& reconOffset, const std::vector<double>& target, const std::vector<double>& targetMask, const std::vector<double>& transCoeffs, const std::vector<size_t>& transIndices);
	//MMProjectionCostFunction_Local(Tensor* pMM, const std::vector<double>& mean, const std::vector<double>& reconOffset, const size_t numTargets, const size_t targetDim, const std::vector<double>& target, const std::vector<double>& targetMask, const std::vector<double>& transCoeffs, const std::vector<size_t>& transIndices);

	MMProjectionCostFunction_Local(Tensor* pMM, const std::vector<double>& mean, const std::vector<double>& reconOffset, const size_t numTargets, const size_t targetDim, const std::vector<double>& target, const std::vector<double>& targetMask
											, const std::vector<double>& transCoeffs, const std::vector<size_t>& transIndices, const double smoothNeighborWeight, const std::vector<std::vector<size_t>>& precomputedSmoothIndices, const std::vector<std::vector<double>>& precomputedSmoothWeights);

	~MMProjectionCostFunction_Local();

	void fixMode(size_t mode, const std::vector<double>& values);
	bool isModeFixed(size_t mode)
	{
		if (mode == 2)
			return m_modeFixed[0];
		if (mode == 3)
			return m_modeFixed[1];
		return false;
	}

	virtual void compute(const vnl_vector<double>& x, double *f, vnl_vector<double>* g);

private:
	bool computeBiLaplacianEnergy(double& neighborSmoothingEnergy, std::vector<double>& neighborSmoothingGradient);
};


#endif //__MMPROJECTIONCOSTFUNCTIONLOCAL_H__

