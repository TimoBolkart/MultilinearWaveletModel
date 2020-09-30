#ifndef MULTILINEARMODELHANLDER_H
#define MULTILINEARMODELHANLDER_H

#include "DataContainer.h"
#include "MultilinearModel.h"

#include <vector>
#include <string>

class ANNkd_tree;

class MultilinearModelHandler
{
public:
	MultilinearModelHandler();

	~MultilinearModelHandler();

	//Clear the data structures
	void clear();

	
	//Reconstruct coordinates of the model by mode-multiplications of the weight vectors with the multilinear model
	//order of weights: weightVector = w2_1 ... w2_m2 w3_1 ... w3_m3 ... wN_1 ... wN_mN
	void reconstructForWeights(const std::vector<double>& weightVector, std::vector<double>& outPoints);

	
	//!Get mean coefficients of mode.
	//! \param mode		specific mode (normally 2 = shape mode, 3 = expression mode)
	void getModeMean(const size_t mode, std::vector<double>& modeMean);

	//! Get upper bound of multilinear model space.
	//! \param mode		specific mode (normally 2 = shape mode, 3 = expression mode)
	//! \param k			edge length of hyperbox for current component (k >= 0.0)
	//! \param index		coefficient index
	//! \return				upper bound
	double getUpperBound(const size_t mode, const double k, const size_t index);

	//! Get lower bound of multilinear model space.
	//! \param mode		specific mode (normally 2 = shape mode, 3 = expression mode)
	//! \param k			edge length of hyperbox for current component (k >= 0.0)
	//! \param index		coefficient index
	//! \return				lower bound
	double getLowerBound(const size_t mode, const double k, const size_t index);
	
	bool importMultilinearModelBinary(std::fstream& input);
	
	bool importMultilinearModelText(std::fstream& input);

	void getModeDimensions(std::vector<size_t>& outVec)
	{
		copyArray(m_modeDimensions, outVec);
	}
	
	void getTruncModeDimensions(std::vector<size_t>& outVec)
	{
		copyArray(m_truncModeDimensions, outVec);
	}
	
	void getMean(std::vector<double>& outVec)
	{
		copyArray(m_mean, outVec);
	}

	Tensor& getTensor()		{ return m_tensor; }

private:
	double getWeightForVariation(const double variation, const double meanWeight, const size_t mode)
	{
		if((mode<1) || (m_modeDimensions.size() < mode-1))
		{
			return 0.0;
		}

		const double tmp = sqrt(static_cast<double>(m_modeDimensions[mode-1]));
		return variation/tmp+meanWeight;
	}

	double getVariationForWeight(const double weight, const double meanWeight, const size_t mode)
	{
		if((mode<1) || (m_modeDimensions.size() < mode-1))
		{
			return 0.0;
		}

		const double tmp = sqrt(static_cast<double>(m_modeDimensions[mode-1]));
		return (weight-meanWeight)*tmp;
	}

	template<typename T>
	void copyArray(const std::vector<T>& inVec, std::vector<T>& outVec)
	{
		outVec.clear();
		
		const size_t numElements = inVec.size();
		outVec.reserve(numElements);

		for(size_t i = 0; i < numElements; ++i)
		{
			outVec.push_back(inVec[i]);
		}
	}

	std::vector<size_t> m_modeDimensions;
	std::vector<size_t> m_truncModeDimensions;
	std::vector<double> m_mean;
	std::vector<double> m_meanWeights;
	std::vector<double> m_singularValues;

	Tensor m_tensor;
};
#endif