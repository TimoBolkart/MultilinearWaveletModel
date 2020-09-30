#include "MultilinearModelHandler.h"	
#include "FileLoader.h"

namespace clapack
{
	extern "C"
	{
		#include "blaswrap.h"
		#include "f2c.h"
		extern int dgesdd_(char *jobz, integer *m, integer *n, doublereal *a
								, integer *lda, doublereal *s, doublereal *u, integer *ldu
								, doublereal *vt, integer *ldvt, doublereal *work, integer *lwork
								, integer *iwork, integer *info);
		extern int dgels_(char *trans, integer *m, integer *n, integer *nrhs, doublereal *a
								, integer *lda, doublereal *b, integer *ldb, doublereal *work
								, integer *lwork, integer *info);
	}
}

MultilinearModelHandler::MultilinearModelHandler()
{

}

MultilinearModelHandler::~MultilinearModelHandler()
{
	
}

void MultilinearModelHandler::clear()
{
	m_truncModeDimensions.clear();
	m_modeDimensions.clear();
	m_mean.clear();
	m_meanWeights.clear();
	m_singularValues.clear();
	m_tensor.clear();
}

void MultilinearModelHandler::reconstructForWeights(const std::vector<double>& weightVector, std::vector<double>& outPoints)
{
	if(m_truncModeDimensions.size()<3)
	{
		return;
	}

	std::vector<double> w2;
	w2.reserve(m_truncModeDimensions[1]);

	for(size_t i = 0; i < m_truncModeDimensions[1]; ++i)
	{
		w2.push_back(weightVector[i]);
	}

	Tensor t2;
	m_tensor.modeMultiply(w2, "T", m_truncModeDimensions[1], 1, 2, t2);

	std::vector<double> w3;
	w3.reserve(m_truncModeDimensions[2]);

	for(size_t i = 0; i < m_truncModeDimensions[2]; ++i)
	{
		const size_t index = m_truncModeDimensions[1]+i;
		w3.push_back(weightVector[index]);
	}

	Tensor result;
	t2.modeMultiply(w3, "T", m_truncModeDimensions[2], 1, 3, result);

	outPoints.clear();
	outPoints.reserve(m_truncModeDimensions[0]);

	for(size_t i = 0; i < m_truncModeDimensions[0]; ++i)
	{
		const double value = result.getElement(i, 0, 0, 0);
		outPoints.push_back(value);
	}

	for(size_t i = 0; i < m_mean.size(); ++i)
	{
		outPoints[i] += m_mean[i];
	}
}

void MultilinearModelHandler::getModeMean(const size_t mode, std::vector<double>& modeMean)
{
	const size_t m2 = m_truncModeDimensions[1];
	const size_t m3 = m_truncModeDimensions[2];

	if(m_meanWeights.size() != m2+m3)
	{
		return;
	}

	modeMean.clear();

	if(mode == 2)
	{
		modeMean.reserve(m2);

		for(size_t i = 0; i < m2; ++i)
		{
			modeMean.push_back(m_meanWeights[i]);
		}
	}
	else if(mode == 3)
	{
		modeMean.reserve(m3);

		for(size_t i = 0; i < m3; ++i)
		{
			modeMean.push_back(m_meanWeights[m2+i]);
		}
	}
}

double MultilinearModelHandler::getUpperBound(const size_t mode, const double k, const size_t index)
{
	if(mode != 2 && mode != 3)
	{
#ifdef DEBUG_OUTPUT
		std::cout << "getUpperBound(...) - wrong mode input" << std::endl;
#endif
		return 0.0;
	}

	const size_t m2 = m_truncModeDimensions[1];
	const size_t currIndex = index + (mode == 3 ? m2 : 0);
	return getWeightForVariation(k, m_meanWeights[currIndex], mode);
}

double MultilinearModelHandler::getLowerBound(const size_t mode, const double k, const size_t index)
{
	if(mode != 2 && mode != 3)
	{
#ifdef DEBUG_OUTPUT
		std::cout << "getLowerBound(...) - wrong mode input" << std::endl;
#endif
		return 0.0;
	}

	return getUpperBound(mode, -k, index);
}

bool MultilinearModelHandler::importMultilinearModelBinary(std::fstream& input)
{
	std::vector<size_t> modeDims;
	std::vector<size_t> truncModeDims;
	std::vector<double> multModel;
	std::vector<double> sVectors; 
	std::vector<double> mean;
	std::vector<double> meanWeights;

	FileLoader loader;
	if(!loader.loadMultilinearModelBinary(input, modeDims, truncModeDims, multModel, sVectors, mean, meanWeights))
	{
		return false;
	}


	clear();

	m_modeDimensions = modeDims;
	m_truncModeDimensions = truncModeDims;
	m_mean = mean;
	m_singularValues = sVectors;
	m_meanWeights = meanWeights;

	const size_t d1 = modeDims[0];
	const size_t m2 = truncModeDims[1];
	const size_t m3 = truncModeDims[2];

	m_tensor.init(multModel, d1, m2, m3);

	return true;
}

bool MultilinearModelHandler::importMultilinearModelText(std::fstream& input)
{
	std::vector<size_t> modeDims;
	std::vector<size_t> truncModeDims;
	std::vector<double> multModel;
	std::vector<double> sVectors; 
	std::vector<double> mean;
	std::vector<double> meanWeights;

	FileLoader loader;
	if(!loader.loadMultilinearModelText(input, modeDims, truncModeDims, multModel, sVectors, mean, meanWeights))
	{
		return false;
	}


	clear();

	m_modeDimensions = modeDims;
	m_truncModeDimensions = truncModeDims;
	m_mean = mean;
	m_singularValues = sVectors;
	m_meanWeights = meanWeights;

	const size_t d1 = modeDims[0];
	const size_t m2 = truncModeDims[1];
	const size_t m3 = truncModeDims[2];

	m_tensor.init(multModel, d1, m2, m3);

	return true;
}
