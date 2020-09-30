///////////////////////////////////////////////////////////////////////////////
//
//	WaveletShapeMultiLinearOptimizer.cpp
//
//	Source file for the CWaveletShapeMultiLinearOptimizer class
//
//	Alan Brunton, January 2013
//
///////////////////////////////////////////////////////////////////////////////


#include <vnl/vnl_cost_function.h>
#include <vnl/algo/vnl_lbfgsb.h>

#include <set>
#include "WaveletShapeMultiLinearOptimizer.h"
#include "MMProjectionCostFunction_Local.h"


void CWaveletShapeMultiLinearOptimizer::initActiveWaveletModel()
{
	size_t mode;

	std::cout << "\tallocate active wavelet models for each mode...\n"; std::cout.flush();
	m_activeWaveletModel.resize(m_nModes);
	for (mode = 0; mode < m_nModes; mode++)
	{
		m_activeWaveletModel[mode] = (CBSplineGridWavelet<C3Vectorf>*)CBSGWFactory::createBSGW(eBSGWFloat3, (int)m_nBaseWidth, (int)m_nBaseHeight, (int)m_nLevels);
	}

	m_pIntermediateModel = (CBSplineGridWavelet<C3Vectorf>*)CBSGWFactory::createBSGW(eBSGWFloat3, (int)m_nBaseWidth, (int)m_nBaseHeight, (int)m_nLevels);

	const size_t numCoeffs = (size_t)m_activeWaveletModel[0]->getNumCoefficients();

	m_pCurrentWeights = new std::vector<double>[numCoeffs];

	m_pVertexMask = new float[numCoeffs];
	for (size_t i = 0; i < numCoeffs; i++)
	{
		m_pVertexMask[i] = 1.f;
		m_pReconVerts[i].set(0.f, 0.f, 0.f);
	}

	for (size_t i = 0; i < numCoeffs; i++)
	{
		std::vector<size_t> truncModeDims;
		m_pMultilinearModels[i].getTruncModeDimensions(truncModeDims);
		const size_t m2 = truncModeDims[1];
		const size_t m3 = truncModeDims[2];
		const size_t numWeights = m2 + m3;

		std::vector<double> w2(m2), w3(m3), weights(numWeights, 0.0), coeff;

		m_pCurrentWeights[i].resize(numWeights, 0.0);

		m_pMultilinearModels[i].getModeMean(2, w2);
		m_pMultilinearModels[i].getModeMean(3, w3);
		for (size_t j = 0; j < m2; j++)
		{
			weights[j] = w2[j];
			m_pCurrentWeights[i][j] = w2[j];
		}
		for (size_t j = 0; j < m3; j++)
		{
			weights[j + m2] = w3[j];
			m_pCurrentWeights[i][j + m2] = w3[j];
		}

		m_pMultilinearModels[i].reconstructForWeights(weights, coeff);

		abutil::C3Vectorf vcoeff(coeff[0], coeff[1], coeff[2]);

		const size_t numIndices = m_pTransformIndices[i].size();
		for (size_t j = 0; j < numIndices; j++)
		{
			const size_t idx = m_pTransformIndices[i][j];
			m_pReconVerts[idx] += vcoeff * (float)m_pTransformCoeffs[i][j];
		}

		m_pIntermediateModel->setCoefficient(i, vcoeff);
	}

	m_fixMode.resize(2, false);
}

void CWaveletShapeMultiLinearOptimizer::clear()
{
	size_t mode;

	for (mode = 0; mode < m_nModes; mode++)
	{
		delete m_activeWaveletModel[mode];
		m_activeWaveletModel[mode] = NULL;
	}
	m_activeWaveletModel.clear();

	delete m_pIntermediateModel;
	m_pIntermediateModel = NULL;

	delete [] m_pCurrentWeights;
	m_pCurrentWeights = NULL;

	delete [] m_pMultilinearModels;
	m_pMultilinearModels = NULL;
	delete [] m_pTransformCoeffs;
	m_pTransformCoeffs = NULL;
	delete [] m_pTransformIndices;
	m_pTransformIndices = NULL;

	delete [] m_pReconVerts;
	m_pReconVerts = NULL;
	delete [] m_pTransformedRecon;
	m_pTransformedRecon = NULL;
	delete [] m_pNNDists;
	m_pNNDists = NULL;
	delete [] m_pNNIdx;
	m_pNNIdx = NULL;

	delete [] m_pObservedPoints;
	m_pObservedPoints = NULL;
	delete [] m_pObservedNormals;
	m_pObservedNormals = NULL;
	m_nObservedPoints = 0;

	delete [] m_pVertexMask;
	m_pVertexMask = NULL;

	delete [] m_pModelLandmarks;
	m_pModelLandmarks = NULL;
	delete [] m_pDataLandmarks;
	m_pDataLandmarks = NULL;
	delete [] m_pLandmarkVerts;
	m_pLandmarkVerts = NULL;
	m_nLandmarks = 0;

	m_nModes = 0;
}

void CWaveletShapeMultiLinearOptimizer::loadModels(std::string strModelFile)
{
	m_profLoadModel.start();

	size_t i;
	size_t modelBaseWidth, modelBaseHeight, modelLevels, modelWidth, modelHeight;

	clear();

	std::fstream modelfile;
	modelfile.open(strModelFile, std::ios_base::in | std::ios_base::binary);

	if (modelfile.bad() || !modelfile.is_open())
	{
		std::cout << "CWaveletShapeMultiLinearOptimizer::loadModels(" << strModelFile << ")\n";
		std::cout << "\tunable to open file wavelet.lmlm\n";
		std::cout.flush();
	}

	modelfile.read((char*)&modelBaseWidth, sizeof(size_t));
	if (modelfile.bad())
	{
		std::cout << "CWaveletShapeMultiLinearOptimizer::loadModels(" << strModelFile << ")\n";
		std::cout << "\terror reading model dimensions\n";
	}
	modelfile.read((char*)&modelBaseHeight, sizeof(size_t));
	if (modelfile.bad())
	{
		std::cout << "CWaveletShapeMultiLinearOptimizer::loadModels(" << strModelFile << ")\n";
		std::cout << "\terror reading model dimensions\n";
	}
	modelfile.read((char*)&modelLevels, sizeof(size_t));
	if (modelfile.bad())
	{
		std::cout << "CWaveletShapeMultiLinearOptimizer::loadModels(" << strModelFile << ")\n";
		std::cout << "\terror reading model dimensions\n";
	}
	modelfile.read((char*)&modelWidth, sizeof(size_t));
	if (modelfile.bad())
	{
		std::cout << "CWaveletShapeMultiLinearOptimizer::loadModels(" << strModelFile << ")\n";
		std::cout << "\terror reading model dimensions\n";
	}
	modelfile.read((char*)&modelHeight, sizeof(size_t));
	if (modelfile.bad())
	{
		std::cout << "CWaveletShapeMultiLinearOptimizer::loadModels(" << strModelFile << ")\n";
		std::cout << "\terror reading model dimensions\n";
	}

	size_t numModels = modelHeight * modelWidth;

	m_pMultilinearModels = new MultilinearModelHandler[numModels];

	m_pTransformCoeffs = new std::vector<double>[numModels];
	m_pTransformIndices = new std::vector<size_t>[numModels];

	for (i = 0; i < numModels; i++)
	{
		size_t numActive;
		modelfile.read((char*)&numActive, sizeof(size_t));
		if (modelfile.bad() || modelfile.eof())
		{
			std::cout << "CWaveletShapeMultiLinearOptimizer::loadModels(" << strModelFile << ")\n";
			std::cout << "\terror reading sparse transform data for model " << i << " from file\n";
			break;
		}

		double* pTransCoeffs = new double[numActive];
		size_t* pTransInds = new size_t[numActive];
		
		modelfile.read((char*)pTransCoeffs, sizeof(double) * numActive);
		if (modelfile.bad() || modelfile.eof())
		{
			std::cout << "CWaveletShapeMultiLinearOptimizer::loadModels(" << strModelFile << ")\n";
			std::cout << "\terror reading sparse transform data for model " << i << " from file\n";
			break;
		}

		modelfile.read((char*)pTransInds, sizeof(size_t) * numActive);
		if (modelfile.bad() || modelfile.eof())
		{
			std::cout << "CWaveletShapeMultiLinearOptimizer::loadModels(" << strModelFile << ")\n";
			std::cout << "\terror reading sparse transform data for model " << i << " from file\n";
			break;
		}

		for (size_t j = 0; j < numActive; j++)
		{
			m_pTransformCoeffs[i].push_back(pTransCoeffs[j]);
			m_pTransformIndices[i].push_back(pTransInds[j]);
		}
		m_pTransformCoeffs[i].shrink_to_fit();
		m_pTransformIndices[i].shrink_to_fit();
	}

	for (i = 0; i < numModels; i++)
	{
		if (!m_pMultilinearModels[i].importMultilinearModelBinary(modelfile))
		{
			std::cout << "CWaveletShapeMultiLinearOptimizer::loadModels(" << strModelFile << ")\n";
			std::cout << "\terror importing multilinear model " << i << " from file\n";
			break;
		}

		std::vector<size_t> modeDims;

		if (i == 0)
		{
			m_pMultilinearModels[i].getModeDimensions(modeDims);
			m_nModes = modeDims.size();
		}
		else
		{
			m_pMultilinearModels[i].getModeDimensions(modeDims);
			if (modeDims.size() != m_nModes)
			{
				std::cout << "CWaveletShapeMultiLinearOptimizer::loadModels(" << strModelFile << ")\n";
				std::cout << "\tmodel " << i << " has incorrect number of modes (" << modeDims.size() << ")\n";
			}
		}
	}

	m_nBaseWidth = modelBaseWidth;
	m_nBaseHeight = modelBaseHeight;
	m_nLevels = modelLevels;
	m_nFullWidth = modelWidth;
	m_nFullHeight = modelHeight;

	size_t numCoeffs = m_nFullHeight * m_nFullWidth;
	m_pReconVerts = new abutil::C3Vectorf[numCoeffs];
	m_pTransformedRecon = new abutil::C3Vectorf[numCoeffs];
	m_pNNDists = new NNAReal[numCoeffs];
	m_pNNIdx = new NNAIndex[numCoeffs];

	std::cout << "initializing active wavelet model...\n"; std::cout.flush();
	initActiveWaveletModel();

	m_profLoadModel.stop();
}

void CWaveletShapeMultiLinearOptimizer::loadModelsText(std::string strModelFile)
{
	m_profLoadModel.start();

	size_t i;
	size_t modelBaseWidth(0), modelBaseHeight(0), modelLevels(0), modelWidth(0), modelHeight(0);

	clear();

	std::fstream modelfile;
	modelfile.open(strModelFile, std::ios_base::in);

	if (modelfile.bad() || !modelfile.is_open())
	{
		std::cout << "CWaveletShapeMultiLinearOptimizer::loadModels(" << strModelFile << ")\n";
		std::cout << "\tunable to open file wavelet.lmm\n";
		std::cout.flush();
	}

	//modelfile.read((char*)&modelBaseWidth, sizeof(size_t));
	modelfile >> modelBaseWidth;
	if (modelfile.bad())
	{
		std::cout << "CWaveletShapeMultiLinearOptimizer::loadModels(" << strModelFile << ")\n";
		std::cout << "\terror reading model dimensions\n";
	}
	//modelfile.read((char*)&modelBaseHeight, sizeof(size_t));
	modelfile >> modelBaseHeight;
	if (modelfile.bad())
	{
		std::cout << "CWaveletShapeMultiLinearOptimizer::loadModels(" << strModelFile << ")\n";
		std::cout << "\terror reading model dimensions\n";
	}
	//modelfile.read((char*)&modelLevels, sizeof(size_t));
	modelfile >> modelLevels;
	if (modelfile.bad())
	{
		std::cout << "CWaveletShapeMultiLinearOptimizer::loadModels(" << strModelFile << ")\n";
		std::cout << "\terror reading model dimensions\n";
	}
	//modelfile.read((char*)&modelWidth, sizeof(size_t));
	modelfile >> modelWidth;
	if (modelfile.bad())
	{
		std::cout << "CWaveletShapeMultiLinearOptimizer::loadModels(" << strModelFile << ")\n";
		std::cout << "\terror reading model dimensions\n";
	}
	//modelfile.read((char*)&modelHeight, sizeof(size_t));
	modelfile >> modelHeight;
	if (modelfile.bad())
	{
		std::cout << "CWaveletShapeMultiLinearOptimizer::loadModels(" << strModelFile << ")\n";
		std::cout << "\terror reading model dimensions\n";
	}

//std::cout << "Num Models: " << modelHeight * modelWidth << std::endl;
//std::cout << "modelBaseWidth: " << modelBaseWidth << std::endl;
//std::cout << "modelBaseHeight: " << modelBaseHeight << std::endl;
//std::cout << "modelLevels: " << modelLevels << std::endl;
//std::cout << "modelWidth: " << modelWidth << std::endl;
//std::cout << "modelHeight: " << modelHeight << std::endl;

	size_t numModels = modelHeight * modelWidth;

	m_pMultilinearModels = new MultilinearModelHandler[numModels];

	m_pTransformCoeffs = new std::vector<double>[numModels];
	m_pTransformIndices = new std::vector<size_t>[numModels];

	for (i = 0; i < numModels; i++)
	{
		size_t numActive;
		//modelfile.read((char*)&numActive, sizeof(size_t));
		modelfile >> numActive;
		if (modelfile.bad() || modelfile.eof())
		{
			std::cout << "CWaveletShapeMultiLinearOptimizer::loadModels(" << strModelFile << ")\n";
			std::cout << "\terror reading sparse transform data for model " << i << " from file\n";
			break;
		}

		double* pTransCoeffs = new double[numActive];
		size_t* pTransInds = new size_t[numActive];
		
		//modelfile.read((char*)pTransCoeffs, sizeof(double) * numActive);
		for(size_t j = 0; j < numActive; ++j)
		{
			modelfile >> pTransCoeffs[j];
			if (modelfile.bad() || modelfile.eof())
			{
				std::cout << "CWaveletShapeMultiLinearOptimizer::loadModels(" << strModelFile << ")\n";
				std::cout << "\terror reading sparse transform data for model " << i << " from file\n";
				break;
			}
		}

		//modelfile.read((char*)pTransInds, sizeof(size_t) * numActive);
		for(size_t j = 0; j < numActive; ++j)
		{
			modelfile >> pTransInds[j];
			if (modelfile.bad() || modelfile.eof())
			{
				std::cout << "CWaveletShapeMultiLinearOptimizer::loadModels(" << strModelFile << ")\n";
				std::cout << "\terror reading sparse transform data for model " << i << " from file\n";
				break;
			}
		}

		for (size_t j = 0; j < numActive; j++)
		{
			m_pTransformCoeffs[i].push_back(pTransCoeffs[j]);
			m_pTransformIndices[i].push_back(pTransInds[j]);
		}
		m_pTransformCoeffs[i].shrink_to_fit();
		m_pTransformIndices[i].shrink_to_fit();
	}

	for (i = 0; i < numModels; i++)
	{
		if (!m_pMultilinearModels[i].importMultilinearModelText(modelfile))
		{
			std::cout << "CWaveletShapeMultiLinearOptimizer::loadModels(" << strModelFile << ")\n";
			std::cout << "\terror importing multilinear model " << i << " from file\n";
			break;
		}

		std::vector<size_t> modeDims;

		if (i == 0)
		{
			m_pMultilinearModels[i].getModeDimensions(modeDims);
			m_nModes = modeDims.size();
		}
		else
		{
			m_pMultilinearModels[i].getModeDimensions(modeDims);
			if (modeDims.size() != m_nModes)
			{
				std::cout << "CWaveletShapeMultiLinearOptimizer::loadModels(" << strModelFile << ")\n";
				std::cout << "\tmodel " << i << " has incorrect number of modes (" << modeDims.size() << ")\n";
			}
		}
	}

	m_nBaseWidth = modelBaseWidth;
	m_nBaseHeight = modelBaseHeight;
	m_nLevels = modelLevels;
	m_nFullWidth = modelWidth;
	m_nFullHeight = modelHeight;

	size_t numCoeffs = m_nFullHeight * m_nFullWidth;
	m_pReconVerts = new abutil::C3Vectorf[numCoeffs];
	m_pTransformedRecon = new abutil::C3Vectorf[numCoeffs];
	m_pNNDists = new NNAReal[numCoeffs];
	m_pNNIdx = new NNAIndex[numCoeffs];

	std::cout << "initializing active wavelet model...\n"; std::cout.flush();
	initActiveWaveletModel();

	modelfile.close();
	m_profLoadModel.stop();
}


void CWaveletShapeMultiLinearOptimizer::updateObservation(std::vector<double>& obsData, std::vector<double>& obsMask)
{
	m_profObsUpdate.start();

	const size_t numCoeffs = m_activeWaveletModel[0]->getFullResHeight() * m_activeWaveletModel[0]->getFullResWidth();
	size_t i, j;

	std::vector<double> obsNorm(numCoeffs * 3, 0.f);

	if (m_bFitLandmarksOnly)
	{
		if (m_level == 0)
		{
			float* pLMDists = new float[m_nLandmarks];

			for (j = 0; j < m_nLandmarks; j++)
			{
				pLMDists[j] = 1e30f;
				if (m_pLandmarkVerts[j] > numCoeffs)
				{
					m_pLandmarkVerts[j] = numCoeffs;
					for (i = 0; i < numCoeffs; i++)
					{
						if (m_pVertexMask[i] > 0.f)
						{
							const float distLM = (m_pReconVerts[i] - m_pModelLandmarks[j]).length();
							if (distLM < pLMDists[j])
							{
								pLMDists[j] = distLM;
								m_pLandmarkVerts[j] = i;
							}
						}
					}
				}
			}

			for (i = 0; i < numCoeffs; i++)
				obsMask[i] = 0.f;


			for (j = 0; j < m_nLandmarks; j++)
			{
				i = m_pLandmarkVerts[j];
				if (i < numCoeffs)
				{
					obsMask[i] = m_landmarkWeights[j];
					const abutil::C3Vectorf modelSpaceNN = m_matDataToModel * m_pDataLandmarks[j];
					obsData[i * 3 + 0] = modelSpaceNN.x;
					obsData[i * 3 + 1] = modelSpaceNN.y;
					obsData[i * 3 + 2] = modelSpaceNN.z;
				}
			}

			delete [] pLMDists;
		}
		else
		{
			for (i = 0; i < numCoeffs; i++)
				obsMask[i] = 0.f;

			for (j = 0; j < m_nLandmarks; j++)
			{
				i = m_pLandmarkVerts[j];
				if (i < numCoeffs)
				{
					obsMask[i] = m_landmarkWeights[j];
					const abutil::C3Vectorf modelSpaceNN = m_matDataToModel * m_pDataLandmarks[j];
					obsData[i * 3 + 0] = modelSpaceNN.x;
					obsData[i * 3 + 1] = modelSpaceNN.y;
					obsData[i * 3 + 2] = modelSpaceNN.z;
				}
			}
		}
	}
	else
	{
		for (i = 0; i < numCoeffs; i++)
			obsMask[i] = m_pVertexMask[i];

		for (i = 0; i < numCoeffs; i++)
		{
			m_pTransformedRecon[i] = m_matModelToData * m_pReconVerts[i];
		}

		m_pNNA->setQueryPoints(numCoeffs, (NNAReal*)m_pTransformedRecon, m_pNNDists, m_pNNIdx);
		m_pNNA->compute();

		if (m_pObservedNormals == NULL)
		{
			for (i = 0; i < numCoeffs; i++)
			{
				if (obsMask[i] == 0.0)	continue;
				j = i * 3;
				const abutil::C3Vectorf obsPoint(m_pObservedPoints + m_pNNIdx[i] * 3);
				const abutil::C3Vectorf modelSpaceNN = m_matDataToModel * obsPoint;
				obsData[j + 0] = modelSpaceNN.x;
				obsData[j + 1] = modelSpaceNN.y;
				obsData[j + 2] = modelSpaceNN.z;
			}
		}
		else
		{
			//we have normals, use point-to-plane distance
			for (i = 0; i < numCoeffs; i++)
			{
				if (obsMask[i] == 0.0)	continue;
				j = i * 3;
				const abutil::C3Vectorf obsPoint(m_pObservedPoints + m_pNNIdx[i] * 3);
				const abutil::C3Vectorf obsNormal = m_pObservedNormals[m_pNNIdx[i]];
				const abutil::C3Vectorf modelSpaceNN = m_matDataToModel * obsPoint;

				if (obsNormal.length() == 0.f)
				{
					obsData[j + 0] = modelSpaceNN.x;
					obsData[j + 1] = modelSpaceNN.y;
					obsData[j + 2] = modelSpaceNN.z;
					continue;
				}

				//do this to only rotate the normal, not translate
				const abutil::C4Vectorf modelSpaceNNNormH = m_matDataToModel * abutil::C4Vectorf(obsNormal.x, obsNormal.y, obsNormal.z, 0.f);
				const abutil::C3Vectorf modelSpaceNNNorm = abutil::C3Vectorf(modelSpaceNNNormH.x, modelSpaceNNNormH.y, modelSpaceNNNormH.z).normalize();
				obsNorm[j + 0] = obsNormal.x;
				obsNorm[j + 1] = obsNormal.y;
				obsNorm[j + 2] = obsNormal.z;

				//get plane offset rel. origin
				const float distToOrigin = modelSpaceNN.dot(modelSpaceNNNorm);

				//project model point onto plane
				const abutil::C3Vectorf projPoint = m_pReconVerts[i] - modelSpaceNNNorm * (m_pReconVerts[i].dot(modelSpaceNNNorm) - distToOrigin);
				
				//clip projection to be no more than m_inPlaneThresh away from nearest neighbor
				abutil::C3Vectorf vNNtoProj = projPoint - modelSpaceNN;
				const float inPlaneDist = vNNtoProj.length();
				vNNtoProj.normalize();

				const abutil::C3Vectorf clipProjPoint = modelSpaceNN + vNNtoProj * min(m_inPlaneThresh, inPlaneDist);

				obsData[j + 0] = clipProjPoint.x;
				obsData[j + 1] = clipProjPoint.y;
				obsData[j + 2] = clipProjPoint.z;

				//m_pNNDists[i] = (m_pReconVerts[i] - projPoint).length();
			}
		}

		float totalValidNN = 0.f;
		for (i = 0; i < numCoeffs; i++)
		{
			if (m_pNNDists[i] > m_nnthresh)
				obsMask[i] = 0.0;
			totalValidNN += obsMask[i];
		}

		float ratio = 1.f;
		if (m_nLandmarks > 0)
		{
			ratio = totalValidNN / (float)m_nLandmarks;//(float)m_nLandmarks / totalValidNN;//

			if(m_neighborSmoothWeight > 0.0)
			{
				ratio += (static_cast<float>(numCoeffs)/static_cast<float>(m_nLandmarks));
			}
		}

		for (j = 0; j < m_nLandmarks; j++)
		{
			i = m_pLandmarkVerts[j];
			if (i < numCoeffs)
			{
				obsMask[i] = m_landmarkWeights[j] * m_landmarkRelativeWeight * ratio;
				const abutil::C3Vectorf modelSpaceNN = m_matDataToModel * m_pDataLandmarks[j];
				obsData[i * 3 + 0] = modelSpaceNN.x;
				obsData[i * 3 + 1] = modelSpaceNN.y;
				obsData[i * 3 + 2] = modelSpaceNN.z;
			}
		}
	}

	m_profObsUpdate.stop();
}

void CWaveletShapeMultiLinearOptimizer::updateObservations(size_t& numObs, std::vector<double>& obsWeights, std::vector<double>& obsData, std::vector<double>& obsMask)
{
	m_profObsUpdate.start();

	const size_t numCoeffs = (size_t)m_pIntermediateModel->getNumCoefficients();
	const size_t numValues = numCoeffs * 3;

	const size_t numPreds = m_predictions.size();

	size_t curObs = 0;

	numObs = 0;

	if (m_pObservedPoints != NULL)
		numObs++;
	//if (m_pPrediction != NULL)
	//	numObs++;
	numObs += numPreds;

	obsWeights.resize(numObs, 0.0);
	obsData.resize(numObs * numValues, 0.0);
	obsMask.resize(numObs * numCoeffs, 0.0);

	if (m_pObservedPoints != NULL)
	{
		updateObservation(obsData, obsMask);
		obsWeights[curObs] = 1.0;
		curObs++;
	}

	//if (m_pPrediction != NULL)
	//{
	//	const size_t predOffset = curObs * numValues;
	//	const size_t maskOffset = curObs * numCoeffs;

	//	const double predWeight = (double)m_predWeight;

	//	for (size_t i = 0; i < numCoeffs; i++)
	//	{
	//		const size_t i3 = i * 3 + predOffset;
	//		const size_t j = i + maskOffset;

	//		obsMask[j]		= predWeight * (double)m_pVertexMask[i];
	//		obsData[i3 + 0]	= (double)m_pPrediction[i].x;
	//		obsData[i3 + 1]	= (double)m_pPrediction[i].y;
	//		obsData[i3 + 2]	= (double)m_pPrediction[i].z;
	//	}

	//	obsWeights[curObs] = predWeight;
	//	curObs++;
	//}

	for (size_t p = 0; p < numPreds; p++)
	{
		if (m_predictions[p] != NULL)
		{
			const size_t predOffset = curObs * numValues;
			const size_t maskOffset = curObs * numCoeffs;

			const double predWeight = (double)m_predWeights[p];

			for (size_t i = 0; i < numCoeffs; i++)
			{
				const size_t i3 = i * 3 + predOffset;
				const size_t j = i + maskOffset;

				obsMask[j]		= predWeight * (double)m_pVertexMask[i];
				obsData[i3 + 0]	= (double)m_predictions[p][i].x;
				obsData[i3 + 1]	= (double)m_predictions[p][i].y;
				obsData[i3 + 2]	= (double)m_predictions[p][i].z;
			}

			obsWeights[curObs] = predWeight;
			curObs++;
		}
	}

	m_profObsUpdate.stop();
}

void CWaveletShapeMultiLinearOptimizer::optimizeActiveWaveletModel()
{
//	_ASSERT(m_priors.size() > 0);

	if (m_level >= m_nLevels)
	{
		std::cout << "CWaveletShapeMultiLinearOptimizer::optimizeActiveWaveletModel(): invalid level, cannot optimize model.\n";
		return;
	}

	const size_t numCoeffs = m_activeWaveletModel[0]->getFullResHeight() * m_activeWaveletModel[0]->getFullResWidth();
	const size_t numCoeffsLevel = m_activeWaveletModel[0]->getLevelHeight(m_level) * m_activeWaveletModel[0]->getLevelWidth(m_level);
	const size_t numValues = numCoeffs * 3;

	size_t i, j;

	std::vector<double> fitData, fitMask, runningRecon;
	std::vector<double> obsWeights;

	size_t numObs;

	fitData.resize(numValues, 0.f);
	fitMask.resize(numCoeffs, 0.f);
	runningRecon.resize(numValues, 0.f);

	//std::cout << "\tupdating observation...\n"; std::cout.flush();
	//updateObservation(fitData, fitMask);
	updateObservations(numObs, obsWeights, fitData, fitMask);

	//std::cout << "\tinitializing reconstruction...\n"; std::cout.flush();
	for (i = 0; i < numCoeffs; i++)
	{
		const size_t i3 = i * 3;
		runningRecon[i3 + 0] = m_pReconVerts[i].x;
		runningRecon[i3 + 1] = m_pReconVerts[i].y;
		runningRecon[i3 + 2] = m_pReconVerts[i].z;
	}

	//std::cout << "\toptimizing coefficients " << m_iStartCoefficient << " through " << numCoeffsLevel - 1 << "...\n"; std::cout.flush();
	m_profLevelOpt.start();

	for (size_t iter = 0; iter < m_nInnerIters; iter++)
	{
		for (i = m_iStartCoefficient; i < numCoeffsLevel; i++)
		{
			//m_profCoeffOpt.start();

			const int icoeff = m_activeWaveletModel[0]->getSerialIndex((int)i);

			//std::cout << "\tcoefficient " << icoeff << "\n"; std::cout.flush();
			std::vector<size_t> truncModeDims;

			m_pMultilinearModels[icoeff].getTruncModeDimensions(truncModeDims);
			const size_t m2 = truncModeDims[1];
			const size_t m3 = truncModeDims[2];
			const size_t numWeights = m2 + m3;

			vnl_vector<long> boundSelection(numWeights, 2);
			vnl_vector<double> lowerBounds(numWeights, 0.0);
			vnl_vector<double> upperBounds(numWeights, 0.0);
			vnl_vector<double> x(numWeights, 0.0);

			std::vector<double> w2(numWeights), w3(numWeights);

			const bool bM2Fix = m_fixMode[0];
			const bool bM3Fix = m_fixMode[1];
			const size_t offsetM3 = bM2Fix ? 0 : m2;

			for (j = 0; j < numWeights; j++)
				x[j] = m_pCurrentWeights[icoeff][j];
			for (j = 0; j < m2; j++)
				w2[j] = m_pCurrentWeights[icoeff][j];
			for (j = 0; j < m3; j++)
				w3[j] = m_pCurrentWeights[icoeff][j + m2];

			const size_t numIndices = m_pTransformIndices[icoeff].size();

			std::vector<double> initCoeff;
			m_pMultilinearModels[icoeff].reconstructForWeights(m_pCurrentWeights[icoeff], initCoeff);
			for (j = 0; j < numIndices; j++)
			{
				const size_t idx = m_pTransformIndices[icoeff][j];
				const size_t k = idx * 3;
				const double tc = m_pTransformCoeffs[icoeff][j];
				runningRecon[k + 0]		-= tc * initCoeff[0];
				runningRecon[k + 1]		-= tc * initCoeff[1];
				runningRecon[k + 2]		-= tc * initCoeff[2];
			}

			if (bM2Fix)
			{
				x.set_size(m3);
				for (j = 0; j < m3; j++)
					x[j] = m_pCurrentWeights[icoeff][j + m2];
			}
			else if (bM3Fix)
			{
				x.set_size(m2);
				for (j = 0; j < m2; j++)
					x[j] = m_pCurrentWeights[icoeff][j];
			}

			if (!bM2Fix)
			{
				for (j = 0; j < m2; j++)
				{
					const double lowerBound = m_pMultilinearModels[icoeff].getLowerBound(2, (double)m_allowDeviation, j);
					const double upperBound = m_pMultilinearModels[icoeff].getUpperBound(2, (double)m_allowDeviation, j);		

					//lowerBounds[j] = lowerBound;
					//upperBounds[j] = upperBound;
					lowerBounds[j] = min(lowerBound, m_pCurrentWeights[icoeff][j]);
					upperBounds[j] = max(upperBound, m_pCurrentWeights[icoeff][j]);
				}
			}

			if (!bM3Fix)
			{
				for (j = 0; j < m3; j++)
				{
					const size_t idx = j + offsetM3;
					const size_t idx2 = j + m2;

					const double lowerBound = m_pMultilinearModels[icoeff].getLowerBound(3, (double)m_allowDeviation, j);
					const double upperBound = m_pMultilinearModels[icoeff].getUpperBound(3, (double)m_allowDeviation, j);		

					//lowerBounds[idx] = lowerBound;
					//upperBounds[idx] = upperBound;
					lowerBounds[idx] = min(lowerBound, m_pCurrentWeights[icoeff][idx2]);
					upperBounds[idx] = max(upperBound, m_pCurrentWeights[icoeff][idx2]);
				}
			}

			std::vector<double> mean;
			m_pMultilinearModels[icoeff].getMean(mean);
			//MMProjectionCostFunction_Local fkt(&(m_pMultilinearModels[icoeff].getTensor()), mean, runningRecon, fitData, fitMask, m_pTransformCoeffs[icoeff], m_pTransformIndices[icoeff]);
			//MMProjectionCostFunction_Local fkt(&(m_pMultilinearModels[icoeff].getTensor()), mean, runningRecon, numObs, numCoeffs, fitData, fitMask, m_pTransformCoeffs[icoeff], m_pTransformIndices[icoeff]);

			MMProjectionCostFunction_Local fkt(&(m_pMultilinearModels[icoeff].getTensor()), mean, runningRecon, numObs, numCoeffs, fitData, fitMask, m_pTransformCoeffs[icoeff], m_pTransformIndices[icoeff], m_neighborSmoothWeight
															, m_precomputedSmoothIndices, m_precomputedSmoothWeights);


			if (bM2Fix)
				fkt.fixMode(2, w2);
			else if (bM3Fix)
				fkt.fixMode(3, w3);

			//std::cout << "\toptimizing...\n"; std::cout.flush();
			vnl_lbfgsb minimizer(fkt);
			minimizer.set_cost_function_convergence_factor(10000000);
			//minimizer.set_projected_gradient_tolerance(0.000001);
			minimizer.set_projected_gradient_tolerance(0.00000001);
			minimizer.set_max_function_evals(1000);
			//minimizer.set_max_function_evals(10000);
			minimizer.set_bound_selection(boundSelection);
			minimizer.set_lower_bound(lowerBounds);
			minimizer.set_upper_bound(upperBounds);
			//minimizer.set_trace(true);
			minimizer.minimize(x);

			if (bM2Fix)
			{
				for (j = 0; j < m3; j++)
					m_pCurrentWeights[icoeff][j + m2] = x[j];
			}
			else if (bM3Fix)
			{
				for (j = 0; j < m2; j++)
					m_pCurrentWeights[icoeff][j] = x[j];
			}
			else
			{
				for (j = 0; j < numWeights; j++)
					m_pCurrentWeights[icoeff][j] = x[j];
			}

			//std::cout << "\tsaving result...\n"; std::cout.flush();
			std::vector<double> resCoeff;
			m_pMultilinearModels[icoeff].reconstructForWeights(m_pCurrentWeights[icoeff], resCoeff);

			//std::cout << "\tupdating recon...\n"; std::cout.flush();
			for (j = 0; j < numIndices; j++)
			{
				const size_t idx = m_pTransformIndices[icoeff][j];
				const size_t k = idx * 3;
				const double tc = m_pTransformCoeffs[icoeff][j];

				runningRecon[k + 0]		+= tc * resCoeff[0];
				runningRecon[k + 1]		+= tc * resCoeff[1];
				runningRecon[k + 2]		+= tc * resCoeff[2];
				m_pReconVerts[idx].x	= (float)runningRecon[k + 0];
				m_pReconVerts[idx].y	= (float)runningRecon[k + 1];
				m_pReconVerts[idx].z	= (float)runningRecon[k + 2];
			}

			m_pIntermediateModel->setCoefficient(icoeff, abutil::C3Vectorf(resCoeff[0], resCoeff[1], resCoeff[2]));

			//m_activeWaveletModel[0]->setCoefficient(icoeff, abutil::C3Vectorf(x[0], x[1], x[2]));
			//m_activeWaveletModel[1]->setCoefficient(icoeff, abutil::C3Vectorf(x[3], x[4], x[5]));

			//m_profCoeffOpt.stop();
		}
	}

	m_iStartCoefficient = numCoeffsLevel;
	m_level++;

	m_profLevelOpt.stop();
}

void CWaveletShapeMultiLinearOptimizer::refineSurface()
{
	refineSurface(m_pReconVerts);
}

void CWaveletShapeMultiLinearOptimizer::refineSurface(abutil::C3Vectorf* pRefinedVerts)
{
	const size_t numCoeffs = (size_t)m_pIntermediateModel->getNumCoefficients();
	const size_t numValues = numCoeffs * 3;

	std::vector<double> fitData, fitMask, totalWeight, totalDataWeight, deformation, defNext, defData, defSmooth;

	fitData.resize(numValues, 0.0);
	fitMask.resize(numCoeffs, 0.0);
	deformation.resize(numValues, 0.0);
	defNext.resize(numValues, 0.0);
	defData.resize(numValues, 0.0);
	defSmooth.resize(numValues, 0.0);
	totalWeight.resize(numCoeffs, 0.0);
	totalDataWeight.resize(numCoeffs, 0.0);

	//std::cout << "\tupdating observation...\n"; std::cout.flush();
	updateObservation(fitData, fitMask);

	m_profRefine.start();

	const size_t fhm1 = m_nFullHeight - 1;
	const size_t fwm1 = m_nFullWidth - 1;

	const int innerOffset[8] = { -m_nFullWidth - 1, -m_nFullWidth, -m_nFullWidth + 1, -1, 1, m_nFullWidth - 1, m_nFullWidth, m_nFullWidth + 1 };

	size_t x, y;

	const double smoothWeight = (double)m_smoothRelativeWeight;
	const double geomSmoothWeight = 0.25 * smoothWeight;

	///////////////////////
	//compute weight totals

	for (y = 0; y < m_nFullHeight; y++)
	{
		for (x = 0; x < m_nFullWidth; x++)
		{
			const int i = (int)(y * m_nFullWidth + x);
			totalDataWeight[i] = fitMask[i];
		}
	}

	//upper left neighbor
	for (y = 1; y < m_nFullHeight; y++)
	{
		for (x = 1; x < m_nFullWidth; x++)
		{
			const int i = (int)(y * m_nFullWidth + x);
			const int i3 = i * 3;
			const int in = i + innerOffset[0];
			const double w = (double)m_pVertexMask[in];
			totalWeight[i] += w;
			totalDataWeight[i] += fitMask[in];
			defSmooth[i3 + 0] += w * (double)m_pReconVerts[in].x;
			defSmooth[i3 + 1] += w * (double)m_pReconVerts[in].y;
			defSmooth[i3 + 2] += w * (double)m_pReconVerts[in].z;
		}
	}

	//top neighbor
	for (y = 1; y < m_nFullHeight; y++)
	{
		for (x = 0; x < m_nFullWidth; x++)
		{
			const int i = (int)(y * m_nFullWidth + x);
			const int i3 = i * 3;
			const int in = i + innerOffset[1];
			const double w = (double)m_pVertexMask[in];
			totalWeight[i] += w;
			totalDataWeight[i] += fitMask[in];
			defSmooth[i3 + 0] += w * (double)m_pReconVerts[in].x;
			defSmooth[i3 + 1] += w * (double)m_pReconVerts[in].y;
			defSmooth[i3 + 2] += w * (double)m_pReconVerts[in].z;
		}
	}

	//upper right neighbor
	for (y = 1; y < m_nFullHeight; y++)
	{
		for (x = 0; x < fwm1; x++)
		{
			const int i = (int)(y * m_nFullWidth + x);
			const int i3 = i * 3;
			const int in = i + innerOffset[2];
			const double w = (double)m_pVertexMask[in];
			totalWeight[i] += w;
			totalDataWeight[i] += fitMask[in];
			defSmooth[i3 + 0] += w * (double)m_pReconVerts[in].x;
			defSmooth[i3 + 1] += w * (double)m_pReconVerts[in].y;
			defSmooth[i3 + 2] += w * (double)m_pReconVerts[in].z;
		}
	}

	//left neighbor
	for (y = 0; y < m_nFullHeight; y++)
	{
		for (x = 1; x < m_nFullWidth; x++)
		{
			const int i = (int)(y * m_nFullWidth + x);
			const int i3 = i * 3;
			const int in = i + innerOffset[3];
			const double w = (double)m_pVertexMask[in];
			totalWeight[i] += w;
			totalDataWeight[i] += fitMask[in];
			defSmooth[i3 + 0] += w * (double)m_pReconVerts[in].x;
			defSmooth[i3 + 1] += w * (double)m_pReconVerts[in].y;
			defSmooth[i3 + 2] += w * (double)m_pReconVerts[in].z;
		}
	}

	//right neighbor
	for (y = 0; y < m_nFullHeight; y++)
	{
		for (x = 0; x < fwm1; x++)
		{
			const int i = (int)(y * m_nFullWidth + x);
			const int i3 = i * 3;
			const int in = i + innerOffset[4];
			const double w = (double)m_pVertexMask[in];
			totalWeight[i] += w;
			totalDataWeight[i] += fitMask[in];
			defSmooth[i3 + 0] += w * (double)m_pReconVerts[in].x;
			defSmooth[i3 + 1] += w * (double)m_pReconVerts[in].y;
			defSmooth[i3 + 2] += w * (double)m_pReconVerts[in].z;
		}
	}

	//lower left neighbor
	for (y = 0; y < fhm1; y++)
	{
		for (x = 1; x < m_nFullWidth; x++)
		{
			const int i = (int)(y * m_nFullWidth + x);
			const int i3 = i * 3;
			const int in = i + innerOffset[5];
			const double w = (double)m_pVertexMask[in];
			totalWeight[i] += w;
			totalDataWeight[i] += fitMask[in];
			defSmooth[i3 + 0] += w * (double)m_pReconVerts[in].x;
			defSmooth[i3 + 1] += w * (double)m_pReconVerts[in].y;
			defSmooth[i3 + 2] += w * (double)m_pReconVerts[in].z;
		}
	}

	//bottom neighbor
	for (y = 0; y < fhm1; y++)
	{
		for (x = 0; x < m_nFullWidth; x++)
		{
			const int i = (int)(y * m_nFullWidth + x);
			const int i3 = i * 3;
			const int in = i + innerOffset[6];
			const double w = (double)m_pVertexMask[in];
			totalWeight[i] += w;
			totalDataWeight[i] += fitMask[in];
			defSmooth[i3 + 0] += w * (double)m_pReconVerts[in].x;
			defSmooth[i3 + 1] += w * (double)m_pReconVerts[in].y;
			defSmooth[i3 + 2] += w * (double)m_pReconVerts[in].z;
		}
	}

	//lower right neighbor
	for (y = 0; y < fhm1; y++)
	{
		for (x = 0; x < fwm1; x++)
		{
			const int i = (int)(y * m_nFullWidth + x);
			const int i3 = i * 3;
			const int in = i + innerOffset[7];
			const double w = (double)m_pVertexMask[in];
			totalWeight[i] += w;
			totalDataWeight[i] += fitMask[in];
			defSmooth[i3 + 0] += w * (double)m_pReconVerts[in].x;
			defSmooth[i3 + 1] += w * (double)m_pReconVerts[in].y;
			defSmooth[i3 + 2] += w * (double)m_pReconVerts[in].z;
		}
	}

	for (y = 0; y < m_nFullHeight; y++)
	{
		for (x = 0; x < m_nFullWidth; x++)
		{
			const int i = (int)(y * m_nFullWidth + x);

			if (totalWeight[i] > 0.0)
				totalWeight[i] = 1.f / totalWeight[i];
			else
				totalWeight[i] = 0.0;

			if (totalDataWeight[i] > 0.0)
				totalDataWeight[i] = 1.f / totalDataWeight[i];
			else
				totalDataWeight[i] = 0.0;
		}
	}

	for (y = 0; y < m_nFullHeight; y++)
	{
		for (x = 0; x < m_nFullWidth; x++)
		{
			const int i = (int)(y * m_nFullWidth + x);
			const int i3 = i * 3;

			//deformation to data
			defData[i3 + 0] = fitMask[i] * (fitData[i3 + 0] - (double)m_pReconVerts[i].x);
			defData[i3 + 1] = fitMask[i] * (fitData[i3 + 1] - (double)m_pReconVerts[i].y);
			defData[i3 + 2] = fitMask[i] * (fitData[i3 + 2] - (double)m_pReconVerts[i].z);

			//deformation to neighbors
			defSmooth[i3 + 0] *= totalWeight[i];
			defSmooth[i3 + 1] *= totalWeight[i];
			defSmooth[i3 + 2] *= totalWeight[i];
			defSmooth[i3 + 0] -= (double)m_pReconVerts[i].x;
			defSmooth[i3 + 1] -= (double)m_pReconVerts[i].y;
			defSmooth[i3 + 2] -= (double)m_pReconVerts[i].z;
		}
	}

	//iteratively optimize deformation field
	for (size_t iter = 0; iter < m_nRefineInnerIters; iter++)
	{
		//////////////////////////////
		//compute updated deformations

		//initialize to deformation to zero
		for (size_t v = 0; v < numValues; v++)
			defNext[v] = 0.0;

		//upper left neighbor
		for (y = 1; y < m_nFullHeight; y++)
		{
			for (x = 1; x < m_nFullWidth; x++)
			{
				const int i = (int)(y * m_nFullWidth + x);
				const int i3 = i * 3;
				const int in = i + innerOffset[0];
				const int in3 = in * 3;

				//deformation smoothing
				defNext[i3 + 0] += smoothWeight * totalDataWeight[i] * fitMask[in] * deformation[in3 + 0];
				defNext[i3 + 1] += smoothWeight * totalDataWeight[i] * fitMask[in] * deformation[in3 + 1];
				defNext[i3 + 2] += smoothWeight * totalDataWeight[i] * fitMask[in] * deformation[in3 + 2];

				//geometric smoothing
				defNext[i3 + 0] += geomSmoothWeight * totalWeight[i] * (double)m_pVertexMask[in] * deformation[in3 + 0];
				defNext[i3 + 1] += geomSmoothWeight * totalWeight[i] * (double)m_pVertexMask[in] * deformation[in3 + 1];
				defNext[i3 + 2] += geomSmoothWeight * totalWeight[i] * (double)m_pVertexMask[in] * deformation[in3 + 2];
			}
		}

		//top neighbor
		for (y = 1; y < m_nFullHeight; y++)
		{
			for (x = 0; x < m_nFullWidth; x++)
			{
				const int i = (int)(y * m_nFullWidth + x);
				const int i3 = i * 3;
				const int in = i + innerOffset[1];
				const int in3 = in * 3;

				//deformation smoothing
				defNext[i3 + 0] += smoothWeight * totalDataWeight[i] * fitMask[in] * deformation[in3 + 0];
				defNext[i3 + 1] += smoothWeight * totalDataWeight[i] * fitMask[in] * deformation[in3 + 1];
				defNext[i3 + 2] += smoothWeight * totalDataWeight[i] * fitMask[in] * deformation[in3 + 2];

				//geometric smoothing
				defNext[i3 + 0] += geomSmoothWeight * totalWeight[i] * (double)m_pVertexMask[in] * deformation[in3 + 0];
				defNext[i3 + 1] += geomSmoothWeight * totalWeight[i] * (double)m_pVertexMask[in] * deformation[in3 + 1];
				defNext[i3 + 2] += geomSmoothWeight * totalWeight[i] * (double)m_pVertexMask[in] * deformation[in3 + 2];
			}
		}

		//upper right neighbor
		for (y = 1; y < m_nFullHeight; y++)
		{
			for (x = 0; x < fwm1; x++)
			{
				const int i = (int)(y * m_nFullWidth + x);
				const int i3 = i * 3;
				const int in = i + innerOffset[2];
				const int in3 = in * 3;

				//deformation smoothing
				defNext[i3 + 0] += smoothWeight * totalDataWeight[i] * fitMask[in] * deformation[in3 + 0];
				defNext[i3 + 1] += smoothWeight * totalDataWeight[i] * fitMask[in] * deformation[in3 + 1];
				defNext[i3 + 2] += smoothWeight * totalDataWeight[i] * fitMask[in] * deformation[in3 + 2];

				//geometric smoothing
				defNext[i3 + 0] += geomSmoothWeight * totalWeight[i] * (double)m_pVertexMask[in] * deformation[in3 + 0];
				defNext[i3 + 1] += geomSmoothWeight * totalWeight[i] * (double)m_pVertexMask[in] * deformation[in3 + 1];
				defNext[i3 + 2] += geomSmoothWeight * totalWeight[i] * (double)m_pVertexMask[in] * deformation[in3 + 2];
			}
		}

		//left neighbor
		for (y = 0; y < m_nFullHeight; y++)
		{
			for (x = 1; x < m_nFullWidth; x++)
			{
				const int i = (int)(y * m_nFullWidth + x);
				const int i3 = i * 3;
				const int in = i + innerOffset[3];
				const int in3 = in * 3;

				//deformation smoothing
				defNext[i3 + 0] += smoothWeight * totalDataWeight[i] * fitMask[in] * deformation[in3 + 0];
				defNext[i3 + 1] += smoothWeight * totalDataWeight[i] * fitMask[in] * deformation[in3 + 1];
				defNext[i3 + 2] += smoothWeight * totalDataWeight[i] * fitMask[in] * deformation[in3 + 2];

				//geometric smoothing
				defNext[i3 + 0] += geomSmoothWeight * totalWeight[i] * (double)m_pVertexMask[in] * deformation[in3 + 0];
				defNext[i3 + 1] += geomSmoothWeight * totalWeight[i] * (double)m_pVertexMask[in] * deformation[in3 + 1];
				defNext[i3 + 2] += geomSmoothWeight * totalWeight[i] * (double)m_pVertexMask[in] * deformation[in3 + 2];
			}
		}

		//right neighbor
		for (y = 0; y < m_nFullHeight; y++)
		{
			for (x = 0; x < fwm1; x++)
			{
				const int i = (int)(y * m_nFullWidth + x);
				const int i3 = i * 3;
				const int in = i + innerOffset[4];
				const int in3 = in * 3;

				//deformation smoothing
				defNext[i3 + 0] += smoothWeight * totalDataWeight[i] * fitMask[in] * deformation[in3 + 0];
				defNext[i3 + 1] += smoothWeight * totalDataWeight[i] * fitMask[in] * deformation[in3 + 1];
				defNext[i3 + 2] += smoothWeight * totalDataWeight[i] * fitMask[in] * deformation[in3 + 2];

				//geometric smoothing
				defNext[i3 + 0] += geomSmoothWeight * totalWeight[i] * (double)m_pVertexMask[in] * deformation[in3 + 0];
				defNext[i3 + 1] += geomSmoothWeight * totalWeight[i] * (double)m_pVertexMask[in] * deformation[in3 + 1];
				defNext[i3 + 2] += geomSmoothWeight * totalWeight[i] * (double)m_pVertexMask[in] * deformation[in3 + 2];
			}
		}

		//lower left neighbor
		for (y = 0; y < fhm1; y++)
		{
			for (x = 1; x < m_nFullWidth; x++)
			{
				const int i = (int)(y * m_nFullWidth + x);
				const int i3 = i * 3;
				const int in = i + innerOffset[5];
				const int in3 = in * 3;

				//deformation smoothing
				defNext[i3 + 0] += smoothWeight * totalDataWeight[i] * fitMask[in] * deformation[in3 + 0];
				defNext[i3 + 1] += smoothWeight * totalDataWeight[i] * fitMask[in] * deformation[in3 + 1];
				defNext[i3 + 2] += smoothWeight * totalDataWeight[i] * fitMask[in] * deformation[in3 + 2];

				//geometric smoothing
				defNext[i3 + 0] += geomSmoothWeight * totalWeight[i] * (double)m_pVertexMask[in] * deformation[in3 + 0];
				defNext[i3 + 1] += geomSmoothWeight * totalWeight[i] * (double)m_pVertexMask[in] * deformation[in3 + 1];
				defNext[i3 + 2] += geomSmoothWeight * totalWeight[i] * (double)m_pVertexMask[in] * deformation[in3 + 2];
			}
		}

		//bottom neighbor
		for (y = 0; y < fhm1; y++)
		{
			for (x = 0; x < m_nFullWidth; x++)
			{
				const int i = (int)(y * m_nFullWidth + x);
				const int i3 = i * 3;
				const int in = i + innerOffset[6];
				const int in3 = in * 3;

				//deformation smoothing
				defNext[i3 + 0] += smoothWeight * totalDataWeight[i] * fitMask[in] * deformation[in3 + 0];
				defNext[i3 + 1] += smoothWeight * totalDataWeight[i] * fitMask[in] * deformation[in3 + 1];
				defNext[i3 + 2] += smoothWeight * totalDataWeight[i] * fitMask[in] * deformation[in3 + 2];

				//geometric smoothing
				defNext[i3 + 0] += geomSmoothWeight * totalWeight[i] * (double)m_pVertexMask[in] * deformation[in3 + 0];
				defNext[i3 + 1] += geomSmoothWeight * totalWeight[i] * (double)m_pVertexMask[in] * deformation[in3 + 1];
				defNext[i3 + 2] += geomSmoothWeight * totalWeight[i] * (double)m_pVertexMask[in] * deformation[in3 + 2];
			}
		}

		//lower right neighbor
		for (y = 0; y < fhm1; y++)
		{
			for (x = 0; x < fwm1; x++)
			{
				const int i = (int)(y * m_nFullWidth + x);
				const int i3 = i * 3;
				const int in = i + innerOffset[7];
				const int in3 = in * 3;

				//deformation smoothing
				defNext[i3 + 0] += smoothWeight * totalDataWeight[i] * fitMask[in] * deformation[in3 + 0];
				defNext[i3 + 1] += smoothWeight * totalDataWeight[i] * fitMask[in] * deformation[in3 + 1];
				defNext[i3 + 2] += smoothWeight * totalDataWeight[i] * fitMask[in] * deformation[in3 + 2];

				//geometric smoothing
				defNext[i3 + 0] += geomSmoothWeight * totalWeight[i] * (double)m_pVertexMask[in] * deformation[in3 + 0];
				defNext[i3 + 1] += geomSmoothWeight * totalWeight[i] * (double)m_pVertexMask[in] * deformation[in3 + 1];
				defNext[i3 + 2] += geomSmoothWeight * totalWeight[i] * (double)m_pVertexMask[in] * deformation[in3 + 2];
			}
		}

		//offset by constant deformations (data + geometric smoothing terms)
		//for (size_t v = 0; v < numValues; v++)
		//	defNext[v] += defData[v] + geomSmoothWeight * defSmooth[v];
		for (size_t c = 0; c < numCoeffs; c++)
		{
			const size_t c3 = c * 3;
			const double normFactor = 1.0 / (smoothWeight + geomSmoothWeight + fitMask[c]);
			defNext[c3 + 0] += defData[c3 + 0] + geomSmoothWeight * defSmooth[c3 + 0];
			defNext[c3 + 1] += defData[c3 + 1] + geomSmoothWeight * defSmooth[c3 + 1];
			defNext[c3 + 2] += defData[c3 + 2] + geomSmoothWeight * defSmooth[c3 + 2];
			defNext[c3 + 0] *= normFactor;
			defNext[c3 + 1] *= normFactor;
			defNext[c3 + 2] *= normFactor;
		}

		//update deformation
		for (size_t v = 0; v < numValues; v++)
			deformation[v] = defNext[v];
	}

	//now the deformation field to the geometry
	for (y = 0; y < m_nFullHeight; y++)
	{
		for (x = 0; x < m_nFullWidth; x++)
		{
			const int i = (int)(y * m_nFullWidth + x);
			const int i3 = i * 3;
			pRefinedVerts[i].x = m_pReconVerts[i].x + (float)deformation[i3 + 0];
			pRefinedVerts[i].y = m_pReconVerts[i].y + (float)deformation[i3 + 1];
			pRefinedVerts[i].z = m_pReconVerts[i].z + (float)deformation[i3 + 2];
		}
	}

	m_profRefine.stop();
}

void CWaveletShapeMultiLinearOptimizer::highlightCoefficient(const abutil::C3Vectorf& baseColor, const abutil::C3Vectorf& highlightColor, const size_t icoeff, abutil::C3Vectorf* pResult)
{
	if (m_pIntermediateModel == NULL)
		return;

	const size_t numCoeffs = m_pIntermediateModel->getNumCoefficients();

	for (size_t c = 0; c < numCoeffs; c++)
		pResult[c] = baseColor;

	const size_t numIndices = m_pTransformIndices[icoeff].size();
	for (size_t i = 0; i < numIndices; i++)
	{
		const size_t ic = m_pTransformIndices[icoeff][i];
		const double tc = m_pTransformCoeffs[icoeff][i];
		const float t = (float)fabs(tc);
		pResult[ic] = highlightColor;// * t + baseColor * (1.f - t);
	}
}

void CWaveletShapeMultiLinearOptimizer::reportProfiling(FILE* pfOutput)
{
	m_profLoadModel.computeStatistics();
	m_profLevelOpt.computeStatistics();
	m_profCoeffOpt.computeStatistics();
	m_profObsUpdate.computeStatistics();
	m_profRefine.computeStatistics();
	m_profLoadModel.reportStatistics(pfOutput);
	m_profLevelOpt.reportStatistics(pfOutput);
	m_profCoeffOpt.reportStatistics(pfOutput);
	m_profObsUpdate.reportStatistics(pfOutput);
	m_profRefine.reportStatistics(pfOutput);
}

void CWaveletShapeMultiLinearOptimizer::computeGridNeighbors(std::vector<std::vector<size_t>>& neighbors)
{
	neighbors.resize(m_nFullHeight*m_nFullWidth);
	m_invalidSmoothVertices.resize(m_nFullHeight*m_nFullWidth, false);

	float maskThreshold(0.5);

	//Collect all indices of vertices that have a boundary in the 1-ring-neighborhood
	std::set<size_t> oneRingBoundaryNeighbors;

	for(size_t y = 0; y < m_nFullHeight; ++y)
	{
		for(size_t x = 0; x < m_nFullWidth; ++x)
		{
			const size_t i = y*m_nFullWidth+x;
			if(m_pVertexMask[i] <= maskThreshold)
			{
				continue;
			}		

			if(x == 0 || x == m_nFullWidth-1 || y == 0 || y == m_nFullHeight-1)
			{
				oneRingBoundaryNeighbors.insert(i);
				continue;
			}

			//Check if one of the neighbors is invalid
			if(m_pVertexMask[i-1-m_nFullWidth] <= maskThreshold
				|| m_pVertexMask[i-m_nFullWidth] <= maskThreshold
				|| m_pVertexMask[i-m_nFullWidth+1] <= maskThreshold
				|| m_pVertexMask[i-1] <= maskThreshold
				|| m_pVertexMask[i+1] <= maskThreshold
				|| m_pVertexMask[i-1+m_nFullWidth] <= maskThreshold
				|| m_pVertexMask[i+m_nFullWidth] <= maskThreshold
				|| m_pVertexMask[i+m_nFullWidth+1] <= maskThreshold)
			{
				oneRingBoundaryNeighbors.insert(i);
			}
		}
	}

	//Collect all indices of vertices that have a boundary in the 2-ring-neighborhood
	std::set<size_t> twoRingBoundaryNeighbors;

	for(size_t y = 1; y < m_nFullHeight-1; ++y)
	{
		for(size_t x = 1; x < m_nFullWidth-1; ++x)
		{
			const size_t i = y*m_nFullWidth+x;
			if(m_pVertexMask[i] <= maskThreshold)
			{
				continue;
			}		

			//Check if one of the neighbors has an invalid neighbor
			if(oneRingBoundaryNeighbors.find(i-1-m_nFullWidth) != oneRingBoundaryNeighbors.end()
				|| oneRingBoundaryNeighbors.find(i-m_nFullWidth) != oneRingBoundaryNeighbors.end()
				|| oneRingBoundaryNeighbors.find(i-m_nFullWidth+1) != oneRingBoundaryNeighbors.end()
				|| oneRingBoundaryNeighbors.find(i-1) != oneRingBoundaryNeighbors.end()
				|| oneRingBoundaryNeighbors.find(i+1) != oneRingBoundaryNeighbors.end()
				|| oneRingBoundaryNeighbors.find(i-1+m_nFullWidth) != oneRingBoundaryNeighbors.end()
				|| oneRingBoundaryNeighbors.find(i+m_nFullWidth) != oneRingBoundaryNeighbors.end()
				|| oneRingBoundaryNeighbors.find(i+m_nFullWidth+1) != oneRingBoundaryNeighbors.end())
			{
				twoRingBoundaryNeighbors.insert(i);
			}		
		}
	}

	for(size_t y = 0; y < m_nFullHeight; ++y)
	{
		for(size_t x = 0; x < m_nFullWidth; ++x)
		{
			//Current index
			const size_t i = y*m_nFullWidth+x;
			if(m_pVertexMask[i] <= maskThreshold)
			{
				//Point has mask <= 0.5 and therefore is not valid.
				m_invalidSmoothVertices[i] = true;
				continue;
			}

			if(oneRingBoundaryNeighbors.find(i) != oneRingBoundaryNeighbors.end())
			{
				//At least one of the neighbors of current vertex is not valid.
				//Invalid points in the 1-ring neighborhood.
				m_invalidSmoothVertices[i] = true;
				continue;
			}

			if(twoRingBoundaryNeighbors.find(i) != twoRingBoundaryNeighbors.end())
			{
				//At least one of the neighbors of the neighbors of current vertex is not valid.
				//Invalid points in the 2-ring neighborhood.
				m_invalidSmoothVertices[i] = true;
				//continue;
			}

			std::vector<size_t> currNeighbors;
			currNeighbors.push_back(i-1-m_nFullWidth);
			currNeighbors.push_back(i-m_nFullWidth);
			currNeighbors.push_back(i-m_nFullWidth+1);
			currNeighbors.push_back(i-1);
			currNeighbors.push_back(i+1);
			currNeighbors.push_back(i-1+m_nFullWidth);
			currNeighbors.push_back(i+m_nFullWidth);	
			currNeighbors.push_back(i+m_nFullWidth+1);

			neighbors[i] = currNeighbors;		
		}
	}
}

void CWaveletShapeMultiLinearOptimizer::precomputeSmoothingWeights()
{
	const size_t numGridNeighbors = m_gridNeighbors.size();
	m_precomputedSmoothIndices.resize(numGridNeighbors);
	m_precomputedSmoothWeights.resize(numGridNeighbors);
	
	std::map<size_t, double>::iterator mapIter;
	std::map<size_t, double>::iterator endMapIter;

	for(size_t pIndex = 0; pIndex < numGridNeighbors; ++pIndex)
	{
		if(m_invalidSmoothVertices[pIndex])
		{
			continue;
		}

		const std::vector<size_t>& pNeighbors = m_gridNeighbors[pIndex];
		if(pNeighbors.empty())
		{
			continue;
		}

		std::map<size_t, double> indexWeightMap;
		indexWeightMap.insert(std::make_pair(pIndex, 1.0));

		const size_t numPNeighbors = pNeighbors.size();	
		const double factor1 = -(2.0/static_cast<double>(numPNeighbors));

		for(size_t j = 0; j < numPNeighbors; ++j)
		{
			const size_t qIndex = pNeighbors[j];
			const std::vector<size_t>& qNeighbors = m_gridNeighbors[qIndex];

			if(pNeighbors.empty())
			{
				break;
			}

			mapIter = indexWeightMap.find(qIndex);
			endMapIter = indexWeightMap.end();
			if(mapIter == endMapIter)
			{
				indexWeightMap.insert(std::make_pair(qIndex, factor1));
			}
			else
			{
				mapIter->second += factor1;
			}

			const size_t numQNeighbors = qNeighbors.size();
			const double factor2 = (1.0/static_cast<double>(numPNeighbors*numQNeighbors));

			for(size_t k = 0; k < numQNeighbors; ++k)
			{
				const size_t kIndex = qNeighbors[k];

				mapIter = indexWeightMap.find(kIndex);
				endMapIter = indexWeightMap.end();
				if(mapIter == endMapIter)
				{
					indexWeightMap.insert(std::make_pair(kIndex, factor2));
				}
				else
				{
					mapIter->second += factor2;				
				}
			}
		}

		const size_t numElements = indexWeightMap.size();

		std::vector<size_t>& indices = m_precomputedSmoothIndices[pIndex];
		indices.reserve(numElements);

		std::vector<double>& weights = m_precomputedSmoothWeights[pIndex];
		weights.reserve(numElements);

		mapIter=indexWeightMap.begin();
		endMapIter = indexWeightMap.end();
		for(; mapIter != endMapIter; ++mapIter)
		{
			indices.push_back(mapIter->first);
			weights.push_back(mapIter->second);
		}
	}
}