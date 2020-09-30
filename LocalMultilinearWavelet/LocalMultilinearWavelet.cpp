#include "StereoFace.h"
#include "BSplineGridWavelet.h"
#include "BSGWFactory.h"
#include "NearestNeighborAssistant.h"
#include "WaveletShapeMultiLinearOptimizer.h"
#include "FaceData.h"
#include "GlobalImageMatch.h"

#include "FileLoader.h"
#include "MultilinearModelHandler.h"

#include <iostream>
#include <fstream>
#include <time.h>
#include <string>
#include <stdio.h>

//Model specific constants used for the projection plane
#define SF_EMPERICAL_VIEW_X					-0.6704
#define SF_EMPERICAL_VIEW_Y					-0.7420
#define SF_EMPERICAL_VIEW_Z					0.0

//constants for indices of some landmarks on mean face of neutral expression
//for high-res scaled align version of parameterized BUFE3D database
//used to align template with 2D grid for wavelet decomposition
#define SF_BUFE3D_MEAN_NE_NOSETIP			4303
#define SF_BUFE3D_MEAN_NE_NOSEBRIDGE		3692
#define SF_BUFE3D_MEAN_NE_LEFTEYEOC			4736
#define SF_BUFE3D_MEAN_NE_RIGHTEYEOC		2358

//Width and height of the initial base mesh used for subdivision
#define SF_MODEL_BASE_MESH_WIDTH			5
#define SF_MODEL_BASE_MESH_HEIGHT		7
//Number of subdivision levels
#define SF_MODEL_LEVELS						6

const double INIT_SMOOTHING_NEIGHBOR_WEIGHT	= 1000.0;
const double SMOOTHING_NEIGHBOR_WEIGHT			= 100.0;

//Outputs a smoothed version of the fitting result
//#define OUTPUT_REFINED_MESH

//Define the number of fitting refinement steps (needs to be > 0 and <= SF_MODEL_LEVELS)
#define SF_REFINE_GEOMETRY_LEVEL			6

CProfile g_profInitNN;

//Help functions:

//Define plane for 2D sterographic projection and compute the projection.
GLuint initTemplateGeometryMap(CBSplineGridWavelet<C3Vectorf> *& wavelet, GLfloat *& pGeomMapMask, double *& pTemplate2DPositions,
							   double & templateMapLeft, double & templateMapRight, double & templateMapBottom, double & templateMapTop,
							   GLuint *& pGeomMapIndices, GLuint & ntexTemplateGeometryMap, GLuint & ntexTemplateDepthMap, double *& matTemplateStereoMap, 
							   IplImage *& pTemplate3DMap, abutil::C3Vectorf *& pGeomMapNormals);
IplImage* generateTemplateGeometryMap(CBSplineGridWavelet<C3Vectorf> * wavelet, GLfloat * pGeomMapMask, double * pFaceVertices, GLuint & nprogTemplateGeometryGen, 
									  GLuint & nfb, GLuint & ntexTemplateGeometryMap, GLuint & ntexTemplateDepthMap, double * matTemplateStereoMap,
									  double & templateMapLeft, double & templateMapRight, double & templateMapBottom, double & templateMapTop, double * pTemplate2DPositions, 
									  IplImage * pTemplate3DMap, GLuint * pGeomMapIndices);
void initGeometryMapIndices(int nWidth, int nHeight, GLuint * pIndices);
void maskGeometryMapIndices(int nWidth, int nHeight, GLfloat * pMask, GLuint * pIndices);
void initDisplay();

//Save triangle mesh to file.
void saveGeometryMapAsTriangleMesh(char* szFilename, int nWidth, int nHeight, abutil::C3Vectorf* pVertices, float* pMask, float maskThresh, abutil::C3Vectorf* pGeomMapColors);

void initStereoPoints(std::vector<SStereoPatch> patches, float *& pStereoPoints, float *& pStereoNormals, int & nStereoPoints, CNearestNeighborAssistant *& NNA);

//Compute rigid alignment between model points and observed points.
void computeLSAlignmentModelToData(int nPoints, std::vector<double>& modelPoints, std::vector<double>& observedPoints, double* pmatTransform, double* pmatTransformInv);

//Read triangle mesh from file.
bool readMeshAsPatches(const std::string& strFilename, std::vector<SStereoPatch> &patches, float g_nInDataScale);

///////////////////////////////////////////////////////////////////////////////
//main - program starts here
///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
	double templateMapLeft, templateMapRight, templateMapBottom, templateMapTop;
	double * pFaceVertices = NULL, * pFaceNormals = NULL, * pTemplate2DPositions = NULL, * matTemplateStereoMap = NULL;
	int nFaceVertices;
	GLuint nfb, ntexTemplateGeometryMap, ntexTemplateDepthMap, nprogTemplateGeometryGen;
	CBSplineGridWavelet<C3Vectorf> * wavelet = NULL;
	abutil::C3Vectorf * pGeomMapNormals = NULL, * pVertData = NULL;
	IplImage * pTemplate3DMap = NULL;
	GLfloat * pGeomMapMask = NULL;
	GLuint * pGeomMapIndices = NULL;

	//Create a Window (needed for parameterization):
	glutInit(&argc, argv);
	glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH );
	glutInitWindowPosition( 0, 0 );
	glutInitWindowSize( 1024, 768 );
	int windowId = glutCreateWindow( "StereoFace" );

	g_profInitNN.setName("init nearest neighbor look-up");

	std::cout << "initializing GL stuff...\n"; std::cout.flush();
	//Setup:
	initDisplay();
	//initialize the framebuffer
	fbInit(nfb);

	if(argc <= 7)
	{
		printf("Wrong number of command line arguments.\n");
		return 1;
	}

	int index(0);
	char* szModelFile = argv[++index];
	std::cout << "Statistical model: " << szModelFile << "\n"; std::cout.flush();

	const char* cstrTemplateFileName = argv[++index];
	std::cout << "Template geometry file: " << cstrTemplateFileName << "\n"; std::cout.flush();

	const char* cstrModelLandmarksFileName = argv[++index];
	std::cout << "Model landmarks file: " << cstrModelLandmarksFileName << "\n"; std::cout.flush();

	std::cout << "Initializing face data...\n"; std::cout.flush();
	if(!CFaceData::init(cstrTemplateFileName))
	{
		printf("Problem initializing face data\n");
		return 1;
	}

	nFaceVertices = CFaceData::m_nVertices;
	pFaceVertices = new double[nFaceVertices * 3];
	_ASSERT(pFaceVertices != NULL);
	pFaceNormals = new double[nFaceVertices * 3];
	_ASSERT(pFaceNormals != NULL);

	wavelet = (CBSplineGridWavelet<C3Vectorf>*)CBSGWFactory::createBSGW(eBSGWFloat3, SF_MODEL_BASE_MESH_WIDTH, SF_MODEL_BASE_MESH_HEIGHT, SF_MODEL_LEVELS);
	_ASSERT(wavelet != NULL);
	
	pTemplate3DMap = NULL;
	pGeomMapIndices = NULL;
	pTemplate2DPositions = NULL;
	pGeomMapMask = NULL;

	matTemplateStereoMap = new double[16];

	std::cout << "Initializing template geometry resampling...\n"; std::cout.flush();
	nprogTemplateGeometryGen = initTemplateGeometryMap(wavelet, pGeomMapMask, pTemplate2DPositions, 
		templateMapLeft, templateMapRight, templateMapBottom, templateMapTop, pGeomMapIndices, ntexTemplateGeometryMap, ntexTemplateDepthMap, 
		matTemplateStereoMap, pTemplate3DMap, pGeomMapNormals);

	std::cout << "Generating template geometry map...\n"; std::cout.flush();
	generateTemplateGeometryMap(wavelet, pGeomMapMask, CFaceData::m_pAvgVertices, nprogTemplateGeometryGen, nfb, ntexTemplateGeometryMap, ntexTemplateDepthMap, 
								matTemplateStereoMap, templateMapLeft, templateMapRight, templateMapBottom, templateMapTop, pTemplate2DPositions, pTemplate3DMap, pGeomMapIndices);

	if(nprogTemplateGeometryGen == -1)
	{
		printf("Problem in initTemplateGeometryMap\n");
		return 1;
	}

	//Read a 3D point cloud (possibly generated using stereo data)
	const char* cstrTargetFileName = argv[++index];;

	std::cout << "reading input data: " << cstrTargetFileName << "...\n"; std::cout.flush();
	std::vector<SStereoPatch> patches;
	if (!readMeshAsPatches(cstrTargetFileName, patches, 1.0))
	{
		std::cout << "unable to load input file " << cstrTargetFileName << "\n";
		return 1;
	}

	FileLoader lmkLoader;

	std::vector<double> modelLandmarks;
	std::vector<bool> modelLandmarksLoaded;
	if(!lmkLoader.loadLandmarks(std::string(cstrModelLandmarksFileName), modelLandmarks, modelLandmarksLoaded))
	{
		printf("Problems reading model landmarks %s.\n", cstrModelLandmarksFileName);
		return 1;
	}

	const char* cstrTargetLandmarksFileName = argv[++index];

	std::vector<double> targetDataLandmarks;
	std::vector<bool> targetDataLandmarksLoaded;
	if(!lmkLoader.loadLandmarks(std::string(cstrTargetLandmarksFileName), targetDataLandmarks, targetDataLandmarksLoaded))
	{
		printf("Problems reading target landmarks %s.\n", cstrTargetLandmarksFileName);
		return 1;
	}

	std::vector<double> cleanModelLandmarks;
	std::vector<double> cleanDataLandmarks;

	const size_t maxNumLandmarks = std::min<size_t>(modelLandmarksLoaded.size(), targetDataLandmarksLoaded.size());

	int numRigidLmks(0), numFitLmks(0);
	for(size_t iLmk = 0; iLmk < maxNumLandmarks; ++iLmk)
	{
		if(!modelLandmarksLoaded[iLmk] || !targetDataLandmarksLoaded[iLmk])
		{
			continue;
		}

		cleanModelLandmarks.push_back(modelLandmarks[3*iLmk + 0]);
		cleanModelLandmarks.push_back(modelLandmarks[3*iLmk + 1]);
		cleanModelLandmarks.push_back(modelLandmarks[3*iLmk + 2]);

		cleanDataLandmarks.push_back(targetDataLandmarks[3*iLmk + 0]);
		cleanDataLandmarks.push_back(targetDataLandmarks[3*iLmk + 1]);
		cleanDataLandmarks.push_back(targetDataLandmarks[3*iLmk + 2]);

		if(iLmk < 8)
		{
			++numRigidLmks;
		}

		++numFitLmks;
	}

	double matAlignModelToData[16], matAlignDataToModel[16];

	//Read the output filename:
	const char* outputFileName = argv[++index];

	double nnthresh = 10.0;
	if (argc > index)
		nnthresh = atof(argv[++index]);

	double allowDeviation = 1.0;
	if (argc > index)
		allowDeviation = atof(argv[++index]);

	double lmkSmoothingWeight = INIT_SMOOTHING_NEIGHBOR_WEIGHT;
	if (argc > index)
		lmkSmoothingWeight = atof(argv[++index]);

	double nnSmoothingWeight = SMOOTHING_NEIGHBOR_WEIGHT;
	if (argc > index)
		nnSmoothingWeight = atof(argv[++index]);

	std::cout << "initializing nearest neighbor queries...\n"; std::cout.flush();
	//initialize stereo point cloud and nearest neighbor assistant
	CNearestNeighborAssistant* pNNA = new CNearestNeighborAssistant();

	float* pStereoPoints = NULL;
	float* pStereoNormals = NULL;
	int nStereoPoints = 0;
	g_profInitNN.start();
	initStereoPoints(patches, pStereoPoints, pStereoNormals, nStereoPoints, pNNA);
	g_profInitNN.stop();

	std::cout << "initializing shape optimizer...\n"; std::cout.flush();
	std::cout << "loading statistical models from file...\n"; std::cout.flush();
	//load the learned model from file
	CWaveletShapeMultiLinearOptimizer wsmlo;
	
	//Load binary model
	wsmlo.loadModels(szModelFile);
	
	//Load text model
	//wsmlo.loadModelsText(szModelFile);
	wsmlo.setVertexMask(pGeomMapMask);
	wsmlo.computeGridNeighbors();
	wsmlo.precomputeSmoothingWeights();

	uchar* pRawData;
	cvGetRawData(pTemplate3DMap, &pRawData);

	//find landmarks on template, use to get look-up for resampled surface
	std::vector<size_t> modelLmksResIdx;
	modelLmksResIdx.resize(numFitLmks);

	for (size_t ilmk = 0; ilmk < numFitLmks; ilmk++)
	{
		const int m = ilmk * 3;
		const abutil::C3Vectord landmark(cleanModelLandmarks[m + 0], cleanModelLandmarks[m + 1], cleanModelLandmarks[m + 2]);

		int ivert = CFaceData::m_nVertices;
		double distMin = 1e30f;
		for (int i = 0; i < CFaceData::m_nVertices; i++)
		{
			const int k = i * 3;
			const abutil::C3Vectord vert(CFaceData::m_pAvgVertices[k + 0], CFaceData::m_pAvgVertices[k + 1], CFaceData::m_pAvgVertices[k + 2]);
			const double dist = (vert - landmark).length();
			if (dist < distMin)
			{
				ivert = i;
				distMin = dist;
			}
		}

		if (ivert < CFaceData::m_nVertices)
		{
			const int iv2 = ivert * 2;
			const double u = pTemplate2DPositions[iv2 + 0];
			const double v = pTemplate2DPositions[iv2 + 1];
			const double u_scale = ((u - templateMapLeft) / (templateMapRight - templateMapLeft)) * (double)wavelet->getFullResWidth();
			const double v_scale = ((v - templateMapBottom) / (templateMapTop - templateMapBottom)) * (double)wavelet->getFullResHeight();
			int iu, iv;
			if (u_scale - floor(u_scale) > 0.5)
				iu = 1 + (int)u_scale;
			else
				iu = (int)u_scale;
			if (v_scale - floor(v_scale) > 0.5)
				iv = 1 + (int)v_scale;
			else
				iv = (int)v_scale;

			//std::cout << ilmk << ": " << iu << ", " << iv << "\n";
			const size_t iresampleVert = (size_t)(iv * wavelet->getFullResWidth() + iu);
			//std::cout << "\t" << iresampleVert << "\n";
			modelLmksResIdx[ilmk] = iresampleVert;
		}
	}

	//set landmarks
	wsmlo.setLandmarks(numFitLmks, modelLmksResIdx, cleanDataLandmarks);

	//fit model to given landmarks for initialization
	wsmlo.setFitToLandmarksOnly(true);
	wsmlo.setNeighborSmoothWeight(lmkSmoothingWeight);

	//set valid range of model parameter by specifying the size of the prior box
	//prior box of size 1 for landmark fitting, since the landmark positions are reliable
	wsmlo.setAllowedDeviation(1.f);

	size_t lmkrIter = 0;
	std::cout << "optimizing shape model (landmarks only)...\n"; std::cout.flush();
	//optimize/fit model
	for (int j = 0; j < 3/*SF_REFINE_GEOMETRY_LEVEL*/; j++)
	{
		std::cout << "\toptimizing level " << j << "...\n"; std::cout.flush();
		const int levelIters = 4;//SF_REFINE_GEOMETRY_LEVEL;// - j;
		for (int k = 0; k < levelIters; k++)
		{
			//retrieve vertex position associated with landmark indices
			pVertData = (C3Vectorf*)pRawData;
			wsmlo.getInternalReconstruction(pVertData);
			for (size_t ilmk = 0; ilmk < numFitLmks; ilmk++)
			{
				const size_t ilmk3 = ilmk * 3;
				cleanModelLandmarks[ilmk3 + 0]	 = pVertData[modelLmksResIdx[ilmk]].x;
				cleanModelLandmarks[ilmk3 + 1]	 = pVertData[modelLmksResIdx[ilmk]].y;
				cleanModelLandmarks[ilmk3 + 2]	 = pVertData[modelLmksResIdx[ilmk]].z;
			}

			//compute alignment
			if (lmkrIter == 0)
				computeLSAlignmentModelToData(numRigidLmks, cleanModelLandmarks, cleanDataLandmarks, matAlignModelToData, matAlignDataToModel);
			else
				computeLSAlignmentModelToData(numFitLmks, cleanModelLandmarks, cleanDataLandmarks, matAlignModelToData, matAlignDataToModel);

			wsmlo.setNNA(pNNA, (NNAReal*)pStereoPoints, (NNAReal*)pStereoNormals, (size_t)nStereoPoints, matAlignModelToData, matAlignDataToModel, nnthresh, 0.003);
			/*wsmlo.setVertexMask(pGeomMapMask);*/

			wsmlo.optimizeActiveWaveletModel();
			if (k < levelIters - 1)
				wsmlo.resetFitting(j);

			wsmlo.getInternalReconstructionDataCoords(pVertData);

			lmkrIter++;
		}
	}

	wsmlo.resetFitting();

	//fit model to all data vertices
	wsmlo.setFitToLandmarksOnly(false);
	wsmlo.setNeighborSmoothWeight(nnSmoothingWeight);

	wsmlo.setAllowedDeviation(allowDeviation);

	std::cout << "optimizing shape model...\n"; std::cout.flush();
	//optimize/fit model
	for (int j = 0; j < SF_REFINE_GEOMETRY_LEVEL; j++)
	{
		wsmlo.setLandmarkRelativeWeight(0.5);

		const int levelIters = SF_REFINE_GEOMETRY_LEVEL - j;
		for (int k = 0; k < levelIters; k++)
		{
			std::cout << "\toptimizing level " << j << "...\n"; std::cout.flush();
			wsmlo.optimizeActiveWaveletModel();

			if (k < levelIters - 1)
				wsmlo.resetFitting(j);

			wsmlo.setLandmarkRelativeWeight(0.f);
		}
	}

	char szFinalFilename[8192];
	cvGetRawData(pTemplate3DMap, &pRawData);
	pVertData = (C3Vectorf*)pRawData;

	wsmlo.getInternalReconstructionDataCoords(pVertData);
	sprintf(szFinalFilename, "%s.off", outputFileName);
	saveGeometryMapAsTriangleMesh(szFinalFilename, wavelet->getFullResWidth(), wavelet->getFullResHeight(), pVertData, pGeomMapMask, 0.5f, NULL);

#ifdef OUTPUT_REFINED_MESH
	//refine surface
	wsmlo.setSmoothRelativeWeight(20.f);
	//wsmlo.setNumRefineInnerIters(1);
	wsmlo.setNumRefineInnerIters(4);
	wsmlo.refineSurface();

	//will be the same as the final level unless refinement is done
	wsmlo.getInternalReconstructionDataCoords(pVertData);
	sprintf(szFinalFilename, "%s_refined.off", outputFileName);
	saveGeometryMapAsTriangleMesh(szFinalFilename, wavelet->getFullResWidth(), wavelet->getFullResHeight(), pVertData, pGeomMapMask, 0.5f, NULL);
#endif

	wsmlo.reportProfiling(stdout);
	g_profInitNN.computeStatistics();
	g_profInitNN.reportStatistics(stdout);

	glutDestroyWindow(windowId);
	delete [] pFaceVertices;
	delete [] pFaceNormals;
	delete [] matTemplateStereoMap;
	delete wavelet;
	
	//Free stuff allocated in initTemplateGeometryMap
	cvReleaseImage(&pTemplate3DMap);
	delete [] pGeomMapMask;
	delete [] pTemplate2DPositions;
	delete [] pGeomMapIndices;

	return 0;
}


GLuint initTemplateGeometryMap(CBSplineGridWavelet<C3Vectorf> *& wavelet, GLfloat *& pGeomMapMask, double *& pTemplate2DPositions,
							   double & templateMapLeft, double & templateMapRight, double & templateMapBottom, double & templateMapTop,
							   GLuint *& pGeomMapIndices, GLuint & ntexTemplateGeometryMap, GLuint & ntexTemplateDepthMap, double *& matTemplateStereoMap, 
							   IplImage *& pTemplate3DMap, abutil::C3Vectorf *& pGeomMapNormals)
{
	GLint success;
	CvSize mapSize;
	double empNorm[3];

	//init image
	mapSize.width = wavelet->getFullResWidth();
	mapSize.height = wavelet->getFullResHeight();

	pTemplate3DMap = cvCreateImage(mapSize, IPL_DEPTH_32F, 3);
	_ASSERT(pTemplate3DMap != NULL);

	pGeomMapNormals = new C3Vectorf[mapSize.height * mapSize.width];
	_ASSERT(pGeomMapNormals != NULL);

	pGeomMapMask = new GLfloat[mapSize.height * mapSize.width]; 
	_ASSERT(pGeomMapMask != NULL);

	//init 2D mapping
	pTemplate2DPositions = new double[CFaceData::m_nVertices * 2];
	_ASSERT(pTemplate2DPositions != NULL);

	//select emperically derived plane for stereographic projection
	empNorm[0] = SF_EMPERICAL_VIEW_X;
	empNorm[1] = SF_EMPERICAL_VIEW_Y;
	empNorm[2] = SF_EMPERICAL_VIEW_Z;
	normalize(empNorm);

	std::cout << "computing projection plane...\n"; std::cout.flush();
	//use landmark indices to fix horizontal and vertical grid directions
	int landmarkIndices[4] = { SF_BUFE3D_MEAN_NE_NOSETIP, SF_BUFE3D_MEAN_NE_NOSEBRIDGE, SF_BUFE3D_MEAN_NE_LEFTEYEOC, SF_BUFE3D_MEAN_NE_RIGHTEYEOC };
	CFaceData::empericalPlaneFromLandmarkIndices(empNorm, landmarkIndices, 4, matTemplateStereoMap);
	matTemplateStereoMap[12] = 0.0;
	matTemplateStereoMap[13] = 0.0;
	matTemplateStereoMap[14] = 0.0;
	matTemplateStereoMap[15] = 1.0;

	std::cout << "performing stereographic projection...\n"; std::cout.flush();
	//compute 2D stereographic projections of vertices, and extents of projection
	CFaceData::stereographicProject(CFaceData::m_pAvgVertices, CFaceData::m_nVertices, matTemplateStereoMap, pTemplate2DPositions, &templateMapLeft, &templateMapBottom, &templateMapRight, &templateMapTop);
	
	std::cout << "geometry map: [" << templateMapLeft << ", " << templateMapRight << "] x [" << templateMapBottom << ", " << templateMapTop << "]" << std::endl; std::cout.flush();

	//init indices for display
	int nQuadsX = pTemplate3DMap->width - 1;
	int nQuadsY = pTemplate3DMap->height - 1;
	int nQuads = nQuadsX * nQuadsY;
	int nIndices = nQuads * 4;
	pGeomMapIndices = new GLuint[nIndices];
	_ASSERT(pGeomMapIndices != NULL);

	initGeometryMapIndices(pTemplate3DMap->width, pTemplate3DMap->height, pGeomMapIndices);

	//init texture(s)
	ntexTemplateGeometryMap = 0;
	glGenTextures(1, &ntexTemplateGeometryMap);
	glBindTexture(GL_TEXTURE_2D, ntexTemplateGeometryMap);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, 0);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F_ARB, mapSize.width, mapSize.height, 0, GL_RGBA, GL_FLOAT, NULL);

	glBindTexture(GL_TEXTURE_2D, 0);

	ntexTemplateDepthMap = 0;
	glGenTextures(1, &ntexTemplateDepthMap);
	glBindTexture(GL_TEXTURE_2D, ntexTemplateDepthMap);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, 0);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT32_ARB, mapSize.width, mapSize.height, 0, GL_DEPTH_COMPONENT, GL_FLOAT, NULL);
	checkGLError("glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT32_ARB, mapSize.width, mapSize.height, 0, GL_DEPTH_COMPONENT, GL_FLOAT, NULL)");

	glBindTexture(GL_TEXTURE_2D, 0);

	//init shaders
	GLuint nvsTemplateGeometryGen = loadVertexShader("TemplateGeometryProject.vs");
	glCompileShader(nvsTemplateGeometryGen);
	glGetShaderiv(nvsTemplateGeometryGen, GL_COMPILE_STATUS, &success);
	if (!success)
	{
		shaderError("Error compiling TemplateGeometryProject.vs", nvsTemplateGeometryGen);
		return -1;
	}

	GLuint nfsTemplateGeometryGen = loadFragmentShader("TemplateGeometryMap.fs");
	glCompileShader(nfsTemplateGeometryGen);
	glGetShaderiv(nfsTemplateGeometryGen, GL_COMPILE_STATUS, &success);
	if (!success)
	{
		shaderError("Error compiling TemplateGeometryMap.fs", nfsTemplateGeometryGen);
		return -1;
	}

	GLuint nprogTemplateGeometryGen = glCreateProgram();
	glAttachShader(nprogTemplateGeometryGen, nvsTemplateGeometryGen);
	glAttachShader(nprogTemplateGeometryGen, nfsTemplateGeometryGen);
	glLinkProgram(nprogTemplateGeometryGen);
	glGetProgramiv(nprogTemplateGeometryGen, GL_LINK_STATUS, &success);
	if (!success)
	{
		programError("Error linking nprogTemplateGeometryGen", nprogTemplateGeometryGen);
		return -1;
	}
	glValidateProgram(nprogTemplateGeometryGen);
	glGetProgramiv(nprogTemplateGeometryGen, GL_VALIDATE_STATUS, &success);
	if (!success)
	{
		programError("Error validating nprogTemplateGeometryGen", nprogTemplateGeometryGen);
		return -1;
	}

	return nprogTemplateGeometryGen;
}

void initGeometryMapIndices(int nWidth, int nHeight, GLuint* pIndices)
{
	int i, j, k;
	int nQuadsX = nWidth - 1;
	int nQuadsY = nHeight - 1;
	int nQuads = nQuadsX * nQuadsY;
	int nIndices = nQuads * 4;

	k = 0;
	for (i = 0; i < nQuadsY; i++)
	{
		for (j = 0; j < nQuadsX; j++)
		{
			pIndices[k]		= i * nWidth + j;
			pIndices[k + 1]	= (i + 1) * nWidth + j;
			pIndices[k + 2]	= (i + 1) * nWidth + j + 1;
			pIndices[k + 3]	= i * nWidth + j + 1;
			k += 4;
		}
	}
}

IplImage* generateTemplateGeometryMap(CBSplineGridWavelet<C3Vectorf>* wavelet, GLfloat * pGeomMapMask, double* pFaceVertices, GLuint & nprogTemplateGeometryGen, 
									  GLuint & nfb, GLuint & ntexTemplateGeometryMap, GLuint & ntexTemplateDepthMap, double * matTemplateStereoMap,
									  double & templateMapLeft, double & templateMapRight, double & templateMapBottom, double & templateMapTop, double * pTemplate2DPositions,
									  IplImage * pTemplate3DMap, GLuint * pGeomMapIndices)
{
	int i;

	double * matRigidFaceAlign = new double[16];

	for(i = 0; i < 16; i++) matRigidFaceAlign[i] = 0;
	matRigidFaceAlign[0] = matRigidFaceAlign[5] = matRigidFaceAlign[10] = matRigidFaceAlign[15] = 1.0;

	float matRigid[16];

	for(i = 0; i < 16; i++)
		matRigid[i] = (float)matRigidFaceAlign[i];

	fbBind(nfb);
	fbAttachTexture(nfb, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, ntexTemplateGeometryMap);
	fbAttachTexture(nfb, GL_DEPTH_ATTACHMENT_EXT, GL_TEXTURE_2D, ntexTemplateDepthMap);
	glDrawBuffer(GL_COLOR_ATTACHMENT0_EXT);
	glViewport(0, 0, wavelet->getFullResWidth(), wavelet->getFullResHeight());
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glUseProgram(nprogTemplateGeometryGen);

	GLint iloc = glGetUniformLocation(nprogTemplateGeometryGen, "matRigidAlign");
	glUniformMatrix4fv(iloc, 1, true, matRigid);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadTransposeMatrixd(matTemplateStereoMap);
//	glMultTransposeMatrixd(matRigidFaceAlign);
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	gluOrtho2D(templateMapLeft * 1.02, templateMapRight * 1.02, templateMapBottom * 1.02, templateMapTop * 1.02);
	//gluOrtho2D(templateMapLeft, templateMapRight, templateMapBottom, templateMapTop);

	//draw face with texture coordinates
	CFaceData::drawFaceWithTexCoords(pFaceVertices, pTemplate2DPositions);

	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
	glFinish();
	glUseProgram(0);
	fbDetachAll(nfb);
	fbUnbind();

	uchar* pData;
	int nStep;
	cvGetRawData(pTemplate3DMap, &pData, &nStep);

	glBindTexture(GL_TEXTURE_2D, ntexTemplateGeometryMap);
	glGetTexImage(GL_TEXTURE_2D, 0, GL_RGB, GL_FLOAT, pData);
	glGetTexImage(GL_TEXTURE_2D, 0, GL_ALPHA, GL_FLOAT, pGeomMapMask);
	glBindTexture(GL_TEXTURE_2D, 0);

	//flip the image vertically (i.e. about x-axis)
//	cvFlip(pTemplate3DMap);

	float* pGeomMaskCopy = new float[pTemplate3DMap->height * pTemplate3DMap->width];
	_ASSERT(pGeomMaskCopy != NULL);
	memcpy(pGeomMaskCopy, pGeomMapMask, pTemplate3DMap->height * pTemplate3DMap->width * sizeof(float));

	const float mask_blend_eps = 0.2f;
	int x, y, blendedVerts, nNbrs;
	float w, weightTotal;
	abutil::C3Vectorf* pVerts = (abutil::C3Vectorf*)pData;
	abutil::C3Vectorf vHome;
	do
	{
		blendedVerts = 0;
		i = 0;
		for (y = 0; y < pTemplate3DMap->height; y++)
		{
			for (x = 0; x < pTemplate3DMap->width; x++)
			{
				if (pGeomMaskCopy[i] < 1.f - mask_blend_eps)
				{
					weightTotal = 0.f;
					vHome.set(0.f, 0.f, 0.f);
					nNbrs = 0;
					if (y > 0)
					{
						w = pGeomMaskCopy[i - pTemplate3DMap->width];
						vHome += pVerts[i - pTemplate3DMap->width] * w;
						weightTotal += w;
						nNbrs++;
					}
					if (x > 0)
					{
						w = pGeomMaskCopy[i - 1];
						vHome += pVerts[i - 1] * w;
						weightTotal += w;
						nNbrs++;
					}
					if (y < pTemplate3DMap->height - 1)
					{
						w = pGeomMaskCopy[i + pTemplate3DMap->width];
						vHome += pVerts[i + pTemplate3DMap->width] * w;
						weightTotal += w;
						nNbrs++;
					}
					if (x < pTemplate3DMap->width - 1)
					{
						w = pGeomMaskCopy[i + 1];
						vHome += pVerts[i + 1] * w;
						weightTotal += w;
						nNbrs++;
					}

					if (weightTotal > mask_blend_eps)
					{
						vHome /= weightTotal;
						weightTotal /= (float)nNbrs;

						pVerts[i] = vHome;
						pGeomMaskCopy[i] = 1.f;//weightTotal;
					}

					blendedVerts++;
				}

				i++;
			}
		}
	} while (blendedVerts > 0);

	//reset geometry map indices, then apply mask
	initGeometryMapIndices(pTemplate3DMap->width, pTemplate3DMap->height, pGeomMapIndices);
	maskGeometryMapIndices(pTemplate3DMap->width, pTemplate3DMap->height, pGeomMapMask, pGeomMapIndices);

	//set full resolution data of the wavelet
	wavelet->setFullResolutionData((C3Vectorf*)pData);

	delete [] matRigidFaceAlign;
	delete [] pGeomMaskCopy;

	return pTemplate3DMap;
}

void maskGeometryMapIndices(int nWidth, int nHeight, GLfloat* pMask, GLuint* pIndices)
{
	int i, j, k;
	int iv0, iv1, iv2, iv3;
	int nQuadsX = nWidth - 1;
	int nQuadsY = nHeight - 1;
	int nQuads = nQuadsX * nQuadsY;
	int nIndices = nQuads * 4;
//	GLuint iMaskValue = (GLuint)(-1);
	GLuint iMaskValue = CGlobalImageMatch::gim_indexMaskValue;

	k = 0;
	for (i = 0; i < nQuadsY; i++)
	{
		for (j = 0; j < nQuadsX; j++)
		{
			iv0 = pIndices[k];
			iv1 = pIndices[k + 1];
			iv2 = pIndices[k + 2];
			iv3 = pIndices[k + 3];
			if (pMask[iv0] < 0.5 || pMask[iv1] < 0.5 || pMask[iv2] < 0.5f || pMask[iv3] < 0.5f)
			{
				pIndices[k]		= iMaskValue;
				pIndices[k + 1]	= iMaskValue;
				pIndices[k + 2]	= iMaskValue;
				pIndices[k + 3] = iMaskValue;
			}
			k += 4;
		}
	}
}

void initDisplay()
{
	glewInit();

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_TEXTURE_2D);
	glDepthFunc(GL_LEQUAL);
	glShadeModel(GL_SMOOTH);
	glClearDepth(1.0);

	checkGLError("initDisplay pre-shaders");

	if (!glewIsSupported( "GL_VERSION_2_0 GL_ARB_fragment_program GL_ARB_vertex_program GL_EXT_framebuffer_object" ))
	{
		printf("Minimal extensions not supported\n");
	}
}

void saveGeometryMapAsTriangleMesh(char* szFilename, int nWidth, int nHeight, abutil::C3Vectorf* pVertices, float* pMask, float maskThresh, abutil::C3Vectorf* pGeomMapColors)
{
	_ASSERT(szFilename != NULL);
	_ASSERT(pVertices != NULL);
	_ASSERT(pMask != NULL);

	int x, y, i;

	GLuint idx;
	abutil::CXArray<GLuint> indices;
	abutil::CXArray<abutil::C3Vectorf> vertices;

	i = 0;
	for (y = 0; y < nHeight; y++)
	{
		for (x = 0; x < nWidth; x++)
		{
			if (pMask[i] > maskThresh)
				vertices.pushEnd(pVertices[i]);
			i++;
		}
	}

	i = nWidth;
	for (y = 1; y < nHeight; y++)
	{
		for (x = 0; x < nWidth; x++)
		{
			if (pMask[i] > maskThresh)
			{
				if (x < nWidth - 1)
				{
					if (pMask[i - nWidth] > maskThresh && pMask[i - nWidth + 1] > maskThresh)
					{
						idx = (GLuint)i;
						indices.pushEnd(idx);
						idx = (GLuint)(i - nWidth);
						indices.pushEnd(idx);
						idx = (GLuint)(i - nWidth + 1);
						indices.pushEnd(idx);
					}
				}

				if (x > 0)
				{
					if (pMask[i - 1] > maskThresh && pMask[i - nWidth] > maskThresh)
					{
						idx = (GLuint)(i - nWidth);
						indices.pushEnd(idx);
						idx = (GLuint)i;
						indices.pushEnd(idx);
						idx = (GLuint)(i - 1);
						indices.pushEnd(idx);
					}
					else if (pMask[i - nWidth] > maskThresh && pMask[i - nWidth - 1] > maskThresh)
					{
						idx = (GLuint)i;
						indices.pushEnd(idx);
						idx = (GLuint)(i - nWidth - 1);
						indices.pushEnd(idx);
						idx = (GLuint)(i - nWidth);
						indices.pushEnd(idx);
					}
					else if (pMask[i - 1] > maskThresh && pMask[i - nWidth - 1] > maskThresh)
					{
						idx = (GLuint)(i - 1);
						indices.pushEnd(idx);
						idx = (GLuint)(i - nWidth - 1);
						indices.pushEnd(idx);
						idx = (GLuint)i;
						indices.pushEnd(idx);
					}
				}
			}

			i++;
		}
	}

	DataContainer mesh;
	FileLoader loader;

	for (y = 0; y < nHeight; y++)
	{
		for (x = 0; x < nWidth; x++)
		{
			i = y * nWidth + x;
			mesh.getVertexList().push_back(new Vec3d(pVertices[i].x, pVertices[i].y, pVertices[i].z));

			if (pGeomMapColors != NULL)
			{
				mesh.getVertexColorList().push_back(new Vec3d(pGeomMapColors[i].x, pGeomMapColors[i].y, pGeomMapColors[i].z));
			}
		}
	}

	for (i = 0; i < indices.length(); i+=3)
	{
		mesh.getVertexIndexList().push_back(new Vec3i(indices[i + 0], indices[i + 1], indices[i + 2]));
	}

	if (!loader.saveFile(std::string(szFilename), mesh))
	{
		std::cout << "error saving mesh " << szFilename << "\n";
	}
}

void initStereoPoints(std::vector<SStereoPatch> patches, float *& pStereoPoints, float *& pStereoNormals, int & nStereoPoints, CNearestNeighborAssistant *& NNA)
{
	int i;

	if (pStereoPoints != NULL)
	{
		delete [] pStereoPoints;
		pStereoPoints = NULL;
		nStereoPoints = 0;
	}
	if (pStereoNormals != NULL)
	{
		delete [] pStereoNormals;
		pStereoNormals = NULL;
	}

	//this is where the check needs to go for what method was used to
	//get the initial stereo results, to decide where to initialize from
	//for now just assume Furukawa's method and initialize from patch set

	nStereoPoints = (int)patches.size();
	pStereoPoints = new float[nStereoPoints * 3];
	_ASSERT(pStereoPoints != NULL);
	pStereoNormals = new float[nStereoPoints * 3];
	_ASSERT(pStereoNormals != NULL);

	for (i = 0; i < nStereoPoints; i++)
	{
		pStereoPoints[i * 3 + 0]	= patches[i].px;
		pStereoPoints[i * 3 + 1]	= patches[i].py;
		pStereoPoints[i * 3 + 2]	= patches[i].pz;

		pStereoNormals[i * 3 + 0]	= patches[i].nx;
		pStereoNormals[i * 3 + 1]	= patches[i].ny;
		pStereoNormals[i * 3 + 2]	= patches[i].nz;
	}

	//now reset nearest neighbor assistant to use these points as reference points
	NNA->setReferencePoints(nStereoPoints, pStereoPoints);
}

void computeLSAlignmentModelToData(int nPoints, std::vector<double>& modelPoints, std::vector<double>& observedPoints, double* pmatTransform, double* pmatTransformInv)
{
	const double SQRT3 = sqrt(3.0);

	int i, j, nCoords = nPoints * 3;
	long int dim, num, nrhs, nWork, info;
	char job, trans;
	double alpha, beta;

	job = 'A';
	trans = 'N';
	alpha = 1.0;
	beta = 0.0;
	num = nPoints;
	dim = 4;
	nrhs = 3;
	nWork = 5 * nPoints;//this is enough for either dgels or dgesvd

	double* pNormModelPoints = new double[3 * nPoints];
	_ASSERT(pNormModelPoints != NULL);
	double* pNormObservedPoints = new double[3 * nPoints];
	_ASSERT(pNormObservedPoints != NULL);

	double* pA = new double[4 * nPoints];
	_ASSERT(pA != NULL);
	double* pB = new double[3 * nPoints];
	_ASSERT(pB != NULL);
	double* pWork = new double[nWork];
	_ASSERT(pWork != NULL);

	double R[9], S[3], U[9], VT[9];
	double matSimilModel[16], matSimilObservedInv[16], matSimilModelInv[16], matSimilObserved[16];
	double matTemp1[16], matTemp2[16], matTemp3[16], matTemp4[16];
	double x, y, z, cx, cy, cz, rmsdist, rmsscale;
	double rcpPoints = 1.0 / (double)nPoints;

	//normalize data points s.t. centroids are at origin and RMS distance (from origin) is sqrt(3)
	//save similarity transforms that do this
	//first model points
	cx = cy = cz = rmsdist = 0.0;

	//compute centroid
	for (i = 0; i < nCoords; i+=3)
	{
		cx += modelPoints[i];
		cy += modelPoints[i + 1];
		cz += modelPoints[i + 2];
	}
	cx *= rcpPoints;
	cy *= rcpPoints;
	cz *= rcpPoints;

	//subtract centroid and compute RMS distance
	for (i = 0; i < nCoords; i+=3)
	{
		pNormModelPoints[i]		= x = modelPoints[i] - cx;
		pNormModelPoints[i + 1]	= y = modelPoints[i + 1] - cy;
		pNormModelPoints[i + 2]	= z = modelPoints[i + 2] - cz;
		rmsdist += x*x + y*y + z*z;
	}
	rmsdist *= rcpPoints;
	rmsdist = sqrt(rmsdist);
	rmsscale = SQRT3 / rmsdist;

	//scale points
	for (i = 0; i < nCoords; i+=3)
	{
		pNormModelPoints[i]		*= rmsscale;
		pNormModelPoints[i + 1]	*= rmsscale;
		pNormModelPoints[i + 2]	*= rmsscale;
	}

	//"record" similarity transform
	matSimilModel[0]	= rmsscale;
	matSimilModel[1]	= 0.0;
	matSimilModel[2]	= 0.0;
	matSimilModel[3]	= 0.0;
	matSimilModel[4]	= 0.0;
	matSimilModel[5]	= rmsscale;
	matSimilModel[6]	= 0.0;
	matSimilModel[7]	= 0.0;
	matSimilModel[8]	= 0.0;
	matSimilModel[9]	= 0.0;
	matSimilModel[10]	= rmsscale;
	matSimilModel[11]	= 0.0;
	matSimilModel[12]	= -rmsscale * cx;
	matSimilModel[13]	= -rmsscale * cy;
	matSimilModel[14]	= -rmsscale * cz;
	matSimilModel[15]	= 1.0;

	matSimilModelInv[0]		= 1.0 / rmsscale;
	matSimilModelInv[1]		= 0.0;
	matSimilModelInv[2]		= 0.0;
	matSimilModelInv[3]		= 0.0;
	matSimilModelInv[4]		= 0.0;
	matSimilModelInv[5]		= 1.0 / rmsscale;
	matSimilModelInv[6]		= 0.0;
	matSimilModelInv[7]		= 0.0;
	matSimilModelInv[8]		= 0.0;
	matSimilModelInv[9]		= 0.0;
	matSimilModelInv[10]	= 1.0 / rmsscale;
	matSimilModelInv[11]	= 0.0;
	matSimilModelInv[12]	= cx;
	matSimilModelInv[13]	= cy;
	matSimilModelInv[14]	= cz;
	matSimilModelInv[15]	= 1.0;

	//now observed points
	cx = cy = cz = rmsdist = 0.0;

	//compute centroid
	for (i = 0; i < nCoords; i+=3)
	{
		cx += observedPoints[i];
		cy += observedPoints[i + 1];
		cz += observedPoints[i + 2];
	}
	cx *= rcpPoints;
	cy *= rcpPoints;
	cz *= rcpPoints;

	//subtract centroid and compute RMS distance
	for (i = 0; i < nCoords; i+=3)
	{
		pNormObservedPoints[i]		= x = observedPoints[i] - cx;
		pNormObservedPoints[i + 1]	= y = observedPoints[i + 1] - cy;
		pNormObservedPoints[i + 2]	= z = observedPoints[i + 2] - cz;
		rmsdist += x*x + y*y + z*z;
	}
	rmsdist *= rcpPoints;
	rmsdist = sqrt(rmsdist);
	rmsscale = SQRT3 / rmsdist;

	//scale points
	for (i = 0; i < nCoords; i+=3)
	{
		pNormObservedPoints[i]		*= rmsscale;
		pNormObservedPoints[i + 1]	*= rmsscale;
		pNormObservedPoints[i + 2]	*= rmsscale;
	}

	//"record" inverse similarity transform
	matSimilObservedInv[0]	= 1.0 / rmsscale;
	matSimilObservedInv[1]	= 0.0;
	matSimilObservedInv[2]	= 0.0;
	matSimilObservedInv[3]	= 0.0;
	matSimilObservedInv[4]	= 0.0;
	matSimilObservedInv[5]	= 1.0 / rmsscale;
	matSimilObservedInv[6]	= 0.0;
	matSimilObservedInv[7]	= 0.0;
	matSimilObservedInv[8]	= 0.0;
	matSimilObservedInv[9]	= 0.0;
	matSimilObservedInv[10]	= 1.0 / rmsscale;
	matSimilObservedInv[11]	= 0.0;
	matSimilObservedInv[12]	= cx;
	matSimilObservedInv[13]	= cy;
	matSimilObservedInv[14]	= cz;
	matSimilObservedInv[15]	= 1.0;

	matSimilObserved[0]		= rmsscale;
	matSimilObserved[1]		= 0.0;
	matSimilObserved[2]		= 0.0;
	matSimilObserved[3]		= 0.0;
	matSimilObserved[4]		= 0.0;
	matSimilObserved[5]		= rmsscale;
	matSimilObserved[6]		= 0.0;
	matSimilObserved[7]		= 0.0;
	matSimilObserved[8]		= 0.0;
	matSimilObserved[9]		= 0.0;
	matSimilObserved[10]	= rmsscale;
	matSimilObserved[11]	= 0.0;
	matSimilObserved[12]	= -rmsscale * cx;
	matSimilObserved[13]	= -rmsscale * cy;
	matSimilObserved[14]	= -rmsscale * cz;
	matSimilObserved[15]	= 1.0;

	//fill matrices for least-squares estimation
	for (i = 0; i < nPoints; i++)
	{
		j = i * 3;
		pA[i]				= pNormModelPoints[j];
		pA[nPoints + i]		= pNormModelPoints[j + 1];
		pA[2 * nPoints + i]	= pNormModelPoints[j + 2];
		pB[i]				= pNormObservedPoints[j];
		pB[nPoints + i]		= pNormObservedPoints[j + 1];
		pB[2 * nPoints + i]	= pNormObservedPoints[j + 2];
	}

	dim = 3;
	clapack::dgels_(&trans, &num, &dim, &nrhs, pA, &num, pB, &num, pWork, &nWork, &info);
	if (info != 0)
	{
		printf("computeLSAlignmentModelToData(...): clapack::dgels_(...) returned %i\n", info);
		goto computeLSAlignmentModelToData_EXIT;
	}

	//copy rotation separately for SVD
	R[0] = pB[0];			R[3] = pB[1];				R[6] = pB[2];
	R[1] = pB[nPoints];		R[4] = pB[nPoints + 1];		R[7] = pB[nPoints + 2];
	R[2] = pB[2 * nPoints];	R[5] = pB[2 * nPoints + 1];	R[8] = pB[2 * nPoints + 2];

	clapack::dgesvd_(&job, &job, &nrhs, &nrhs, R, &nrhs, S, U, &nrhs, VT, &nrhs, pWork, &nWork, &info);
	if (info != 0)
	{
		printf("computeLSAlignmentModelToData(...): clapack::dgesvd_(...) returned %i\n", info);
		goto computeLSAlignmentModelToData_EXIT;
	}

	//pure rotation R = U*VT
	clapack::dgemm_(&trans, &trans, &nrhs, &nrhs, &nrhs, &alpha, U, &nrhs, VT, &nrhs, &beta, R, &nrhs);

	matTemp1[0]		= R[0];
	matTemp1[1]		= R[1];
	matTemp1[2]		= R[2];
	matTemp1[3]		= 0.0;
	matTemp1[4]		= R[3];
	matTemp1[5]		= R[4];
	matTemp1[6]		= R[5];
	matTemp1[7]		= 0.0;
	matTemp1[8]		= R[6];
	matTemp1[9]		= R[7];
	matTemp1[10]	= R[8];
	matTemp1[11]	= 0.0;
	matTemp1[12]	= 0.0;
	matTemp1[13]	= 0.0;
	matTemp1[14]	= 0.0;
	matTemp1[15]	= 1.0;

	if (pmatTransformInv != NULL)
	{
		//transpose of rotation
		matTemp3[0]		= R[0];
		matTemp3[1]		= R[3];
		matTemp3[2]		= R[6];
		matTemp3[3]		= 0.0;
		matTemp3[4]		= R[1];
		matTemp3[5]		= R[4];
		matTemp3[6]		= R[7];
		matTemp3[7]		= 0.0;
		matTemp3[8]		= R[2];
		matTemp3[9]		= R[5];
		matTemp3[10]	= R[8];
		matTemp3[11]	= 0.0;
		matTemp3[12]	= 0.0;
		matTemp3[13]	= 0.0;
		matTemp3[14]	= 0.0;
		matTemp3[15]	= 1.0;
	}

	//right multiply by matSimilModel: matTemp2 = matTemp1 * matSimilModel
	clapack::integer p, r, q;
	p = 4;
	r = 4;
	q = 4;
	clapack::dgemm_(&trans, &trans, &p, &r, &q, &alpha, matTemp1, &q, matSimilModel, &q, &beta, matTemp2, &p);

	if (pmatTransformInv != NULL)
	{
		//right multiply transpose by matSimilObserved: matTemp4 = matTemp3 * matSimilObserved
		clapack::dgemm_(&trans, &trans, &p, &r, &q, &alpha, matTemp3, &q, matSimilObserved, &q, &beta, matTemp4, &p);
	}

	//left multiply by matSimilObservedInv: matTemp1 = matSimilObservedInv * matTemp2
	clapack::dgemm_(&trans, &trans, &p, &r, &q, &alpha, matSimilObservedInv, &p, matTemp2, &q, &beta, matTemp1, &p);

	if (pmatTransformInv != NULL)
	{
		//left multiply by matSimilModelInv: matTemp3 = matSimilModelInv * matTemp4
		clapack::dgemm_(&trans, &trans, &p, &r, &q, &alpha, matSimilModelInv, &p, matTemp4, &q, &beta, matTemp3, &p);
	}

	//put into row-major order (not sure why I use row-major when OpenGL and CLAPACK use column-major, but...)
	pmatTransform[0]	= matTemp1[0];
	pmatTransform[1]	= matTemp1[4];
	pmatTransform[2]	= matTemp1[8];
	pmatTransform[3]	= matTemp1[12];
	pmatTransform[4]	= matTemp1[1];
	pmatTransform[5]	= matTemp1[5];
	pmatTransform[6]	= matTemp1[9];
	pmatTransform[7]	= matTemp1[13];
	pmatTransform[8]	= matTemp1[2];
	pmatTransform[9]	= matTemp1[6];
	pmatTransform[10]	= matTemp1[10];
	pmatTransform[11]	= matTemp1[14];
	pmatTransform[12]	= matTemp1[3];
	pmatTransform[13]	= matTemp1[7];
	pmatTransform[14]	= matTemp1[11];
	pmatTransform[15]	= matTemp1[15];
	if (pmatTransformInv != NULL)
	{
		pmatTransformInv[0]		= matTemp3[0];
		pmatTransformInv[1]		= matTemp3[4];
		pmatTransformInv[2]		= matTemp3[8];
		pmatTransformInv[3]		= matTemp3[12];
		pmatTransformInv[4]		= matTemp3[1];
		pmatTransformInv[5]		= matTemp3[5];
		pmatTransformInv[6]		= matTemp3[9];
		pmatTransformInv[7]		= matTemp3[13];
		pmatTransformInv[8]		= matTemp3[2];
		pmatTransformInv[9]		= matTemp3[6];
		pmatTransformInv[10]	= matTemp3[10];
		pmatTransformInv[11]	= matTemp3[14];
		pmatTransformInv[12]	= matTemp3[3];
		pmatTransformInv[13]	= matTemp3[7];
		pmatTransformInv[14]	= matTemp3[11];
		pmatTransformInv[15]	= matTemp3[15];
	}

	if (pmatTransformInv != NULL)
	{
		trans = 'T';
		clapack::dgemm_(&trans, &trans, &p, &r, &q, &alpha, pmatTransform, &p, pmatTransformInv, &q, &beta, matTemp3, &p);
	}

computeLSAlignmentModelToData_EXIT:
	delete [] pNormModelPoints;
	delete [] pNormObservedPoints;
	delete [] pA;
	delete [] pB;
	delete [] pWork;
}

bool readMeshAsPatches(const std::string& strFilename, std::vector<SStereoPatch> &patches, float g_nInDataScale)
{
	patches.clear();

	FileLoader loader;
	DataContainer data;
	
	if (!loader.loadFile(strFilename, data))
	{
		std::cout << "error loading patches from " << strFilename << "\n";
		return false;
	}

	const size_t numVerts = data.getNumVertices();
	patches.resize(numVerts);

	for (size_t v = 0; v < numVerts; v++)
	{
		patches[v].px = (*(data.getVertexList()[v]))[0] * g_nInDataScale;
		patches[v].py = (*(data.getVertexList()[v]))[1] * g_nInDataScale;
		patches[v].pz = (*(data.getVertexList()[v]))[2] * g_nInDataScale;
		patches[v].nx = patches[v].ny = patches[v].nz = 0.0;
		patches[v].tx = patches[v].ty = patches[v].tz = 0.0;
		patches[v].bx = patches[v].by = patches[v].bz = 0.0;
		patches[v].u = patches[v].v = 0.0;
		patches[v].r = patches[v].g = patches[v].b = 0.f;
	}

	const size_t numTris = data.getVertexIndexList().size();
	for (size_t t = 0; t < numTris; t++)
	{
		const size_t v1 = (*(data.getVertexIndexList()[t]))[0];
		const size_t v2 = (*(data.getVertexIndexList()[t]))[1];
		const size_t v3 = (*(data.getVertexIndexList()[t]))[2];

		Vec3d e12, e13, e23, e21, e32, e31;

		e12[0] = (*(data.getVertexList()[v2]))[0] - (*(data.getVertexList()[v1]))[0];
		e12[1] = (*(data.getVertexList()[v2]))[1] - (*(data.getVertexList()[v1]))[1];
		e12[2] = (*(data.getVertexList()[v2]))[2] - (*(data.getVertexList()[v1]))[2];

		e13[0] = (*(data.getVertexList()[v3]))[0] - (*(data.getVertexList()[v1]))[0];
		e13[1] = (*(data.getVertexList()[v3]))[1] - (*(data.getVertexList()[v1]))[1];
		e13[2] = (*(data.getVertexList()[v3]))[2] - (*(data.getVertexList()[v1]))[2];

		e23[0] = (*(data.getVertexList()[v3]))[0] - (*(data.getVertexList()[v2]))[0];
		e23[1] = (*(data.getVertexList()[v3]))[1] - (*(data.getVertexList()[v2]))[1];
		e23[2] = (*(data.getVertexList()[v3]))[2] - (*(data.getVertexList()[v2]))[2];

		e21[0] = -e12[0];
		e21[1] = -e12[1];
		e21[2] = -e12[2];
		e31[0] = -e13[0];
		e31[1] = -e13[1];
		e31[2] = -e13[2];
		e32[0] = -e23[0];
		e32[1] = -e23[1];
		e32[2] = -e23[2];

//		const double area = 0.5;

		Vec3d norm;

		e12.crossProduct(e13, norm);
		const double area = norm.length() * 0.5;
		norm.normalize();
		patches[v1].nx += norm[0] * area;
		patches[v1].ny += norm[1] * area;
		patches[v1].nz += norm[2] * area;

		e23.crossProduct(e21, norm);
		norm.normalize();
		patches[v2].nx += norm[0] * area;
		patches[v2].ny += norm[1] * area;
		patches[v2].nz += norm[2] * area;

		e31.crossProduct(e32, norm);
		norm.normalize();
		patches[v3].nx += norm[0] * area;
		patches[v3].ny += norm[1] * area;
		patches[v3].nz += norm[2] * area;
	}

	for (size_t v = 0; v < numVerts; v++)
	{
		const double normlen = sqrt(patches[v].nx*patches[v].nx + patches[v].ny*patches[v].ny + patches[v].nz*patches[v].nz);
		if (normlen == 0.0)
		{
			//std::cout << "bad normal! " << v << "\n";
			patches[v].nx = 0.0;
			patches[v].ny = 0.0;
			patches[v].nz = 0.0;
		}
		else
		{
			const double rcpNormLen = 1.0 / normlen;
			patches[v].nx *= rcpNormLen;
			patches[v].ny *= rcpNormLen;
			patches[v].nz *= rcpNormLen;
		}
	}

	return true;
}
