///////////////////////////////////////////////////////////////////////////////
//
//	FaceData.cpp
//
//	Source file for the CFaceData class
//
//	Alan Brunton
//
///////////////////////////////////////////////////////////////////////////////


#include "FaceData.h"
#include "abutil.h"

#include <map>

///////////////////////////////////////////////////////////////////////////////
//file-scope function prototypes
///////////////////////////////////////////////////////////////////////////////
void addFaceIndexHelper(const int vertexIndex, const int faceIndex, std::map<int, std::vector<int>>& vertexFaceIndices);
void getFaceIndicesForVertices(const std::vector<std::vector<int>>& faceIndices, std::map<int, std::vector<int>>& faceIndicesForVertices);
bool readOFF_geometry(const char *szFilename, std::vector<double>& coords, std::vector<std::vector<int>>& faceIndices);


///////////////////////////////////////////////////////////////////////////////
//CFaceData static members
///////////////////////////////////////////////////////////////////////////////

int* CFaceData::m_pIndices					= NULL;
int CFaceData::m_nIndices					= 0;
double* CFaceData::m_pAvgVertices			= NULL;
int CFaceData::m_nVertices					= 0;
double* CFaceData::m_pAvgNormals			= NULL;
double CFaceData::m_centroidX				= 0.0;
double CFaceData::m_centroidY				= 0.0;
double CFaceData::m_centroidZ				= 0.0;

bool CFaceData::init(const char* cstrFaceFileName)
{
	std::vector<double> coords;
	std::vector<std::vector<int>> faceIndices;
	if(!readOFF_geometry(cstrFaceFileName, coords, faceIndices))
	{
		return false;
	}

	std::map<int, std::vector<int>> faceIndicesForVertices;
	getFaceIndicesForVertices(faceIndices, faceIndicesForVertices);

	m_nVertices = static_cast<int>(coords.size())/3;
	m_pAvgVertices = new double[3 * m_nVertices];
	_ASSERT(m_pAvgVertices != NULL);

	m_pAvgNormals = new double[3 * m_nVertices];
	_ASSERT(m_pAvgNormals != NULL);

	m_nIndices = (int) faceIndices.size() * 3;
	m_pIndices = new int[m_nIndices];
	_ASSERT(m_pIndices != NULL);

	int nElems = (int) faceIndices.size();

	for (int nId = 0; nId < m_nVertices; ++nId)
	{
		abutil::C3Vectord vp;
		vp.set(coords[3*nId], coords[3*nId+1], coords[3*nId+2]);

		m_pAvgVertices[3*nId]	= vp.x;
		m_pAvgVertices[3*nId+1]	= vp.y;
		m_pAvgVertices[3*nId+2]	= vp.z;

		m_centroidX += vp.x;
		m_centroidY += vp.y;
		m_centroidZ += vp.z;

		abutil::C3Vectord vnorm;
		vnorm.set(0.0, 0.0, 0.0);

		std::map<int, std::vector<int>>::const_iterator findIter = faceIndicesForVertices.find(nId);
		if(findIter != faceIndicesForVertices.end())
		{
			const std::vector<int>& vertexFaceIndices = findIter->second;
			const int nNodeElem = static_cast<int>(vertexFaceIndices.size());
			for (int faceId = 0; faceId < nNodeElem; ++faceId)
			{
				int n1(0);
				int n2(0);

				const int currFaceIndex = vertexFaceIndices[faceId];
				if (nId == faceIndices[currFaceIndex][0])
				{
					n1 = faceIndices[currFaceIndex][1];
					n2 = faceIndices[currFaceIndex][2];
				}
				else if (nId == faceIndices[currFaceIndex][1])
				{
					n1 = faceIndices[currFaceIndex][2];
					n2 = faceIndices[currFaceIndex][0];
				}
				else
				{
					n1 = faceIndices[currFaceIndex][0];
					n2 = faceIndices[currFaceIndex][1];
				}

				abutil::C3Vectord ve1;
				ve1.set(coords[3*n1], coords[3*n1+1], coords[3*n1+2]);
				ve1 -= vp;

				abutil::C3Vectord ve2;
				ve2.set(coords[3*n2], coords[3*n2+1], coords[3*n2+2]);
				ve2 -= vp;

				abutil::C3Vectord vn = ve1.cross(ve2);
				vn.normalize();
				vnorm += vn;
			}

			vnorm.normalize();
		}

		m_pAvgNormals[3*nId]	= vnorm.x;
		m_pAvgNormals[3*nId+1]	= vnorm.y;
		m_pAvgNormals[3*nId+2]	= vnorm.z;
	}

	m_centroidX /= (double) m_nVertices;
	m_centroidY /= (double) m_nVertices;
	m_centroidZ /= (double) m_nVertices;

	for (int i = 0; i < nElems; ++i)
	{
		const std::vector<int>& currFaceIndices = faceIndices[i];
		m_pIndices[3*i] = currFaceIndices[0];
		m_pIndices[3*i+1] = currFaceIndices[1];
		m_pIndices[3*i+2] = currFaceIndices[2];
	}

	return true;
}

void CFaceData::drawFaceWithTexCoords(double* pVertices, double* pTexCoord, double* pNormals)
{
	int i, iv, itc;

	if (pNormals == NULL)
		pNormals = m_pAvgNormals;

	glBegin(GL_TRIANGLES);

	glColor3f(0.8, 1.0, 0.8);

	for (i = 0; i < m_nIndices; i += 3)
	{
		iv = m_pIndices[i] * 3;
		itc = m_pIndices[i] * 2;
		glNormal3d(pNormals[iv], pNormals[iv + 1], pNormals[iv + 2]);
		glTexCoord2d(pTexCoord[itc], pTexCoord[itc + 1]);
		glVertex3d(pVertices[iv], pVertices[iv + 1], pVertices[iv + 2]);

		iv = m_pIndices[i + 1] * 3;
		itc = m_pIndices[i + 1] * 2;
		glNormal3d(pNormals[iv], pNormals[iv + 1], pNormals[iv + 2]);
		glTexCoord2d(pTexCoord[itc], pTexCoord[itc + 1]);
		glVertex3d(pVertices[iv], pVertices[iv + 1], pVertices[iv + 2]);

		iv = m_pIndices[i + 2] * 3;
		itc = m_pIndices[i + 2] * 2;
		glNormal3d(pNormals[iv], pNormals[iv + 1], pNormals[iv + 2]);
		glTexCoord2d(pTexCoord[itc], pTexCoord[itc + 1]);
		glVertex3d(pVertices[iv], pVertices[iv + 1], pVertices[iv + 2]);
	}

	glEnd();
}

void CFaceData::empericalPlaneFromLandmarkIndices(double* pNormal, int* pLandmarkIndices, int numIndices, double* pmatTransform)
{
	if (numIndices < m_nTemplateAlignMinIndices)
	{
		for (int i = 0; i < 16; i++)
			pmatTransform[i] = 0.0;
		return;
	}

	double distNoseTipOrigin;
	abutil::C3Vectord centroid, origin;
	abutil::C3Vectord outEyeLeft, outEyeRight, noseTip, noseBridge;
	abutil::C3Vectord xaxis, yaxis, zaxis, ytemp, ztemp;

	centroid.set(m_centroidX, m_centroidY, m_centroidZ);

	outEyeLeft.set(		m_pAvgVertices[pLandmarkIndices[m_iTemplateAlignOutsideLeftEyeCorner] * 3 + 0], 
						m_pAvgVertices[pLandmarkIndices[m_iTemplateAlignOutsideLeftEyeCorner] * 3 + 1], 
						m_pAvgVertices[pLandmarkIndices[m_iTemplateAlignOutsideLeftEyeCorner] * 3 + 2]);
	outEyeRight.set(	m_pAvgVertices[pLandmarkIndices[m_iTemplateAlignOutsideRightEyeCorner] * 3 + 0], 
						m_pAvgVertices[pLandmarkIndices[m_iTemplateAlignOutsideRightEyeCorner] * 3 + 1], 
						m_pAvgVertices[pLandmarkIndices[m_iTemplateAlignOutsideRightEyeCorner] * 3 + 2]);
	noseTip.set(		m_pAvgVertices[pLandmarkIndices[m_iTemplateAlignNoseTip] * 3 + 0], 
						m_pAvgVertices[pLandmarkIndices[m_iTemplateAlignNoseTip] * 3 + 1], 
						m_pAvgVertices[pLandmarkIndices[m_iTemplateAlignNoseTip] * 3 + 2]);
	noseBridge.set(		m_pAvgVertices[pLandmarkIndices[m_iTemplateAlignNoseBridge] * 3 + 0], 
						m_pAvgVertices[pLandmarkIndices[m_iTemplateAlignNoseBridge] * 3 + 1], 
						m_pAvgVertices[pLandmarkIndices[m_iTemplateAlignNoseBridge] * 3 + 2]);

	xaxis = outEyeLeft - outEyeRight;
	xaxis.normalize();
	yaxis = noseBridge - noseTip;
	yaxis.normalize();
	zaxis.set(pNormal[0], pNormal[1], pNormal[2]);
	zaxis.normalize();

	ytemp = zaxis.cross(xaxis);
	ytemp.normalize();

	ztemp = xaxis.cross(ytemp);
	ztemp.normalize();

	distNoseTipOrigin = (noseTip - centroid).length();
	origin = centroid;

	if (ytemp.dot(yaxis) > 0.5 && ztemp.dot(zaxis) > 0.5)
	{
		yaxis = ytemp;
		zaxis = ztemp;
	}
	else
	{
		std::cout << "Template alignment: axes messed up\n";
	}

	origin = noseTip - zaxis * distNoseTipOrigin;

	pmatTransform[0]	= xaxis.x;
	pmatTransform[1]	= xaxis.y;
	pmatTransform[2]	= xaxis.z;
	pmatTransform[3]	= -xaxis.dot(origin);
	pmatTransform[4]	= yaxis.x;
	pmatTransform[5]	= yaxis.y;
	pmatTransform[6]	= yaxis.z;
	pmatTransform[7]	= -yaxis.dot(origin);
	pmatTransform[8]	= zaxis.x;
	pmatTransform[9]	= zaxis.y;
	pmatTransform[10]	= zaxis.z;
	pmatTransform[11]	= -zaxis.dot(origin);
}

void CFaceData::stereographicProject(double *pVertices, int nVertices, double *pmatTransform, double *pVertStereo, double* pXmin, double* pYmin, double* pXmax, double* pYmax)
{
	int i;
	double xw, yw, zw, xf, yf, zf, xs, ys, rcpZ;
	double xmin, ymin, xmax, ymax;

	xmin = ymin = xmax = ymax = 0.0;

	//debug:
	double zmin = 1e30, zmax = -1e30;

	for (i = 0; i < nVertices; i++)
	{
		xw = pVertices[i * 3];
		yw = pVertices[i * 3 + 1];
		zw = pVertices[i * 3 + 2];

		xf = pmatTransform[0] * xw + pmatTransform[1] * yw + pmatTransform[2] * zw + pmatTransform[3];
		yf = pmatTransform[4] * xw + pmatTransform[5] * yw + pmatTransform[6] * zw + pmatTransform[7];
		zf = pmatTransform[8] * xw + pmatTransform[9] * yw + pmatTransform[10] * zw + pmatTransform[11];
		//debug:
		if (zf < zmin)
			zmin = zf;
		if (zf > zmax)
			zmax = zf;
		rcpZ = 1.0 / (zf + 1.0);
//		xs = xf * rcpZ;
//		ys = yf * rcpZ;
		xs = 2.0 * xf * rcpZ;
		ys = 2.0 * yf * rcpZ;

		pVertStereo[i * 2]		= xs;
		pVertStereo[i * 2 + 1]	= ys;

		if (xs < xmin)
			xmin = xs;
		if (xs > xmax)
			xmax = xs;
		if (ys < ymin)
			ymin = ys;
		if (ys > ymax)
			ymax = ys;
	}

	std::cout << "CFaceData::stereographicProject(...): z-range [" << zmin << ", " << zmax << "]\n";

	*pXmin = xmin;
	*pYmin = ymin;
	*pXmax = xmax;
	*pYmax = ymax;
}

///////////////////////////////////////////////////////////////////////////////
//file-scope function definitions
///////////////////////////////////////////////////////////////////////////////
void addFaceIndexHelper(const int vertexIndex, const int faceIndex, std::map<int, std::vector<int>>& vertexFaceIndices)
{
	std::map<int, std::vector<int>>::iterator findIter = vertexFaceIndices.find(vertexIndex);
	if(findIter != vertexFaceIndices.end())
	{
		std::vector<int>& currVertexFaceIndices = findIter->second;
		currVertexFaceIndices.push_back(faceIndex);
	}
	else
	{
		std::vector<int> currVertexFaceIndices;
		currVertexFaceIndices.push_back(faceIndex);
		vertexFaceIndices.insert(std::make_pair(vertexIndex, currVertexFaceIndices));
	}
}

void getFaceIndicesForVertices(const std::vector<std::vector<int>>& faceIndices, std::map<int, std::vector<int>>& faceIndicesForVertices)
{
	for(size_t faceId = 0; faceId < faceIndices.size(); ++faceId)
	{
		const std::vector<int>& currFaceIndices = faceIndices[faceId];

		for(size_t i = 0; i < currFaceIndices.size(); ++i)
		{
			addFaceIndexHelper(currFaceIndices[i], faceId, faceIndicesForVertices);
		}
	}
}

bool readOFF_geometry(const char *szFilename, std::vector<double>& coords, std::vector<std::vector<int>>& faceIndices)
{
	std::fstream fileStream;
	fileStream.open(szFilename, std::ios::in);

	if(!fileStream.is_open())
	{
		std::cout << "Unable to open file " << szFilename << std::endl;
		return false;
	}

	char tmp[256];
	fileStream.getline(tmp, 256);

	std::string sstrTmp(tmp);
	bool bColor = sstrTmp == "COFF";

	int numVertices(0);
	fileStream >> numVertices;

	int numFaces(0);
	fileStream >> numFaces;

	int numEdges(0);
	fileStream >> numEdges;

	for(int i = 0; i < numVertices; ++i)
	{
		double x(0.0),  y(0.0), z(0.0);
		fileStream >> x;
		fileStream >> y;
		fileStream >> z;

		coords.push_back(x);
		coords.push_back(y);
		coords.push_back(z);

		if(bColor)
		{
			//Color currently ignored
			double r(0.0), g(0.0), b(0.0), alpha(0.0);
			fileStream >> r;
			fileStream >> g;
			fileStream >> b;
			fileStream >> alpha;
		}

		if(fileStream.eof() || fileStream.fail() || fileStream.bad())
		{
			break;
		}
	}

	for(int i = 0; i < numFaces; ++i)
	{
		int n(0), a(0),  b(0), c(0);
		fileStream >> n;
		fileStream >> a;
		fileStream >> b;
		fileStream >> c;

		if(n!=3)
		{
			std::cout << "Only triangles supported";
			continue;
		}

		std::vector<int> currTriangle;
		currTriangle.push_back(a);
		currTriangle.push_back(b);
		currTriangle.push_back(c);

		faceIndices.push_back(currTriangle);

		if(fileStream.eof() || fileStream.fail() || fileStream.bad())
		{
			break;
		}
	}

	fileStream.close();
	return true;
}