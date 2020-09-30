#include "FileLoader.h"

#include <fstream>
#include <string>

//#ifdef DEBUG_OUTPUT
#include <iostream>
//#endif


#include <string>

bool check_fstream_error(std::fstream& strm, std::string strMsg)
{
	if (strm.bad())
	{
		std::cout << "file error: " << strMsg << std::endl;
		return false;
	}
	return true;
}

std::string FileLoader::getFileExtension(const std::string& sstrFileName)
{
	if(sstrFileName.empty())
	{
		return "";
	}

	const size_t pos = sstrFileName.rfind(".");
	if(pos != std::string::npos)
	{
		return sstrFileName.substr(pos+1);
	}

	return "";
}

bool FileLoader::fileExist(const std::string& sstrFileName)
{
	if(sstrFileName.empty())
	{
		return false;
	}

	std::fstream inStream;
	inStream.open(sstrFileName, std::ios::in);

	if(inStream.is_open())
	{
		inStream.close();
		return true;
	}
	
	return false;
}

std::string FileLoader::getFileName(const std::string& sstrFullFileName)
{
	std::string sstrFileName = "";

	if(sstrFullFileName.empty())
	{
		return sstrFileName;
	}

	if(!FileLoader::fileExist(sstrFullFileName))
	{
		return sstrFileName;
	}

	size_t posEnd = sstrFullFileName.rfind(".");
	if(posEnd == std::string::npos)
	{
		posEnd = sstrFullFileName.size();
	}

	size_t posStart = sstrFullFileName.rfind("/");
	if(posStart != std::string::npos)
	{
		sstrFileName = sstrFullFileName.substr(posStart+1, posEnd-posStart-1);
		return sstrFileName;
	}

	posStart = sstrFullFileName.rfind("\\");
	if(posStart != std::string::npos)
	{
		sstrFileName = sstrFullFileName.substr(posStart+1, posEnd-posStart-1);
		return sstrFileName;
	}

	return sstrFileName;
}

std::string FileLoader::getFilePath(const std::string& sstrFileName)
{
	std::string sstrFilePath = "";

	if(sstrFileName.empty())
	{
		return sstrFilePath;
	}

	if(!FileLoader::fileExist(sstrFileName))
	{
		return sstrFilePath;
	}

	size_t pos = sstrFileName.rfind("/");
	if(pos != std::string::npos)
	{
		sstrFilePath = sstrFileName.substr(0, pos);
		return sstrFilePath;
	}

	pos = sstrFileName.rfind("\\");
	if(pos != std::string::npos)
	{
		sstrFilePath = sstrFileName.substr(0, pos);
		return sstrFilePath;
	}

	return sstrFilePath;
}

bool FileLoader::loadFile(const std::string& fileName, DataContainer& outData)
{
	if(fileName.empty())
	{
		return false;
	}

	const std::string suffix = FileLoader::getFileExtension(fileName);
	if(suffix==std::string("off"))
	{
		return loadOFF(fileName, outData);
	}

	return false;
}

bool FileLoader::loadLandmarks(const std::string& sstrFileName, std::vector<double>& landmarks, std::vector<bool>& loaded)
{
	if(!FileLoader::fileExist(sstrFileName))
	{
		return false;
	}

	std::string sstrSuffix = FileLoader::getFileExtension(sstrFileName);
	if(sstrSuffix=="txt")
	{
		return loadSimpleLandmarks(sstrFileName, landmarks, loaded);
	}

	return false;
}

bool FileLoader::loadMultilinearModelBinary(std::fstream& input, std::vector<size_t>& modeDims, std::vector<size_t>& truncModeDims, std::vector<double>& multModel
															, std::vector<double>& sVectors, std::vector<double>& mean, std::vector<double>& meanWeights)
{
	if (!input.good() || !input.is_open())
		return false;

	size_t numModes, modelSize;

	input.read((char*)&numModes, sizeof(size_t));
	if (!check_fstream_error(input, "FileLoader::loadMultilinearModelBindary(...): read numModes"))			return false;
	//std::cout << "FileLoader::loadMultilinearModelBindary(...): read numModes: " << numModes << "\n"; std::cout.flush();

	modeDims.resize(numModes);
	truncModeDims.resize(numModes);

	input.read((char*)&modeDims[0], numModes * sizeof(size_t));
	if (!check_fstream_error(input, "FileLoader::loadMultilinearModelBindary(...): read modeDims"))			return false;

	input.read((char*)&truncModeDims[0], numModes * sizeof(size_t));
	if (!check_fstream_error(input, "FileLoader::loadMultilinearModelBindary(...): read truncModeDims"))	return false;

	input.read((char*)&modelSize, sizeof(size_t));
	if (!check_fstream_error(input, "FileLoader::loadMultilinearModelBindary(...): read modelSize"))		return false;
	//std::cout << "FileLoader::loadMultilinearModelBindary(...): read modelSize: " << modelSize << "\n"; std::cout.flush();

	multModel.resize(modelSize);
	input.read((char*)&multModel[0], modelSize * sizeof(double));
	if (!check_fstream_error(input, "FileLoader::loadMultilinearModelBindary(...): read multModel"))		return false;

	size_t /*uMatSize = 0, */sVecSize = 0, meanWeightVecSize = 0;
	for (size_t i = 1; i < numModes; i++)
	{
		const size_t numEntries = modeDims[i]*truncModeDims[i];
		sVecSize += truncModeDims[i];
		meanWeightVecSize += truncModeDims[i];
	}

	sVectors.resize(sVecSize);

	for(size_t i = 1; i < numModes; ++i)
	{
		const size_t numEntries = truncModeDims[i];
		size_t startIndex = 0;
		for(size_t j = 1; j < i; ++j)
		{
			startIndex += truncModeDims[j];
		}

		input.read((char*)&sVectors[startIndex], numEntries * sizeof(double));
		if (!check_fstream_error(input, "FileLoader::loadMultilinearModelBindary(...): read SVector"))		return false;
	}

	mean.resize(modeDims[0]);
	input.read((char*)&mean[0], modeDims[0] * sizeof(double));
	if (!check_fstream_error(input, "FileLoader::loadMultilinearModelBindary(...): read mean"))				return false;

	meanWeights.resize(meanWeightVecSize);

	for(size_t i = 1; i < numModes; ++i)
	{
		const size_t numEntries = truncModeDims[i];
		size_t startIndex = 0;
		for(size_t j = 1; j < i; ++j)
		{
			startIndex += truncModeDims[j];
		}

		input.read((char*)&meanWeights[startIndex], numEntries * sizeof(double));
		if (!check_fstream_error(input, "FileLoader::loadMultilinearModelBindary(...): read meanWeights"))		return false;
	}

	return true;
}

bool FileLoader::loadMultilinearModelText(std::fstream& input, std::vector<size_t>& modeDims, std::vector<size_t>& truncModeDims, std::vector<double>& multModel
															, std::vector<double>& sVectors, std::vector<double>& mean, std::vector<double>& meanWeights)
{
	if (!input.good() || !input.is_open())
		return false;

	size_t numModes, modelSize;

	//input.read((char*)&numModes, sizeof(size_t));
	input >> numModes;
	if (!check_fstream_error(input, "FileLoader::loadMultilinearModelBindary(...): read numModes"))			return false;
	//std::cout << "FileLoader::loadMultilinearModelBindary(...): read numModes: " << numModes << "\n"; std::cout.flush();

	modeDims.resize(numModes);
	truncModeDims.resize(numModes);

	//input.read((char*)&modeDims[0], numModes * sizeof(size_t));
	for(size_t i = 0; i < numModes; ++i)
	{
		input >> modeDims[i];
		if (!check_fstream_error(input, "FileLoader::loadMultilinearModelBindary(...): read modeDims"))			return false;
	}

	//input.read((char*)&truncModeDims[0], numModes * sizeof(size_t));
	for(size_t i = 0; i < numModes; ++i)
	{
		input >> truncModeDims[i];
		if (!check_fstream_error(input, "FileLoader::loadMultilinearModelBindary(...): read truncModeDims"))	return false;
	}

	//input.read((char*)&modelSize, sizeof(size_t));
	input >> modelSize;
	if (!check_fstream_error(input, "FileLoader::loadMultilinearModelBindary(...): read modelSize"))		return false;
	//std::cout << "FileLoader::loadMultilinearModelBindary(...): read modelSize: " << modelSize << "\n"; std::cout.flush();

	multModel.resize(modelSize);
	//input.read((char*)&multModel[0], modelSize * sizeof(double));
	for(size_t i = 0; i < modelSize; ++i)
	{
		input >> multModel[i];
		if (!check_fstream_error(input, "FileLoader::loadMultilinearModelBindary(...): read multModel"))		return false;
	}

	size_t /*uMatSize = 0, */sVecSize = 0, meanWeightVecSize = 0;
	for (size_t i = 1; i < numModes; i++)
	{
		const size_t numEntries = modeDims[i]*truncModeDims[i];
		sVecSize += truncModeDims[i];
		meanWeightVecSize += truncModeDims[i];
	}

	sVectors.resize(sVecSize);

	for(size_t i = 1; i < numModes; ++i)
	{
		const size_t numEntries = truncModeDims[i];
		size_t startIndex = 0;
		for(size_t j = 1; j < i; ++j)
		{
			startIndex += truncModeDims[j];
		}

		//input.read((char*)&sVectors[startIndex], numEntries * sizeof(double));
		for(size_t j = 0; j < numEntries; ++j)
		{
			input >> sVectors[startIndex+j];
			if (!check_fstream_error(input, "FileLoader::loadMultilinearModelBindary(...): read SVector"))		return false;
		}
	}

	mean.resize(modeDims[0]);
	//input.read((char*)&mean[0], modeDims[0] * sizeof(double));
	for(size_t i = 0; i < modeDims[0]; ++i)
	{
		input >> mean[i];
		if (!check_fstream_error(input, "FileLoader::loadMultilinearModelBindary(...): read mean"))				return false;
	}

	meanWeights.resize(meanWeightVecSize);

	for(size_t i = 1; i < numModes; ++i)
	{
		const size_t numEntries = truncModeDims[i];
		size_t startIndex = 0;
		for(size_t j = 1; j < i; ++j)
		{
			startIndex += truncModeDims[j];
		}

		//input.read((char*)&meanWeights[startIndex], numEntries * sizeof(double));
		for(size_t j = 0; j < numEntries; ++j)
		{
			input >> meanWeights[startIndex+j];
			if (!check_fstream_error(input, "FileLoader::loadMultilinearModelBindary(...): read meanWeights"))		return false;
		}
	}

	return true;
}

bool FileLoader::saveFile(const std::string& sstrFileName, const DataContainer& data)
{
	if(sstrFileName.empty())
	{
		return false;
	}

	const std::string suffix = FileLoader::getFileExtension(sstrFileName);
	if(suffix==std::string("off"))
	{
		return writeOff(sstrFileName, data);
	}

	return false;
}

bool FileLoader::loadSimpleLandmarks(const std::string& sstrFileName, std::vector<double>& landmarks, std::vector<bool>& loaded)
{
	const char* cstrFileName = sstrFileName.c_str();
	if(cstrFileName==NULL)
	{
		return false;
	}

	FILE* pFile(NULL);
	errno_t err = fopen_s(&pFile, cstrFileName, "r");
	if(pFile==NULL || err!=0)
	{
		return false;
	}

	std::vector<double> tmpLandmarks;
	if(!loadDataFile(sstrFileName, tmpLandmarks))
	{
		std::cout << "Loading simple landmarks failed " << sstrFileName << std::endl;
		return false;
	}

	if(tmpLandmarks.size() % 3 != 0)
	{
		std::cout << "Wrong landmarks dimension" << std::endl;
		return false;
	}

	const size_t numLandmarks = tmpLandmarks.size()/3;

	loaded.clear();
	loaded.resize(numLandmarks);
	for(int i = 0; i < numLandmarks; ++i)
	{
		loaded[i] = true;
	}

	landmarks = tmpLandmarks;

	return fclose(pFile)==0;
}

bool FileLoader::loadOFF(const std::string& sstrFileName, DataContainer& outData)
{
	const char* cstrFileName = sstrFileName.c_str();
	if(cstrFileName==NULL)
	{
		return false;
	}

	FILE* pFile(NULL);
	errno_t err = fopen_s(&pFile, cstrFileName, "r");
	if(pFile==NULL || err!=0)
	{
		return false;
	}

	char output[100];
	readNextNode(pFile, output);

	bool bColorOff(false);
	if(strcmp(output, "OFF") == 0)
	{
		bColorOff = false;
	}
	else if(strcmp(output, "COFF") == 0)
	{
		bColorOff = true;
	}
	else
	{
		return false;
	}

	int numVertices(0);
	readNextNumber(pFile, INTEGER, output, numVertices);

	int numFaces(0);
	readNextNumber(pFile, INTEGER, output, numFaces);

	int numEdges(0);
	readNextNumber(pFile, INTEGER, output, numEdges);

	if(numVertices < 1)// || numFaces < 1)
	{
		return false;
	}

	std::vector<Vec3d*>& vertexList = outData.getVertexList();
	std::vector<Vec3i*>& vertexIndexList = outData.getVertexIndexList();
	std::vector<Vec3d*>& vertexColors = outData.getVertexColorList();


	bool bEnd(false);
	for(int vertex = 0; vertex < numVertices; ++vertex)
	{
		char tmpOutput[100];
	
		Vec3d* pVertex = new Vec3d(0.0,0.0,0.0);
		for(int i = 0; i < 3; ++i)
		{
			double number(0.0);		
			if(!readNextNumber(pFile, FLOATING_POINT, tmpOutput, number))
			{
				bEnd = true;
				break;
			}

			(*pVertex)[i] = number;
		}

		if(bEnd)
		{
			delete pVertex;
			return false;
		}
		else
		{
			vertexList.push_back(pVertex);
		}

		if(bColorOff)
		{
			Vec3d* pColor = new Vec3d(0.0,0.0,0.0);
			for(int i = 0; i < 3; ++i)
			{
				double colorValue(0);		
				if(!readNextNumber(pFile, FLOATING_POINT, tmpOutput, colorValue))
				{
					bEnd = true;
					break;
				}

				(*pColor)[i] = colorValue;
			}

			if(bEnd)
			{
				delete pColor;
			return false;
			}
			else
			{
				vertexColors.push_back(pColor);
			}

			double alphaValue(0);		
			if(!readNextNumber(pFile, FLOATING_POINT, tmpOutput, alphaValue))
			{
				return false;
			}
		}
	}

	for(int face = 0; face < numFaces; ++face)
	{
		char tmpOutput[100];

		int numPolyPoints(0);		
		if(!readNextNumber(pFile, INTEGER, tmpOutput, numPolyPoints))
		{
			break;
		}

		if(numPolyPoints!=3)
		{
#ifdef DEBUG_OUTPUT
			std::cout << "loadOff() - only triangles supported" << std::endl; 
#endif
			return false;
		}

		Vec3i* pPoly = new Vec3i(0,0,0);
		for(int i = 0; i < 3; ++i)
		{
			int vertexIndex(0);		
			if(!readNextNumber(pFile, INTEGER, tmpOutput, vertexIndex))
			{
				bEnd = true;
				break;
			}

			(*pPoly)[i] = vertexIndex;
		}

		if(bEnd)
		{
			delete pPoly;
			break;
		}
		else
		{
			vertexIndexList.push_back(pPoly);
		}
	}

	return fclose(pFile)==0;
}

bool FileLoader::writeOff(const std::string& sstrFileName, const DataContainer& data)
{
	const size_t numVertices = data.getNumVertices();
	const size_t numFaces = data.getVertexIndexList().size();
	const size_t numEdges = 0;

	const char* cstrFileName = sstrFileName.c_str();
	if(cstrFileName==NULL)
	{
		return false;
	}

	std::fstream output;
	output.open(cstrFileName, std::ios::out);
	output.precision(7);
	output.setf(std::ios::fixed,std::ios::floatfield);

	const std::vector<Vec3d*>& vertices = data.getVertexList();
	const std::vector<Vec3d*>& vertexColors = data.getVertexColorList();

	bool bValidColors = vertices.size() == vertexColors.size();

	if(bValidColors)
	{
		output << "COFF" << std::endl;
	}
	else
	{
		output << "OFF" << std::endl;
	}

	output << numVertices << " " << numFaces << " " << numEdges << std::endl;

	for(size_t i = 0; i < numVertices; ++i)
	{
		const Vec3d& coord = *(vertices[i]);
		output << coord[0] << " " << coord[1] << " " << coord[2] << " ";
		if(bValidColors)
		{
			const Vec3d& color = *(vertexColors[i]);
			output << color[0] << " " << color[1] << " " << color[2] << " " << "1.0";
		}

		output << std::endl;
	}

	const std::vector<Vec3i*>& vertexIndexList = data.getVertexIndexList();
	std::vector<Vec3i*>::const_iterator currFaceIter = vertexIndexList.begin();
	std::vector<Vec3i*>::const_iterator endFaceIter = vertexIndexList.end();
	for( ; currFaceIter != endFaceIter; ++currFaceIter)
	{
		const Vec3i* pFace = *currFaceIter;
		output << "3" << " " << (*pFace)[0] << " " << (*pFace)[1]  << " " << (*pFace)[2]  << std::endl;
	}

	output.close();
	return true;
}

bool FileLoader::loadDataFile(const std::string& sstrDataFileName, std::vector<double>& data)
{
	const char* cstrFileName = sstrDataFileName.c_str();
	if(cstrFileName==NULL)
	{
		return false;
	}

	FILE* pFile(NULL);
	errno_t err = fopen_s(&pFile, cstrFileName, "r");
	if(pFile==NULL || err!=0)
	{
		return false;
	}

	while(true)
	{
		double val(0.0);
		char strOutput[1000];
		if(!readNextNumber(pFile, FLOATING_POINT, strOutput, val))
		{
			break;
		}

		data.push_back(val);
	}

	return fclose(pFile)==0;
}

bool FileLoader::readNextNode(FILE* pFile, char* cstrOutput)
{
	//Read string of characters until a whitespace (blank, newline, tab) is found.
	int ret = fscanf(pFile, "%100s", cstrOutput);

	//Ignore comments marked by # at the beginning of the line
	while(ret!=EOF && cstrOutput[0]=='#')
	{
		//Read until the end of the line is reached
		do
		{
			ret=fgetc(pFile);
		}
		while(ret!=EOF && ret!='\n');

		ret=fscanf(pFile, "%100s", cstrOutput);
	}

	// separate "{", "}", "[" and "]"
	size_t len = strlen(cstrOutput); // returns length of string s
	if (len > 1)
	{
		char currChar = cstrOutput[len-1];

		if(currChar == '{' || currChar == '}'
			|| currChar == '[' || currChar == ']')
		{
			ungetc(cstrOutput[len-1], pFile);
			cstrOutput[len-1] = 0;
		}
	}

	return ret==EOF;
}

bool FileLoader::processShapeNode(FILE* pFile, char* cstrOutput, DataContainer& outPoly)
{
	// Vertex data
	std::vector<Vec3d*>& vertexList = outPoly.getVertexList();
	std::vector<Vec3i*>& vertexIndexList = outPoly.getVertexIndexList();

	// Texture data
	std::vector<Vec2d*>& textureList = outPoly.getTextureList();
	std::vector<Vec3i*>& textureIndexList = outPoly.getTextureIndexList();

	//const std::string& textureName = outPoly.getTextureName();

	bool bEnd(false);
	while(!bEnd)
	{
		bEnd = readNextNode(pFile, cstrOutput);

		if(bEnd)
		{
			break;
		}

		if(!bEnd && isEqual(cstrOutput, "ImageTexture"))
		{
			char cstrTextureName[200];
			bEnd &= processImageTexture(pFile, cstrOutput, cstrTextureName);
			
			std::string sstrTextureName(cstrTextureName);
			if(!sstrTextureName.empty())
			{
				outPoly.setTextureName(sstrTextureName);
			}
		}
		
		if(!bEnd && isEqual(cstrOutput, "coord"))
		{
			processCoordinates(pFile, cstrOutput, vertexList);
		}
		
		if(!bEnd && isEqual(cstrOutput, "coordIndex"))
		{
			processIndices(pFile, cstrOutput, vertexIndexList);
		}
		
		if(!bEnd && isEqual(cstrOutput, "TextureCoordinate"))
		{
			processCoordinates(pFile, cstrOutput, textureList);
		}
		
		if(!bEnd && isEqual(cstrOutput, "texCoordIndex"))
		{
			processIndices(pFile, cstrOutput, textureIndexList);
		}
	}

	return bEnd;
}

bool FileLoader::processImageTexture(FILE* pFile, char* cstrOutput, char* cstrTextureName)
{
	bool bEnd = readNodeBlocksUntil(pFile, "url", cstrOutput);

	assert(!bEnd);
	assert(!isEqual(cstrOutput, "}"));
	assert(isEqual(cstrOutput, "url"));

	bEnd = readNextNode(pFile, cstrOutput);
	
	if(isEqual(cstrOutput, "["))
	{
		bEnd = readNextNode(pFile, cstrOutput);
	}

	// Remove "/" and "\" from the beginning of the name
	int nFirstLetter(0);
	for(size_t i = 0; i < 100; ++i)
	{
		if(cstrOutput[i] != '\"' && cstrOutput[i] != '/')
		{
			break;
		}

		++nFirstLetter;
	}

	strcpy(cstrTextureName, cstrOutput+nFirstLetter);

	return bEnd;
}

void FileLoader::processBlockStructure(FILE* pFile, char* strOutput, std::vector<double>& values)
{
	bool bEnd = readNodeBlocksUntil(pFile, "{", strOutput);;
	while(!bEnd)
	{
		double number(0.0);
		if(!readNextNumber(pFile, FLOATING_POINT, strOutput, number))
		{
			break;
		}

		values.push_back(number);
	}
}