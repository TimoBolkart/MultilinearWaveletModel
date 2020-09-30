#ifndef FILELOADER_H
#define FILELOADER_H

#include "DataContainer.h"


class FileLoader
{
	enum Type
	{
		FLOATING_POINT,
		INTEGER
	};

public:
	//! Get the file extension of a full file name.
	//! e.g. ...\...\Test.txt returns txt
	//! \param sstrFile			full file name
	//! \return file extension
	static std::string getFileExtension(const std::string& sstrFileName);

	//! Checks if a file name exists.
	//! \param sstrFile			full file name
	//! \return true if file exists
	static bool fileExist(const std::string& sstrFileName);

	//! Get the file name of a full file name.
	//! e.g. ...\...\Test.txt returns Test
	//! \param sstrFile			full file name
	//! \return file name if file exists
	static std::string getFileName(const std::string& sstrFullFileName);

	//! Get the file name of a full file name.
	//! e.g. ...\...\Test.txt returns ...\...
	//! \param sstrFile			full file name
	//! \return file path if file exists
	static std::string getFilePath(const std::string& sstrFileName);

	
	bool loadFile(const std::string& fileName, DataContainer& outData);

	//! Load landmarks. Extend this to load other landmark file formats.
	//! \param sstrFileName			full landmark file name
	//! \param landmarks				vertices of loaded landmarks
	//! \param loaded					flags that indicate which vertices are successfully loaded (could be used in the presence of occlusions)
	//! \return true if successful
	bool loadLandmarks(const std::string& sstrFileName, std::vector<double>& landmarks, std::vector<bool>& loaded);

	//! Load the multilinear model.
	bool loadMultilinearModelBinary(std::fstream& input, std::vector<size_t>& modeDims, std::vector<size_t>& truncModeDims, std::vector<double>& multModel
												  , std::vector<double>& sVectors, std::vector<double>& mean, std::vector<double>& meanWeights);

	bool loadMultilinearModelText(std::fstream& input, std::vector<size_t>& modeDims, std::vector<size_t>& truncModeDims, std::vector<double>& multModel
												  , std::vector<double>& sVectors, std::vector<double>& mean, std::vector<double>& meanWeights);

	//! Save geometry file.
	bool saveFile(const std::string& sstrFileName, const DataContainer& data);

private:
	//! Specific loader for off files. Use similar structure to consider other file formats.	
	bool loadOFF(const std::string& sstrFileName, DataContainer& outData);
	
	//! Load ordered set of landmarks from file.
	//! \param sstrFileName
	//! \param landmarks
	//! \param loaded
	//! \return true if successful
	bool loadSimpleLandmarks(const std::string& sstrFileName, std::vector<double>& landmarks, std::vector<bool>& loaded);

	//! Specific writer for off files.	
	bool writeOff(const std::string& sstrFileName, const DataContainer& data);

	//! Load file of double values.
	//! \param sstrDataFileName	full file name
	//! \param data					loaded data
	//! \return true if successful
	bool loadDataFile(const std::string& sstrDataFileName, std::vector<double>& data);

	bool readNextNode(FILE* pFile, char* cstrOutput);

	template<typename Unit>
	bool readNextNumber(FILE* pFile, const Type type, char* strOutput, Unit& number)
	{
		bool bEnd = readNextNode(pFile, strOutput);
		if(bEnd)
		{
			return false;
		}

		return convertToNumber(strOutput, type, number);
	}

	bool processShapeNode(FILE* pFile, char* cstrOutput, DataContainer& outPoly);

	bool processImageTexture(FILE* pFile, char* cstrOutput, char* cstrTextureName);

	void processBlockStructure(FILE* pFile, char* strOutput, std::vector<double>& values);

	template<typename Unit>
	bool convertToNumber(const char* strIn, const Type& type, Unit& out)
	{
		const char* cstrFormat = type==FLOATING_POINT ? "%lf%*c" : "%ld%*c";
		return sscanf(strIn, cstrFormat, &out)==1; 
	}

	template<typename Unit, size_t DIM>
	void processCoordinates(FILE* pFile, char* cstrOutput, std::vector<VecNX<Unit, DIM>*>& vertexList)
	{
		bool bEnd = readNodeBlocksUntil(pFile, "[", cstrOutput);
		assert(!bEnd && !isEqual(cstrOutput, "}"));

		do
		{
			VecNX<double, DIM>* pVec = new VecNX<double, DIM>();
			if(!readDataBlock(pFile, FLOATING_POINT, cstrOutput, *pVec))
			{
				delete pVec;
				break;
			}
			
			vertexList.push_back(pVec);
		}
		while(!isEqual(cstrOutput, "}"));
	}

	template<typename Unit, size_t DIM>
	void processIndices(FILE* pFile, char* /*cstrOutput*/, std::vector<VecNX<Unit, DIM>*>& indexList)
	{
		char endChar[3] = { '[', ']', '}' };

		char tmpOutput[100];
		bool bEnd = readBlockUntil(pFile, endChar, tmpOutput);
		assert(!bEnd && !isEqual(tmpOutput, "}"));

		while(!bEnd)
		{
			std::vector<int> values;
			bEnd = readIndexBlock(pFile, values);
			if(bEnd || values.size() != DIM)
			{
				break;
			}

			VecNX<Unit, DIM>* pVec = new VecNX<Unit, DIM>();
			for(size_t i = 0; i < DIM; ++i)
			{
				(*pVec)[i] = values[i];
			}

			indexList.push_back(pVec);
		}
	}

	template<typename Unit, size_t DIM>
	bool readDataBlock(FILE* pFile, const Type& type, char* cstrOutput, VecNX<Unit, DIM>& vec)
	{
		for(size_t i = 0; i < DIM; ++i)
		{
			Unit value = static_cast<Unit>(0.0);
			if(readNextNode(pFile, cstrOutput) || !convertToNumber(cstrOutput, type, value)
				|| isEqual(cstrOutput, "}") || isEqual(cstrOutput, "]"))
			{
				return false;
			}

			vec[i] = value;
		}

		return true;
	}

	bool readIndexBlock(FILE* pFile, std::vector<int>& values)
	{
		bool bEnd(false);
		while(!bEnd)
		{
			char cstrTmpOutput[20] = {0};

			bEnd = getIndexNumberString(pFile, cstrTmpOutput);

			size_t len = strlen(cstrTmpOutput);
			if(len>0)
			{
				int value(0);
				if(!convertToNumber(cstrTmpOutput, INTEGER, value))
				{
					return true;
				}

				if(value==-1)
				{
					return false;
				}

				values.push_back(value);
			}
			
			if(bEnd)
			{
				return true;
			}
		}

		return bEnd;
	}

	bool getIndexNumberString(FILE* pFile, char* cstrOutput)
	{
		char cstrTmpOutput[20] = {0};

		bool bEnd(false);
		while(!bEnd)
		{
			size_t len = strlen(cstrTmpOutput);

			char currChar = fgetc(pFile);
			if(currChar == ',' || currChar == ' ' || currChar == '\n')
			{
				if(len == 0)
				{
					continue;
				}
				else
				{
					break;
				}
			}
			else if(currChar == '}' || currChar == ']' || currChar == EOF)
			{
				bEnd = true;
				break;
			}

			cstrTmpOutput[len] = currChar;
		}

		size_t len = strlen(cstrTmpOutput);
		for(size_t i = 0; i < len; ++i)
		{
			cstrOutput[i] = cstrTmpOutput[i];
		}

		return bEnd;
	}

	//! Reads next block until it reaches the specified character.
	//! Stops reading if the end of file is reached.
	//! Stops reading if "}" or "]" is reached.
	bool readNodeBlocksUntil(FILE* pFile, const char* cstrEnd, char* cstrOutput)
	{
		bool bEnd(false);
		do
		{
			bEnd = readNextNode(pFile, cstrOutput);
		}
		while(!bEnd && !isEqual(cstrOutput, cstrEnd) 
				&& !isEqual(cstrOutput, "}") && !isEqual(cstrOutput, "]"));

		return bEnd;
	}

	bool readBlockUntil(FILE* pFile, const char* cstrEnd, char* cstrOutput)
	{
		const size_t numStopChar = strlen(cstrEnd);
		if(numStopChar == 0)
		{
			return false;
		}

		char cstrTmpOutput[100] = {0};

		bool bEnd(false);

		bool bStop(false);
		while(!bStop)
		{
			const size_t len = strlen(cstrTmpOutput);

			char currChar = fgetc(pFile);
			for(size_t i = 0; i < numStopChar; ++i)
			{
				if(currChar == cstrEnd[i])
				{
					bStop = true;
					break;
				}
			}

			if(currChar!= EOF)
			{
				cstrTmpOutput[len] = currChar;
			}
			else
			{
				bEnd = true;
			}
		}

		size_t len = strlen(cstrTmpOutput);
		for(size_t i = 0; i < len; ++i)
		{
			cstrOutput[i] = cstrTmpOutput[i];
		}

		return bEnd;
	}

	bool isEqual(const char* cstrS1, const char* cstrS2)
	{
		return strcmp(cstrS1, cstrS2) == 0;
	}
};

#endif