#ifndef DATACONTAINER_H
#define DATACONTAINER_H

#include <vector>
#include "VectorNX.h"

//CoordsVec			vector type of coordinates
//						e.g. Vec2d for 2d double vectors, Vec3d for 3d double vectors,...
//TexVec				vector type of texture coordinates
//						e.g. Vec2d for 2d double vectors, Vec3d for 3d double vectors,...
//IndexDim			number of indices which define a face
//						e.g. 3 for triangle, 4 for quad,...

template<typename CoordsVec, typename TexVec, size_t IndexDim>
class ImportMesh
{
	//Indices always integer values
	typedef VecNX<int, IndexDim> IndexVec;

public:
	ImportMesh()
	: m_textureName("")
	{

	}

	ImportMesh(const ImportMesh& mesh)
	{
		m_textureName = mesh.getTextureName();
		copyPtrArray(mesh.getVertexList(), m_vertexList);
		copyPtrArray(mesh.getVertexIndexList(), m_vertexIndexList);
		copyPtrArray(mesh.getVertexNormalList(), m_vertexNormalList);
		copyPtrArray(mesh.getVertexColorList(), m_vertexColorList);
		copyPtrArray(mesh.getTextureList(), m_textureList);
		copyPtrArray(mesh.getTextureIndexList(), m_textureIndexList);
	}

	~ImportMesh()
	{
		clear();
	}

	ImportMesh& operator=(const ImportMesh& mesh)
	{
		clear();

		m_textureName = mesh.getTextureName();
		copyPtrArray(mesh.getVertexList(), m_vertexList);
		copyPtrArray(mesh.getVertexIndexList(), m_vertexIndexList);
		copyPtrArray(mesh.getVertexNormalList(), m_vertexNormalList);
		copyPtrArray(mesh.getVertexColorList(), m_vertexColorList);
		copyPtrArray(mesh.getTextureList(), m_textureList);
		copyPtrArray(mesh.getTextureIndexList(), m_textureIndexList);

		return *this;
	}

	void clear()
	{
		clearPtrArray(m_vertexList);
		clearPtrArray(m_vertexIndexList);
		clearPtrArray(m_vertexNormalList);
		clearPtrArray(m_vertexColorList);
		clearPtrArray(m_textureList);
		clearPtrArray(m_textureIndexList);
	}

	size_t getNumVertices()
	{
		return m_vertexList.size();
	}

	size_t getNumVertices() const
	{
		return m_vertexList.size();
	}

	std::vector<CoordsVec*>& getVertexList()
	{
		return m_vertexList;
	}

	const std::vector<CoordsVec*>& getVertexList() const
	{
		return m_vertexList;
	}

	std::vector<IndexVec*>& getVertexIndexList()
	{
		return m_vertexIndexList;
	}

	const std::vector<IndexVec*>& getVertexIndexList() const
	{
		return m_vertexIndexList;
	}

	void clearVertexColorList()
	{
		clearPtrArray(m_vertexColorList);
	}

	std::vector<Vec3d*>& getVertexColorList()
	{
		return m_vertexColorList;
	}

	void clearVertexNormalList()
	{
		clearPtrArray(m_vertexNormalList);
	}

	std::vector<Vec3d*>& getVertexNormalList()
	{
		return m_vertexNormalList;
	}

	const std::vector<Vec3d*>& getVertexNormalList() const
	{
		return m_vertexNormalList;
	}

	const std::vector<Vec3d*>& getVertexColorList() const
	{
		return m_vertexColorList;
	}

	std::vector<TexVec*>& getTextureList()
	{
		return m_textureList;
	}

	const std::vector<TexVec*>& getTextureList() const
	{
		return m_textureList;
	}

	std::vector<IndexVec*>& getTextureIndexList()
	{
		return m_textureIndexList;
	}

	const std::vector<IndexVec*>& getTextureIndexList() const
	{
		return m_textureIndexList;
	}

	void setTextureName(std::string sstrTextureName)
	{
		m_textureName = sstrTextureName;
	}

	const std::string& getTextureName()
	{
		return m_textureName;
	}

	const std::string& getTextureName() const 
	{
		return m_textureName;
	}

private:
	template<typename ArrayType>
	void clearPtrArray(std::vector<ArrayType*>& vec)
	{
		std::vector<ArrayType*>::iterator currIter = vec.begin();
		std::vector<ArrayType*>::iterator endIter = vec.end();
		for(; currIter!=endIter; ++currIter)
		{
			if(*currIter != NULL)
			{
				delete *currIter;
			}
		}

		vec.clear();
	}

	template<typename ArrayType>
	void copyPtrArray(const std::vector<ArrayType*>& inVec, std::vector<ArrayType*>& outVec)
	{
		outVec.reserve(inVec.size());

		std::vector<ArrayType*>::const_iterator currIter = inVec.begin();
		std::vector<ArrayType*>::const_iterator endIter = inVec.end();
		for(; currIter!=endIter; ++currIter)
		{
			ArrayType* pNewType = new ArrayType(**currIter);
			outVec.push_back(pNewType);
		}
	}

	//ImportMesh(const ImportMesh&);

	//ImportMesh &operator=(const ImportMesh&);	

	// Vertex data
	std::vector<CoordsVec*> m_vertexList;
	std::vector<IndexVec*> m_vertexIndexList;
	std::vector<Vec3d*> m_vertexNormalList;
	std::vector<Vec3d*> m_vertexColorList;

	// Texture data
	std::vector<TexVec*> m_textureList;
	std::vector<IndexVec*> m_textureIndexList;

	std::string m_textureName;
};

typedef ImportMesh<Vec3d,Vec2d,3> DataContainer;

#endif