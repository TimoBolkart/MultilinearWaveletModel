#ifndef VECTORNX_H
#define VECTORNX_H

#include <math.h>
#include <assert.h>

#define VEC_PI       3.14159265358979323846

template<typename T, size_t DIM>
class VecNX
{
public:
	VecNX()
	{
		for(size_t i = 0; i < DIM; ++i)
		{
			m_values[i] = static_cast<T>(0.0);
		}
	}

	VecNX(T x, T y)
	{
		assert(DIM==2);
		m_values[0] = x;
		m_values[1] = y;
	}

	VecNX(T x, T y, T z)
	{
		assert(DIM==3);
		m_values[0] = x;
		m_values[1] = y;
		m_values[2] = z;
	}

	~VecNX()
	{}

	int getDim() const { return DIM; }

	T& operator[](size_t i)
	{
		assert(i<DIM);
		return m_values[i];
	}

	const T& operator[](size_t i) const
	{
		assert(i<DIM);
		return m_values[i];
	}

	bool normalize()
	{
		const T norm = length();
		if(norm == static_cast<T>(0.0))
		{
			return false;
		}
		
		scalarMult(1.0 / norm);
		return true;
	}

	T sqrLength() const
	{
		return dotProduct(*this);
	}

	T length() const
	{
		return sqrt(sqrLength());	
	}

	void scalarMult(T mult)
	{
		for(size_t i = 0; i < DIM; ++i)
		{
			m_values[i]*=mult;
		}
	}

	T dotProduct(const VecNX<T,DIM>& vec) const
	{
		T val = static_cast<T>(0.0);
		for(size_t i = 0; i < DIM; ++i)
		{
			val += m_values[i]*vec[i];
		}

		return val;
	}

	double angle(const VecNX<T,DIM>& vec) const
	{
		T product = dotProduct(vec);
		T lengthA = length();
		T lengthB = vec.length();
		const double angle = acos(static_cast<double>(product/(lengthA*lengthB)));
		return angle*180.0/VEC_PI;
	}

	void crossProduct(const VecNX<T,3>& vec, VecNX<T,3>& outVec) const
	{
		assert(DIM==3);
		outVec[0] = m_values[1]*vec[2]-m_values[2]*vec[1];
		outVec[1] = m_values[2]*vec[0]-m_values[0]*vec[2];
		outVec[2] = m_values[0]*vec[1]-m_values[1]*vec[0];
	}

private:
	T m_values[DIM];
};

typedef VecNX<double,2> Vec2d;
typedef VecNX<double,3> Vec3d;

typedef VecNX<int,2> Vec2i;
typedef VecNX<int,3> Vec3i;
typedef VecNX<int,4> Vec4i;

#endif