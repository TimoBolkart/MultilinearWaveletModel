///////////////////////////////////////////////////////////////////////////////
//
//	Matrix.h
//
//	Header file for the C4x4Matrix and C3x3Matrix classes
//	updated to use templates
//
//	Alan Brunton 2009
//
///////////////////////////////////////////////////////////////////////////////


#ifndef __AB_MATRIX_H__
#define __AB_MATRIX_H__


#include "Vector.h"
#include "Quaternion.h"


///////////////////////////////////////////////////////////////////////////////
//C3x3Matrix
///////////////////////////////////////////////////////////////////////////////

template<typename scalar>
class C3x3Matrix
{
public:

	scalar						_11, _12, _13;
	scalar						_21, _22, _23;
	scalar						_31, _32, _33;

	__AB_CUDIFY__ C3x3Matrix()
	{
	}

	__AB_CUDIFY__ void setIdentity()
	{
		::memset(&_11, 0, sizeof(scalar) * 9);
		_11 = _22 = _33 = 1.0;
	}
};

typedef C3x3Matrix<float>	C3x3Matrixf;
typedef C3x3Matrix<double>	C3x3Matrixd;


///////////////////////////////////////////////////////////////////////////////
//C4x4Matrix
///////////////////////////////////////////////////////////////////////////////

template<typename scalar>
class C4x4Matrix
{
public:

	scalar						_11, _12, _13, _14;
	scalar						_21, _22, _23, _24;
	scalar						_31, _32, _33, _34;
	scalar						_41, _42, _43, _44;

	__AB_CUDIFY__ C4x4Matrix()
	{
#if defined(_AB_CUDA_)
#else
		memset(&_11, 0, sizeof(scalar) * 16);
#endif //defined(_AB_CUDA_)
	}
	__AB_CUDIFY__ C4x4Matrix(const C4x4Matrix& m)
	{
		copy(m);
	}

	__AB_CUDIFY__ void copy(const C4x4Matrix& m)
	{
//		_11 = m._11;	_12 = m._12;	_13 = m._13;	_14 = m._14;
//		_21 = m._21;	_22 = m._22;	_23 = m._23;	_24 = m._24;
//		_31 = m._31;	_32 = m._32;	_33 = m._33;	_24 = m._34;
//		_41 = m._41;	_42 = m._42;	_43 = m._43;	_44 = m._44;
		::memcpy(&_11, &m._11, sizeof(scalar) * 16);
	}

	__AB_CUDIFY__ scalar* asArray()
	{
		return &_11;
	}
	__AB_CUDIFY__ void toArray(scalar* pArray)
	{
		memcpy(pArray, &_11, sizeof(scalar) * 16);
	}

	__AB_CUDIFY__ void setIdentity()
	{
		::memset(&_11, 0, sizeof(scalar) * 16);
		_11 = _22 = _33 = _44 = 1.0;
	}

	__AB_CUDIFY__ C4x4Matrix& operator=(const C4x4Matrix& m)
	{
		copy(m);
		return *this;
	}

	__AB_CUDIFY__ C4x4Matrix& operator+=(const C4x4Matrix& m)
	{
		_11 += m._11;	_12 += m._12;	_13 += m._13;	_14 += m._14;
		_21 += m._21;	_22 += m._22;	_23 += m._23;	_24 += m._24;
		_31 += m._31;	_32 += m._32;	_33 += m._33;	_24 += m._34;
		_41 += m._41;	_42 += m._42;	_43 += m._43;	_44 += m._44;
		return *this;
	}
	__AB_CUDIFY__ C4x4Matrix operator+(const C4x4Matrix& m)
	{
		C4x4Matrix result(*this);
		return result;
	}

	__AB_CUDIFY__ C3Vector<scalar> operator*(const C3Vector<scalar>& v)
	{
		return C3Vector<scalar>(_11 * v.x + _12 * v.y + _13 * v.z + _14, 
								_21 * v.x + _22 * v.y + _23 * v.z + _24, 
								_31 * v.x + _32 * v.y + _33 * v.z + _34);
	}

	__AB_CUDIFY__ C4Vector<scalar> operator*(const C4Vector<scalar>& v)
	{
		return C4Vector<scalar>(_11 * v.x + _12 * v.y + _13 * v.z + _14 * v.w, 
								_21 * v.x + _22 * v.y + _23 * v.z + _24 * v.w, 
								_31 * v.x + _32 * v.y + _33 * v.z + _34 * v.w, 
								_41 * v.x + _42 * v.y + _43 * v.z + _44 * v.w);
	}

	__AB_CUDIFY__ C4x4Matrix operator*(const C4x4Matrix& m)
	{
		C4x4Matrix c;

		//first row
		c._11 = _11 * m._11 + _12 * m._21 + _13 * m._31 + _14 * m._41;
		c._12 = _11 * m._12 + _12 * m._22 + _13 * m._32 + _14 * m._42;
		c._13 = _11 * m._13 + _12 * m._23 + _13 * m._33 + _14 * m._43;
		c._14 = _11 * m._14 + _12 * m._24 + _13 * m._34 + _14 * m._44;

		//second row
		c._21 = _21 * m._11 + _22 * m._21 + _23 * m._31 + _24 * m._41;
		c._22 = _21 * m._12 + _22 * m._22 + _23 * m._32 + _24 * m._42;
		c._23 = _21 * m._13 + _22 * m._23 + _23 * m._33 + _24 * m._43;
		c._24 = _21 * m._14 + _22 * m._24 + _23 * m._34 + _24 * m._44;

		//third row
		c._31 = _31 * m._11 + _32 * m._21 + _33 * m._31 + _34 * m._41;
		c._32 = _31 * m._12 + _32 * m._22 + _33 * m._32 + _34 * m._42;
		c._33 = _31 * m._13 + _32 * m._23 + _33 * m._33 + _34 * m._43;
		c._34 = _31 * m._14 + _32 * m._24 + _33 * m._34 + _34 * m._44;

		//fourth row
		c._41 = _41 * m._11 + _42 * m._21 + _43 * m._31 + _44 * m._41;
		c._42 = _41 * m._12 + _42 * m._22 + _43 * m._32 + _44 * m._42;
		c._43 = _41 * m._13 + _42 * m._23 + _43 * m._33 + _44 * m._43;
		c._44 = _41 * m._14 + _42 * m._24 + _43 * m._34 + _44 * m._44;

		return c;
	}
	__AB_CUDIFY__ C4x4Matrix& operator*=(const C4x4Matrix& m)
	{
		copy( (*this) * m );
		return *this;
	}

	__AB_CUDIFY__ C4x4Matrix& transpose()
	{
		scalar s;
		s = _12;
		_12 = _21;
		_21 = s;
		
		s = _13;
		_13 = _31;
		_31 = s;

		s = _14;
		_14 = _41;
		_41 = s;

		s = _23;
		_23 = _32;
		_32 = s;
		
		s = _24;
		_24 = _42;
		_42 = s;

		s = _34;
		_34 = _43;
		_43 = s;

		return *this;
	}

	__AB_CUDIFY__ void setTranslation(const C3Vector<scalar>& v)
	{
		_14 = v.x;
		_24 = v.y;
		_34 = v.z;
	}
	__AB_CUDIFY__ void setTranslationMatrix(const C3Vector<scalar>& v)
	{
		setIdentity();
		_14 = v.x;
		_24 = v.y;
		_34 = v.z;
	}
	__AB_CUDIFY__ void setScaleMatrix(scalar s)
	{
		setIdentity();
		_11 = _22 = _33 = s;
	}
	__AB_CUDIFY__ void setScaleMatrix(const C3Vector<scalar>& v)
	{
		setIdentity();
		_11 = v.x;
		_22 = v.y;
		_33 = v.z;
	}
	__AB_CUDIFY__ void setRotationAboutX(scalar radians)
	{
		setIdentity();
		_22 = cos(radians);
		_23 = -sin(radians);
		_32 = -_23;
		_33 = _22;
	}
	__AB_CUDIFY__ void setRotationAboutY(scalar radians)
	{
		setIdentity();
		_11 = cos(radians);
		_13 = sin(radians);
		_31 = -_13;
		_33 = _11;
	}
	__AB_CUDIFY__ void setRotationAboutZ(scalar radians)
	{
		setIdentity();
		_11 = cos(radians);
		_12 = -sin(radians);
		_21 = -_12;
		_22 = _11;
	}

	__AB_CUDIFY__ void setRotation(const CQuaternion<scalar>& q)
	{
		scalar w2 = q.w*q.w, x2 = q.x*q.x, y2 = q.y*q.y, z2 = q.z*q.z;
		scalar wx = q.w * q.x, wy = q.w * q.y, wz = q.w * q.z;
		scalar xy = q.x * q.y, xz = q.x * q.z;
		scalar yz = q.y * q.z;

		_11 = w2 + x2 - y2 - z2;
		_12 = 2.0 * (xy - wz);
		_13 = 2.0 * (wy + xz);
		_21 = 2.0 * (wz + xy);
		_22 = w2 - x2 + y2 - z2;
		_23 = 2.0 * (yz - wx);
		_31 = 2.0 * (xz - wy);
		_32 = 2.0 * (wx + yz);
		_33 = w2 - x2 - y2 + z2;
	}
	__AB_CUDIFY__ void setRotationMatrix(const CQuaternion<scalar>& q)
	{
		setIdentity();
		setRotation(q);
	}

	__AB_CUDIFY__ void getRotation(CQuaternion<scalar>& q)
	{
#if 0
		scalar trace = _11 + _22 + _33 + 1.0;
		scalar sqTrace;
		scalar factor;

		if (trace > 0.0)
		{
			sqTrace = (scalar) ::sqrt((double) trace);
			factor = 0.5 / sqTrace;
			q.w = sqTrace * 0.5;
			q.x = (_32 - _23) * factor;
			q.y = (_13 - _31) * factor;
			q.z = (_21 - _12) * factor;
		}
		else
		{
			if (_11 > _22 && _11 > _33)
			{
				sqTrace = 2.0 * (scalar) ::sqrt((double) (1.0 + _11 - _22 - _33));
				factor = 1.0 / sqTrace;
				q.w = (_23 - _32) * factor;
				q.x = sqTrace * 0.25;
				q.y = (_12 + _21) * factor;
				q.z = (_13 + _31) * factor;
			}
			else if (_22 > _33)
			{
				sqTrace = 2.0 * (scalar) ::sqrt((double) (1.0 + _22 - _11 - _33));
				factor = 1.0 / sqTrace;
				q.w = (_13 - _31) * factor;
				q.x = (_12 + _21) * factor;
				q.y = sqTrace * 0.25;
				q.z = (_23 + _32) * factor;
			}
			else
			{
				sqTrace = 2.0 * (scalar) ::sqrt((double) (1.0 + _33 - _11 - _22));
				factor = 1.0 / sqTrace;
				q.w = (_12 - _21) * factor;
				q.x = (_13 + _31) * factor;
				q.y = (_23 + _32) * factor;
				q.z = sqTrace * 0.25;
			}
		}
#else
		int i, imax;
		scalar temp[4];
		scalar mag = (scalar)-1;
		scalar factor;
		
		temp[0] = 1.0 + _11 + _22 + _33;
		temp[1] = 1.0 + _11 - _22 - _33;
		temp[2] = 1.0 - _11 + _22 - _33;
		temp[3] = 1.0 - _11 - _22 + _33;

		for (i = 0; i < 4; i++)
		{
			if (temp[i] > mag)
			{
				mag = temp[i];
				imax = i;
			}
		}

		switch (imax)
		{
		case 0:
			q.w = sqrt(temp[0]) * 0.5;
			factor = 1.0 / (4.0 * q.w);
			q.x = (_32 - _23) * factor;
			q.y = (_13 - _31) * factor;
			q.z = (_21 - _12) * factor;
			break;
		case 1:
			q.x = sqrt(temp[1]) * 0.5;
			factor = 1.0 / (4.0 * q.x);
			q.w = (_32 - _23) * factor;
			q.y = (_21 + _12) * factor;
			q.z = (_13 + _31) * factor;
			break;
		case 2:
			q.y = sqrt(temp[2]) * 0.5;
			factor = 1.0 / (4.0 * q.y);
			q.w = (_13 - _31) * factor;
			q.x = (_21 + _12) * factor;
			q.z = (_32 + _23) * factor;
			break;
		case 3:
			q.z = sqrt(temp[3]) * 0.5;
			factor = 1.0 / (4.0 * q.z);
			q.w = (_21 - _12) * factor;
			q.x = (_13 + _31) * factor;
			q.y = (_32 + _23) * factor;
			break;
		}
#endif

		q.normalize();
	}

	__AB_CUDIFY__ C4Vector<scalar> column(int nCol)
	{
		switch (nCol)
		{
		case 1:
			return C4Vector<scalar>(_11, _21, _31, _41);
		case 2:
			return C4Vector<scalar>(_12, _22, _32, _42);
		case 3:
			return C4Vector<scalar>(_13, _23, _33, _43);
		case 4:
			return C4Vector<scalar>(_14, _24, _34, _44);
		}
		return C4Vector<scalar>(0.0, 0.0, 0.0, 0.0);
	}
	__AB_CUDIFY__ C4Vector<scalar> row(int nRow)
	{
		switch (nRow)
		{
		case 1:
			return C4Vector<scalar>(_11, _12, _13, _14);
		case 2:
			return C4Vector<scalar>(_21, _22, _23, _24);
		case 3:
			return C4Vector<scalar>(_31, _32, _33, _34);
		case 4:
			return C4Vector<scalar>(_41, _42, _43, _44);
		}
		return C4Vector<scalar>(0.0, 0.0, 0.0, 0.0);
	}

	__AB_CUDIFY__ void getNormalTransformNoScale(C4x4Matrix& matNorm)
	{
		matNorm = *this;
		matNorm._14 = matNorm._24 = matNorm._34 = 0.0;
	}
	__AB_CUDIFY__ void getNormalTransformUniScale(C4x4Matrix& matNorm)
	{
		C3Vector<scalar> c1(_11, _21, _31), c1o;
		C3Vector<scalar> c2(_12, _22, _32), c2o;
		C3Vector<scalar> c3(_13, _23, _33), c3o;
		matNorm.setIdentity();
		c2.crossFast(c3, c1o);
		c3.crossFast(c1, c2o);
		c1.crossFast(c2, c3o);
		c1o.normalize();
		c2o.normalize();
		c3o.normalize();
		matNorm._11 = c1o.x;	matNorm._12 = c2o.x;	matNorm._13 = c3o.x;
		matNorm._21 = c1o.y;	matNorm._22 = c2o.y;	matNorm._23 = c3o.y;
		matNorm._31 = c1o.z;	matNorm._32 = c2o.z;	matNorm._33 = c3o.z;
	}
	__AB_CUDIFY__ void getNormalTransformGeneral(C4x4Matrix& matNorm)
	{
		C3Vector<scalar> c1(_11, _21, _31), c1o;
		C3Vector<scalar> c2(_12, _22, _32), c2o;
		C3Vector<scalar> c3(_13, _23, _33), c3o;
		matNorm.setIdentity();
		c2.crossFast(c3, c1o);
		c3.crossFast(c1, c2o);
		c1.crossFast(c2, c3o);
		c1o.normalize();
		c2o.normalize();
		c3o.normalize();
		c1o /= c1.length();
		c2o /= c2.length();
		c3o /= c3.length();
		matNorm._11 = c1o.x;	matNorm._12 = c2o.x;	matNorm._13 = c3o.x;
		matNorm._21 = c1o.y;	matNorm._22 = c2o.y;	matNorm._23 = c3o.y;
		matNorm._31 = c1o.z;	matNorm._32 = c2o.z;	matNorm._33 = c3o.z;
	}

	__AB_CUDIFY__ void setLookAt(C3Vector<scalar>& eye, C3Vector<scalar>& at, C3Vector<scalar>& up)
	{
		C3Vector<scalar> zaxis = eye - at;
		zaxis.normalize();

		C3Vector<scalar> xaxis = up.cross(zaxis);
		xaxis.normalize();

		C3Vector<scalar> yaxis = zaxis.cross(xaxis);
		yaxis.normalize();

		setIdentity();

		_11 = xaxis.x;	_12 = xaxis.y;	_13 = xaxis.z;	_14 = -xaxis.dot(eye);
		_21 = yaxis.x;	_22 = yaxis.y;	_23 = yaxis.z;	_24 = -yaxis.dot(eye);
		_31 = zaxis.x;	_32 = zaxis.y;	_33 = zaxis.z;	_34 = -zaxis.dot(eye);
	}

	__AB_CUDIFY__ void setPerspectiveProjection(scalar width, scalar height, scalar centerX, scalar centerY, scalar nearZ, scalar farZ)
	{
		scalar rcpDiff;
		scalar rcpW, rcpH;

		rcpDiff = 1.0 / (nearZ - farZ);
		rcpW = 1.0 / width;
		rcpH = 1.0 / height;

		setIdentity();
		_11 = 2.0 * nearZ * rcpW;
		_13 = 2.0 * centerX * rcpW;
		_22 = 2.0 * nearZ * rcpH;
		_23 = 2.0 * centerY * rcpH;
		_33 = farZ * rcpDiff;
		_34 = nearZ * farZ * rcpDiff;
		_43 = -1.0;
		_44 = 0.0;
	}

	__AB_CUDIFY__ void setOrthographicProjection(scalar width, scalar height, scalar centerX, scalar centerY, scalar nearZ, scalar farZ)
	{
		scalar rcpDiff;
		scalar rcpW, rcpH;

		rcpDiff = 1.0 / (nearZ - farZ);
		rcpW = 1.0 / width;
		rcpH = 1.0 / height;

		setIdentity();
		_11 = 2.0 * rcpW;
		_14 = 2.0 * centerX * -rcpW;
		_22 = 2.0 * rcpH;
		_24 = 2.0 * centerY * -rcpH;
		_33 = rcpDiff;
		_34 = nearZ * rcpDiff;
	}

	__AB_CUDIFY__ void setPerspectiveProjection(scalar fieldOfView, scalar aspectRatio, scalar nearZ, scalar farZ)
	{
		scalar width, height;

		width = 2.0 * nearZ * ::tan(fieldOfView * 0.5);
		height = width / aspectRatio;
		setPerspectiveProjection(width, height, 0.0, 0.0, nearZ, farZ);
	}

	void print()
	{
		::printf("%f %f %f %f\n", _11, _12, _13, _14);
		::printf("%f %f %f %f\n", _21, _22, _23, _24);
		::printf("%f %f %f %f\n", _31, _32, _33, _34);
		::printf("%f %f %f %f\n", _41, _42, _43, _44);
	}
};

typedef C4x4Matrix<float>	C4x4Matrixf;
typedef C4x4Matrix<double>	C4x4Matrixd;


#endif //__AB_MATRIX_H__


