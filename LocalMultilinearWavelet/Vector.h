///////////////////////////////////////////////////////////////////////////////
//
//	Vector.h
//
//	Header file for the C3Vector and C4Vector classes
//	updated to use templates
//
//	Alan Brunton 2009
//
///////////////////////////////////////////////////////////////////////////////


#ifndef __AB_VECTOR_H__
#define __AB_VECTOR_H__


///////////////////////////////////////////////////////////////////////////////
//C3Vector
///////////////////////////////////////////////////////////////////////////////
template<typename scalar>
class C3Vector
{
public:

	scalar			x, y, z;

	__AB_CUDIFY__ C3Vector(): x(0.0), y(0.0), z(0.0) { }
	__AB_CUDIFY__ C3Vector(scalar sx, scalar sy, scalar sz): x(sx), y(sy), z(sz) { }
	__AB_CUDIFY__ C3Vector(scalar sx, scalar sy): x(sx), y(sy), z(0.0) { }
	__AB_CUDIFY__ C3Vector(scalar s): x(s), y(s), z(s) { }
	__AB_CUDIFY__ C3Vector(const C3Vector& v): x(v.x), y(v.y), z(v.z) { }
	__AB_CUDIFY__ C3Vector(scalar* pv): x(pv[0]), y(pv[1]), z(pv[2]) { }

	__AB_CUDIFY__ C3Vector& set(scalar sx, scalar sy, scalar sz)
	{
		x = sx;
		y = sy;
		z = sz;
		return *this;
	}
	__AB_CUDIFY__ C3Vector& set(scalar* pv)
	{
		x = pv[0];
		y = pv[1];
		z = pv[2];
		return *this;
	}

	__AB_CUDIFY__ C3Vector& operator=(const C3Vector& v)
	{
		x = v.x;
		y = v.y;
		z = v.z;
		return *this;
	}
	__AB_CUDIFY__ C3Vector& operator*=(const scalar s)
	{
		x *= s;
		y *= s;
		z *= s;
		return *this;
	}
	__AB_CUDIFY__ C3Vector& operator/=(const scalar s)
	{
		scalar rcpS = 0.0;
//		if (fabs(s) > g_epsilon)
		if (fabs(s) > (scalar)ABUTIL_EPSILON)
			rcpS = 1.0 / s;
		x *= rcpS;
		y *= rcpS;
		z *= rcpS;
		return *this;
	}
	__AB_CUDIFY__ C3Vector& operator+=(const C3Vector& v)
	{
		x += v.x;
		y += v.y;
		z += v.z;
		return *this;
	}
	__AB_CUDIFY__ C3Vector& operator-=(const C3Vector& v)
	{
		x -= v.x;
		y -= v.y;
		z -= v.z;
		return *this;
	}

	__AB_CUDIFY__ C3Vector operator*(const scalar s) const
	{
		C3Vector v(x, y, z);
		v *= s;
		return v;
	}
	__AB_CUDIFY__ C3Vector operator/(const scalar s) const
	{
		C3Vector v(x, y, z);
		v /= s;
		return v;
	}
	__AB_CUDIFY__ C3Vector operator+(const C3Vector& v) const
	{
		C3Vector result(x, y, z);
		result += v;
		return result;
	}
	__AB_CUDIFY__ C3Vector operator-(const C3Vector& v) const
	{
		C3Vector result(x, y, z);
		result -= v;
		return result;
	}
	__AB_CUDIFY__ C3Vector operator-() const
	{
		return C3Vector(-x, -y, -z);
	}

	__AB_CUDIFY__ scalar dot(const C3Vector& v) const
	{
		return x*v.x + y*v.y + z*v.z;
	}
	__AB_CUDIFY__ scalar lengthSquared() const
	{
		return dot(*this);
	}
	__AB_CUDIFY__ scalar length() const
	{
		return (scalar) ::sqrt(dot(*this));
	}

	__AB_CUDIFY__ C3Vector& normalize()
	{
		scalar rcpMag = 1.0 / length();
		x *= rcpMag;
		y *= rcpMag;
		z *= rcpMag;
		return *this;
	}

	__AB_CUDIFY__ C3Vector cross(const C3Vector& v) const
	{
		C3Vector result;
		result.x = y * v.z - z * v.y;
		result.y = z * v.x - x * v.z;
		result.z = x * v.y - y * v.x;
		return result;
	}
	__AB_CUDIFY__ void crossFast(const C3Vector& v, C3Vector& result) const
	{
		result.x = y * v.z - z * v.y;
		result.y = z * v.x - x * v.z;
		result.z = x * v.y - y * v.x;
	}

	__AB_CUDIFY__ bool operator==(const C3Vector& v)
	{
		C3Vector vdiff(x, y, z);
		vdiff -= v;
		return ::fabs(vdiff.x) < g_epsilon && ::fabs(vdiff.y) < g_epsilon && ::fabs(vdiff.z) < g_epsilon;
	}
	__AB_CUDIFY__ bool operator!=(const C3Vector& v)
	{
		C3Vector vdiff(x, y, z);
		vdiff -= v;
		return fabs(vdiff.x) >= g_epsilon || fabs(vdiff.y) >= g_epsilon || fabs(vdiff.z) >= g_epsilon;
	}
};

typedef C3Vector<float>		C3Vectorf;
typedef C3Vector<double>	C3Vectord;


///////////////////////////////////////////////////////////////////////////////
//C4Vector
///////////////////////////////////////////////////////////////////////////////

template<typename scalar>
class C4Vector
{
public:

	scalar			x, y, z, w;

	__AB_CUDIFY__ C4Vector(): x(0.0), y(0.0), z(0.0), w(1.0) { }
	__AB_CUDIFY__ C4Vector(scalar sx, scalar sy, scalar sz, scalar sw): x(sx), y(sy), z(sz), w(sw) { }
	__AB_CUDIFY__ C4Vector(scalar sx, scalar sy, scalar sz): x(sx), y(sy), z(sz), w(1.0) { }
	__AB_CUDIFY__ C4Vector(scalar sx, scalar sy): x(sx), y(sy), z(0.0), w(1.0) { }
	__AB_CUDIFY__ C4Vector(const C4Vector<scalar>& v): x(v.x), y(v.y), z(v.z), w(v.w) { }
	__AB_CUDIFY__ C4Vector(const C3Vector<scalar>& v): x(v.x), y(v.y), z(v.z), w(1.0) { }

	__AB_CUDIFY__ C4Vector& set(scalar sx, scalar sy, scalar sz, scalar sw)
	{
		x = sx;
		y = sy;
		z = sz;
		w = sw;
		return *this;
	}
	__AB_CUDIFY__ C4Vector& set(scalar sx, scalar sy, scalar sz)
	{
		x = sx;
		y = sy;
		z = sz;
		return *this;
	}

	__AB_CUDIFY__ C4Vector& operator=(const C4Vector& v)
	{
		x = v.x;
		y = v.y;
		z = v.z;
		w = v.w;
		return *this;
	}
	__AB_CUDIFY__ C4Vector& operator=(const C3Vector<scalar>& v)
	{
		x = v.x;
		y = v.y;
		z = v.z;
		w = 1.0;
		return *this;
	}
	__AB_CUDIFY__ C4Vector& operator*=(const scalar s)
	{
		x *= s;
		y *= s;
		z *= s;
		w *= s;
		return *this;
	}
	__AB_CUDIFY__ C4Vector& operator/=(const scalar s)
	{
		scalar rcpS = 1.0 / s;
		x *= rcpS;
		y *= rcpS;
		z *= rcpS;
		w *= rcpS;
		return *this;
	}
	__AB_CUDIFY__ C4Vector& operator+=(const C4Vector& v)
	{
		x += v.x;
		y += v.y;
		z += v.z;
		w += v.w;
		return *this;
	}
	__AB_CUDIFY__ C4Vector& operator-=(const C4Vector& v)
	{
		x -= v.x;
		y -= v.y;
		z -= v.z;
		w -= v.w;
		return *this;
	}

	__AB_CUDIFY__ C4Vector operator*(const scalar s)
	{
		C4Vector v(*this);
		v *= s;
		return v;
	}
	__AB_CUDIFY__ C4Vector operator/(const scalar s)
	{
		C4Vector v(*this);
		v /= s;
		return v;
	}
	__AB_CUDIFY__ C4Vector operator+(const C4Vector& v)
	{
		C4Vector result(*this);
		result += v;
		return result;
	}
	__AB_CUDIFY__ C4Vector operator-(const C4Vector& v)
	{
		C4Vector result(*this);
		result -= v;
		return result;
	}
	__AB_CUDIFY__ C4Vector operator-()
	{
		return C4Vector(-x, -y, -z, -w);
	}

	__AB_CUDIFY__ scalar dot(const C4Vector& v)
	{
		return x*v.x + y*v.y + z*v.z + w*v.w;
	}
	__AB_CUDIFY__ scalar lengthSquared()
	{
		return dot(*this);
	}
	__AB_CUDIFY__ scalar length()
	{
		return ::sqrt(dot(*this));
	}

	__AB_CUDIFY__ C4Vector& normalize()
	{
		scalar rcpMag = 1.0 / length();
		x *= rcpMag;
		y *= rcpMag;
		z *= rcpMag;
		w *= rcpMag;
		return *this;
	}

	__AB_CUDIFY__ void unitW()
	{
		scalar rcpW = 0.0;
		if (::abs(w) > g_epsilon)
			rcpW = 1.0 / w;
		x *= rcpW;
		y *= rcpW;
		z *= rcpW;
		w *= rcpW;
	}

	__AB_CUDIFY__ bool operator==(const C4Vector& v)
	{
		C4Vector vdiff(*this);
		vdiff -= v;
		return ::abs(vdiff.x) < g_epsilon && ::abs(vdiff.y) < g_epsilon && ::abs(vdiff.z) < g_epsilon && ::abs(vdiff.w) < g_epsilon;
	}
};

typedef C4Vector<float>		C4Vectorf;
typedef C4Vector<double>	C4Vectord;


#endif //__AB_VECTOR_H__


