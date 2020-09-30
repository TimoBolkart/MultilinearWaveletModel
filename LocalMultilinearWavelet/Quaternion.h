///////////////////////////////////////////////////////////////////////////////
//
//	Quaternion.h
//
//	Header file for the CQuaternion class
//	updated to use template
//
//	Alan Brunton 2009
//
///////////////////////////////////////////////////////////////////////////////


#ifndef __AB_QUATERNION_H__
#define __AB_QUATERNION_H__


#include "Vector.h"


template<typename scalar>
class CQuaternion
{
public:
	 
	scalar								w, x, y, z;

	__AB_CUDIFY__ CQuaternion(): w(1.0), x(0.0), y(0.0), z(0.0)
	{
	}
	__AB_CUDIFY__ CQuaternion(scalar sw, scalar sx, scalar sy, scalar sz): w(sw), x(sx), y(sy), z(sz)
	{
	}
	__AB_CUDIFY__ CQuaternion(const C3Vector<scalar>& v): w(0.0), x(v.x), y(v.y), z(v.z)
	{
	}
	__AB_CUDIFY__ CQuaternion(const CQuaternion& q): w(q.w), x(q.x), y(q.y), z(q.z)
	{
	}
	__AB_CUDIFY__ CQuaternion(scalar radians, C3Vector<scalar>& axis)
	{
		scalar half = radians * 0.5;
		scalar c = (scalar) ::cos((double) half);
		scalar s = (scalar) ::sin((double) half);
		scalar rcpA = 1.0 / axis.length();
		w = c;
		x = axis.x * rcpA * s;
		y = axis.y * rcpA * s;
		z = axis.z * rcpA * s;
	}
	__AB_CUDIFY__ CQuaternion(scalar azimuth, scalar altitude)
	{
		CQuaternion<scalar> q1(azimuth, C3Vector<scalar>(0.0, 1.0, 0.0));
		CQuaternion<scalar> q2(altitude, C3Vector<scalar>(1.0, 0.0, 0.0));
		*this = q1 * q2;
	}
	__AB_CUDIFY__ CQuaternion(scalar* pq): w(pq[0]), x(pq[1]), y(pq[2]), z(pq[3])
	{
	}

	__AB_CUDIFY__ CQuaternion& set(scalar qw, scalar qx, scalar qy, scalar qz)
	{
		w = qw;
		x = qx;
		y = qy;
		z = qz;
		return *this;
	}
	__AB_CUDIFY__ CQuaternion& set(scalar* pq)
	{
		w = pq[0];
		x = pq[1];
		y = pq[2];
		z = pq[3];
		return *this;
	}

	__AB_CUDIFY__ void get(scalar* pq)
	{
		pq[0] = w;
		pq[1] = x;
		pq[2] = y;
		pq[3] = z;
	}

	__AB_CUDIFY__ C3Vector<scalar> vectorPart() { return C3Vector<scalar>(x, y, z); }

	__AB_CUDIFY__ CQuaternion& operator=(const CQuaternion& q)
	{
		w = q.w;
		x = q.x;
		y = q.y;
		z = q.z;
		return *this;
	}

	__AB_CUDIFY__ CQuaternion& operator+=(const CQuaternion& q)
	{
		w += q.w;
		x += q.x;
		y += q.y;
		z += q.z;
		return *this;
	}
	__AB_CUDIFY__ CQuaternion operator+(const CQuaternion& q)
	{
		CQuaternion result(*this);
		result += q;
		return result;
	}

	__AB_CUDIFY__ CQuaternion& operator*=(scalar s)
	{
		w *= s;
		x *= s;
		y *= s;
		z *= s;
		return *this;
	}
	__AB_CUDIFY__ CQuaternion operator*(scalar s)
	{
		CQuaternion q(*this);
		q *= s;
		return q;
	}

	__AB_CUDIFY__ void mulGrassman(const CQuaternion& q, CQuaternion& qresult)
	{
		//these two pieces of code are equivalent
#if 1
		qresult.w = w * q.w - (x * q.x + y * q.y + z * q.z);
		qresult.x = w * q.x + q.w * x + y * q.z - z * q.y;
		qresult.y = w * q.y + q.w * y + z * q.x - x * q.z;
		qresult.z = w * q.z + q.w * z + x * q.y - y * q.x;
#else
		C3Vector<scalar> v1(x, y, z), v2(q.x, q.y, q.z), vresult;
		qresult.w = w * q.w - v1.dot(v2);
		vresult = v2 * w + v1 * q.w + v1.cross(v2);
		qresult.x = vresult.x;
		qresult.y = vresult.y;
		qresult.z = vresult.z;
#endif
	}
	__AB_CUDIFY__ CQuaternion operator*(const CQuaternion& q)
	{
		CQuaternion qresult;
		mulGrassman(q, qresult);
		return qresult;
	}
	__AB_CUDIFY__ CQuaternion& operator*=(const CQuaternion& q)
	{
		CQuaternion qresult(*this);
		qresult.mulGrassman(q, *this);
		return *this;
	}
	__AB_CUDIFY__ CQuaternion operator-()
	{
		CQuaternion qres(*this);
		qres.w = -qres.w;
		qres.x = -qres.x;
		qres.y = -qres.y;
		qres.z = -qres.z;
		return qres;
	}

	__AB_CUDIFY__ scalar dot(const CQuaternion& q) { return w * q.w + x * q.x + y * q.y + z * q.z; }
	__AB_CUDIFY__ scalar lengthSquared() { return w*w + x*x + y*y + z*z; }
	__AB_CUDIFY__ scalar length() { return (scalar) sqrt(w*w + x*x + y*y + z*z); }

	__AB_CUDIFY__ CQuaternion& normalize()
	{
		scalar scale = 1.0 / length();
		w *= scale;
		x *= scale;
		y *= scale;
		z *= scale;
		return *this;
	}

	__AB_CUDIFY__ CQuaternion conjugate() { return CQuaternion(w, -x, -y, -z); }
	__AB_CUDIFY__ CQuaternion inverse()
	{
		scalar rcp2 = 1.0 / lengthSquared();
		return conjugate() *= rcp2;
	}

	__AB_CUDIFY__ CQuaternion& slerp(CQuaternion& q1, CQuaternion& q2, scalar t)
	{
		const scalar eps = (scalar)1e-10;
		const scalar one_minus_eps = 1.0 - eps;

		scalar c;
		scalar angle, rcpS, f1, f2;
		CQuaternion q1p(q1), q2p(q2);				//q1' and q2' are mutable versions of q1 and q2

		c = q1.dot(q2);								//cosine of the angle between them
		c /= (q1.length() * q2.length());

		if (c > one_minus_eps)
		{
			//angle ~0, use linear interpolation
			q1p *= (1.0 - t);
			q2p *= t;
			*this = q1p += q2p;
		}
		else
		{
			if (c < 0.0)
			{
				q2p = -q2p;
				c = q1p.dot(q2p);						//cosine of the angle between them
				c /= (q1p.length() * q2p.length());
			}

			angle = acos(c);
			rcpS = 1.0 / sin(angle);

			f1 = rcpS * sin((1.0 - t) * angle);
			f2 = rcpS * sin(t * angle);

			q1p *= f1;
			q2p *= f2;
			*this = q1p += q2p;
		}

		return *this;
	}
};

typedef CQuaternion<float>		CQuaternionf;
typedef CQuaternion<double>		CQuaterniond;


#endif //__AB_QUATERNION_H__


