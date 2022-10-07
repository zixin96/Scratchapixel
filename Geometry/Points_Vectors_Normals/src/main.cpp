#include <iostream>
using namespace std;

template <typename T>
class Vec3
{
public:
	// 3 most basic ways of initializing a vector

	//!? Notice T(0) here
	Vec3() : x(T(0)), y(T(0)), z(T(0))
	{
	}

	Vec3(const T& xx) : x(xx), y(xx), z(xx)
	{
	}

	Vec3(T xx, T yy, T zz) : x(xx), y(yy), z(zz)
	{
	}

	// length can be a method from the class...
	T length()
	{
		return sqrt(x * x + y * y + z * z);
	}

	// as a method of the class Vec3
	Vec3<T>& normalize()
	{
		T len = length();

		//!? could use dot product to compute len2
		// T len2 = dot(*this);

		if (len > 0)
		{
			T invLen = 1 / len;
			x *= invLen, y *= invLen, z *= invLen;
		}

		return *this;
	}

	T dot(const Vec3<T>& v) const
	{
		return x * v.x + y * v.y + z * v.z;
	}

	// as a method of the class...
	Vec3<T> cross(const Vec3<T>& v) const
	{
		return Vec3<T>(
		               y * v.z - z * v.y,
		               z * v.x - x * v.z,
		               x * v.y - y * v.x);
	}

	Vec3<T> operator +(const Vec3<T>& v) const
	{
		return Vec3<T>(x + v.x, y + v.y, z + v.z);
	}

	Vec3<T> operator -(const Vec3<T>& v) const
	{
		return Vec3<T>(x - v.x, y - v.y, z - v.z);
	}

	Vec3<T> operator *(const T& r) const
	{
		return Vec3<T>(x * r, y * r, z * r);
	}

	T x, y, z;
};

// ... or you can also compute the length in a function which is not part of the class
template <typename T>
T length(const Vec3<T>& v)
{
	return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

// or as a utility function
template <typename T>
void normalize(Vec3<T>& v)
{
	T len2 = v.x * v.x + v.y * v.y + v.z * v.z;
	// avoid division by 0
	if (len2 > 0)
	{
		T invLen = 1 / sqrt(len2);
		v.x *= invLen, v.y *= invLen, v.z *= invLen;
	}
}

template <typename T>
T dot(const Vec3<T>& a, const Vec3<T>& b)
{
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

// or as an utility function
template <typename T>
Vec3<T> cross(const Vec3<T>& a, const Vec3<T>& b)
{
	return Vec3<T>(
	               a.y * b.z - a.z * b.y,
	               a.z * b.x - a.x * b.z,
	               a.x * b.y - a.y * b.x);
}

typedef Vec3<float> Vec3f;


template <typename T>
class Matrix44
{
public:
	Matrix44()
	{
	}

	const T* operator [](uint8_t i) const { return m[i]; }

	T* operator [](uint8_t i) { return m[i]; }

	Matrix44 operator *(const Matrix44& rhs) const
	{
		Matrix44 mult;
		for (uint8_t i = 0; i < 4; ++i)
		{
			for (uint8_t j = 0; j < 4; ++j)
			{
				mult[i][j] = m[i][0] * rhs[0][j] +
				             m[i][1] * rhs[1][j] +
				             m[i][2] * rhs[2][j] +
				             m[i][3] * rhs[3][j];
			}
		}

		return mult;
	}

	//!? transforming points
	void multVecMatrix(const Vec3<T>& src, Vec3<T>& dst) const
	{
		dst.x = src.x * m[0][0] + src.y * m[1][0] + src.z * m[2][0] + m[3][0];
		dst.y = src.x * m[0][1] + src.y * m[1][1] + src.z * m[2][1] + m[3][1];
		dst.z = src.x * m[0][2] + src.y * m[1][2] + src.z * m[2][2] + m[3][2];
		T w   = src.x * m[0][3] + src.y * m[1][3] + src.z * m[2][3] + m[3][3];
		if (w != 1 && w != 0)
		{
			dst.x /= w;
			dst.y /= w;
			dst.z /= w;
		}
	}

	//!? transforming vectors
	void multDirMatrix(const Vec3<T>& src, Vec3<T>& dst) const
	{
		dst.x = src.x * m[0][0] + src.y * m[1][0] + src.z * m[2][0];
		dst.y = src.x * m[0][1] + src.y * m[1][1] + src.z * m[2][1];
		dst.z = src.x * m[0][2] + src.y * m[1][2] + src.z * m[2][2];
	}

	// initialize the coefficients of the matrix with the coefficients of the identity matrix
	T m[4][4] = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
};

typedef Matrix44<float> Matrix44f;


// void func1(const int& a)
// {
// 	cout << "A" << endl;
// }

void func1(int& a)
{
	cout << "B" << endl;
}

void func1(int&& a)
{
	cout << "C" << endl;
}

int main()
{
	int a = 3;
	func1(a);
	// Vec3<float> a(2.f, 1.f, 1.f);
	// a.normalize();
	//
	// Vec3f b(1.f, 2.f, 3.f);
	// normalize(b);
	//
	// float lengthA = a.length();
	//
	// cout << lengthA << endl;
	//
	// cout << length(a) << endl;
	//
	// cout << length(b) << endl;

	// int a = 3;
	// cin >> a;
	// int array1[a];

	return 0;
}
