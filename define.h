#ifndef _DEFINE_H
#define _DEFINE_H

#include <vector>
#include <cerrno>
#include <limits>
#include <string>

using namespace std;

#define GSL_RANGE_CHECK_OFF

typedef double					DATA;
typedef vector<DATA>			NODE;

extern const DATA DATA_NAN;
extern const DATA PI;
extern const DATA PI_2;
extern const DATA R2D;
extern const DATA D2R;
extern const DATA LOG2;
extern const DATA INV_LOG2;
extern const DATA INV_255;
extern const DATA INV_360;

extern const string WORK_DIR;
extern const string TEST_DATA_DIR;
extern const string OUTPUT_DIR;

template<class T>
class Vect2D
{
public:
	Vect2D();
	Vect2D(T x, T y);
	Vect2D(const Vect2D<T> &vFrom);
	~Vect2D();
	Vect2D<T>& operator=(const Vect2D<T> &v)
	{
		m_x = v.m_x;
		m_y = v.m_y;
		return *this;
	}

	T  m_x;		T  m_y;

private:
};

template <class T>
class Vect3D
{
public:
	Vect3D(T x, T y, T z); //restrict(cpu, amp);
	Vect3D(const Vect3D<T> &vFrom); //restrict(cpu, amp);
	Vect3D(); //restrict(cpu, amp);
	~Vect3D(); //restrict(cpu, amp);
	Vect3D<T>& operator=(const Vect3D<T> &v) //restrict(cpu, amp)
	{
		m_x = v.m_x;
		m_y = v.m_y;
		m_z = v.m_z;
		return *this;
	}

	T m_x;		T m_y;		T  m_z;

private:
};

template <class T>
class Vect4D
{
public:
	Vect4D();
	Vect4D(T r, T g, T b, T a);
	Vect4D(const Vect4D<T> &vFrom);
	~Vect4D();
	Vect4D<T>& operator=(const Vect4D<T> &v)
	{
		m_r = v.m_r;
		m_g = v.m_g;
		m_b = v.m_b;
		m_a = v.m_a;
		return *this;
	}

	union{T m_r;	T m_x;};
	union{T m_g;	T m_y;};	
	union{T m_b;	T m_z;};		
	union{T m_a;	T m_c;};	

private:
};

template<class T> T FindMin(Vect3D<T> &in);
template<class T> T FindMax(Vect3D<T> &in);

template<class T> void Mult(Vect3D<T> &in, DATA scl);

//*************************************************************************************************

template<class T>
Vect2D<T>::Vect2D()
	:m_x(0), m_y(0) {}

template<class T> 
Vect2D<T>::Vect2D(T x, T y)
	:m_x(x), m_y(y) {}

template<class T>
Vect2D<T>::Vect2D(const Vect2D &vFrom)
	:m_x(vFrom.m_x), m_y(vFrom.m_y) {}

template<class T>
Vect2D<T>::~Vect2D() {}

template<class T>
Vect3D<T>::Vect3D(T x, T y, T z) //restrict(cpu, amp)
	:m_x(x), m_y(y), m_z(z) {}

template<class T>
Vect3D<T>::Vect3D(const Vect3D<T> &vFrom) //restrict(cpu, amp)
	:m_x(vFrom.m_x), m_y(vFrom.m_y), m_z(vFrom.m_z) {}

template<class T>
Vect3D<T>::Vect3D() //restrict(cpu, amp)
	:m_x(0), m_y(0), m_z(0) {}

template<class T>
Vect3D<T>::~Vect3D() //restrict(cpu, amp)
{}

template<class T>
Vect4D<T>::Vect4D()
	:m_r(0), m_g(0), m_b(0), m_a(0) {}

template<class T>
Vect4D<T>::Vect4D(T r, T g, T b, T a)
	:m_r(r), m_g(g), m_b(b), m_a(a) {}

template<class T>
Vect4D<T>::Vect4D(const Vect4D<T> &vFrom)
	:m_r(vFrom.m_r), m_g(vFrom.m_g), m_b(vFrom.m_b), m_a(vFrom.m_a) {}

template<class T>
Vect4D<T>::~Vect4D() {}

template<class T>
T FindMin(Vect3D<T> &in)
{
	T min = (in.m_x<in.m_y)? in.m_x: in.m_y;
	if(min > in.m_z)	min = in.m_z;
	else {}
	return min;
}

template<class T>
T FindMax(Vect3D<T> &in)
{
	T max = (in.m_x>in.m_y)? in.m_x: in.m_y;
	if(max < in.m_z)	max = in.m_z;
	else {}
	return max;
}

template<class T> 
void Mult(Vect3D<T> &in, DATA scl)
{
	in.m_x *= scl;
	in.m_y *= scl;
	in.m_z *= scl;
}

//*************************************************************************************************

void _MyAssert_(bool st, unsigned ln, string fn);
#define MyAssert(st) (_MyAssert_((st), __LINE__, __FILE__))

void PrintCount(unsigned c, unsigned lim);
void PrintLine(unsigned c, unsigned lineC);

#endif