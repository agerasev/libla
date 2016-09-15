#pragma once

#include<type_traits>

template <typename T, int N>
struct vec;

/* Template tools */

namespace __vec_tools__
{

template <typename T, typename S, int N>
void memcopy(T *dst, const S *src, int d = 1)
{
	for(int i = 0; i < N; ++i)
	{
		dst[i] = static_cast<T>(src[d*i]);
	}
}

template <typename T, int Left>
struct Unroller
{
	
	template <typename S, typename ... Args>
	static void unroll(T *ptr, S arg, Args ... args)
	{
		*ptr = static_cast<T>(arg);
		Unroller<T,Left-1>::unroll(ptr+1,args...);
	}
	
	template <typename S, int _N, typename ... Args>
	static void unroll(T *ptr, vec<S,_N> arg, Args ... args)
	{
		static_assert(_N <= Left, "arg size larger than needed");
		__vec_tools__::memcopy<T,S,_N>(ptr,arg._data);
		Unroller<T,Left-_N>::unroll(ptr+_N,args...);
	}
	
};

template <typename T>
struct Unroller<T,0>
{
	
	static void unroll(T *)
	{
		
	}
	
};

}

template <typename T, int N>
struct vec
{
	static_assert(std::is_arithmetic<T>::value, "vector must contain arithmetic type");
	static_assert(N > 0, "vector size must be greater than zero");
	
	/* Types */
	
	typedef T type;
	static const int size = N;
	
	/* Data */
	
	T _data[N];
	
	template <typename S>
	void memcopy(const S *src, int d = 1)
	{
		__vec_tools__::memcopy<T,S,N>(_data,src,d);
	}
	
	template <typename ... Args>
	void unroll(Args ... args)
	{
		__vec_tools__::Unroller<T,N>::unroll(_data,args...);
	}

	/* Constructors */
	
	vec()
	{
		
	}
	
	template <typename ... Args>
	vec(Args ... args)
	{
		/* TODO: Initialize with vectors with less dimensions */
		unroll(args...);
	}
	
	/* Copy constructors */
	
	template <typename S>
	vec(const vec<S,N> &v)
	{
		memcopy(v._data);
	}
	
	/* Assign operators */
	
	template <typename S>
	vec &operator = (const vec<S,N> &v)
	{
		memcopy(v._data);
		return *this;
	}
	
	/* Access operators */
	
	T *data() 
	{
		return _data;
	}
	
	const T *data() const 
	{
		return _data;
	}
	
	T &operator [](int n)
	{
		return _data[n];
	}
	
	T operator [](int n) const 
	{
		return _data[n];
	}
	
	T &operator ()(int n) 
	{
		return _data[n];
	}
	
	T operator ()(int n) const 
	{
		return _data[n];
	}
	
	template <int _N>
	vec<T,_N> sub(int n = 0) const 
	{
		vec<T,_N> r;
		r.memcopy(_data + n);
		return r;
	}
	
		/* Index access */
	
	T x() const
	{
		return _data[0];
	}
	
	T &x()
	{
		return _data[0];
	}
	
	T y() const
	{
		return _data[1];
	}
	
	T &y()
	{
		return _data[1];
	}
	
	T z() const
	{
		return _data[2];
	}
	
	T &z()
	{
		return _data[2];
	}
	
	T w() const
	{
		return _data[3];
	}
	
	T &w()
	{
		return _data[3];
	}
};

/* Addition */

template<typename T, typename S, int N>
vec<typename std::common_type<T,S>::type,N> operator +(const vec<T,N> &a, const vec<S,N> &b) 
{
	vec<typename std::common_type<T,S>::type,N> c;
	for(int i = 0; i < N; ++i) 
	{
		c._data[i] = a._data[i] + b._data[i];
	}
	return c;
}

/* Component product */

template<typename T, typename S, int N>
vec<typename std::common_type<T,S>::type,N> operator *(const vec<T,N> &a, const vec<S,N> &b) 
{
	vec<typename std::common_type<T,S>::type,N> c;
	for(int i = 0; i < N; ++i) 
	{
		c._data[i] = a._data[i]*b._data[i];
	}
	return c;
}

/* Multiplication by constant */

template<typename T, typename S, int N>
typename std::enable_if<
	std::is_arithmetic<S>::value, 
	vec<typename std::common_type<T,S>::type,N>
>::type operator *(S a, const vec<T,N> &b) 
{
	vec<typename std::common_type<T,S>::type,N> c;
	for(int i = 0; i < N; ++i) 
	{
		c._data[i] = a*b._data[i];
	}
	return c;
}

/* Dot product */

template<typename T, typename S, int N>
typename std::common_type<T,S>::type dot(const vec<T,N> &a, const vec<S,N> &b) 
{
	T c = static_cast<typename std::common_type<T,S>::type>(0);
	for(int i = 0; i < N; ++i) 
	{
		c += a._data[i]*b._data[i];
	}
	return c;
}

/* Pseudo cross product */

template <typename T, typename S>
typename std::common_type<T,S>::type cross(const vec<T,2> &a, const vec<S,2> &b)
{
	return a[0]*b[1] - a[1]*b[0];
}

template <typename T, typename S>
typename std::common_type<T,S>::type operator ^(const vec<T,2> &a, const vec<S,2> &b)
{
	return cross(a, b);
}

/* Cross product */

template <typename T, typename S>
vec<typename std::common_type<T,S>::type,3> cross(const vec<T,3> &a, const vec<S,3> &b) {
	return vec<typename std::common_type<T,S>::type,3>(
	  a[1]*b[2] - b[1]*a[2],
	  a[2]*b[0] - b[2]*a[0],
	  a[0]*b[1] - b[0]*a[1]
	);
}

template <typename T, typename S>
vec<typename std::common_type<T,S>::type,3> operator ^(const vec<T,3> &a, const vec<S,3> &b) {
	return cross(a, b);
}

template <typename T, typename S>
vec<typename std::common_type<T,S>::type,4> cross(const vec<T,4> &a, const vec<S,4> &b) {
	return vec<typename std::common_type<T,S>::type,4>(
	  a[1]*b[2] - b[1]*a[2],
	  a[2]*b[0] - b[2]*a[0],
	  a[0]*b[1] - b[0]*a[1],
	  a[3]*b[3]
	);
}

template <typename T, typename S>
vec<typename std::common_type<T,S>::type,4> operator ^(const vec<T,4> &a, const vec<S,4> &b) {
	return cross(a, b);
}

/* Division */

template<typename T, typename S, int N>
vec<typename std::common_type<T,S>::type,N> operator /(const vec<T,N> &a, const vec<S,N> &b) 
{
	vec<typename std::common_type<T,S>::type, N> r;
	for(int i = 0; i < N; ++i) 
	{
		r._data[i] = a._data[i]/b._data[i];
	}
	return r;
}

/* Derivative operations */

template<typename T, typename S, int N>
typename std::enable_if<
	std::is_arithmetic<S>::value, 
	vec<typename std::common_type<T,S>::type,N>
>::type operator *(const vec<T,N> &b, S a)
{
	return a*b;
}

template<typename T, int N>
vec<T,N> operator +(const vec<T,N> &a) 
{
	return a;
}

template<typename T, int N>
vec<T,N> operator -(const vec<T,N> &a) 
{
	return static_cast<T>(-1)*a;
}

template<typename T, typename S, int N>
vec<typename std::common_type<T,S>::type,N> operator -(const vec<T,N> &a, const vec<S,N> &b) 
{
	return a+(-b);
}

template<typename T, typename S, int N>
typename std::enable_if<
  std::is_floating_point<typename std::common_type<T,S>::type>::value, 
  vec<typename std::common_type<T,S>::type,N>
  >::type
 operator /(const vec<T,N> &b, S a) 
{
	return b*(static_cast<typename std::common_type<T,S>::type>(1)/static_cast<typename std::common_type<T,S>::type>(a));
}

template<typename T, typename S, int N>
typename std::enable_if<
  !std::is_floating_point<typename std::common_type<T,S>::type>::value, 
  vec<typename std::common_type<T,S>::type,N>
  >::type
 operator /(const vec<T,N> &b, S a) 
{
	vec<typename std::common_type<T,S>::type,N> c;
	for(int i = 0; i < N; ++i) 
	{
		c._data[i] = b._data[i]/a;
	}
	return c;
}

/* Assign operations */

template<typename T, typename S, int N>
vec<T,N> &operator +=(vec<T,N> &a, const vec<S,N> &b) 
{
	return a = a + b;
}

template<typename T, typename S, int N>
vec<T,N> &operator -=(vec<T,N> &a, const vec<S,N> &b) 
{
	return a = a - b;
}

template<typename T, typename S, int N>
vec<T,N> &operator *=(vec<T,N> &a, S b) 
{
	return a = a*b;
}

template<typename T, typename S, int N>
vec<T,N> &operator /=(vec<T,N> &a, S b) 
{
	return a = a/b;
}

/* Math */

#include<cmath>

template<typename T, int N>
inline T sqr(const vec<T,N> &v) 
{
	return dot(v, v);
}

template<typename T, int N>
inline vec<T,N> abs(const vec<T,N> &v) 
{
	vec<T,N> ret;
	for(int i = 0; i < N; ++i) 
	{
		ret._data[i] = v._data[i] > 0 ? v._data[i] : -v._data[i];
	}
	return ret;
}

template<typename T, int N>
inline T length(const vec<T,N> &v) 
{
	return sqrt(sqr(v));
}

template<typename T, int N>
inline vec<T,N> normalize(const vec<T,N> &v) 
{
	return v/length(v);
}

/* Comparison */
template<typename T, int N>
inline vec<bool,N> operator ==(const vec<T,N> &a, const vec<T,N> &b) 
{
	vec<bool,N> ret;
	for(int i = 0; i < N; ++i) 
	{
		ret._data[i] = (a._data[i] == b._data[i]);
	}
	return ret;
}

template<typename T, int N>
inline vec<bool,N> operator !=(const vec<T,N> &a, const vec<T,N> &b) 
{
	vec<bool,N> ret;
	for(int i = 0; i < N; ++i) 
	{
		ret._data[i] = (a._data[i] != b._data[i]);
	}
	return ret;
}

template<typename T, int N>
inline vec<bool,N> operator <(const vec<T,N> &a, const vec<T,N> &b) 
{
	vec<bool,N> ret;
	for(int i = 0; i < N; ++i) 
	{
		ret._data[i] = (a._data[i] < b._data[i]);
	}
	return ret;
}

template<typename T, int N>
inline vec<bool,N> operator >(const vec<T,N> &a, const vec<T,N> &b) 
{
	vec<bool,N> ret;
	for(int i = 0; i < N; ++i) 
	{
		ret._data[i] = (a._data[i] > b._data[i]);
	}
	return ret;
}

template<typename T, int N>
inline vec<bool,N> operator <=(const vec<T,N> &a, const vec<T,N> &b) 
{
	vec<bool,N> ret;
	for(int i = 0; i < N; ++i) 
	{
		ret._data[i] = (a._data[i] <= b._data[i]);
	}
	return ret;
}

template<typename T, int N>
inline vec<bool,N> operator >=(const vec<T,N> &a, const vec<T,N> &b) 
{
	vec<bool,N> ret;
	for(int i = 0; i < N; ++i) 
	{
		ret._data[i] = (a._data[i] >= b._data[i]);
	}
	return ret;
}


/* Bool vector functions */

template <typename T, int N>
bool all(const vec<T, N> &v) {
	bool c = true;
	for(int i = 0; i < N; ++i) {
		c = c && v[i];
	}
	return c;
}

template <typename T, int N>
bool any(const vec<T, N> &v) {
	bool c = false;
	for(int i = 0; i < N; ++i) {
		c = c || v[i];
	}
	return c;
}

/* Max/min */
template <typename T, int N>
T max(const vec<T, N> &v) {
	T m = v[0];
	for(int i = 1; i < N; ++i) {
		m = v[i] > m ? v[i] : m;
	}
	return m;
}

template <typename T, int N>
T min(const vec<T, N> &v) {
	T m = v[0];
	for(int i = 1; i < N; ++i) {
		m = v[i] < m ? v[i] : m;
	}
	return m;
}

/* Type aliases and constants */

typedef vec<double,2> dvec2;
typedef vec<float,2>  fvec2;
typedef vec<int,2>    ivec2;
typedef vec<bool,2>   bvec2;

const dvec2 nulldvec2(0.0,0.0);
const fvec2 nullfvec2(0.0f,0.0f);
const ivec2 nullivec2(0,0);

typedef dvec2 vec2;
const vec2 nullvec2 = nulldvec2;

typedef vec<double,3> dvec3;
typedef vec<float,3>  fvec3;
typedef vec<int,3>    ivec3;
typedef vec<bool,3>   bvec3;

const dvec3 nulldvec3(0.0,0.0,0.0);
const fvec3 nullfvec3(0.0f,0.0f,0.0f);
const ivec3 nullivec3(0,0,0);

typedef dvec3 vec3;
const vec3 nullvec3 = nulldvec3;

typedef vec<double,4> dvec4;
typedef vec<float,4>  fvec4;
typedef vec<int,4>    ivec4;
typedef vec<bool,4>   bvec4;

const dvec4 nulldvec4(0.0,0.0,0.0,0.0);
const fvec4 nullfvec4(0.0f,0.0f,0.0f,0.0f);
const ivec4 nullivec4(0,0,0,0);

typedef dvec4 vec4;
const vec4 nullvec4 = nulldvec4;
