#pragma once

#include<type_traits>

#include"vec.hpp"

namespace __mat_tools__
{

template <typename T, typename S, int N, int M>
void memcopy(T *dst, const S *src, int dy = 0, int dx = 1)
{
	for(int i = 0; i < M; ++i) {
		for(int j = 0; j < N; ++j) {
			dst[i*N + j] = static_cast<T>(src[(dx*N + dy)*i + dx*j]);
		}
	}
}

template <typename T, int N, int Left>
struct Unroller;

template <typename T, int N, int Left, typename S>
struct ArgumentAssigner
{
	template <typename ... Args>
	static void assign(T *ptr, S arg, Args ... args)
	{
		*ptr = static_cast<T>(arg);
		Unroller<T,N,Left-1>::unroll(ptr+1,args...);
	}
};

template <typename T, int N, int Left, typename S>
struct ArgumentAssigner<T,N,Left,vec<S,N>>
{
	template <typename ... Args>
	static void assign(T *ptr, const vec<S,N> &arg, Args ... args)
	{
		static_assert(!(Left%N), "vectors are not aligned with rows");
		__vec_tools__::memcopy<T,S,N>(ptr,static_cast<vec<T,N>>(arg).data);
		Unroller<T,N,Left-N>::unroll(ptr+N,args...);
	}
};

template <typename T, int N, int Left>
struct Unroller
{
	template <typename S, typename ... Args>
	static void unroll(T *ptr, S arg, Args ... args)
	{
		ArgumentAssigner<T,N,Left,S>::assign(ptr,arg,args...);
	}
};

template <typename T, int N>
struct Unroller<T,N,0>
{
	static void unroll(T *)
	{
		
	}
};

}

/* Matrix stricture */
/* N - columns, M - rows */

template <typename T, int N, int M>
struct mat
{
	/* Types */
	
	typedef T type;
	static const int row_size = N;
	static const int column_size = M;
	
	/* Data */
	
	T data[N*M];
	
	template <typename S>
	void memcopy(const S *src, int dy = 0, int dx = 1)
	{
		__mat_tools__::memcopy<T,S,N,M>(data,src,dy,dx);
	}
	
	template <typename ... Args>
	void unroll(Args ... args)
	{
		__mat_tools__::Unroller<T,N,N*M>::unroll(data,args...);
	}
	
	/* Constructors */
	
	mat()
	{
		
	}
	
	template <typename ... Args>
	mat(Args ... args)
	{
		unroll(args...);
	}
	
	/* Copy constructors */
	
	template <typename S>
	mat(const mat<S,N,M> &m)
	{
		memcopy(m.data);
	}
	
	/* Assign operators */

	template <typename S>
	mat<T,N,M> &operator =(const mat<S,N,M> &m)
	{
		memcopy(m.data);
		return *this;
	}

	/* Access operators */
	
	T &get(int x, int y)
	{
		return data[y*N + x];
	}
	
	T get(int x, int y) const
	{
		return data[y*N + x];
	}
	
	T &operator ()(int x, int y)
	{
		return get(x,y);
	}
	
	T operator ()(int x, int y) const
	{
		return get(x,y);
	}
	
	T &operator ()(const vec<int,2> &v)
	{
		return get(v[0],v[1]);
	}
	
	T operator ()(const vec<int,2> &v) const
	{
		return get(v[0],v[1]);
	}
	
	vec<T,N> row(int n) const
	{
		vec<T,N> ret;
		ret.memcopy(data + n*N);
		return ret;
	}
	
	vec<T,M> col(int n) const
	{
		vec<T,N> ret;
		ret.memcopy(data + n,N);
		return ret;
	}
	
	T &operator [](const vec<int,2> &v)
	{
		return get(v[0],v[1]);
	}
	
	T operator [](const vec<int,2> &v) const
	{
		return get(v[0],v[1]);
	}
	
	T *operator[](int n)
	{
		return data + n*N;
	}
	
	const T *operator[](int n) const
	{
		return data + n*N;
	}
	
	/* Submatrix */
	
	mat<T,N-1,M-1> sub(int x, int y) const
	{
		mat<T,N-1,M-1> c;
		for(int ix = 0, jx = 0; ix < N; ++ix,++jx)
		{
			if(ix == x) 
			{
				--jx;
				continue;
			}
			for(int iy = 0, jy = 0; iy < M; ++iy,++jy)
			{
				if(iy == y) 
				{
					--jy;
					continue;
				}
				c(jx,jy) = (*this)(ix,iy);
			}
		}
		return c;
	}
};

namespace __mat_tools__
{

template <typename T, int N>
struct Determinator
{
	static T det(const mat<T,N,N> &m)
	{
		T c = static_cast<T>(0);
		const int rc = 0;
		for(int i = 0; i < N; ++i)
		{
			c += m(i,rc)*cofactor(m,i,rc);
		}
		return c;
	}
};

template <typename T>
struct Determinator<T,1>
{
	static T det(const mat<T,1,1> &m)
	{
		return m.data[0];
	}
};

}

/* Transpose */

template <typename T, int N, int M>
mat<T,M,N> transpose(const mat<T,N,M> &m)
{
	return mat<T,M,N>().memcopy(m.data,-M*N+1,N);
}

/* Determinant */

template <typename T, int N>
T cofactor(const mat<T,N,N> &m, int x, int y)
{
	return (1 - 2*((x+y)%2))*det(m.sub(x,y));
}

template <typename T, int N>
T det(const mat<T,N,N> &m)
{
	return __mat_tools__::Determinator<T,N>::det(m);
}

/* Adjugate matrix */

template <typename T, int N>
mat<T,N,N> adj(const mat<T,N,N> &m)
{
	mat<T,N,N> ret;
	for(int i = 0; i < N; ++i)
	{
		for(int j = 0; j < N; ++j)
		{
			ret(j,i) = cofactor(m,i,j);
		}
	}
	return ret;
}

/* Inverse matrix */
template <typename T, int N>
mat<T,N,N> invert(const mat<T,N,N> &m)
{
	return adj(m)/det(m);
}

/* Basic operations */

/* Addition */

template <typename T, typename S, int N, int M>
mat<typename std::common_type<T,S>::type,N,M> operator +(const mat<T,N,M> &a, const mat<S,N,M> &b) 
{
	mat<typename std::common_type<T,S>::type,N,M> c;
	for(int i = 0; i < N*M; ++i) 
	{
		c.data[i] = a.data[i] + b.data[i];
	}
	return c;
}

/* Multiplication by constant */

template <typename T, typename S, int N, int M>
mat<typename std::common_type<T,S>::type,N,M> operator *(const mat<T,N,M> &m, S s) 
{
	mat<typename std::common_type<T,S>::type,N,M> c;
	for(int i = 0; i < N*M; ++i) 
	{
		c.data[i] = s*m.data[i];
	}
	return c;
}

/* Multiplication by vector */

template <typename T, typename S, int N, int M>
vec<typename std::common_type<T,S>::type,M> operator *(const mat<T,N,M> &m, const vec<S,N> &v)
{
	vec<typename std::common_type<T,S>::type,M> c;
	for(int i = 0; i < M; ++i) 
	{
		c(i) = m.row(i)*v;
	}
	return c;
}

/* Product of matrices */

template <typename T, typename S, int N, int M, int L>
mat<typename std::common_type<T,S>::type,N,M> operator *(const mat<T,L,M> &a, const mat<S,N,L> &b) 
{
  mat<typename std::common_type<T,S>::type,N,M> c;
  for(int i = 0; i < M; ++i)
  {
	  for(int j = 0; j < N; ++j)
	  {
		  c(j,i) = a.row(i)*b.col(j);
	  }
  }
  return c;
}

/* Derivative operations */

template <typename T, int N, int M>
mat<T,N,M> operator +(const mat<T,N,M> &a) 
{
	return a;
}

template <typename T, int N, int M>
mat<T,N,M> operator -(const mat<T,N,M> &a) 
{
	return static_cast<T>(-1)*a;
}

template <typename T, typename S, int N, int M>
mat<typename std::common_type<T,S>::type,N,M> operator -(const mat<T,N,M> &a, const mat<S,N,M> &b) 
{
	return a + (-b);
}

template <typename T, typename S, int N, int M>
mat<typename std::common_type<T,S>::type,N,M> operator *(S s, const mat<T,N,M> &m) 
{
	return m*s;
}

template <typename T, typename S, int N, int M>
mat<typename std::common_type<T,S>::type,N,M> operator /(const mat<T,N,M> &m, S s) 
{
	return m*(static_cast<typename std::common_type<T,S>::type>(1)/static_cast<typename std::common_type<T,S>::type>(s));
}

/* Assign operators */

template <typename T, int N, int M>
mat<T,N,M> &operator +=(mat<T,N,M> &a, const mat<T,N,M> &b) 
{
	return a = a + b;
}

template <typename T, int N, int M>
mat<T,N,M> &operator -=(mat<T,N,M> &a, const mat<T,N,M> &b) 
{
	return a = a + b;
}

template <typename T, int N, int M>
mat<T,N,M> &operator *=(mat<T,N,M> &m, T s) 
{
	return m = m*s;
}

template <typename T, int N, int M>
mat<T,N,M> &operator /=(mat<T,N,M> &m, T s) 
{
	return m = m*s;
}

/* Type aliases and constants */

typedef mat<double,2,2> dmat2;
typedef mat<float,2,2>  fmat2;
typedef mat<int,2,2>	   imat2;
typedef mat<double,3,3> dmat3;
typedef mat<float,3,3>  fmat3;
typedef mat<int,3,3>	   imat3;
typedef mat<double,4,4> dmat4;
typedef mat<float,4,4>  fmat4;
typedef mat<int,4,4>	   imat4;
typedef dmat2 mat2;
typedef dmat3 mat3;
typedef dmat4 mat4;

const imat2 nullimat2(0,0,0,0);
const dmat2 nulldmat2 = nullimat2;
const fmat2 nullfmat2 = nullimat2;
const imat3 nullimat3(0,0,0,0,0,0,0,0,0);
const dmat3 nulldmat3 = nullimat3;
const fmat3 nullfmat3 = nullimat3;
const imat4 nullimat4(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
const dmat4 nulldmat4 = nullimat4;
const fmat4 nullfmat4 = nullimat4;
const mat2 nullmat2 = nulldmat2;
const mat3 nullmat3 = nulldmat3;
const mat4 nullmat4 = nulldmat4;

const imat2 uniimat2(1,0,0,1);
const dmat2 unidmat2 = uniimat2;
const fmat2 unifmat2 = uniimat2;
const imat3 uniimat3(1,0,0,0,1,0,0,0,1);
const dmat3 unidmat3 = uniimat3;
const fmat3 unifmat3 = uniimat3;
const imat4 uniimat4(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1);
const dmat4 unidmat4 = uniimat4;
const fmat4 unifmat4 = uniimat4;
const mat2 unimat2 = unidmat2;
const mat3 unimat3 = unidmat3;
const mat4 unimat4 = unidmat4;
