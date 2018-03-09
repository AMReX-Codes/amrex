#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// These are about ten times faster than the STL analogues, since here the
// number of elements is known at compile time.  These functions are as
// fast as hand-unrolled loops (but you need -O3 or better on gcc).
//
// For examples, see lib/src/BoxTools/BaseFabImplem.H in the "template" branch
// of the Chombo_prev module under ~tdsternberg/my-cvsroot.
//
// Author: Ted
//

#ifndef _METAPROGRAMS_H_
#define _METAPROGRAMS_H_

#include <functional>

namespace Metaprograms
{

template<typename T> struct Identity
{
    T const& operator()( T const& t )
    {
      return t;
    }
};

//-------------------------------------------------------------------
template<int N, typename T, typename RT, typename PlusT, typename TimesT>
struct InnerProduct
{
    //
    // Binary mode with null-constructed functors.
    //
    RT operator()( T const* v1, T const* v2 )
    {
        return PlusT()( InnerProduct<N-1,T,RT,PlusT,TimesT>()(v1,v2),
                        TimesT()(v1[N-1],v2[N-1]) );
    }
    RT operator()( T const* v, T const& x )
    {
        return PlusT()( InnerProduct<N-1,T,RT,PlusT,TimesT>()(v,x),
                        TimesT()(v[N-1],x) );

    }
};

template<typename PlusT, typename T, typename RT, typename TimesT>
struct InnerProduct<1,T,RT,PlusT,TimesT>
{
    RT operator()( T const* v1, T const* v2 )
    {
        return TimesT()( v1[0], v2[0] );
    }
    RT operator()( T const* v, T const& x )
    {
        return TimesT()( v[0], x );
    }
};
//-------------------------------------------------------------------

//-------------------------------------------------------------------
template<int N, typename T, typename CompareT>
inline bool pointwiseCompare( T const* v1, T const* v2 )
{
    return InnerProduct< N,T,bool,
                         std::logical_and<bool>,
                         CompareT >()( v1, v2 );
}

template<int N, typename T, typename CompareT>
inline bool pointwiseCompare( T const* v, T const& x )
{
    return InnerProduct< N,T,bool,
                         std::logical_and<bool>,
                         CompareT >()( v, x );
}
//-------------------------------------------------------------------

//-------------------------------------------------------------------
template<int N, typename T, typename ReduceT, typename TransformT=Identity<T> >
struct Accum
{
    T operator()( T const * v ) const
    {
        return ReduceT()( TransformT()(v[N-1]), Accum<N-1,T,ReduceT,TransformT>()(v) );
    }

    // Stateful TransformT
    T operator()( T const * v, TransformT const& xform ) const
    {
        return ReduceT()( xform(v,N-1), Accum<N-1,T,ReduceT,TransformT>()(v,xform) );
    }
};

template<typename T, typename ReduceT, typename TransformT >
struct Accum<1,T,ReduceT,TransformT>
{
    T operator()( T const * v ) const
    {
        return TransformT()(v[0]);
    }

    // Stateful TransformT
    T operator()( T const * v, TransformT const& xform ) const
    {
        return xform(v,0);
    }
};
//-------------------------------------------------------------------

//-------------------------------------------------------------------
template<int N, typename T> struct LexLT
{
    bool operator()( const T* v1, const T* v2 )
    {
        return doit( v1+N, v2+N );
    }

    static bool doit( const T* v1, const T* v2 )
    {
        if      ( v1[-N] < v2[-N] ) return true;
        else if ( v1[-N] > v2[-N] ) return false;
        else                        return LexLT<N-1,T>::doit(v1,v2);
    }
};

template<typename T> struct LexLT<1,T>
{
    bool operator()( const T* v1, const T* v2 )
    {
        return doit( v1+1, v2+1 );
    }
    static bool doit( const T* v1, const T* v2 )
    {
        if ( v1[-1] < v2[-1] ) return true;
        else                   return false;
    }
};
//-------------------------------------------------------------------

//-------------------------------------------------------------------
template<int N, typename T, typename FunctorT> struct Transform
{
    void operator()( T* v )
    {
        v[N-1] = FunctorT()(v[N-1]);
        Transform<N-1,T,FunctorT>()(v);
    }
    void operator()( T* v, const T& x )
    {
        v[N-1] = FunctorT()(v[N-1],x);
        Transform<N-1,T,FunctorT>()(v,x);
    }
    void operator()( T* v1, const T* v2 )
    {
        v1[N-1] = FunctorT()(v1[N-1],v2[N-1]);
        Transform<N-1,T,FunctorT>()(v1,v2);
    }
};

template<typename T, typename FunctorT> struct Transform<1,T,FunctorT>
{
    void operator()( T* v )
    {
        v[0] = FunctorT()(v[0]);
    }
    void operator()( T* v, const T& x )
    {
        v[0] = FunctorT()(v[0],x);
    }
    void operator()( T* v1, const T* v2 )
    {
        v1[0] = FunctorT()(v1[0],v2[0]);
    }
};
//-------------------------------------------------------------------

//-------------------------------------------------------------------
template<int N, int P> struct Pow
{
    static const int value = N * Pow<N,P-1>::value;
};

template<int N> struct Pow<N,1>
{
    static const int value = N;
};
//-------------------------------------------------------------------

//-------------------------------------------------------------------
/**
  Arbitrarily deep nested loop.
  Example: This code will print all the integers from 0 to B^N-1, in base B:

  template<int N> struct BaseRep
  {
    void operator()( int* index )
    {
      for ( int i=N-1; i>=0; --i ) std::cout << index[i];
      std::cout << '\n';
    }
  };

  main()
  {
  #define N 4
  #define B 3
    int index[N];
    BaseRep<N> op;
    Metaprograms::NestedLoop< N,BaseRep<N> >()( index, 0, B, op );
  }
*/
template<int N, class OP> struct NestedLoop
{
  void operator()( int * index, int lo, int hi, OP& op ) const
  {
      for ( index[N-1]=lo; index[N-1]<hi; ++index[N-1] )
      {
          NestedLoop<N-1,OP>()( index, lo, hi, op );
      }
  }
  void operator()( int * index, const int * lo, const int * hi, OP& op ) const
  {
      for ( index[N-1]=lo[N-1]; index[N-1]<hi[N-1]; ++index[N-1] )
      {
          NestedLoop<N-1,OP>()( index, lo, hi, op );
      }
  }
};

template<class OP> struct NestedLoop<0,OP>
{
    void operator()( int * index, int lo, int hi, OP& op ) const
    {
        op( index );
    }
    void operator()( int * index, const int * lo, const int * hi, OP& op ) const
    {
        op( index );
    }
};
//-------------------------------------------------------------------

//-------------------------------------------------------------------
template<int N, class OP> struct NestedPrestagedLoop
{
  void operator()( int * index, int lo, int hi, OP& op ) const
  {
      for ( index[N-1]=lo; index[N-1]<hi; ++index[N-1] )
      {
          op.prestage(N-1);
          NestedPrestagedLoop<N-1,OP>()( index, lo, hi, op );
      }
  }
  void operator()( int * index, const int * lo, const int * hi, OP& op ) const
  {
      for ( index[N-1]=lo[N-1]; index[N-1]<hi[N-1]; ++index[N-1] )
      {
          op.prestage(N-1);
          NestedPrestagedLoop<N-1,OP>()( index, lo, hi, op );
      }
  }
};

template<class OP> struct NestedPrestagedLoop<0,OP>
{
    void operator()( int * index, int lo, int hi, OP& op ) const
    {
        op( index );
    }
    void operator()( int * index, const int * lo, const int * hi, OP& op ) const
    {
        op( index );
    }
};
//-------------------------------------------------------------------

//-------------------------------------------------------------------
/** Named for defunct D_TERM macro.
 *  Applies an operation once for n = 0,...N-1.
*/
template<int N, class OP> struct dterm
{
    void operator()( OP& op ) const
    {
        dterm<N-1,OP>()( op );
        op( N-1 );
    }
};

template<class OP> struct dterm<1,OP>
{
    void operator()( OP& op ) const
    {
        op( 0 );
    }
};
//-------------------------------------------------------------------

} // namespace Metaprograms

#endif // include guard
