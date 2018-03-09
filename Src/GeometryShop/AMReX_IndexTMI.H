///heh, TMI
#ifndef _INDEXTMI_H_
#define _INDEXTMI_H_

#include <iostream>
using std::ostream;
using std::istream;
using std::ws;

#include "AMReX_IndexTM.H"

#include <cmath>
#include <algorithm>
#include <functional>
#include <numeric>
#include "AMReX_MetaPrograms.H"

namespace amrex
{

  template<typename T, int N> ostream& operator<< (ostream            & a_os,
                                                   const IndexTM<T,N> & a_p)
  {
    a_os << '(';

    for (int i = 0; i < N; ++i)
    {
      a_os << a_p[i];

      if (i < N-1)
      {
        a_os << ',';
      }
    }

    a_os << ')';

    if (a_os.fail())
    {
      MayDay::Abort("operator<<(ostream&,Index&) failed");
    }

    return a_os;
  }

//
// Copied from <Utility.H>
//
#define CH_IGNORE_MAX 100000

  template<typename T, int N> istream& operator>> (istream      & a_is,
                                                   IndexTM<T,N> & a_p)
  {
    a_is >> ws;

    char c;
    a_is >> c;
    a_is.putback(c);

    if (c == '(')
    {
      a_is.ignore(CH_IGNORE_MAX, '(') >> a_p[0];

      for (int i = 1; i < N; ++i)
      {
        a_is.ignore(CH_IGNORE_MAX, ',') >> a_p[i];
      }

      a_is.ignore(CH_IGNORE_MAX, ')');
    }
    else if (c == '<')
    {
      a_is.ignore(CH_IGNORE_MAX, '<') >> a_p[0];

      for (int i = 1; i < N; ++i)
      {
        a_is.ignore(CH_IGNORE_MAX, ',') >> a_p[1];
      }
      a_is.ignore(CH_IGNORE_MAX, '>');
    }
    else
    {
      MayDay::Abort("operator>>(istream&,Index&): expected \'(\' or \'<\'");
    }

    if (a_is.fail())
    {
      MayDay::Abort("operator>>(istream&,Index&) failed");
    }

    return a_is;
  }

  template<typename T, int N> void IndexTM<T,N>::printOn (ostream & a_os) const
  {
    a_os << "Index: " << *this << "\n";
  }

  template<typename T, int N> void IndexTM<T,N>::p() const
  {
    pout() << *this << "\n";
  }

  template<typename T, int N> void IndexTM<T,N>::dumpOn (ostream & a_os) const
  {
    a_os << "Index " << *this << "\n";
  }

//
// Static object initialization.
//
/*
  template<typename T, int N> int IndexTM<T,N>::InitStatics()
  {
  IndexTM<T,N>& pz = const_cast<IndexTM<T,N>&>(IndexTM<T,N>::Zero);

  pz.setAll(0);

  IndexTM<T,N>& pu = const_cast<IndexTM<T,N>&>(IndexTM<T,N>::Unit);

  pu.setAll(1);

  // No danger of Index::Zero and Unit not having been allocated, as ARM section
  // 3.4 says "The initialization of nonlocal static objects in a translation unit
  // is done before the first use of any function...defined in that translation
  // unit."
  //
  // Had to go through the const_cast stuff because it's nice to be able to declare
  // Index::Zero and Index::Unit as const.

  return 0; // arbitrary
  }
*/

//
// Inlines.
//
  template<typename T, int N> inline IndexTM<T,N>::IndexTM(const T *a_a)
    : GenericArithmeticable< T, IndexTM<T,N> >(this)
  {
    memcpy(m_vect, a_a, N*sizeof(T));
  }

  template<typename T, int N> inline IndexTM<T,N>::IndexTM(const char* a_reference) 
    : GenericArithmeticable< T, IndexTM<T,N> >(this)
  {
    if     (a_reference[0] == 48) setAll(0);
    else if (a_reference[0] == 49) setAll(1);
    else
      MayDay::Error("unknown static initialization for IndexTM");
  }

  template<typename T, int N> inline IndexTM<T,N>::IndexTM(const IndexTM<T,N> & a_iv)
    : GenericArithmeticable< T, IndexTM<T,N> >(this)
  {
    memcpy(m_vect, a_iv.m_vect, N*sizeof(T));
  }

  template<typename T, int N> inline IndexTM<T,N>& IndexTM<T,N>::operator=(const IndexTM<T,N> & a_iv)
  {
    memcpy(m_vect, a_iv.m_vect, N*sizeof(T));
    return *this;
  }

  template<typename T, int N> inline T& IndexTM<T,N>::operator[] (int a_i)
  {
    CH_assert(a_i >= 0 && a_i < N);
    return m_vect[a_i];
  }

  template<typename T, int N> inline T IndexTM<T,N>::operator[] (int a_i) const
  {
    CH_assert(a_i >= 0 && a_i < N);
    return m_vect[a_i];
  }

  template<typename T, int N> inline void IndexTM<T,N>::setVal(int a_i,
                                                               T   a_val)
  {
    CH_assert(a_i >= 0 && a_i < N);
    m_vect[a_i] = a_val;
  }

  template<typename T, int N> inline void IndexTM<T,N>::setAll(T a_val)
  {
    for (int i = 0; i < N; i++)
    {
      m_vect[i] = a_val;
    }
  }

  template<typename T, int N> inline const T* IndexTM<T,N>::dataPtr() const
  {
    return m_vect;
  }

  template<typename T, int N> inline T* IndexTM<T,N>::dataPtr()
  {
    return m_vect;
  }

  template<typename T, int N> inline const T* IndexTM<T,N>::getVect() const
  {
    return m_vect;
  }

  template<typename T, int N> inline bool IndexTM<T,N>::operator== (const IndexTM & a_p) const
  {
    return Metaprograms::pointwiseCompare<N,T,std::equal_to<T> >(m_vect,a_p.m_vect);
  }

  template<typename T, int N> inline bool IndexTM<T,N>::operator!= (const IndexTM & a_p) const
  {
    return !(operator==(a_p));
  }

  template<typename T, int N> inline bool IndexTM<T,N>::lexLT(const IndexTM & a_s) const
  {
    return Metaprograms::LexLT<N,T>()(m_vect,a_s.m_vect);
  }

  template<typename T, int N> inline bool IndexTM<T,N>::lexGT(const IndexTM & a_s) const
  {
    return ! Metaprograms::LexLT<N,T>()(m_vect,a_s.m_vect);
  }

  template<typename T, int N> inline IndexTM<T,N> IndexTM<T,N>::operator+ () const
  {
    return IndexTM<T,N>(*this);
  }

  template<typename T, int N> inline IndexTM<T,N> IndexTM<T,N>::operator- () const
  {
    IndexTM<T,N> result(*this);
    for (int i = 0; i < N; ++i)
    {
      result.m_vect[i] *= -1;
    }
    return result;
  }

  template<typename T, int N> inline T IndexTM<T,N>::dotProduct(const IndexTM<T,N> & a_rhs ) const
  {
    return Metaprograms::InnerProduct<N,T,T,
                                      std::plus<T>,
                                      std::multiplies<T> >()(m_vect,
                                                             a_rhs.m_vect);
  }

  template<typename T, int N> inline T IndexTM<T,N>::sum() const
  {
    return Metaprograms::Accum<N,T,std::plus<T> >()(m_vect);
  }

  template<typename T, int N> inline T IndexTM<T,N>::product () const
  {
    return Metaprograms::Accum<N,T,std::multiplies<T> >()(m_vect);
  }

  template<typename T, int N> template<typename OP> bool IndexTM<T,N>::operatorCompare(const IndexTM<T,N> & a_p,
                                                                                       const OP           & a_op) const
  {
    return Metaprograms::pointwiseCompare<N,T,OP>(m_vect,a_p.m_vect);
  }

  template<typename T, int N> template<typename OP> inline IndexTM<T,N>& IndexTM<T,N>::operatorOpEquals(const T  & a_s,
                                                                                                        const OP & a_op)
  {
    Metaprograms::Transform<N,T,OP>()(m_vect,a_s);
    return *this;
  }

  template<typename T, int N> template<typename OP> inline IndexTM<T,N>& IndexTM<T,N>::operatorOpEquals(const IndexTM<T,N> & a_p,
                                                                                                        const OP           & a_op)
  {
    Metaprograms::Transform<N,T,OP>()(m_vect,a_p.m_vect);
    return *this;
  }

  template<typename T, int N> inline IndexTM<T,N>& IndexTM<T,N>::reciprocal()
  {
    std::transform(m_vect, m_vect+N, m_vect, bind1st(std::divides<T>(),T(1)));
    return *this;
  }

  template<typename T> static bool abscompare(const T & a_a,
                                              const T & a_b)
  {
    return Abs(a_a) < Abs(a_b);
  }
  template<typename T, int N> inline int IndexTM<T,N>::minDir(bool a_doAbs) const
  {
    if (a_doAbs)
    {
      return std::min_element(m_vect, m_vect+N, std::ptr_fun(abscompare<T>)) - m_vect;
    }
    else
    {
      return std::min_element(m_vect, m_vect+N) - m_vect;
    }
  }

  template<typename T, int N> inline int IndexTM<T,N>::maxDir(bool a_doAbs) const
  {
    if (a_doAbs)
    {
      return std::max_element(m_vect, m_vect+N, std::ptr_fun(abscompare<T>)) - m_vect;
    }
    else
    {
      return std::max_element(m_vect, m_vect+N) - m_vect;
    }
  }

  template<typename T> static T ourmin(T a_a,
                                       T a_b)
  {
    return ((a_a < a_b) ? a_a : a_b);
  }

  template<typename T> static T ourmax(T a_a,
                                       T a_b)
  {
    return (a_a > a_b ? a_a : a_b);
  }

  template<typename T, int N> inline IndexTM<T,N>& IndexTM<T,N>::min(const IndexTM<T,N> & a_p)
  {
    std::transform(m_vect, m_vect+N, a_p.m_vect, m_vect,
                   std::ptr_fun(ourmin<T>));
    return *this;
  }

  template<typename T, int N> inline IndexTM<T,N>& IndexTM<T,N>::max(const IndexTM<T,N> & a_p)
  {
    std::transform(m_vect, m_vect+N, a_p.m_vect, m_vect,
                   std::ptr_fun(ourmax<T>));
    return *this;
  }

  template<typename T, int N> inline IndexTM<T,N>& IndexTM<T,N>::scale(T a_s)
  {
    return (*this) *= a_s;
  }

  template<typename T, int N> inline IndexTM<T,N>& IndexTM<T,N>::reflect (T   a_refIx,
                                                                          int a_idir)
  {
    CH_assert(a_idir >= 0 && a_idir < N);
    m_vect[a_idir] = -m_vect[a_idir] + 2*a_refIx;
    return *this;
  }

  template<typename T, int N> inline IndexTM<T,N>& IndexTM<T,N>::shift(int a_coord,
                                                                       T   a_s)
  {
    CH_assert(a_coord >= 0 && a_coord < N);
    m_vect[a_coord] += a_s;
    return *this;
  }

  template<typename T, int N> inline IndexTM<T,N>& IndexTM<T,N>::shift(const IndexTM<T,N> & a_iv)
  {
    return (*this) += a_iv;
  }

  template<typename T, int N> inline IndexTM<T,N>& IndexTM<T,N>::diagShift(T a_s)
  {
    return (*this) += a_s;
  }

  template<typename T, int N> inline IndexTM<T,N> scale (const IndexTM<T,N> & a_p,
                                                         T                   a_s)
  {
    return a_p * a_s;
  }

  template<typename T, int N> inline IndexTM<T,N> diagShift(const IndexTM<T,N> & a_p,
                                                            T                    a_s)
  {
    return a_p + a_s;
  }

  template<typename T, int N> inline IndexTM<T,N> min(const IndexTM<T,N> & a_p1,
                                                      const IndexTM<T,N> & a_p2)
  {
    IndexTM<T,N> result(a_p1);
    return result.min(a_p2);
  }

  template<typename T, int N> inline IndexTM<T,N> max(const IndexTM<T,N> & a_p1,
                                                      const IndexTM<T,N> & a_p2)
  {
    IndexTM<T,N> result(a_p1);
    return result.max(a_p2);
  }

  template<typename T, int N> inline IndexTM<T,N> BASISV_TM(int a_dir)
  {
    CH_assert(a_dir >= 0 && a_dir < N);
    IndexTM<T,N> tmp = IndexTM<T,N>::Zero ;
    tmp.dataPtr()[a_dir] = T(1);
    return tmp;
  }

  template<typename T, int N> inline IndexTM<T,N> reflect(const IndexTM<T,N> & a_a,
                                                          T                    a_refIx,
                                                          int                  a_idir)
  {
    IndexTM<T,N> result(a_a);
    return result.reflect(a_refIx, a_idir);
  }

  template<typename T> static T ourcoarsen(T a_a,
                                           T a_b)
  {
    return (a_a < 0 ? T(-Abs(a_a+1))/a_b-1 : a_a/a_b);
  }

  template<typename T, int N> inline IndexTM<T,N> coarsen(const IndexTM<T,N> & a_p,
                                                          T                    a_s)
  {
    IndexTM<T,N> result(a_p);
    return result.coarsen(a_s);
  }

  template<typename T, int N> inline IndexTM<T,N> coarsen(const IndexTM<T,N> & a_p1,
                                                          const IndexTM<T,N> & a_p2)
  {
    IndexTM<T,N> result(a_p1);
    return result.coarsen(a_p2);
  }

/*
  template<typename T, int N> IndexTM<T,N>::operator IndexTM<Real,N>()
  {
  IndexTM<Real,N> result;
  for (int i = 0; i < N; ++i)
  {
  result.dataPtr()[i] = Real(m_vect[i]);
  }
  return result;
  }
*/

  template<typename T, int N> inline IndexTM<T,N>& IndexTM<T,N>::coarsen(T a_s)
  {
    CH_assert(a_s > 0);
    std::transform(m_vect, m_vect+N, m_vect,
                   std::bind2nd(std::ptr_fun(ourcoarsen<T>),a_s));
    return *this;
  }

  template<typename T, int N> inline IndexTM<T,N>& IndexTM<T,N>::coarsen(const IndexTM<T,N> & a_p)
  {
    CH_assert(a_p > (IndexTM<T,N>::Zero));
    std::transform(m_vect, m_vect+N, a_p.m_vect, m_vect, std::ptr_fun(ourcoarsen<T>));
    return *this;
  }

  template<typename T, int N> void IndexTM<T,N>::linearIn(const void* a_inBuf)
  {
    memcpy(m_vect, (T*)a_inBuf, N*sizeof(T));
  }

  template<typename T, int N> void IndexTM<T,N>::linearOut(void* a_outBuf) const
  {
    memcpy((T*)a_outBuf, m_vect, N*sizeof(T));
  }

  template<typename T, int N> IndexTM<T,N>::IndexTM(T a_i)
    : GenericArithmeticable< T, IndexTM<T,N> >(this)
  {
    STATIC_ASSERT(N == 1);
    m_vect[0] = a_i;
  }

  template<typename T, int N> IndexTM<T,N>::IndexTM(T a_i,
                                                    T a_j)
    : GenericArithmeticable< T, IndexTM<T,N> >(this)
  {
    STATIC_ASSERT(N == 2);
    m_vect[0] = a_i;
    m_vect[1] = a_j;
  }

  template<typename T, int N> IndexTM<T,N>::IndexTM(T a_i,
                                                    T a_j,
                                                    T a_k)
    : GenericArithmeticable< T, IndexTM<T,N> >(this)
  {
    STATIC_ASSERT(N == 3);
    m_vect[0] = a_i;
    m_vect[1] = a_j;
    m_vect[2] = a_k;
  }

  template<typename T, int N> IndexTM<T,N>::IndexTM(T a_i,
                                                    T a_j,
                                                    T a_k,
                                                    T a_l)
    : GenericArithmeticable< T, IndexTM<T,N> >(this)
  {
    STATIC_ASSERT(N == 4);
    m_vect[0] = a_i;
    m_vect[1] = a_j;
    m_vect[2] = a_k;
    m_vect[3] = a_l;
  }

  template<typename T, int N> IndexTM<T,N>::IndexTM(T a_i,
                                                    T a_j,
                                                    T a_k,
                                                    T a_l,
                                                    T a_m)
    : GenericArithmeticable< T, IndexTM<T,N> >(this)
  {
    STATIC_ASSERT(N == 5);
    m_vect[0] = a_i;
    m_vect[1] = a_j;
    m_vect[2] = a_k;
    m_vect[3] = a_l;
    m_vect[4] = a_m;
  }

  template<typename T, int N> IndexTM<T,N>::IndexTM(T a_i,
                                                    T a_j,
                                                    T a_k,
                                                    T a_l,
                                                    T a_m,
                                                    T a_n)
    : GenericArithmeticable< T, IndexTM<T,N> >(this)
  {
    STATIC_ASSERT(N == 6);
    m_vect[0] = a_i;
    m_vect[1] = a_j;
    m_vect[2] = a_k;
    m_vect[3] = a_l;
    m_vect[4] = a_m;
    m_vect[5] = a_n;
  }


}

#endif // include guard
