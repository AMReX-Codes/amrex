#ifndef AMREX_GENERICARITHMETICI_H_
#define AMREX_GENERICARITHMETICI_H_

#include <functional>

//
// Comparison operators
//

  template<typename ScalarT, typename SelfT> bool
  GenericArithmeticable<ScalarT,SelfT>::operator<( const SelfT& rhs ) const
  {
    return m_child->operatorCompare( rhs, std::less<ScalarT>() );
  }

  template<typename ScalarT, typename SelfT> bool
  GenericArithmeticable<ScalarT,SelfT>::operator>( const SelfT& rhs ) const
  {
    return m_child->operatorCompare( rhs, std::greater<ScalarT>() );
  }

  template<typename ScalarT, typename SelfT> bool
  GenericArithmeticable<ScalarT,SelfT>::operator<=( const SelfT& rhs ) const
  {
    return m_child->operatorCompare( rhs, std::less_equal<ScalarT>() );
  }

  template<typename ScalarT, typename SelfT> bool
  GenericArithmeticable<ScalarT,SelfT>::operator>=( const SelfT& rhs ) const
  {
    return m_child->operatorCompare( rhs, std::greater_equal<ScalarT>() );
  }

//
// In-place arithmetic operators
//
  template<typename ScalarT, typename SelfT> SelfT&
  GenericArithmeticable<ScalarT,SelfT>::operator+=(const SelfT& rhs)
  {
    return m_child->operatorOpEquals( rhs, std::plus<ScalarT>() );
  }

  template<typename ScalarT, typename SelfT> SelfT&
  GenericArithmeticable<ScalarT,SelfT>::operator+=(const ScalarT& a)
  {
    return m_child->operatorOpEquals( a, std::plus<ScalarT>() );
  }

  template<typename ScalarT, typename SelfT> SelfT&
  GenericArithmeticable<ScalarT,SelfT>::operator-=(const SelfT& rhs)
  {
    return m_child->operatorOpEquals( rhs, std::minus<ScalarT>() );
  }

  template<typename ScalarT, typename SelfT> SelfT&
  GenericArithmeticable<ScalarT,SelfT>::operator-=(const ScalarT& a)
  {
    return m_child->operatorOpEquals( a, std::minus<ScalarT>() );
  }

  template<typename ScalarT, typename SelfT> SelfT&
  GenericArithmeticable<ScalarT,SelfT>::operator*=(const SelfT& rhs)
  {
    return m_child->operatorOpEquals( rhs, std::multiplies<ScalarT>() );
  }

  template<typename ScalarT, typename SelfT> SelfT&
  GenericArithmeticable<ScalarT,SelfT>::operator*=(const ScalarT& a)
  {
    return m_child->operatorOpEquals( a, std::multiplies<ScalarT>() );
  }

  template<typename ScalarT, typename SelfT> SelfT&
  GenericArithmeticable<ScalarT,SelfT>::operator/=(const SelfT& rhs)
  {
    return m_child->operatorOpEquals( rhs, std::divides<ScalarT>() );
  }

  template<typename ScalarT, typename SelfT> SelfT&
  GenericArithmeticable<ScalarT,SelfT>::operator/=(const ScalarT& a)
  {
    return m_child->operatorOpEquals( a, std::divides<ScalarT>() );
  }

  template< class C >
  C operatorOp( const C& c1, const C& c2, C& (C::*op)(const C&) )
  {
    C result( c1 );
    (result.*op)( c2 );
    return result;
  }

// Same, but for scalar second arg.
  template< class C >
  C operatorOpScalar( const C& c1, const typename C::scalar_type& x,
                      C& (C::*op)(const typename C::scalar_type&) )
  {
    C result( c1 );
    (result.*op)( x );
    return result;
  }

//
// Args C and C
//
  template<class C> typename C::self_type operator+( const C& c1, const C& c2 )
  {
    return operatorOp< C >( c1, c2, &C::operator+= );
  }
  template<class C> typename C::self_type operator-( const C& c1, const C& c2 )
  {
    return operatorOp< C >( c1, c2, &C::operator-= );
  }
  template<class C> typename C::self_type operator*( const C& c1, const C& c2 )
  {
    return operatorOp< C >( c1, c2, &C::operator*= );
  }
  template<class C> typename C::self_type operator/( const C& c1, const C& c2 )
  {
    return operatorOp< C >( c1, c2, &C::operator/= );
  }

//
// Args C and const C::scalar_type&
//
  template<class C> typename C::self_type operator+( const C& c, const typename C::scalar_type& x )
  {
    return operatorOpScalar<C>( c, x, &C::operator+= );
  }
  template<class C> typename C::self_type operator+( const typename C::scalar_type& x, const C& c )
  {
    return c + x;
  }

  template<class C> typename C::self_type operator*( const C& c, const typename C::scalar_type& x )
  {
    return operatorOpScalar<C>( c, x, &C::operator*= );
  }
  template<class C> typename C::self_type operator*( const typename C::scalar_type& x, const C& c )
  {
    return c * x;
  }

  template<class C> typename C::self_type operator-( const C& c, const typename C::scalar_type& x )
  {
    return operatorOpScalar<C>( c, x, &C::operator-= );
  }
  template<class C> typename C::self_type operator-( const typename C::scalar_type& x, const C& c )
  {
    C temp( c );
    return (-temp) + x;
  }

  template<class C> typename C::self_type operator/( const C& c, const typename C::scalar_type& x )
  {
    return operatorOpScalar<C>( c, x, &C::operator/= );
  }
  template<class C> typename C::self_type operator/( const typename C::scalar_type& x, const C& c )
  {
    C temp( c );
    return temp.reciprocal() * x;
    // Warning: if x is an int, then temp.reciprocal() will always be 0
    // and the return value from this function will also always be 0.
  }


#endif // include guard

