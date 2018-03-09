//
// For any class that has operator+=, defines a non-member operator+.  Ditto
// operator-=, *=, /=.
// To use these functions, your class must derive from
// class GenericArithmeticable (defined below), and furthermore define
// a unary operator-(), and a reciprocal() function.
//
// All of this is done more elegantly and completely in boost/operators.hpp;
// let's keep that in mind if we ever decide to use Boost.
//
// Author: Ted
//

#ifndef _GENERICARITHMETIC_H_
#define _GENERICARITHMETIC_H_


namespace amrex
{
/** Class that you need to derive from, if you want the *global* operators below
 *  to be available for your class.
 *  For example:
 *  class Foo : public GenericArithmeticable<double,Foo>
 *  {
 *    public:
 *      Foo() : GenericArithmeticable<double,Foo>(this)
 *      {
 *      }
 *      Foo& operator+=( const Foo& );
 *      Foo& operator+=( const double& );
 *      [etc]
 *  };

 *
 * If, in addition, you want the member-function operators, you'll need to
 * implement some generic operator-appliers in your derived class: see, in
 * GenericArithmeticI.H, operatorCompare() and operatorOpEquals().
 * If you implement operator<, operator+= etc in your derived class, then that'll
 * override the ones here (regardless of whether you've defined operatorCompare()
 * and operatorOpEquals() in your derived class, though it would then be weird if
 * you did).
 *
 */
  template<typename ScalarT, typename SelfT>
  struct GenericArithmeticable
  {
    // Intel 8.1 icpc wants this
    virtual ~GenericArithmeticable()
      {
      }
    typedef ScalarT scalar_type;
    typedef SelfT   self_type;

  public:
    GenericArithmeticable( SelfT * s ) : m_child( s )
      {
      }

    bool operator<(const SelfT&) const;
    bool operator>(const SelfT&) const;
    bool operator<=(const SelfT&) const;
    bool operator>=(const SelfT&) const;

    SelfT& operator+=(const SelfT&);
    SelfT& operator+=(const ScalarT&);
    SelfT& operator-=(const SelfT&);
    SelfT& operator-=(const ScalarT&);
    SelfT& operator*=(const SelfT&);
    SelfT& operator*=(const ScalarT&);
    SelfT& operator/=(const SelfT&);
    SelfT& operator/=(const ScalarT&);

  private:
    SelfT * m_child;
  };

// The self_type is there so we can restrict these to classes that define
// that; otherwise, this is just too promiscuous; when a class already has
// an operator+ member function, the compiler considers the existence of
// that one, and the one here, to be an intolerable ambiguity.
// C::self_type should just be C itself.
// C::scalar_type should be the kind of scalar that you hope to add, multiply
// etc against a C.  For Chombo::IntVect, for example, scalar_type is int.
  template<class C> typename C::self_type operator+( const C&, const C& );
  template<class C> typename C::self_type operator+( const C&, const typename C::scalar_type& );
  template<class C> typename C::self_type operator+( const typename C::scalar_type&, const C& );

  template<class C> typename C::self_type operator*( const C&, const C& );
  template<class C> typename C::self_type operator*( const C&, const typename C::scalar_type& );
  template<class C> typename C::self_type operator*( const typename C::scalar_type&, const C& );

  template<class C> typename C::self_type operator-( const C&, const C& );
  template<class C> typename C::self_type operator-( const C&, const typename C::scalar_type& );
  template<class C> typename C::self_type operator-( const typename C::scalar_type&, const C& );

  template<class C> typename C::self_type operator/( const C&, const C& );
  template<class C> typename C::self_type operator/( const C&, const typename C::scalar_type& );
  template<class C> typename C::self_type operator/( const typename C::scalar_type&, const C& );

}

#include "AMReX_GenericArithmeticI.H"

#endif // include guard
