/*
 *      .o.       ooo        ooooo ooooooooo.             ooooooo  ooooo 
 *     .888.      `88.       .888' `888   `Y88.            `8888    d8'  
 *    .8"888.      888b     d'888   888   .d88'  .ooooo.     Y888..8P    
 *   .8' `888.     8 Y88. .P  888   888ooo88P'  d88' `88b     `8888'     
 *  .88ooo8888.    8  `888'   888   888`88b.    888ooo888    .8PY888.    
 * .8'     `888.   8    Y     888   888  `88b.  888    .o   d8'  `888b   
 *o88o     o8888o o8o        o888o o888o  o888o `Y8bod8P' o888o  o88888o 
 *
 */
#ifndef AMREX_REALVECT_H_
#define AMREX_REALVECT_H_

#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <cmath>

#include "AMReX_Box.H"
#include "AMReX_REAL.H"
#include "AMReX_SPACE.H"
#include <vector>
#include <cassert>
#include "AMReX_IntVect.H"


//template<typename T, int n>
//class IndexTM;

/// A Real vector in SpaceDim-dimensional space
/**
  The class RealVect is an implementation of a Real vector in a
  SpaceDim-dimensional space.
  RealVect values are accessed using the operator[] function, as for a normal
  C++ array.  In addition, the basic arithmetic operators have been overloaded
  to implement scaling and translation operations.
*/

namespace amrex
{

class RealVect
{
public:

  /**
     \name Constructors and Accessors
  */
  /*@{*/

  ///
  /**
     Construct a RealVect whose components are zero.
  */
  RealVect ();

  explicit RealVect (const std::vector<Real>&);

  ///
  /**
     Construct a RealVect given the specific values for its
     coordinates.  AMREX_D_DECL is a macro that sets the constructor to
     take AMREX_SPACEDIM arguments.
  */
  RealVect (AMREX_D_DECL(Real i, Real j, Real k));

  ///
  /**
     The copy constructor.
  */
  RealVect (const RealVect& rhs);


  ///
  /**
     Construct a RealVect from an IntVect by coercing each component
     from <tt>int</tt> to Real.
  */
  RealVect (const IntVect & iv)
  {
    for (int d=0 ; d<SpaceDim ; ++d)
    {
      vect[d] = (Real)iv[d];
    }
  }

  ///
  /**
     The assignment operator.
  */
  RealVect& operator= (const RealVect& rhs);

  ///
  /**
     Returns a modifiable lvalue reference to the <i>i</i>'th coordinate of the
     RealVect.
  */
  inline
  Real& operator[] (int i);

  ///
  /**
     Returns the <i>i</i>'th coordinate of the RealVect.
  */
  inline
  const Real& operator[] (int i) const;

  /*@}*/

  /**
     \name Comparison Operators
  */
  /*@{*/

  ///
  /**
     Returns true if this RealVect is equivalent to argument RealVect.  All
     comparisons between analogous components must be satisfied.
  */
  bool operator== (const RealVect& p) const;

  ///
  /**
     Returns true if this RealVect is different from argument RealVect.
     All comparisons between analogous components must be satisfied.
  */
  bool operator!= (const RealVect& p) const;

  ///
  /**
     Returns true if this RealVect is less than argument RealVect.  All
     comparisons between analogous components must be satisfied.  Note
     that, since the comparison is component-wise, it is possible for
     an RealVect to be neither greater than, less than, nor equal to
     another.
  */
  bool operator< (const RealVect& p) const;

  ///
  /**
     Returns true if this RealVect is less than or equal to argument
     RealVect.  All comparisons between analogous components must be
     satisfied.  Note that, since the comparison is component-wise, it
     is possible for an RealVect to be neither greater than or equal
     to, less than or equal to, nor equal to another.
  */
  bool operator<= (const RealVect& p) const;

  ///
  /**
     Returns true if this RealVect is greater than argument RealVect.
     All comparisons between analogous components must be satisfied.
     Note that, since the comparison is component-wise, it is possible
     for an RealVect to be neither greater than, less than, nor equal
     to another.
  */
  bool operator> (const RealVect& p) const;

  ///
  /**
     Returns true if this RealVect is greater than or equal to argument
     RealVect.  All comparisons between analogous components must be
     satisfied.  Note that, since the comparison is component-wise, it
     is possible for an RealVect to be neither greater than or equal
     to, less than or equal to, nor equal to another.
  */

  bool operator>= (const RealVect& p) const;

  /*@}*/

  /**
     \name Arithmetic Operators
  */
  /*@{*/

  ///
  /**
     Modifies this RealVect by addition of a scalar to each component.
  */
  RealVect& operator+= (Real s);

  ///
  /**
     Returns a RealVect that is this RealVect with a scalar s added
     to each component.
  */
  RealVect operator+ (Real s) const;

  ///
  /**
     Modifies this RealVect by component-wise addition by argument.
  */
  RealVect& operator+= (const RealVect& p);

  ///
  /**
     Modifies this RealVect by subtraction of a scalar from each
     component.
  */
  RealVect& operator-= (Real s);

  ///
  /**
     Modifies this RealVect by component-wise subtraction by argument.
  */
  RealVect& operator-= (const RealVect& p);

  ///
  /**
     Returns a RealVect that is this RealVect with a scalar s subtracted
     from each component.
  */
  RealVect operator- (Real s) const;

  ///
  /**
     Modifies this RealVect by multiplying each component by a scalar.
  */
  RealVect& operator*= (Real s);

  ///
  /**

  */
  Real dotProduct(const RealVect& a_rhs) const;

  ///
  /**
     Modifies this RealVect by component-wise multiplication by argument.
  */
  RealVect& operator*= (const RealVect& p);

//XXX  ///
//XXX  /**
//XXX     Returns component-wise product of this RealVect and argument.
//XXX  */
//XXX  RealVect operator* (const RealVect& p) const;

  ///
  /**
     Returns a RealVect that is this RealVect with each component
     multiplied by a scalar.
  */
  RealVect operator* (Real s) const;

  ///
  /**
     Modifies this RealVect by dividing each component by a scalar.
  */
  RealVect& operator/= (Real s);

  ///
  /**
     Modifies this RealVect by component-wise division by argument.
  */
  RealVect& operator/= (const RealVect& p);

//XXX  ///
//XXX  /**
//XXX     Returns component-wise quotient of this RealVect by argument.
//XXX  */
//XXX  RealVect operator/ (const RealVect& p) const;

  ///
  /**
     Returns a RealVect that is this RealVect with each component
     divided by a scalar.
  */
  RealVect operator/ (Real s) const;

  ///
  /**
     Modifies this RealVect by multiplying each component by a scalar.
  */
  RealVect& scale (Real s);

  /*@}*/

  /**
     \name Other arithmetic operators
  */
  /*@{*/

  ///
  /**
     Modifies this RealVect by taking component-wise min with RealVect
     argument.
  */
  RealVect& min (const RealVect& p);

  ///
  /**
     Returns the RealVect that is the component-wise minimum of two
     argument RealVects.
  */
  friend inline RealVect min (const RealVect& p1,
                              const RealVect& p2);

  ///
  /**
     Modifies this RealVect by taking component-wise max with RealVect
     argument.
  */
  RealVect& max (const RealVect& p);

  ///
  /**
     Returns the RealVect that is the component-wise maximum of two
     argument RealVects.
  */
  friend inline RealVect max (const RealVect& p1,
                              const RealVect& p2);

  /*@}*/

  /**
     \name Unary operators
  */
  /*@{*/

  ///
  /**
     Unary plus -- for completeness.
  */
  RealVect operator+ () const;

  ///
  /**
     Unary minus -- negates all components of this RealVect.
  */
  RealVect operator- () const;

  ///
  /**
     Sum of all components of this RealVect.
  */
  Real sum () const;

  ///
  /**
     sqrt(sum squares)
  */
  Real vectorLength() const;

  ///
  /**
     sum squares--no square root
  */
  Real radSquared() const;

  ///
  /**
     Product of all components of this RealVect.
  */
  Real product () const;

  ///
  /**
     Component with the minimum value of this RealVect (returns 0 if they are all the same).
     a_doAbs : if true then take the absolute value before comparing
  */
  int minDir(const bool& a_doAbs) const;

  ///
  /**
     Component with the maximum value of this RealVect (returns 0 if they are all the same).
     a_doAbs : if true then take the absolute value before comparing
  */
  int maxDir(const bool& a_doAbs) const;

  /*@}*/

  /**
     \name Data pointer functions
  */
  /*@{*/

  ///
  /**
     Only for sending stuff to Fortran
   */
  const Real* dataPtr() const;

  ///
  /**
     Only for sending stuff to Fortran
   */
  Real* dataPtr() ;

  /*@}*/

  /**
     \name Constants
  */
  /*@{*/

  ///
  /**
     Returns a basis vector in the given coordinate direction.<br>
     In 2-D:<br>
     BASISREALV(0) == (1.,0.);
     BASISREALV(1) == (0.,1.).<br>
     In 3-D:<br>
     BASISREALV(0) == (1.,0.,0.);
     BASISREALV(1) == (0.,1.,0.);
     BASISREALV(2) == (0.,0.,1.).<br>
     Note that the coordinate directions are based at zero.
  */
  friend RealVect BASISREALV(int dir);

  ///
  /**
     This is a RealVect all of whose components are equal to zero.
  */
  static const RealVect Zero;

  ///
  /**
     This is a RealVect all of whose components are equal to one.
  */
  static const RealVect Unit;

  /*@}*/

  /**
     \name Arithmetic friend functions
  */
  /*@{*/

  ///
  /**
     Returns a RealVect that is a RealVect <i>p</i> with
     a scalar <i>s</i> added to each component.
  */
  friend RealVect operator+ (Real            s,
                             const RealVect& p);

  ///
  /**
     Returns <i>s - p</i>.
  */
  friend RealVect operator- (Real            s,
                             const RealVect& p);

  ///
  /**
     Returns a RealVect that is a RealVect <i>p</i> with each component
     multiplied by a scalar <i>s</i>.
  */
  friend RealVect operator* (Real            s,
                             const RealVect& p);
  ///
  /**
     Returns a RealVect that is a RealVect <i>p</i> with each component
     divided by a scalar <i>s</i>.
  */
  friend RealVect operator/ (Real            s,
                             const RealVect& p);

  ///
  /**
     Returns component-wise sum of RealVects <i>s</i> and <i>p</i>.
  */
  friend RealVect operator+ (const RealVect& s,
                             const RealVect& p);

  ///
  /**
     Returns <i>s - p</i>.
  */
  friend RealVect operator- (const RealVect& s,
                             const RealVect& p);

  ///
  /**
     Returns component-wise product of <i>s</i> and <i>p</i>.
  */
  friend RealVect operator* (const RealVect& s,
                             const RealVect& p);
  ///
  /**
     Returns component-wise quotient <i>p / s</i>.
  */
  friend RealVect operator/ (const RealVect& s,
                             const RealVect& p);

  ///
  /**
     Returns a RealVect obtained by multiplying each of the components
     of the given RealVect by a scalar.
  */
  friend inline RealVect scale (const RealVect& p,
                                Real            s);

  /*@}*/

  ///
  /**
     Print to the given output stream in ASCII.
  */
  friend std::ostream& operator<< (std::ostream& ostr,
                                   const RealVect& p);

  friend std::istream& operator>> (std::istream& is, RealVect& iv);
  friend class HDF5Handle;

  static size_t io_offset;

protected:

  /**
     The individual components of this RealVect.
   */
  Real vect[SpaceDim];

};

inline Real& RealVect::operator[] (int i)
{
  assert(i>=0 && i < SpaceDim);
  return vect[i];
}

inline const Real& RealVect::operator[] (int i) const
{
  assert(i>=0 && i < SpaceDim);
  return vect[i];
}

inline RealVect::RealVect (const RealVect &iv)
{
  AMREX_D_EXPR(vect[0]=iv.vect[0], vect[1]=iv.vect[1], vect[2]=iv.vect[2]);
}

inline
RealVect&
RealVect::operator-= (Real s)
{
  AMREX_D_EXPR(vect[0] -= s, vect[1] -= s, vect[2] -= s);
  return *this;
}

inline
RealVect&
RealVect::operator-= (const RealVect& p)
{
  AMREX_D_EXPR(vect[0] -= p[0], vect[1] -= p[1], vect[2] -= p[2]);

  return *this;
}

inline
RealVect
RealVect::operator+ () const
{
  return RealVect(*this);
}

inline
RealVect
RealVect::operator- () const
{
  return RealVect(AMREX_D_DECL(-vect[0], -vect[1], -vect[2]));
}

inline
RealVect&
RealVect::scale (Real s)
{
  AMREX_D_EXPR(vect[0] *= s, vect[1] *= s, vect[2] *= s);
  return *this;
}

inline
Real
RealVect::sum () const
{
  return AMREX_D_TERM(vect[0], + vect[1], + vect[2]);
}

inline
Real
RealVect::vectorLength () const
{
  Real len = this->radSquared();
  len = std::sqrt(len);

  return len;
}

inline
Real
RealVect::radSquared() const
{
  Real len = 0;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      len = len + vect[idir]*vect[idir];
    }

  return len;
}

inline
Real
RealVect::product () const
{
  return AMREX_D_TERM(vect[0], * vect[1], * vect[2]);
}

inline
RealVect
scale (const RealVect& p,
       Real            s)
{
  return RealVect(AMREX_D_DECL(s * p[0], s * p[1], s * p[2]));
}

inline
bool
RealVect::operator< (const RealVect& p) const
{
  return AMREX_D_TERM(vect[0] < p[0], && vect[1] < p[1], && vect[2] < p[2]);
}

inline
bool
RealVect::operator<= (const RealVect& p) const
{
  return AMREX_D_TERM(vect[0] <= p[0], && vect[1] <= p[1], && vect[2] <= p[2]);
}

inline
bool
RealVect::operator> (const RealVect& p) const
{
  return AMREX_D_TERM(vect[0] > p[0], && vect[1] > p[1], && vect[2] > p[2]);
}

inline
bool
RealVect::operator>= (const RealVect& p) const
{
  return AMREX_D_TERM(vect[0] >= p[0], && vect[1] >= p[1], && vect[2] >= p[2]);
}

inline
RealVect&
RealVect::min (const RealVect& p)
{
  AMREX_D_EXPR(vect[0] = std::min(vect[0], p.vect[0]),
         vect[1] = std::min(vect[1], p.vect[1]),
         vect[2] = std::min(vect[2], p.vect[2]));

  return *this;
}

inline
RealVect&
RealVect::max (const RealVect& p)
{
  AMREX_D_EXPR(vect[0] = std::max(vect[0], p.vect[0]),
         vect[1] = std::max(vect[1], p.vect[1]),
         vect[2] = std::max(vect[2], p.vect[2]));
  return *this;
}

inline
RealVect
min (const RealVect& p1,
     const RealVect& p2)
{
  RealVect p(p1);
  return p.min(p2);
}

inline
RealVect
max (const RealVect& p1,
     const RealVect& p2)
{
  RealVect p(p1);
  return p.max(p2);
}

extern RealVect BASISREALV(int idir);

}

#endif
