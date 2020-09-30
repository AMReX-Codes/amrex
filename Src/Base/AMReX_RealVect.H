#ifndef AMREX_REALVECT_H_
#define AMREX_REALVECT_H_

#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <cmath>

#include <AMReX_Box.H>
#include <AMReX_REAL.H>
#include <AMReX_SPACE.H>
#include <vector>
#include <AMReX_IntVect.H>
#include <AMReX_Utility.H>
#include <AMReX_Math.H>

namespace amrex
{

/**
  \brief A Real vector in SpaceDim-dimensional space

  The class RealVect is an implementation of a Real vector in a
  SpaceDim-dimensional space.
  RealVect values are accessed using the operator[] function, as for a normal
  C++ array.  In addition, the basic arithmetic operators have been overloaded
  to implement scaling and translation operations.
*/

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
  AMREX_GPU_HOST_DEVICE
  constexpr RealVect () noexcept : vect{AMREX_D_DECL(0.,0.,0.)} {}

  explicit RealVect (const std::vector<Real>& vr) noexcept : vect{AMREX_D_DECL(vr[0],vr[1],vr[2])} {
      BL_ASSERT(vr.size() == AMREX_SPACEDIM);
  }

  ///
  /**
     Construct a RealVect given the specific values for its
     coordinates.  AMREX_D_DECL is a macro that sets the constructor to
     take AMREX_SPACEDIM arguments.
  */
#if (AMREX_SPACEDIM > 1)
  AMREX_GPU_HOST_DEVICE
  constexpr RealVect (AMREX_D_DECL(Real i, Real j, Real k)) noexcept : vect{AMREX_D_DECL(i,j,k)} {}
#endif

  AMREX_GPU_HOST_DEVICE
  explicit constexpr RealVect (Real i) noexcept : vect{AMREX_D_DECL(i,i,i)} {}

  ///
  /**
     Construct a RealVect setting the coordinates to the corresponding values in
     the Real array a.
  */
  AMREX_GPU_HOST_DEVICE
  explicit RealVect (const Real* a) noexcept : vect{AMREX_D_DECL(a[0],a[1],a[2])} {}

  ///
  /**
     Construct a RealVect from an IntVect by coercing each component
     from <tt>int</tt> to Real.
  */
  AMREX_GPU_HOST_DEVICE
  RealVect (const IntVect & iv) noexcept : vect{AMREX_D_DECL((Real)iv[0],(Real)iv[1],(Real)iv[2])} {}

  ///
  /**
     Prevent rvalue assignment to the <i>i</i>'th coordinate of the
     RealVect.
  */
  Real& operator[] (int i) && = delete;

  ///
  /**
     Returns a modifiable lvalue reference to the <i>i</i>'th coordinate of the
     RealVect.
  */
  AMREX_GPU_HOST_DEVICE
  inline
  Real& operator[] (int i) & noexcept
  { AMREX_ASSERT(i>=0 && i < SpaceDim); return vect[i]; }

  ///
  /**
     Returns the <i>i</i>'th coordinate of the RealVect.
  */
  AMREX_GPU_HOST_DEVICE
  inline
  const Real& operator[] (int i) const& noexcept
  { AMREX_ASSERT(i>=0 && i < SpaceDim); return vect[i]; }

  /*@}*/

  /**
     \name Iterators
  */
  /*@{*/

  ///
  /**
     Return pointer to the first component.
  */
  AMREX_GPU_HOST_DEVICE
  inline
  Real* begin () noexcept { return vect; }
  
  ///
  /**
     Return pointer to the first component.
  */
  AMREX_GPU_HOST_DEVICE
  inline
  const Real* begin () const noexcept { return vect; }
  
  ///
  /**
     Return pointer to the one past last component
  */
  AMREX_GPU_HOST_DEVICE
  inline
  Real* end () noexcept { return vect + AMREX_SPACEDIM; }
  
  ///
  /**
     Return pointer to the one past last component
  */
  AMREX_GPU_HOST_DEVICE
  inline
  const Real* end () const noexcept { return vect + AMREX_SPACEDIM; }
  
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
  AMREX_GPU_HOST_DEVICE
  bool operator== (const RealVect& p) const noexcept {
      return AMREX_D_TERM(vect[0] == p[0], && vect[1] == p[1], && vect[2] == p[2]);
  }

  ///
  /**
     Returns true if this RealVect is different from argument RealVect.
     All comparisons between analogous components must be satisfied.
  */
  AMREX_GPU_HOST_DEVICE
  bool operator!= (const RealVect& p) const noexcept {
      return AMREX_D_TERM(vect[0] != p[0], || vect[1] != p[1], || vect[2] != p[2]);
  }

  ///
  /**
     Returns true if this RealVect is less than argument RealVect.  All
     comparisons between analogous components must be satisfied.  Note
     that, since the comparison is component-wise, it is possible for
     an RealVect to be neither greater than, less than, nor equal to
     another.
  */
  AMREX_GPU_HOST_DEVICE inline
  bool operator< (const RealVect& p) const noexcept;

  ///
  /**
     Returns true if this RealVect is less than or equal to argument
     RealVect.  All comparisons between analogous components must be
     satisfied.  Note that, since the comparison is component-wise, it
     is possible for an RealVect to be neither greater than or equal
     to, less than or equal to, nor equal to another.
  */
  AMREX_GPU_HOST_DEVICE inline
  bool operator<= (const RealVect& p) const noexcept;

  ///
  /**
     Returns true if this RealVect is greater than argument RealVect.
     All comparisons between analogous components must be satisfied.
     Note that, since the comparison is component-wise, it is possible
     for an RealVect to be neither greater than, less than, nor equal
     to another.
  */
  AMREX_GPU_HOST_DEVICE inline
  bool operator> (const RealVect& p) const noexcept;

  ///
  /**
     Returns true if this RealVect is greater than or equal to argument
     RealVect.  All comparisons between analogous components must be
     satisfied.  Note that, since the comparison is component-wise, it
     is possible for an RealVect to be neither greater than or equal
     to, less than or equal to, nor equal to another.
  */

  AMREX_GPU_HOST_DEVICE inline
  bool operator>= (const RealVect& p) const noexcept;

  /*@}*/

  /**
     \name Arithmetic Operators
  */
  /*@{*/

  ///
  /**
     Modifies this RealVect by addition of a scalar to each component.
  */
  AMREX_GPU_HOST_DEVICE inline
  RealVect& operator+= (Real s) noexcept;

  ///
  /**
     Returns a RealVect that is this RealVect with a scalar s added
     to each component.
  */
  AMREX_GPU_HOST_DEVICE inline
  RealVect operator+ (Real s) const noexcept;

  ///
  /**
     Modifies this RealVect by component-wise addition by argument.
  */
  AMREX_GPU_HOST_DEVICE inline
  RealVect& operator+= (const RealVect& p) noexcept;

  ///
  /**
     Modifies this RealVect by subtraction of a scalar from each
     component.
  */
  AMREX_GPU_HOST_DEVICE inline
  RealVect& operator-= (Real s) noexcept;

  ///
  /**
     Modifies this RealVect by component-wise subtraction by argument.
  */
  AMREX_GPU_HOST_DEVICE inline
  RealVect& operator-= (const RealVect& p) noexcept;

  ///
  /**
     Returns a RealVect that is this RealVect with a scalar s subtracted
     from each component.
  */
  AMREX_GPU_HOST_DEVICE inline
  RealVect operator- (Real s) const noexcept;

  ///
  /**
     Modifies this RealVect by multiplying each component by a scalar.
  */
  AMREX_GPU_HOST_DEVICE inline
  RealVect& operator*= (Real s) noexcept;

  ///
  /**
     Returns the result of the scalar product with another RealVect
  */
  AMREX_GPU_HOST_DEVICE inline
  Real dotProduct(const RealVect& a_rhs) const noexcept;

#if (AMREX_SPACEDIM == 3)
  ///
  /**
     Returns the result of the cross product with another RealVect
  */
  AMREX_GPU_HOST_DEVICE inline
  RealVect crossProduct(const RealVect& a_rhs) const noexcept;
#endif

  ///
  /**
     Modifies this RealVect by component-wise multiplication by argument.
  */
  AMREX_GPU_HOST_DEVICE inline
  RealVect& operator*= (const RealVect& p) noexcept;

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
  AMREX_GPU_HOST_DEVICE inline
  RealVect operator* (Real s) const noexcept;

  ///
  /**
     Modifies this RealVect by dividing each component by a scalar.
  */
  AMREX_GPU_HOST_DEVICE inline
  RealVect& operator/= (Real s) noexcept;

  ///
  /**
     Modifies this RealVect by component-wise division by argument.
  */
  AMREX_GPU_HOST_DEVICE inline
  RealVect& operator/= (const RealVect& p) noexcept;

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
  AMREX_GPU_HOST_DEVICE inline
  RealVect operator/ (Real s) const noexcept;

  ///
  /**
     Modifies this RealVect by multiplying each component by a scalar.
  */
  AMREX_GPU_HOST_DEVICE inline
  RealVect& scale (Real s) noexcept;

  ///
  /**
     Returns an IntVect whose components are the std::floor of the RealVect 
     components.
  */
  AMREX_GPU_HOST_DEVICE inline
  IntVect floor () const noexcept;

  ///
  /**
     Returns an IntVect whose components are the std::floor of the RealVect 
     components.
  */
  AMREX_GPU_HOST_DEVICE inline
  IntVect ceil () const noexcept;

  ///
  /**
     Returns an IntVect whose components are the std::floor of the RealVect 
     components.
  */
  AMREX_GPU_HOST_DEVICE inline
  IntVect round () const noexcept;

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
  AMREX_GPU_HOST_DEVICE inline
  RealVect& min (const RealVect& p) noexcept;

  ///
  /**
     Returns the RealVect that is the component-wise minimum of two
     argument RealVects.
  */
  friend AMREX_GPU_HOST_DEVICE inline
  RealVect min (const RealVect& p1, const RealVect& p2) noexcept;

  ///
  /**
     Modifies this RealVect by taking component-wise max with RealVect
     argument.
  */
  AMREX_GPU_HOST_DEVICE inline
  RealVect& max (const RealVect& p) noexcept;

  ///
  /**
     Returns the RealVect that is the component-wise maximum of two
     argument RealVects.
  */
  friend AMREX_GPU_HOST_DEVICE inline
  RealVect max (const RealVect& p1, const RealVect& p2) noexcept;

  /*@}*/

  /**
     \name Unary operators
  */
  /*@{*/

  ///
  /**
     Unary plus -- for completeness.
  */
  AMREX_GPU_HOST_DEVICE inline
  RealVect operator+ () const noexcept;

  ///
  /**
     Unary minus -- negates all components of this RealVect.
  */
  AMREX_GPU_HOST_DEVICE inline
  RealVect operator- () const noexcept;

  ///
  /**
     Sum of all components of this RealVect.
  */
  AMREX_GPU_HOST_DEVICE inline
  Real sum () const noexcept;

  ///
  /**
     sqrt(sum squares)
  */
  AMREX_GPU_HOST_DEVICE inline
  Real vectorLength() const noexcept;

  ///
  /**
     sum squares--no square root
  */
  AMREX_GPU_HOST_DEVICE inline
  Real radSquared() const noexcept;

  ///
  /**
     Product of all components of this RealVect.
  */
  AMREX_GPU_HOST_DEVICE inline
  Real product () const noexcept;

  ///
  /**
     Component with the minimum value of this RealVect (returns 0 if they are all the same).
     a_doAbs : if true then take the absolute value before comparing
  */
  AMREX_GPU_HOST_DEVICE inline
  int minDir(const bool& a_doAbs) const noexcept;

  ///
  /**
     Component with the maximum value of this RealVect (returns 0 if they are all the same).
     a_doAbs : if true then take the absolute value before comparing
  */
  AMREX_GPU_HOST_DEVICE inline
  int maxDir(const bool& a_doAbs) const noexcept;

  /*@}*/

  /**
     \name Data pointer functions
  */
  /*@{*/

  ///
  /**
     Only for sending stuff to Fortran
   */
  AMREX_GPU_HOST_DEVICE
  const Real* dataPtr() const noexcept { return vect; }

  ///
  /**
     Only for sending stuff to Fortran
   */
  AMREX_GPU_HOST_DEVICE
  Real* dataPtr() noexcept { return vect; }

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
  friend AMREX_GPU_HOST_DEVICE
  RealVect BASISREALV(int dir) noexcept;

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
     \name Static public members functions
  */
  /*@{*/

  ///
  /**
     The zero vector in AMREX_SPACEDIM-dimensional space.
  */
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  static constexpr RealVect TheZeroVector () noexcept {
      return RealVect(0.);
  }

  /**
     The unit vector in AMREX_SPACEDIM-dimensional space.
  */
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  static constexpr RealVect TheUnitVector () noexcept {
      return RealVect(1.);
  }

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
  friend AMREX_GPU_HOST_DEVICE
  RealVect operator+ (Real s, const RealVect& p) noexcept;

  ///
  /**
     Returns <i>s - p</i>.
  */
  friend AMREX_GPU_HOST_DEVICE
  RealVect operator- (Real s, const RealVect& p) noexcept;

  ///
  /**
     Returns a RealVect that is a RealVect <i>p</i> with each component
     multiplied by a scalar <i>s</i>.
  */
  friend AMREX_GPU_HOST_DEVICE
  RealVect operator* (Real s, const RealVect& p) noexcept;
  ///
  /**
     Returns a RealVect that is a RealVect <i>p</i> with each component
     divided by a scalar <i>s</i>.
  */
  friend AMREX_GPU_HOST_DEVICE
  RealVect operator/ (Real s, const RealVect& p) noexcept;

  ///
  /**
     Returns component-wise sum of RealVects <i>s</i> and <i>p</i>.
  */
  friend AMREX_GPU_HOST_DEVICE
  RealVect operator+ (const RealVect& s, const RealVect& p) noexcept;

  ///
  /**
     Returns <i>s - p</i>.
  */
  friend AMREX_GPU_HOST_DEVICE
  RealVect operator- (const RealVect& s, const RealVect& p) noexcept;

  ///
  /**
     Returns component-wise product of <i>s</i> and <i>p</i>.
  */
  friend AMREX_GPU_HOST_DEVICE
  RealVect operator* (const RealVect& s, const RealVect& p) noexcept;
  ///
  /**
     Returns component-wise quotient <i>p / s</i>.
  */
  friend AMREX_GPU_HOST_DEVICE
  RealVect operator/ (const RealVect& s, const RealVect& p) noexcept;

  ///
  /**
     Returns a RealVect obtained by multiplying each of the components
     of the given RealVect by a scalar.
  */
  friend inline AMREX_GPU_HOST_DEVICE
  RealVect scale (const RealVect& p, Real s) noexcept;

  /*@}*/

  ///
  /**
     Print to the given output stream in ASCII.
  */
  friend std::ostream& operator<< (std::ostream& ostr, const RealVect& p);

  friend std::istream& operator>> (std::istream& is, RealVect& iv);

protected:

  /**
     The individual components of this RealVect.
   */
  Real vect[SpaceDim];

};

AMREX_GPU_HOST_DEVICE
inline
RealVect&
RealVect::operator-= (Real s) noexcept
{
  AMREX_D_EXPR(vect[0] -= s, vect[1] -= s, vect[2] -= s);
  return *this;
}

AMREX_GPU_HOST_DEVICE
inline
RealVect&
RealVect::operator*= (Real s) noexcept
{
  AMREX_D_EXPR(vect[0] *= s, vect[1] *= s, vect[2] *= s);
  return *this;
}

AMREX_GPU_HOST_DEVICE
inline
RealVect&
RealVect::operator-= (const RealVect& p) noexcept
{
  AMREX_D_EXPR(vect[0] -= p[0], vect[1] -= p[1], vect[2] -= p[2]);

  return *this;
}

AMREX_GPU_HOST_DEVICE
inline
RealVect
RealVect::operator+ () const noexcept
{
  return RealVect(*this);
}

AMREX_GPU_HOST_DEVICE
inline
RealVect
RealVect::operator- () const noexcept
{
  return RealVect(AMREX_D_DECL(-vect[0], -vect[1], -vect[2]));
}

AMREX_GPU_HOST_DEVICE
inline
RealVect&
RealVect::scale (Real s) noexcept
{
  AMREX_D_EXPR(vect[0] *= s, vect[1] *= s, vect[2] *= s);
  return *this;
}

AMREX_GPU_HOST_DEVICE
inline
IntVect
RealVect::floor () const noexcept
{
  return IntVect(AMREX_D_DECL(static_cast<int>(amrex::Math::floor(vect[0])),
                              static_cast<int>(amrex::Math::floor(vect[1])),
                              static_cast<int>(amrex::Math::floor(vect[2]))));
}

AMREX_GPU_HOST_DEVICE
inline
IntVect
RealVect::ceil () const noexcept
{
  return IntVect(AMREX_D_DECL(static_cast<int>(amrex::Math::ceil(vect[0])),
                              static_cast<int>(amrex::Math::ceil(vect[1])),
                              static_cast<int>(amrex::Math::ceil(vect[2]))));
}

AMREX_GPU_HOST_DEVICE
inline
IntVect
RealVect::round () const noexcept
{
  return IntVect(AMREX_D_DECL(static_cast<int>(amrex::Math::round(vect[0])),
                              static_cast<int>(amrex::Math::round(vect[1])),
                              static_cast<int>(amrex::Math::round(vect[2]))));
}

AMREX_GPU_HOST_DEVICE
inline
Real
RealVect::sum () const noexcept
{
  return AMREX_D_TERM(vect[0], + vect[1], + vect[2]);
}

AMREX_GPU_HOST_DEVICE
inline
Real
RealVect::vectorLength () const noexcept
{
  Real len = this->radSquared();
  len = std::sqrt(len);
  return len;
}

AMREX_GPU_HOST_DEVICE
inline
Real
RealVect::radSquared() const noexcept
{
  Real len = 0;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      len = len + vect[idir]*vect[idir];
    }

  return len;
}

AMREX_GPU_HOST_DEVICE
inline
Real
RealVect::product () const noexcept
{
  return AMREX_D_TERM(vect[0], * vect[1], * vect[2]);
}

AMREX_GPU_HOST_DEVICE
inline
RealVect
scale (const RealVect& p,
       Real            s) noexcept
{
  return RealVect(AMREX_D_DECL(s * p[0], s * p[1], s * p[2]));
}

AMREX_GPU_HOST_DEVICE
inline
bool
RealVect::operator< (const RealVect& p) const noexcept
{
  return AMREX_D_TERM(vect[0] < p[0], && vect[1] < p[1], && vect[2] < p[2]);
}

AMREX_GPU_HOST_DEVICE
inline
bool
RealVect::operator<= (const RealVect& p) const noexcept
{
  return AMREX_D_TERM(vect[0] <= p[0], && vect[1] <= p[1], && vect[2] <= p[2]);
}

AMREX_GPU_HOST_DEVICE
inline
bool
RealVect::operator> (const RealVect& p) const noexcept
{
  return AMREX_D_TERM(vect[0] > p[0], && vect[1] > p[1], && vect[2] > p[2]);
}

AMREX_GPU_HOST_DEVICE
inline
bool
RealVect::operator>= (const RealVect& p) const noexcept
{
  return AMREX_D_TERM(vect[0] >= p[0], && vect[1] >= p[1], && vect[2] >= p[2]);
}

AMREX_GPU_HOST_DEVICE
inline
RealVect&
RealVect::min (const RealVect& p) noexcept
{
    AMREX_D_EXPR(vect[0] = amrex::min(vect[0], p.vect[0]),
                 vect[1] = amrex::min(vect[1], p.vect[1]),
                 vect[2] = amrex::min(vect[2], p.vect[2]));
    return *this;
}

AMREX_GPU_HOST_DEVICE
inline
RealVect&
RealVect::max (const RealVect& p) noexcept
{
    AMREX_D_EXPR(vect[0] = amrex::max(vect[0], p.vect[0]),
                 vect[1] = amrex::max(vect[1], p.vect[1]),
                 vect[2] = amrex::max(vect[2], p.vect[2]));
    return *this;
}

AMREX_GPU_HOST_DEVICE
inline
RealVect
min (const RealVect& p1,
     const RealVect& p2) noexcept
{
  RealVect p(p1);
  return p.min(p2);
}

AMREX_GPU_HOST_DEVICE
inline
RealVect
max (const RealVect& p1,
     const RealVect& p2) noexcept
{
  RealVect p(p1);
  return p.max(p2);
}

AMREX_GPU_HOST_DEVICE inline
Real RealVect::dotProduct(const RealVect& a_rhs) const noexcept
{
    return AMREX_D_TERM(vect[0]*a_rhs.vect[0], +
                        vect[1]*a_rhs.vect[1], +
                        vect[2]*a_rhs.vect[2]);
}

#if (AMREX_SPACEDIM == 3)
AMREX_GPU_HOST_DEVICE inline
  RealVect RealVect::crossProduct(const RealVect& a_rhs) const noexcept
  {
    RealVect tmp(vect[1]*a_rhs[2] - vect[2]*a_rhs[1],
                 vect[2]*a_rhs[0] - vect[0]*a_rhs[2],
                 vect[0]*a_rhs[1] - vect[1]*a_rhs[0]);
    return tmp;
  }
#endif

AMREX_GPU_HOST_DEVICE inline
RealVect&
RealVect::operator+= (Real s) noexcept
{
    AMREX_D_EXPR(vect[0] += s, vect[1] += s, vect[2] += s);
    return *this;
}

AMREX_GPU_HOST_DEVICE inline
RealVect&
RealVect::operator+= (const RealVect& p) noexcept
{
    AMREX_D_EXPR(vect[0] += p[0], vect[1] += p[1], vect[2] += p[2]);
    return *this;
}

AMREX_GPU_HOST_DEVICE inline
RealVect&
RealVect::operator*= (const RealVect &p) noexcept
{
    AMREX_D_EXPR(vect[0] *= p[0], vect[1] *= p[1], vect[2] *= p[2]);
    return *this;
}

AMREX_GPU_HOST_DEVICE inline
RealVect
RealVect::operator* (Real s) const noexcept
{
    RealVect v(AMREX_D_DECL(vect[0]*s, vect[1]*s, vect[2]*s));
    return v;
}

AMREX_GPU_HOST_DEVICE inline
RealVect
RealVect::operator- (Real s) const noexcept
{
    RealVect v(AMREX_D_DECL(vect[0]-s, vect[1]-s, vect[2]-s));
    return v;
}

AMREX_GPU_HOST_DEVICE inline
RealVect
RealVect::operator+ (Real s) const noexcept
{
    RealVect v(AMREX_D_DECL(vect[0]+s, vect[1]+s, vect[2]+s));
    return v;
}

AMREX_GPU_HOST_DEVICE inline
RealVect&
RealVect::operator/= (Real s) noexcept
{
    AMREX_D_EXPR(vect[0] /= s, vect[1] /= s, vect[2] /= s);
    return *this;
}

AMREX_GPU_HOST_DEVICE inline
RealVect&
RealVect::operator/= (const RealVect& p) noexcept
{
    AMREX_D_EXPR(vect[0] /= p[0], vect[1] /= p[1], vect[2] /= p[2]);
    return *this;
}

AMREX_GPU_HOST_DEVICE inline
RealVect
RealVect::operator/ (Real s) const noexcept
{
    RealVect result( AMREX_D_DECL( vect[0] / s, vect[1] / s, vect[2] / s));
    return result ;
}

AMREX_GPU_HOST_DEVICE inline
int
RealVect::minDir(const bool& a_doAbs) const noexcept
{
    int mDir = 0;
    for (int idir=0; idir<SpaceDim; idir++)
    {
        if (a_doAbs)
        {
            if (amrex::Math::abs(vect[idir]) < amrex::Math::abs(vect[mDir]))
            {
                mDir = idir;
            }
        }
        else
        {
            if (vect[idir] < vect[mDir])
            {
                mDir = idir;
            }
        }
    }
    return mDir;
}

AMREX_GPU_HOST_DEVICE inline
int
RealVect::maxDir(const bool& a_doAbs) const noexcept
{
    int mDir = 0;
    for (int idir=0; idir<SpaceDim; idir++)
    {
        if (a_doAbs)
        {
            if (amrex::Math::abs(vect[idir]) > amrex::Math::abs(vect[mDir]))
            {
                mDir = idir;
            }
        }
        else
        {
            if (vect[idir] > vect[mDir])
            {
                mDir = idir;
            }
        }
    }
    return mDir;
}

AMREX_GPU_HOST_DEVICE inline
RealVect
BASISREALV (int dir) noexcept
{
    AMREX_ASSERT(dir >= 0 && dir < SpaceDim);
    RealVect tmp(AMREX_D_DECL(0.,0.,0.));
    tmp.vect[dir] = 1;
    return tmp;
}

AMREX_GPU_HOST_DEVICE inline
RealVect
operator/ (Real            s,
           const RealVect& p) noexcept
{
    return RealVect(AMREX_D_DECL(s/p[0], s/p[1], s/p[2]));
}

AMREX_GPU_HOST_DEVICE inline
RealVect
operator+ (Real            s,
           const RealVect& p) noexcept
{
    return RealVect(AMREX_D_DECL(p[0] + s, p[1] + s, p[2] + s));
}

AMREX_GPU_HOST_DEVICE inline
RealVect
operator- (Real            s,
           const RealVect& p) noexcept
{
    return RealVect(AMREX_D_DECL(s - p[0], s - p[1], s - p[2]));
}

AMREX_GPU_HOST_DEVICE inline
RealVect
operator* (Real            s,
           const RealVect& p) noexcept
{
    return RealVect(AMREX_D_DECL(s * p[0], s * p[1], s * p[2]));
}

AMREX_GPU_HOST_DEVICE inline
RealVect
operator/ (const RealVect& s,
           const RealVect& p) noexcept
{
    return RealVect(AMREX_D_DECL(s[0] / p[0], s[1] /p[1], s[2] / p[2]));
}

AMREX_GPU_HOST_DEVICE inline
RealVect
operator+ (const RealVect& s,
           const RealVect& p) noexcept
{
    return RealVect(AMREX_D_DECL(p[0] + s[0], p[1] +s[1], p[2] + s[2]));
}

AMREX_GPU_HOST_DEVICE inline
RealVect
operator- (const RealVect& s,
           const RealVect& p) noexcept
{
    return RealVect(AMREX_D_DECL(s[0] - p[0], s[1] - p[1], s[2] - p[2]));
}

AMREX_GPU_HOST_DEVICE inline
RealVect
operator* (const RealVect& s,
           const RealVect& p) noexcept
{
    return RealVect(AMREX_D_DECL(p[0] * s[0], p[1] *s[1], p[2] * s[2]));
}

}

#endif
