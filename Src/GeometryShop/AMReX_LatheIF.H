#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _LATHEIF_H_
#define _LATHEIF_H_

#include "MayDay.H"
#include "RealVect.H"
#include "Vector.H"

#include "BaseIF.H"

#include "NamespaceHeader.H"

///
/**
    This implicit function takes one or two implicit functions and uses them
    to produce a generalized surface of revolution.

    In 3D, if one function is given, then it is restricted to x >= 0 in the
    x-y plane (z = 0), revolved around the y axis, and the y axis is rotated
    to the z axis.  In 2D, the function is restricted to x >= 0 and the x axis
    (y = 0) and revolved around the origin.

    In 3D, if two functions are given then the first is restricted to x >= 0
    in the x-y plane (z = 0), revolved around the y axis, and the y axis is
    rotated to the z axis.  At each angle between 0 and 2*Pi, the second
    function restricted to the x axis (y = z = 0) is used to get an angle to
    rotate the first function in the x-y plane with respect to a given point
    before it is rotated about the y axis.  In 2D, the first function is used
    as above and the second function is ignored.

    Note: the second function (if given) should have the same value at x = 0
    and x = 2*Pi (y = z = 0) or a values that differs by some integer multiple
    of 2*Pi.
 */
class LatheIF: public BaseIF
{
public:
  ///
  /**
      Constructor specifying one implicit function to be rotated and whether
      the domain is on the inside (a_inside), i.e. where the function is
      negative.
   */
  LatheIF(const BaseIF& a_impFunc1,
          const bool&   a_inside);

  ///
  /**
      Constructor specifying two implicit functions, one to be rotated
      and the other to cause an additional rotation about a_point, and
      whether the domain is on the inside (a_inside), i.e. where the
      function is negative.
   */
  LatheIF(const BaseIF&   a_impFunc1,
          const BaseIF&   a_impFunc2,
          const RealVect& a_point,
          const bool&     a_inside);

  /// Copy constructor
  LatheIF(const LatheIF& a_inputIF);

  /// Destructor
  virtual ~LatheIF();

  ///
  /**
      Return the value of the function at a_point.
   */
  virtual Real value(const RealVect& a_point) const;

  virtual BaseIF* newImplicitFunction() const;

  virtual bool fastIntersection(const RealVect& a_low,
                                const RealVect& a_high)const
  {
    return false;
// #if CH_SPACEDIM == 3
//     if (m_impFunc2 == NULL)
//       {
//         return true;
//       }
//     else
//       {
//         return false;
//       }
// #elif CH_SPACEDIM == 2
//     return false;
// #else
//     MayDay::Abort("need higher dim in LatheIF\n");
// #endif
  }

  virtual GeometryService::InOut InsideOutside(const RealVect& a_low,
                                               const RealVect& a_high) const ;

  ///
  /**
     Pass this call onto the IFs contained in this IF class.
  */
  virtual void boxLayoutChanged(const DisjointBoxLayout & a_newBoxLayout,
                                const RealVect          & a_dx)
  {
    m_impFunc1->boxLayoutChanged(a_newBoxLayout,a_dx);

    if (m_impFunc2 != NULL)
    {
      m_impFunc2->boxLayoutChanged(a_newBoxLayout,a_dx);
    }
  }

protected:
  BaseIF* m_impFunc1; // implicit function to rotate
  BaseIF* m_impFunc2; // implicit function to use for rotation

  RealVect              m_point;    // point for center of second rotation

  bool                  m_inside;   // inside flag

private:
  LatheIF()
  {
    MayDay::Abort("LatheIF uses strong construction");
  }

  void operator=(const LatheIF& a_inputIF)
  {
    MayDay::Abort("LatheIF doesn't allow assignment");
  }
};

#include "NamespaceFooter.H"
#endif
