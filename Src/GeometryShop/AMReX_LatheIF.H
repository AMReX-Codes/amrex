#ifndef _LATHEIF_H_
#define _LATHEIF_H_

#include "AMReX_RealVect.H"
#include "AMReX_BaseIF.H"

namespace amrex
{


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
    
    /// Copy constructor

    
    /// Destructor
    virtual ~LatheIF();
    
    ///
    /**
       Return the value of the function at a_point.
    */
    virtual Real value(const RealVect& a_point) const;
    
    
    virtual BaseIF* newImplicitFunction() const;
  protected:
    BaseIF* m_impFunc1; // implicit function to rotate
    
    RealVect              m_point;    // point for center of second rotation
    
    bool                  m_inside;   // inside flag
    
  private:
    LatheIF(const LatheIF& a_inputIF);    
    LatheIF();
    void operator=(const LatheIF& a_inputIF);

  };
}
    
#endif
