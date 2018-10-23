#ifndef AMREX_PLANEIF_H_
#define AMREX_PLANEIF_H_

#include "AMReX_RealVect.H"
#include "AMReX_BaseIF.H"

namespace amrex
{
  ///
  /**
     This implicit function specifies a plane.
  */
  class PlaneIF: public BaseIF
  {
  public:
    ///
    /**
       Constructor specifying plane normal (a_normal), a point on the plance
       (a_point), and whether the domain is the same side as the normal
       (a_inside).
    */
    PlaneIF(const RealVect& a_normal,
            const RealVect& a_point,
            const bool&     a_inside);

    /// Destructor
    virtual ~PlaneIF()
    {
    }

    ///return the partial derivative at the point
    virtual Real derivative(const  IntVect& a_deriv,
                            const RealVect& a_point) const;

    ///
    /**
       Return the value of the function at a_point.  When delineating a domain,
       the level set value=0 represents the boundary and value<0 is inside the
       fluid.
    */
    virtual Real value(const RealVect& a_point) const;


    ///
    /**
       Return a newly allocated derived class.  The responsibility
       for deleting the memory is left to the calling function.
    */
    virtual BaseIF* newImplicitFunction() const;

  protected:
    RealVect m_normal;
    RealVect m_point;
    bool     m_inside;
  };
}

#endif
