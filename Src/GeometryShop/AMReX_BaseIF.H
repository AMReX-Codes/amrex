#ifndef _BASEIF_H_
#define _BASEIF_H_

#include "AMReX_RealVect.H"
#include "AMReX_Box.H"
#include "AMReX_REAL.H"
#include "AMReX_SPACE.H"


///
/**
   This is the base class for an implicit function specification of geometry.
   All that is needed is a constructor/destructor, a method to give the value
   of the function at any point in space (in 2D or 3D), and a factory method
   (these will probably all be very similar).
*/
namespace amrex
{

  class BaseIF
  {

  public:

    /// Default constructor
    BaseIF()
    {
    }

    /// Default destructor
    virtual ~BaseIF()
      {
      }
  
    ///return the partial derivative at the point
    virtual Real derivative(const  IntVect& a_deriv,
                            const RealVect& a_point) const
      {
        amrex::Error("derivative function not implemented for this implicit function");
        return 0.;
      }

    ///
    virtual Real value(const IntVect  & a_deriv,
                       const RealVect & a_point) const
      {
        Real retval= 0;
        if (a_deriv == IntVect::Zero)
        {
          retval = value(a_point);
        }
        else
        {
          retval = derivative(a_deriv, a_point);
        }
        return retval;
      }
    ///
    /**
       Return the value of the function at a_point.  When delineating a domain,
       the level set value=0 represents the boundary and value<0 is inside the
       fluid.
    */
    virtual Real value(const RealVect& a_point) const = 0;

    ///
    /**
       Return a newly allocated derived class.  The responsibility
       for deleting the memory is left to the calling function.
    */
    virtual BaseIF* newImplicitFunction() const = 0;


    static void corners(const Box     & a_region, 
                        const RealVect& a_origin, 
                        const Real    & a_dx,
                        RealVect      & a_lo, 
                        RealVect      & a_hi)
    {
      a_lo = a_origin + RealVect(a_region.smallEnd())*a_dx;
      a_hi = a_origin + RealVect(a_region.bigEnd()+IntVect::TheUnitVector())*a_dx;
    }

  };

}


#endif
