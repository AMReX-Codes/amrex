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
