#ifndef _AnisotropicPLANEIF_H_
#define _AnisotropicPLANEIF_H_

#include "AMReX_RealVect.H"
#include "AMReX_BaseIF.H"
#include "AMReX_PlaneIF.H"

namespace amrex
{
  ///
  /**
     This implicit function specifies a plane.
  */
  class AnisotropicDxPlaneIF: public PlaneIF
  {
  public:
    ///
    /**
       Constructor specifying plane normal (a_normal), a point on the plance
       (a_point), and whether the domain is the same side as the normal
       (a_inside).
    */
    AnisotropicDxPlaneIF(const RealVect& a_normal,
                       const RealVect& a_point,
                       const bool&     a_inside,
                       const RealVect& a_dxVector)
      :PlaneIF(a_normal, a_point, a_inside)
      {
        m_dxVector = a_dxVector;
      }

    /// Destructor
    virtual ~AnisotropicDxPlaneIF()
    {
    }

    ///
    /**
       Return the value of the function at a_point.  When delineating a domain,
       the level set value=0 represents the boundary and value<0 is inside the
       fluid.
    */
    virtual Real value(const RealVect& a_point) const
      {
        RealVect scaledPt;
        for(int idir = 0; idir < SpaceDim; idir++)
        {
          Real ratio = m_dxVector[idir]/m_dxVector[0];
//          Real ratio = m_dxVector[0]/m_dxVector[idir];
          scaledPt[idir] = ratio*a_point[idir];
        }
        Real retval = PlaneIF::value(scaledPt);
        return retval;
      }


  BaseIF* 
  newImplicitFunction() const
  {
    AnisotropicDxPlaneIF* newPtr = new AnisotropicDxPlaneIF(m_normal,
                                                            m_point,
                                                            m_inside,
                                                            m_dxVector);


    return static_cast<BaseIF*>(newPtr);
  }
  protected:
    RealVect m_dxVector;
  };
}

#endif
